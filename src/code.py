#!/usr/bin/env python
#
# This file is a part of pindel.
# Copyright (C) 2013 DNAnexus, Inc.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

import os, subprocess, time, datetime, re, multiprocessing
import dxpy

def ExportVCF(kwargs, output_path, ref_fn):
    ref_name_version = dxpy.describe(kwargs["reference_fasta"])["name"]
    ref_name_version = ref_name_version.rstrip(".fa")
    vcf_out_fn = kwargs["output_prefix"] + '.pindel.vcf'
    
    command_args = ["pindel2vcf"]
    command_args.append("-r {input}".format(input=ref_fn))
    command_args.append("-P {input}".format(input=output_path))
    command_args.append("-v {input}".format(input=vcf_out_fn))
    if kwargs["vcf_gatk_compatible"]:
        command_args.append("-G")  

    if "export_vcf_advanced_options" in kwargs: 
        command_args.append(kwargs["export_vcf_advanced_options"])
    else: 
        ref_date = str(datetime.date.today())
        command_args.append("-R {input}".format(input=ref_name_version))
        command_args.append("-d ''")

    try:
        vcf_command = " ".join(command_args)
        print "Executing: " + vcf_command
        print subprocess.check_output(vcf_command, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError, e: 
        print e
        print e.output
        raise dxpy.AppError("APP ERROR: App was not able to convert pindel to vcf. Please check pindel2vcf inputs")

    vcf_dxlink = dxpy.dxlink(dxpy.upload_local_file(vcf_out_fn))
    return vcf_dxlink

def RenumberMergedOutput(filename, out_fn):
    print "\nRenumbering variants in merged file"
    awk_command = "awk 'BEGIN {n=0; OFS= \"\\t\"; FS=\"\\t\"} /^[0-9]+\\t[A-Z]+/ {$1 = n; n++} {print $0}' %s > %s"%(filename, out_fn)
    print awk_command
    subprocess.check_call(awk_command, shell=True)
    return out_fn

def SplitBamForSubjobs(kwargs, bam_names, bam_config_fn=None):
    num_threads = kwargs["num_threads_per_instance"]
    print "\nSplitting bam for subjobs"
    
    # Assuming that all bam files have the same chromosomes (is this safe?)
    subprocess.check_output("samtools view -H {input_bam} > header.txt".format(input_bam=bam_names[0]),
                                    shell=True)
    with open('header.txt') as fh: 
        header = [line.rstrip('\n') for line in fh]
    print "Input header: "
    for line in header: 
        print line
    
    print "Save unmapped reads as bam files to merge into subjob files"
    unmapped = {}
    for bam in bam_names:
        fn = bam.rstrip('.bam')+'_unmapped'
        command = "samtools view -@ {n} -u -b -f 4 {bam} > {unmapped}".format(n=num_threads, bam=bam, unmapped=fn)
        print command
        subprocess.check_call(command, shell=True)
        unmapped[bam] = fn
    
    groups = SplitGenomeFromSam(header, kwargs["num_instances"])
    subjobs = []
    subjob_no = 0
    for group in groups:
        group = " ".join(group)
        subjob_bam_fn = [] 

        for bam in bam_names:
            start_time = time.time()
            print "\nMerging {bam} with unmapped reads for pindel subjobs".format(bam=bam)
            out_fn = bam.rstrip('.bam') + '_' + str(subjob_no) + '.bam'
            
            command = "samtools view -@ {n} -bh {bam} {group} > tmp.bam".format(n=num_threads, bam=bam, group=group)
            subprocess.check_call(command, shell=True)
            
            split_command = "samtools merge -@ {n} {out} {unmapped} tmp.bam ".format(n=num_threads, out=out_fn, unmapped=unmapped[bam])
            print split_command
            subprocess.check_call(split_command, shell=True)
            
            print "Samtools view and merge ran in: {min} minutes".format(min=float((time.time()-start_time)/60))
            subjob_bam_fn.append(out_fn)

        subjob_kwargs = kwargs.copy()
        subjob_bam_fn, subjob_bam_idx_fn = IndexBams(bam_names=subjob_bam_fn)
      
        print "Uploading split bam files: " + str(subjob_bam_fn)
        subjob_bam_ids = [dxpy.dxlink(dxpy.upload_local_file(bam)) for bam in subjob_bam_fn]
        print "Uploading split bam index files: " + str(subjob_bam_idx_fn)
        subjob_bam_idx_ids = [dxpy.dxlink(dxpy.upload_local_file(idx)) for idx in subjob_bam_idx_fn]
        
        subjob_kwargs["mappings_files"] = subjob_bam_ids
        subjob_kwargs["bam_index_files"] = subjob_bam_idx_ids
        
        print "Updating bam config file for subjob"
        if bam_config_fn:
            new_config_fn = "subjob_config_" + str(subjob_no) + '.txt'
            with open(bam_config_fn, 'r') as config_fh, open(new_config_fn, 'w') as write_fh:
                for line in config_fh: 
                    line = line.split('\t')
                    bam_name = line[0]
                    out_fn = bam_name.rstrip('.bam') + '_' + str(subjob_no) + '.bam'
                    write_fh.write(out_fn + '\t' + "\t".join(line[1:]) + '\n')      
            
            print "Uploading new config file: " + str(new_config_fn)
            subjob_kwargs["bam_config_file"] = dxpy.dxlink(dxpy.upload_local_file(new_config_fn))
    
        job = dxpy.new_dxjob(subjob_kwargs, "process")
        print "Started subjob #{n}: {job_id}".format(n=subjob_no, job_id=job.get_id())
        subjobs.append(job)
        subjob_no += 1
    
    return subjobs

def SplitGenomeFromSam(header, splits):
    samples = []
    chromosomes = []
    groups = [[]]
    sizes = {}
    totalSize = 0
    
    for line in header:
        if line[:3] == "@SQ":
            parse = re.findall("SN:([^\t\n]*).*LN:([^\t\n]*)", line.strip())
            sizes[parse[0][0]] = int(parse[0][1])
            totalSize += int(parse[0][1])
            chromosomes.append(parse[0][0])
    
    maxSize = totalSize/float(splits)
    currentGroup = 0
    currentSize = 0
    for x in chromosomes:
        if len(groups[currentGroup]) != 0 and sizes[x] + currentSize > maxSize and splits > 0:
            currentGroup += 1
            splits -= 1
            groups.append([])
            maxSize = totalSize/float(splits)
            currentSize = 0
        groups[currentGroup].append(x)
        currentSize += sizes[x]
        totalSize -= sizes[x]
    
    print "\n" + str(groups)
    return groups

def DownloadRefFasta(reference_fasta_id, ref_fn="reference_fasta"):
    dxpy.download_dxfile(reference_fasta_id, ref_fn)
    UnpackInput(ref_fn)
    return ref_fn

def FindFileType(path):
    file_output = subprocess.check_output(["file", path])
    if re.search("gzip compressed data", file_output) != None:  
        return "gzip" 
    elif re.search("bzip2 compressed data", file_output) != None:        
        return "bzip2"    
    elif re.search("XZ compressed data", file_output) != None:        
        return "xz"    
    elif re.search("tar archive", file_output) != None:
        return "tar"
    elif re.search("ASCII", file_output) != None:
        return "text"
    else: 
        return "other"

def UnpackInput(ref_fn="reference_fasta"): 
    while True: 
        if FindFileType(ref_fn) == "text":             
            subprocess.check_call("mv -v {ref_fn} {ref_fn}.fa".format(ref_fn=ref_fn),  shell = True)
            break
        elif FindFileType(ref_fn) == "gzip":
            try: 
                subprocess.check_call("mv -v {ref_fn} {ref_fn}.gz && gunzip -v {ref_fn}.gz".format(ref_fn=ref_fn), shell = True)
            except subprocess.CalledProcessError as e: 
                # code is 2 then it's only a warning and we should continue
                if e.returncode != 2:
                    raise e
        elif FindFileType(ref_fn) == "bzip2":            
            subprocess.check_call("mv -v {ref_fn} {ref_fn}.bz2 && bunzip2 -v {ref_fn}.bz2".format(ref_fn=ref_fn), shell = True)
        elif FindFileType(ref_fn) == "xz":
            subprocess.check_call("mv -v {ref_fn} {ref_fn}.xz && unxz -v {ref_fn}.xz".format(ref_fn=ref_fn), shell = True)
        elif FindFileType(ref_fn) == "tar":        
            subprocess.check_call("mv -v {ref_fn} {ref_fn}.tar && tar xvf {ref_fn}.tar".format(ref_fn=ref_fn) , shell = True)
            break
        else:    
            raise Exception("Unrecognized file type: " + FindFileType(ref_fn))
    subprocess.check_call("cat *.fa > {fa}".format(fa=ref_fn), shell = True)
    print "Sucessfully unpacked FASTA genome into " + ref_fn
    return ref_fn

def DownloadFilesFromArray(input_ids):
    print "\nDownloading {n} files".format(n=len(input_ids))
    if len(input_ids) < 1:
        raise dxpy.AppInternalError("No files were given as input")
    filenames = []
    
    start_time = time.time()
    for id in input_ids:
        fn = dxpy.describe(id)["name"]
        filenames.append(fn)
        dxpy.download_dxfile(dxid=id, filename=fn)
    print "Downloaded {files} in {min} minutes".format(files=sorted(filenames), min=float((time.time()-start_time)/60))
    return sorted(filenames)

def SortBams(bam_names, num_threads):
    print "\nSorting bams: " + str(bam_names)
    sorted_bam_names = []
    for bam in bam_names: 
        sorted_prefix = bam.rstrip('.bam') + '_sorted_by_app' 
        sorted_name = sorted_prefix + '.bam'
        start_time = time.time()
        command = "samtools sort -@ {n} {in_bam} {out_prefix}".format(n=num_threads, in_bam=bam, out_prefix=sorted_prefix)
        print command
        subprocess.check_call(command, shell=True)
        command = "mv {sorted_name} {orig_name}".format(sorted_name=sorted_name, orig_name=bam)
        print command
        subprocess.check_call(command, shell=True)
        print "Sorted {bam} in {min} minutes".format(bam=bam, min=float((time.time()-start_time)/60))
        sorted_bam_names.append(sorted_name)
    return sorted(bam_names)

def IndexBams(bam_names): 
    print "\nIndexing " + str(bam_names)
    bam_idx_names = []
    for bam in bam_names: 
        start_time = time.time()
        command = "samtools index {in_bam}".format(in_bam=bam)
        print command
        subprocess.check_call(command, shell=True)
        print "Indexed {bam} in {min} minutes".format(bam=bam, min=float((time.time()-start_time)/60))
        bam_idx_names.append(bam + '.bai')
    return sorted(bam_names), sorted(bam_idx_names)
        
def CheckBamIdxMatch(bam_names, idx_names):
    print "\nChecking if BAM names and index names match"
    if len(bam_names) != len(idx_names):
        return False
    for i in range(len(bam_names)):
        if bam_names[i].rstrip('.bam') != idx_names[i].rstrip('.bam.bai'):
            print "\tBAM names and index names do not match"
            return False
    print "\tBAM names and index names match"
    return True

def ValidateBamConfig(bam_config_fn, bam_name_array):
    print "\nValidating bam config file"
    with open(bam_config_fn) as config_fh:
        for line in config_fh:
            name = line.split()[0]
            if name not in bam_name_array:
                raise dxpy.AppError("Bam config file contains filenames which do not match input bam files")
    print "\tBam config file is valid"
    return True

def WriteConfigFile(mappings_names, fn, insert_size=0, is_pindel=False):
    print "\nWriting BAM/Pindel config file"
    with open(fn, 'w') as config_fh: 
        for name in mappings_names: 
            if is_pindel: 
                config_fh.write("{fn}\n".format(fn=name))
            else:  
                config_fh.write("{fn}\t{insert}\t{name}\n".format(fn=name, insert=insert_size, name=name.rstrip('.bam')))    
    return fn

def RunSam2Pindel(bam_names, insert_size, seq_platform, num_threads, config_fn):
    print "\nBAM files were not created by BWA. Converting BAM files to Pindel Input Format"
    pindel_files = []
    for bam in bam_names:
        output_name = bam.rstrip('.bam')+'_pindel.txt'
        command = "samtools view -@ {n} {input} | sam2pindel - {output} {insert} {tag} 0 {seq_platform}".format(n=num_threads,
                                                                                                                 input=bam,
                                                                                                          output=output_name,
                                                                                                          insert=insert_size,
                                                                                                          tag=bam.rstrip('.bam'),
                                                                                                          seq_platform=seq_platform)
        print command
        start_time = time.time()
        subprocess.check_call(command, shell=True)
        pindel_files.append(output_name)
        print "Conversion of {bam} took {min} minutes".format(bam=bam, min=float((time.time()-start_time)/60))
    
    pindel_config_fn = WriteConfigFile(mappings_names=pindel_files, fn=config_fn, is_pindel=True)
    return pindel_config_fn

def BuildPindelCommand(kwargs, chrom, input_fn, is_pindel_input_type=False):
    # Download Reference FASTA
    reference_fasta_id = kwargs["reference_fasta"]
    ref_fn = DownloadRefFasta(reference_fasta_id)
    
    create_index = True
    if "fasta_index" in kwargs:
        fasta_idx_id = kwargs["fasta_index"]
        if dxpy.describe(fasta_idx_id)["name"].rstrip(".fa.fai") ==  dxpy.describe(reference_fasta_id)["name"].rstrip(".fa"):
            dxpy.download_dxfile(fasta_idx_id, ref_fn+".fai")
            create_index = False
    if create_index: 
        print "No FASTA index was provided as input. Making one now."
        samtools_command = "samtools faidx {fasta}".format(fasta=ref_fn)
        subprocess.check_call(samtools_command, shell=True)
    
    print "\nBuilding pindel command from app inputs"
    command_args = ["pindel"]    
    output_path = "output/" + kwargs["output_prefix"]
    command_args.append("-o {output_path}".format(output_path=output_path))
    command_args.append("-f {fa}".format(fa=ref_fn))
    
    if is_pindel_input_type: 
        command_args.append("-P {pindel_config}".format(pindel_config=input_fn))
    else:
        command_args.append("-i {bam_config}".format(bam_config=input_fn))
        
    command_args.append("-T {option}".format(option=kwargs["num_threads_per_instance"]))
    command_args.append("-c {chrom}".format(chrom=chrom))
    
    if kwargs["report_only_close_mapped_reads"]:
        command_args.append("-S {option}".format(option=kwargs["report_only_close_mapped_reads"]))
    else:
        command_args.append("-r {option}".format(option=kwargs["report_inversions"]))
        command_args.append("-t {option}".format(option=kwargs["report_duplications"]))
        command_args.append("-l {option}".format(option=kwargs["report_long_insertions"]))
        command_args.append("-k {option}".format(option=kwargs["report_breakpoints"]))
        command_args.append("-s {option}".format(option=kwargs["report_close_mapped_reads"]))     

    if "breakdancer_calls_file" in kwargs:
        breakdancer_fn = DownloadFilesFromArray([kwargs["breakdancer_calls_file"]["$dnanexus_link"]])[0]
        print breakdancer_fn
        command_args.append("-b {option}".format(option=breakdancer_fn))
            
    if "pindel_command_line" in kwargs:
        advanced_command = kwargs["pindel_command_line"]
        if advanced_command.startswith("pindel"):
            advanced_command = advanced_command.replace("pindel", "")
        command_args.append(advanced_command)
 
    command = " ".join(command_args)
    print command
    return command, output_path

def RunPindel(kwargs, pindel_command, output_path):   
    folder = output_path.split("/")[0]
    print "Making folder for output: " + str(folder)
    os.mkdir(folder)
    
    print "Running pindel with: " 
    print '\t' + str(pindel_command)
    start_time = time.time()
    try:
        p = subprocess.check_output(pindel_command, stderr=subprocess.STDOUT, shell=True)
        print p 
        tot_time = time.time() - start_time
        hours = int(tot_time/3600)
        mins = int(float(tot_time%3600)/60)
        secs = tot_time%60    
        print "Pindel ran in: {hrs}h {mins}m {secs}s".format(hrs=hours, mins=mins, secs=secs)    
    except subprocess.CalledProcessError, e:
        print "\n" + str(e.output)
        raise dxpy.AppError("Pindel failed to run. Please check job logs for pindel output. If error is a segmentation fault " +
                            "raised as pindel begins to run, check that reference FASTA file is the same reference used to produce the mappings")
        
    return output_path

def UploadPindelOutputs(kwargs, output_path):
    prefix = kwargs["output_prefix"]
    output_folder = output_path.split("/")[0]
    
    print "\nPindel outputs: "
    print os.listdir(output_folder)
    print subprocess.check_output("wc -l {folder}/*".format(folder=output_folder), shell=True)
    
    suffix = kwargs["variant_suffixes"]
    
    print "Uploading Pindel detected deletions file"
    deletion_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["deletions"], name=prefix+"_"+suffix["deletions"])
    
    print "Uploading Pindel detected short insertions file"
    short_insert_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["short_inserts"], name=prefix+"_"+suffix["short_inserts"])    
    
    print "Uploading Pindel detected inversions file"
    inversion_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["inversions"], name=prefix+"_"+suffix["inversions"])
    
    print "Uploading Pindel detected tandem duplications file"
    tandem_duplication_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["tandem_duplications"], name=prefix+"_"+suffix["tandem_duplications"])
    
    print "Uploading Pindel detected large insertions file"
    large_insert_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["large_inserts"], name=prefix+"_"+suffix["large_inserts"])

    print "Uploading Pindel detected unassigned breakpoints file"
    breakpoint_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["breakpoints"], name=prefix+"_"+suffix["breakpoints"])
    
    app_outputs = {"deletions": dxpy.dxlink(deletion_fh), 
                   "short_inserts": dxpy.dxlink(short_insert_fh),
                   "inversions" : dxpy.dxlink(inversion_fh), 
                   "tandem_duplications" : dxpy.dxlink(tandem_duplication_fh), 
                   "large_inserts" : dxpy.dxlink(large_insert_fh),
                   "breakpoints" : dxpy.dxlink(breakpoint_fh) 
                   }    
    if kwargs["report_close_mapped_reads"] or kwargs["report_only_close_mapped_reads"]:
        print "Uploading Pindel detected Close End Mapped file"
        close_mapped_reads_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["close_mapped_reads"], name=prefix+"_"+suffix["close_mapped_reads"])
        app_outputs["close_mapped_reads"] = dxpy.dxlink(close_mapped_reads_fh)
    
    # Currently "pindel -Q breakdancer_output_fn" option causes pindel to segfault
    #if "breakdancer_calls_file" in kwargs:
    #    print "Uploading Confirmed Breakdancer Outputs file"
    #    breakdancer_fh = dxpy.upload_local_file(filename=output_path+"_"+suffix["breakdancer_outputs"], name=prefix+"_"+suffix["breakdancer_outputs"])
    #    app_outputs["breakdancer_outputs"] = dxpy.dxlink(breakdancer_fh)
        
    return app_outputs

@dxpy.entry_point("main")
def main(**kwargs):        
    mappings_ids = kwargs["mappings_files"]
    mappings_names = sorted([dxpy.describe(id)["name"] for id in mappings_ids])
    
    if "num_threads_per_instance" not in kwargs:
        kwargs["num_threads_per_instance"] = multiprocessing.cpu_count()
    if "num_instances" not in kwargs:
        kwargs["num_instances"] = 1

    # Set output prefix here
    if "output_prefix" not in kwargs:
        kwargs["output_prefix"] = mappings_names[0].rstrip('.bam').rstrip('.txt')
    # Set output suffixes (for consistency through app) 
    kwargs["variant_suffixes"] = {"deletions" : 'D',
                                  "short_inserts" : 'SI', 
                                  "tandem_duplications" : 'TD',
                                  "large_inserts" :'LI',
                                  "inversions" : 'INV',
                                  "breakpoints" : 'BP',
                                  #"breakdancer_outputs": 'BD',
                                  "close_mapped_reads": 'CloseEndMapped'} 
    """
    if kwargs["export_vcf"]:
        print "\nTESTING pindel2vcf command line inputs on dummy inputs"
        ExportVCF(kwargs, output_path="/usr/test_vcf/dummy", ref_fn="/usr/test_vcf/dummy.fa")
    """
    # Check if input files have .bam extension
    if mappings_names[0][-4:] == ".bam":
        for name in mappings_names:
            if name[-4:] != ".bam":
                raise dxpy.AppError("Input mappings files are not all bam files with .bam extensions")
        app_outputs = RunWithBamInput(kwargs=kwargs) 
    else:
        app_outputs = RunWithPindelInput(kwargs=kwargs, sam2pindel=False)
        
    return app_outputs

def RunWithPindelInput(kwargs, sam2pindel=False):
    pindel_config_fn="pindel_config.txt"
    mappings_names = DownloadFilesFromArray(input_ids=kwargs["mappings_files"])
    
    if sam2pindel:
        print "\nInput was not produced by BWA/MOSAIK, running sam2pindel script on input BAM files"
        if "sequence_platform" not in kwargs:
            raise dxpy.AppError("If BAM files were not produced by BWA, must ALSO specify which sequence platform was used to produce the mappings")
        pindel_config_fn = RunSam2Pindel(bam_names=mappings_names, insert_size=kwargs["insert_size"], 
                                         seq_platform=kwargs["sequence_platform"], num_threads=kwargs["num_threads_per_instance"], 
                                         config_fn=pindel_config_fn)
    else:
        print "\nInput is pindel input. Making pindel configuration file" 
        pindel_config_fn = WriteConfigFile(mappings_names=mappings_names, fn=pindel_config_fn, is_pindel=True)   

    chrom = kwargs["chromosome"] if "chromosome" in kwargs else "ALL"
    command, output_path = BuildPindelCommand(kwargs=kwargs, chrom=chrom, input_fn=pindel_config_fn, is_pindel_input_type=True)
    output_path = RunPindel(kwargs=kwargs, pindel_command=command, output_path=output_path)

    app_outputs = UploadPindelOutputs(kwargs=kwargs, output_path=output_path)
    if kwargs["export_vcf"]:
        app_outputs["vcf"] = ExportVCF(kwargs=kwargs, output_path=output_path, ref_fn="reference_fasta")   
    
    return app_outputs
    
def RunWithBamInput(kwargs):  
    mappings_ids = kwargs["mappings_files"]
    mappings_names = sorted([dxpy.describe(id)["name"] for id in mappings_ids])
    num_threads = kwargs["num_threads_per_instance"]
    bam_config_fn = "bam_config.txt"
    
    if "bam_config_file" in kwargs:
        print "\nInput has a BAM config file. Need to download and validate bam config file"
        dxpy.download_dxfile(kwargs["bam_config_file"], bam_config_fn)
        ValidateBamConfig(bam_config_fn=bam_config_fn, bam_name_array=mappings_names)
    else:
        if "insert_size" not in kwargs:
            raise dxpy.AppError("Input files are bam files but neither a bam configuration file, nor an insert size was given as an app input.")
        if kwargs["bam_not_produced_by_bwa"]:
            return RunWithPindelInput(kwargs, sam2pindel=True)
        else:
            bam_config_fn = WriteConfigFile(mappings_names=mappings_names, fn=bam_config_fn,  insert_size=kwargs["insert_size"])   
        
    need_to_index=True
    if "bam_index_files" in kwargs:
        bam_idx_ids = kwargs["bam_index_files"]
        idx_names = sorted([dxpy.describe(id)["name"] for id in bam_idx_ids])
        if CheckBamIdxMatch(bam_names=mappings_names, idx_names=idx_names):
            need_to_index = False
            mappings_names = DownloadFilesFromArray(mappings_ids)
            bam_idx_names = DownloadFilesFromArray(bam_idx_ids)
            
    if need_to_index:
        mappings_names = DownloadFilesFromArray(mappings_ids)
        if not kwargs["assume_sorted"]:
            mappings_names = SortBams(bam_names=mappings_names, num_threads=num_threads)
        mappings_names, bam_idx_names = IndexBams(mappings_names) 
    
    chrom = kwargs["chromosome"] if "chromosome" in kwargs else "ALL"
    
    if "chromosome" in kwargs or kwargs["num_instances"] == 1:        
        command, output_path = BuildPindelCommand(kwargs=kwargs, chrom=chrom, input_fn=bam_config_fn, is_pindel_input_type=False)
        output_path = RunPindel(kwargs=kwargs, pindel_command=command, output_path=output_path)    
        app_outputs = UploadPindelOutputs(kwargs=kwargs, output_path=output_path)
        if kwargs["export_vcf"]:
            app_outputs["vcf"] = ExportVCF(kwargs=kwargs, output_path=output_path, ref_fn="reference_fasta")   

    else: 
        subjob_ids = SplitBamForSubjobs(kwargs, mappings_names, bam_config_fn)
        postprocess_inputs = {"subjob_outputs": [job.get_output_ref("subjob_output") for job in subjob_ids], "kwargs": kwargs}
        postprocess_job = dxpy.new_dxjob(fn_input = postprocess_inputs, fn_name = "postprocess")
        
        app_outputs = {"deletions" : {"job": postprocess_job.get_id(), "field": "deletions"},
                       "short_inserts" : {"job": postprocess_job.get_id(), "field": "short_inserts"}, 
                       "tandem_duplications" : {"job": postprocess_job.get_id(), "field": "tandem_duplications"},
                       "large_inserts" : {"job": postprocess_job.get_id(), "field": "large_inserts"},
                       "inversions" : {"job": postprocess_job.get_id(), "field": "inversions"},
                       "breakpoints" : {"job": postprocess_job.get_id(), "field": "breakpoints"}
                       }
        if kwargs["report_close_mapped_reads"] or kwargs["report_only_close_mapped_reads"]:
            app_outputs["close_mapped_reads"] = {"job": postprocess_job.get_id(), "field": "close_mapped_reads"}
        if kwargs["export_vcf"]:
            app_outputs["vcf"] = {"job": postprocess_job.get_id(), "field": "vcf"}   
        #if "breakdancer_calls_file" in kwargs:
        #    app_outputs["breakdancer_outputs"] = {"job": postprocess_job.get_id(), "field": "breakdancer_outputs"}
    
    dxlinks = []
    if need_to_index:
        if not kwargs["assume_sorted"]:
            for bam in mappings_names:
                uploaded_bam = dxpy.upload_local_file(bam, name=bam.rstrip('.bam')+"_sorted.bam")
                dxlinks.append(dxpy.dxlink(uploaded_bam))
        for idx in bam_idx_names:
            uploaded_idx = dxpy.upload_local_file(idx, name=idx.rstrip('.bam.bai')+"_sorted.bam.bai")
            dxlinks.append(dxpy.dxlink(uploaded_idx))
        app_outputs["sortedbam_and_index_files"] = dxlinks
        
    return app_outputs

@dxpy.entry_point("process")
def process(**kwargs):
    bam_names = DownloadFilesFromArray(kwargs["mappings_files"])
    bam_idx_names = DownloadFilesFromArray(kwargs["bam_index_files"])
    
    # Change name of subjob outputs for easy debugging
    kwargs["output_prefix"] = bam_names[0].rstrip('.bam')
    
    bam_config_fn = "bam_config.txt"
    if "bam_config_file" in kwargs:
       dxpy.download_dxfile(kwargs["bam_config_file"], bam_config_fn)


    # Run Pindel on all chromosomes in BAM and upload outputs
    chrom = "ALL"
    command, output_path = BuildPindelCommand(kwargs=kwargs, chrom=chrom, input_fn=bam_config_fn, is_pindel_input_type=False)
    output_path = RunPindel(kwargs=kwargs, pindel_command=command, output_path=output_path)
    app_outputs = UploadPindelOutputs(kwargs, output_path)
    
    return { "subjob_output": app_outputs }

@dxpy.entry_point("postprocess")
def postprocess(**inputs):
    kwargs = inputs["kwargs"]
    subjob_outputs = inputs["subjob_outputs"] 
    print "\nMerging outputs from {n} subjobs".format(n=len(subjob_outputs))

    output_prefix = kwargs["output_prefix"]
    variant_suffixes = kwargs["variant_suffixes"]
    
    app_output_fn = {}
    for subjob_output in subjob_outputs:
        for type, id in subjob_output.iteritems():
            file_id = id["$dnanexus_link"]
            filename = output_prefix + "_" + variant_suffixes[type]
            
            print "Downloading " + str(file_id) + " into " + filename
            dxpy.download_dxfile(dxid=file_id, filename=filename, append=True)
            app_output_fn[type] = filename

    postprocess_outputs = {}
    need_to_renumber = ["deletions", "short_inserts", "tandem_duplications", "inversions", "large_inserts"]
    for type, fn in app_output_fn.iteritems():
        out_fn = fn
        if type in need_to_renumber:
            out_fn = RenumberMergedOutput(fn, fn+"_renumbered") 
        print "\nUploading {file} as {fn}".format(file=out_fn, fn=fn)
        postprocess_outputs[type] = dxpy.dxlink(dxpy.upload_local_file(out_fn, name=fn))

    if kwargs["export_vcf"]:
        DownloadRefFasta(kwargs["reference_fasta"])
        postprocess_outputs["vcf"] = ExportVCF(kwargs=kwargs, output_path=output_prefix, ref_fn="reference_fasta") 
        
    return postprocess_outputs

print dxpy.run()
