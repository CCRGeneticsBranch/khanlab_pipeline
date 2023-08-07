import os
import io
import sys
import traceback
import json

SAMPLES = []
TARGETS = []
samples = config['samples']
data_dir = config['data_dir']
work_dir = config['work_dir']
pipeline_type = config['type']
pipeline_home = config['pipeline_home']
chgrp_group = config["chgrp_group"]
shell.prefix("""
    set -e -o pipefail
    module purge
    sleep 20s
    MEM=`echo "${{SLURM_MEM_PER_NODE}} / 1024 "|bc`
    LOCAL="/lscratch/${{SLURM_JOBID}}/"
    THREADS=${{SLURM_CPUS_ON_NODE}}
    """)
configfile: pipeline_home +"/config/common.yml"
pipeline_version = config["pipeline_version"]
emails = config["emails"]
config["pipeline_home"] = pipeline_home
config["work_dir"] = work_dir
suffix_R1 = config["FASTQ_suffix_R1"]
suffix_R2 = config["FASTQ_suffix_R2"]
suffix_SE = config["FASTQ_suffix_SE"]

configfile: pipeline_home +"/config/" + pipeline_type + ".yml"

include: "rules/utility.smk"
include: "rules/pipeline." + pipeline_type + ".smk"

onerror:
    shell("echo ' {pipeline_type} pipeline version {pipeline_version} failed on Biowulf. Samples: {SAMPLES}. Working Dir:  {work_dir}' |mutt -e 'my_hdr From:chouh@nih.gov' -s 'Khanlab {pipeline_type} Pipeline Status' `whoami`@mail.nih.gov {emails} ")
onstart:
    shell("echo ' {pipeline_type} pipeline version {pipeline_version} started on Biowulf. Samples: {SAMPLES}. Working Dir:  {work_dir}' |mutt -e 'my_hdr From:chouh@nih.gov' -s 'Khanlab {pipeline_type} Pipeline Status' `whoami`@mail.nih.gov {emails} ")
onsuccess:
    shell("echo ' {pipeline_type} pipeline version {pipeline_version} finished on Biowulf. Samples: {SAMPLES}. Working Dir:  {work_dir}' |mutt -e 'my_hdr From:chouh@nih.gov' -s 'Khanlab {pipeline_type} Pipeline Status' `whoami`@mail.nih.gov {emails} ")
    shell("for s in {SAMPLES};do touch {work_dir}/${{s}}/successful.txt;chgrp -R {chgrp_group} {work_dir}/${{s}};done")
    print("Workflow finished, no error")

rule all:
    input: TARGETS
