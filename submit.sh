#!/usr/bin/bash

module load snakemake/5.24.1

snakemake --directory $work_dir --snakefile $snakefile --configfile $sheet --config type=$type \
	pipeline_home=$pipeline_home work_dir=$work_dir data_dir=$data_dir genome=$genome now=$now \
	--jobname {params.rulename}.{jobid} --nolock  --ri -k -p -r -j 1000 --cores 150 --jobscript $pipeline_home/scripts/jobscript.sh --cluster "sbatch -o {params.log_dir}/{params.rulename}.%j.o -e {params.log_dir}/{params.rulename}.%j.e {params.batch}" \