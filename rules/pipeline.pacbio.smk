from snakemake.utils import R

try:
    #prepare targets
    TARGETS= []    
    genome = config["genome"]
    for sample_id, sample in samples.items():
        if not "SampleFiles" in sample:
            raise Exception('SampleFiles not found in sample sheet')
        SAMPLES.append(sample_id)
        sample_file = data_dir + "/" + sample["SampleFiles"]
        sample["FLBam"] = sample_id + "/" + sample_id + ".flnc.bam"
        if 'Skip_lima' in sample and sample['Skip_lima']:
            sample["FLBam"] = data_dir + "/" + sample["SampleFiles"]
            
        # check source FASTQ files exist
        if not os.path.exists(sample_file):
            raise Exception(sample_file + ' not found')
        TARGETS.append(sample_id + "/" + sample_id + ".sam")
        TARGETS.append(sample_id + "/" + sample_id + ".sorted.sam")
        TARGETS.append(sample_id + "/" + sample_id + ".flnc.mapped.sorted.bam")
        TARGETS.append(sample_id + "/pb.collapsed.group.txt")
        TARGETS.append(sample_id + "/pb.collapsed.rep.fa")
        TARGETS.append(sample_id + "/pb.collapsed.abundance.txt")
        TARGETS.append(sample_id + "/fusion/fusion.annotated.filtered.with_symbol.txt")
        TARGETS.append(sample_id + "/" + sample_id + ".sqanti_sqanti_report.pdf")
        TARGETS.append(sample_id + "/" + sample_id + ".sqanti_classification.filtered_lite.fasta")
        TARGETS.append(sample_id + "/" + sample_id + ".sqanti_corrected.sam")
        TARGETS.append(sample_id + "/" + sample_id + ".sqanti_corrected_sorted.bam.bai")
        #TARGETS.append(sample_id + "/" + sample_id + ".sqanti_classification_genename.txt")
        sample["Sam"] = sample_id + "/" + sample_id + ".sorted.sam"
        sample["TargetSam"] = sample_id + "/" + sample_id + ".sorted.sam"
        sample["Type"] = "transcriptome"
        if 'Ref' in sample:
            config[genome]['ref'] = sample['Ref']
            config[genome]['gtf_gencode'] = sample['GTF']
        #default primers
        if 'Primers' not in sample:
            sample["Primers"] = config["pipeline_home"] + '/config/pacbio/primers.fasta'
        if 'Target' in sample and sample["Target"] != "":
            targeted_sam = sample_id + "/" + sample_id + ".sorted.targeted.sam"
            sample["TargetSam"] = targeted_sam
            sample["Type"] = "targeted"
            TARGETS.append(targeted_sam)
        if not os.path.exists(work_dir + '/' + sample_id):
            os.makedirs(work_dir + '/' + sample_id)
        if not os.path.exists(work_dir + '/' + sample_id + '/log'):
            os.makedirs(work_dir + '/' + sample_id + '/log')
        
except Exception as err:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    output = io.StringIO()
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=output)
    contents = output.getvalue()
    output.close()
    print(contents)    
    shell("echo 'Pacbio pipeline has exception: reason " + contents + ". Working Dir:  {work_dir}' |mutt -e 'my_hdr From:chouh@nih.gov' -s 'Khanlab Pacbio Pipeline Status' `whoami`@mail.nih.gov {emails} ")
    sys.exit()


localrules: all


rule lima:
    input: lambda wildcards: data_dir + '/' + samples[wildcards.sample]["SampleFiles"]
    #output: "{sample}/{sample}.demux.NEB_5p--NEB_Clontech_3p.bam"
    output: "{sample}/{sample}.demux.primer_5p--primer_3p.bam"
    version:
            config["version"]["smrtanalysis"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_lima'],
            primers = lambda wildcards: samples[wildcards.sample]["Primers"],
            prefix = "demux.bam",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "lima"
    shell: 
            """
            module load smrtanalysis/{version}
            mkdir -p {wildcards.sample}
            mkdir -p {wildcards.sample}/log
            lima {input} {params.primers} {wildcards.sample}/{wildcards.sample}.demux.bam --isoseq --peek-guess 2>{wildcards.sample}/lima.err
            """

rule refine:
    #input: "{sample}/{sample}.demux.NEB_5p--NEB_Clontech_3p.bam",
    input: "{sample}/{sample}.demux.primer_5p--primer_3p.bam",
    output: "{sample}/{sample}.flnc.bam"
    version:
            config["version"]["smrtanalysis"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_refine'],
            primers = lambda wildcards: samples[wildcards.sample]["Primers"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "refine"
    shell: 
            """
            module load smrtanalysis/{version}
            isoseq3 refine {input} {params.primers} {output} -j ${{THREADS}}  --require-polya 2>{wildcards.sample}/refine.err
            """

rule bam2fastq:
    input: "{sample}/{sample}.flnc.bam"
    output: "{sample}/{sample}.flnc.fastq"
    version:
            config["version"]["samtools"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_refine'],            
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "bam2fastq"
    shell: 
            """
            module load samtools/{version}
            samtools fastq {input} > {output}
            """

rule map_flnc:
    input: 
            "{sample}/{sample}.flnc.fastq"
    output: 
            sam="{sample}/{sample}.flnc.mapped.sam",
            bam="{sample}/{sample}.flnc.mapped.bam",
            sorted_bam="{sample}/{sample}.flnc.mapped.sorted.bam"
    version:
            config["version"]["minimap2"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_minimap2'],
            bed = lambda wildcards: config[genome]["bed_gencode"],
            ref = lambda wildcards: config[genome]["ref"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "map_flnc"
    shell: 
            """
            module load minimap2/{version}
            minimap2 -ax splice -t ${{THREADS}} -uf --secondary=no -C5 --junc-bed {params.bed} {params.ref} {input} > {output.sam}
            module load samtools
            samtools view -bS {output.sam} > {output.bam}
            samtools sort -o {output.sorted_bam} -@ ${{THREADS}} {output.bam}
            grep -v '^@' {output.sam} | cut -f3 | sort | uniq -c > {wildcards.sample}/{wildcards.sample}.flnc.summary_by_chr.tsv
            samtools index {output.sorted_bam}
            """
            
rule cluster:
    input: lambda wildcards: samples[wildcards.sample]["FLBam"]
    output: 
            fasta = "{sample}/{sample}.clustered.hq.fasta", 
            bam = "{sample}/{sample}.clustered.bam",
            report = "{sample}/{sample}.clustered.cluster_report.csv"
    version:
            config["version"]["smrtanalysis"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_cluster'],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "cluster"
    shell: 
            """
            module load smrtanalysis/{version}
            mkdir -p {wildcards.sample}
            mkdir -p {wildcards.sample}/log
            isoseq3 cluster {input} {output.bam} -j ${{THREADS}} --use-qvs 2>{wildcards.sample}/cluster.err
            gzip -d {output.fasta}.gz
            """
    
rule minimap2:
    input: 
            "{sample}/{sample}.clustered.hq.fasta"
    output: 
            sam="{sample}/{sample}.sam",
            sorted_sam="{sample}/{sample}.sorted.sam"
    version:
            config["version"]["minimap2"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_minimap2'],
            bed = lambda wildcards: config[genome]["bed_gencode"],
            ref = lambda wildcards: config[genome]["ref"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "minimap2"
    shell: 
            """
            module load minimap2/{version}
            minimap2 -ax splice -t ${{THREADS}} -uf --secondary=no -C5 --junc-bed {params.bed} {params.ref} {input} > {output.sam}
            module load samtools
            samtools sort -o {output.sorted_sam} -@ ${{THREADS}} {output.sam}
            grep -v '^@' {output.sorted_sam} | cut -f3 | sort | uniq -c > {wildcards.sample}/{wildcards.sample}.summary_by_chr.tsv
            """

rule targeted:
    input: "{sample}/{sample}.sorted.sam"
    output: "{sample}/{sample}.sorted.targeted.sam"
    version:
            config["version"]["samtools"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_short'],
            target = lambda wildcards: samples[wildcards.sample]["Target"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "targeted"
    shell:
            """
            module load samtools/{version}
            samtools view -bS {input} > {wildcards.sample}/{wildcards.sample}.sorted.bam
            samtools index {wildcards.sample}/{wildcards.sample}.sorted.bam
            samtools view -h {wildcards.sample}/{wildcards.sample}.sorted.bam {params.target} > {output}
            """

rule collapse:
    input:
            fasta = "{sample}/{sample}.clustered.hq.fasta",
            sam = lambda wildcards: samples[wildcards.sample]["TargetSam"]
    output:
            "{sample}/pb.collapsed.group.txt",
            "{sample}/pb.collapsed.rep.fa"
    params:
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_cdna_cupcake'],
            pipeline_home = config["pipeline_home"],
            cdna_cupcake_sif = config["apps"]["cdna_cupcake_sif"],
            work_dir = config["work_dir"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "collapse"
    shell: 
            """
            module load singularity
            export SINGULARITY_BINDPATH="{params.work_dir},/scratch,/lscratch"
            singularity exec {params.cdna_cupcake_sif} collapse_isoforms_by_sam --input {params.work_dir}/{input.fasta} -s {params.work_dir}/{input.sam} --dun-merge-5-shorter -o {params.work_dir}/{wildcards.sample}/pb
            """

rule find_fusion:
    input:
            fasta = "{sample}/{sample}.clustered.hq.fasta",
            sam = lambda wildcards: samples[wildcards.sample]["Sam"],
            csv = "{sample}/{sample}.clustered.cluster_report.csv"
    output:
            gff = "{sample}/fusion/fusion.gff",
            fa = "{sample}/fusion/fusion.rep.fa",
            abundance = "{sample}/fusion/fusion.abundance.txt"
    params:
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_cdna_cupcake'],
            pipeline_home = config["pipeline_home"],
            cdna_cupcake_sif = config["apps"]["cdna_cupcake_sif"],
            conda_home = config["apps"]["conda_home"],
            conda_env = config["apps"]["conda_env"],
            work_dir = config["work_dir"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "find_fusion"
    shell: 
            """
            #module load singularity
            #export SINGULARITY_BINDPATH="{params.work_dir},/scratch,/lscratch"
            #singularity exec {params.cdna_cupcake_sif} fusion_finder --input {params.work_dir}/{input.fasta} -s {params.work_dir}/{input.sam} --cluster_report {params.work_dir}/{input.csv} --dun-merge-5-shorter -o {params.work_dir}/{wildcards.sample}/fusion --min_locus_coverage_bp 500 -d 1000000
            source {params.conda_home}/etc/profile.d/conda.sh
            conda activate {params.conda_env}
            #module load python/3.7

            #conda init bash
            #export PYTHONPATH=$PYTHONPATH:/data/khanlab/projects/hsienchao/apps/cDNA_Cupcake/sequence/
            #export PATH=$PATH:/data/khanlab/projects/hsienchao/apps/utils
            fusion_finder.py --input {params.work_dir}/{input.fasta} -s {params.work_dir}/{input.sam} --cluster_report {params.work_dir}/{input.csv} --dun-merge-5-shorter -o {params.work_dir}/{wildcards.sample}/fusion --min_locus_coverage_bp 500 -d 1000000
            mkdir -p {wildcards.sample}/fusion
            mv {wildcards.sample}/fusion.* {wildcards.sample}/fusion/
            grep -v '^#' {output.abundance} > {output.abundance}.tmp
            mv {output.abundance}.tmp {output.abundance}
            """
            
rule fusion_sqanti:
    input:
            abundance = "{sample}/fusion/fusion.abundance.txt",
            pb_abundance = "{sample}/pb.collapsed.abundance.txt",
            gff = "{sample}/fusion/fusion.gff"
    output:
            annotation="{sample}/fusion/fusion.annotated.filtered.txt",
            annotation_with_symbol="{sample}/fusion/fusion.annotated.filtered.with_symbol.txt",
            classification="{sample}/fusion/sqanti_classification.txt"
    params:
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_cdna_cupcake'],
            sqanti3_path = config["apps"]["sqanti3_path"],
            pipeline_home = config["pipeline_home"],
            gtf = config[genome]["gtf_gencode"],
            ref = config[genome]["ref"],
            ref_path = config[genome]["ref_path"],
            conda_home = config["apps"]["conda_home"],
            conda_env = config["apps"]["conda_env"],
            work_dir = config["work_dir"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "fusion_sqanti"
    shell: 
            """
            rm -f {wildcards.sample}/fusion/sqanti*
            module load singularity
            export SINGULARITY_CACHEDIR={params.sqanti3_path}/.singularity
            export SINGULARITY_BINDPATH="{params.pipeline_home},{params.work_dir},{params.sqanti3_path},/scratch,/lscratch,{params.ref_path}"
            total_fl=`grep '# Total Number of FL reads: ' {input.pb_abundance} | sed 's/# Total Number of FL reads: //'`
            {params.pipeline_home}/scripts/make_paired_abundance.sh {input.abundance} > {wildcards.sample}/fusion/fusion.pairs.abundance.txt
            singularity exec {params.sqanti3_path}/sqanti3_latest.sif bash -c "export R_LIBS={params.sqanti3_path} && sqanti3_qc.py --gtf {params.work_dir}/{input.gff} {params.gtf} {params.ref} --is_fusion -fl {params.work_dir}/{wildcards.sample}/fusion/fusion.pairs.abundance.txt -d {params.work_dir}/{wildcards.sample}/fusion -o sqanti"
            source {params.conda_home}/etc/profile.d/conda.sh
            conda activate {params.conda_env}
            fusion_collate_info.py {wildcards.sample}/fusion/fusion {output.classification} {params.gtf} --genome {params.ref} --total_fl_count ${{total_fl}}
            {params.pipeline_home}/scripts/filter_fusion_annotation.sh {wildcards.sample}/fusion/fusion.annotated.txt > {output.annotation}
            dos2unix {output.annotation}
            {params.pipeline_home}/scripts/add_annotation.pl -a {output.annotation} -i 3 > {output.annotation_with_symbol}.tmp
            {params.pipeline_home}/scripts/add_annotation.pl -a {output.annotation_with_symbol}.tmp -i 7 > {output.annotation_with_symbol}
            rm {output.annotation_with_symbol}.tmp
            """

rule get_abundance:
    input: 
            csv = "{sample}/{sample}.clustered.cluster_report.csv",
            group = "{sample}/pb.collapsed.group.txt",
            fa = "{sample}/pb.collapsed.rep.fa"
    output: 
            "{sample}/pb.collapsed.abundance.txt"
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_cdna_cupcake'],
            pipeline_home = config["pipeline_home"],
            cdna_cupcake_sif = config["apps"]["cdna_cupcake_sif"],
            work_dir = config["work_dir"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "get_abundance"
    shell:
            """
            module load singularity
            export SINGULARITY_BINDPATH="{params.work_dir},/scratch,/lscratch"
            singularity exec {params.cdna_cupcake_sif} get_abundance_post_collapse {params.work_dir}/{wildcards.sample}/pb.collapsed {params.work_dir}/{input.csv}
            """

rule sqanti:
    input: 
            collapse = "{sample}/pb.collapsed.rep.fa", 
            count = "{sample}/pb.collapsed.abundance.txt"
    output: 
            pdf = "{sample}/{sample}.sqanti_sqanti_report.pdf",  
            classification = "{sample}/{sample}.sqanti_classification.txt", 
            faa = "{sample}/{sample}.sqanti_corrected.faa",  
            sgtf = "{sample}/{sample}.sqanti_corrected.gtf",
            sam = "{sample}/{sample}.sqanti_corrected.sam",
            filtered_fa = "{sample}/{sample}.sqanti_classification.filtered_lite.fasta", 
            filtered_classification = "{sample}/{sample}.sqanti_classification.filtered_lite_classification.txt",
            filtered_with_symbol = "{sample}/{sample}.sqanti_classification.filtered_lite_classification_with_symbol.txt",
            bai = "{sample}/{sample}.sqanti_corrected_sorted.bam.bai", 
            bam = "{sample}/{sample}.sqanti_corrected.bam", 
            sort = "{sample}/{sample}.sqanti_corrected_sorted.bam", 
            log = "{sample}/{sample}.sqanti_corrected.bam.log"
            
    version:
            config["version"]["python"]
    params: 
            batch = lambda wildcards: config['cluster'][samples[wildcards.sample]["Type"]]['job_sqanti3'],
            pipeline_home = config["pipeline_home"],
            sqanti3_path = config["apps"]["sqanti3_path"],
            work_dir = config["work_dir"],
            polya = config["apps"]["polya"],
            gtf = config[genome]["gtf_gencode"],
            ref = config[genome]["ref"],
            ref_path = config[genome]["ref_path"],
            cage_peak = config[genome]["cage_peak"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            rulename = "sqanti"
    shell: 
            """
            rm -f {wildcards.sample}/{wildcards.sample}.sqanti*
            module load singularity
            export SINGULARITY_CACHEDIR={params.sqanti3_path}/.singularity
            export SINGULARITY_BINDPATH="{params.work_dir},{params.sqanti3_path},/scratch,/lscratch,{params.ref_path},{params.pipeline_home}"
            singularity exec {params.sqanti3_path}/sqanti3_latest.sif bash -c "export R_LIBS={params.sqanti3_path} && sqanti3_qc.py {params.work_dir}/{input.collapse} {params.gtf} {params.ref} --aligner_choice=minimap2 -t ${{THREADS}} -d {params.work_dir}/{wildcards.sample} -o {wildcards.sample}.sqanti -fl {params.work_dir}/{input.count} --cage_peak {params.cage_peak} --polyA_motif_list {params.polya} --isoAnnotLite"
            singularity exec {params.sqanti3_path}/sqanti3_latest.sif bash -c "export R_LIBS={params.sqanti3_path} && sqanti3_RulesFilter.py {params.work_dir}/{output.classification} {params.work_dir}/{output.faa} {params.work_dir}/{output.sgtf}"
            module load samtools
            samtools view -bS {output.sam} > {output.bam} 2>{output.log}
            samtools sort -@ ${{THREADS}} {output.bam} -o {output.sort}
            samtools index -@ ${{THREADS}} {output.sort}
            dos2unix {output.filtered_classification}
            {params.pipeline_home}/scripts/add_annotation.pl -a {output.filtered_classification} -i 7 > {output.filtered_with_symbol}
            """

rule genebed:
    input: "{sample}.sqanti_classification.txt"
    output: "{sample}.sqanti_classification_genename.txt"
    params: batch = "--nodes=1 --ntasks=8"    
    shell: "export PATH=/mnt/projects/CCR-SF/active/Software/tools/Anaconda/3.7/bin:$PATH; python {genebed} {gtf} {input} {output}"


