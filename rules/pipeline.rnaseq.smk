try:
    #prepare targets
    FASTQS = {}
    BAMS = []
    RSEM= []
    STRINGTIE=[]
    HLA= []

    for sample_id, sample in samples.items():
        sample["Genome"] = config["genome"]
        SAMPLES.append(sample_id)
        FASTQS[sample_id] = []
        #add FASTQ targets
        if not "SampleFiles" in sample:
            raise Exception('SampleFiles not found in sample sheet')
        sample_file = sample["SampleFiles"]
        # check source FASTQ files exist and determine if this is paired_end
        if not os.path.exists(data_dir + "/" + sample_file + "/" + sample_file + suffix_R1):
            raise Exception(data_dir + "/" + sample_file + "/" + sample_file + suffix_R1 + ' not found')
        samples[sample_id]["PE"] = False
        if os.path.exists(data_dir + "/" + sample_file + "/" + sample_file + suffix_R2):
            config['samples'][sample_id]["PE"] = True
        if not "Genome" in sample:
            raise Exception('Genome not found in sample sheet')        
        if "Xenograft" in sample and sample["Xenograft"] == "yes":
            #only add script targets once
            if not "XenograftGenome" in sample:
                raise Exception('XenograftGenome not found in sample sheet')
            #BAMS.append("mapping/xenofilter." + sample_id + "." + sample["Genome"] + "." + sample["XenograftGenome"] + ".ens.genome/Filtered_bams/" + sample_id + "." + sample["Genome"] + "." + sample["XenograftGenome"] + ".ens.genome_Filtered.bam");            
            if (samples[sample_id]["PE"]):
                FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + ".filtered" + suffix_R1)
                FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + ".filtered" + suffix_R2)
            else:
                FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + ".filtered" + suffix_SE)
            #for annotation in config["annotations"]:
                #BAMS.append("mapping/" + sample_id + "." + sample["XenograftGenome"] + "." + annotation + ".genome.bam")
                #BAMS.append("mapping/" + sample_id + "." + sample["XenograftGenome"] + "." + annotation + ".genome.bam.bai")
                #BAMS.append("mapping/" + sample_id + "." + sample["XenograftGenome"] + "." + annotation + ".trans.bam")
                #BAMS.append("mapping/xenofilter." + sample_id + "." + sample["Genome"] + "." + sample["XenograftGenome"] + "." + annotation + ".trans/Filtered_bams/" + sample_id + "." + sample["Genome"] + "." + sample["XenograftGenome"] + "." + annotation + ".trans_Filtered.bam")
                #BAMS.append("mapping/xenofilter." + sample_id + "." + sample["Genome"] + "." + sample["XenograftGenome"] + "." + annotation + ".genome/Filtered_bams/" + sample_id + "." + sample["Genome"] + "." + sample["XenograftGenome"] + "." + annotation + ".genome_Filtered.bam")
            
        else:            
            if (samples[sample_id]["PE"]):
                FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + suffix_R1)
                FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + suffix_R2)
            else:
                FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + suffix_SE)
        for annotation in config["annotations"]:
            BAMS.append(sample_id + "/STAR_" + sample["Genome"] + "_" + annotation + "/" + sample_id + "." + sample["Genome"] + "." + annotation + ".star.genome.bam")
            BAMS.append(sample_id + "/STAR_" + sample["Genome"] + "_" + annotation + "/" + sample_id + "." + sample["Genome"] + "." + annotation + ".star.genome.bam.bai")
            BAMS.append(sample_id + "/STAR_" + sample["Genome"] + "_" + annotation + "/" + sample_id + "." + sample["Genome"] + "." + annotation + ".star.trans.bam")
            RSEM.append(sample_id + "/RSEM_" + sample["Genome"] + "_" + annotation + "/" + sample_id + "." + sample["Genome"] + "." + annotation + ".genes.results")
            RSEM.append(sample_id + "/RSEM_" + sample["Genome"] + "_" + annotation + "/" + sample_id + "." + sample["Genome"] + "." + annotation + ".isoforms.results")
            STRINGTIE.append(sample_id + "/StringTie_" + sample["Genome"] + "_" + annotation + "/" + sample_id + "." + sample["Genome"] + "." + annotation + ".stringtie.gtf")
            STRINGTIE.append(sample_id + "/StringTie_" + sample["Genome"] + "_" + annotation + "/" + sample_id + "." + sample["Genome"] + "." + annotation + ".stringtie.tab")
            f="{sample}/StringTie_{genome}_{annotation}/{sample}.{genome}.{annotation}.stringtie.gtf",
        HLA.append(sample_id + "/HLA/" + sample_id + ".Calls.txt")
        
except Exception as err:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    output = io.StringIO()
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=output)
    contents = output.getvalue()
    output.close()
    print(contents)    
    shell("echo 'RNAseq pipeline has exception: reason " + contents + ". Working Dir:  {work_dir}' |mutt -e 'my_hdr From:chouh@nih.gov' -s 'Khanlab RNAseq Pipeline Status' `whoami`@mail.nih.gov {emails} ")
    sys.exit()
    
TARGETS = BAMS + HLA + RSEM + STRINGTIE

localrules: all, prepareFASTQ, MergeHLA

rule RSEM:
    input:
            bam="{sample}/STAR_{genome}_{annotation}/{sample}.{genome}.{annotation}.star.trans.bam"
    output:
            "{sample}/RSEM_{genome}_{annotation}/{sample}.{genome}.{annotation}.genes.results",
            "{sample}/RSEM_{genome}_{annotation}/{sample}.{genome}.{annotation}.isoforms.results"
    version:
            config["version"]["rsem"]
    params:
            rulename = "RSEM",
            ref = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["rsem_ref_" + wildcards.annotation],
            batch = config['cluster']['job_rsem'],
            library = lambda wildcards: "--strandedness reverse" if samples[wildcards.sample]['Enrichment step'] != "polya" else "",
            paired_end = lambda wildcards: "--paired-end" if samples[wildcards.sample]["PE"] else "",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            work_dir = work_dir
    benchmark:
            "{sample}/benchmark/rsem.{genome}.{annotation}.benchmark.txt"
    shell:
            """
            module load rsem/{version}
            cd ${{LOCAL}}/
            rsem-calculate-expression --no-bam-output {params.paired_end} -p ${{THREADS}} {params.library} --estimate-rspd \
                --bam {params.work_dir}/{input.bam} {params.ref} {wildcards.sample}_{wildcards.annotation}
            mkdir -p {params.work_dir}/{wildcards.sample}/RSEM_{wildcards.genome}_{wildcards.annotation}
            mv -f {wildcards.sample}_{wildcards.annotation}.genes.results {params.work_dir}/{wildcards.sample}/RSEM_{wildcards.genome}_{wildcards.annotation}/{wildcards.sample}.{wildcards.genome}.{wildcards.annotation}.genes.results
            mv -f {wildcards.sample}_{wildcards.annotation}.isoforms.results {params.work_dir}/{wildcards.sample}/RSEM_{wildcards.genome}_{wildcards.annotation}/{wildcards.sample}.{wildcards.genome}.{wildcards.annotation}.isoforms.results
            """

rule XenofilteR:
    input:
            ref_genome = lambda wildcards: wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + "." + wildcards.annotation + "." + wildcards.ref_type + ".bam",
            xeno_genome = lambda wildcards: wildcards.sample +  "/" + wildcards.sample + "." + wildcards.xeno_genome + "." + wildcards.annotation + "." + wildcards.ref_type + ".bam"
    output:
            bam="{sample}/xenofilter.{sample}.{genome}.{xeno_genome}.{annotation}.{ref_type}/Filtered_bams/{sample}.{genome}.{xeno_genome}.{annotation}.{ref_type}_Filtered.bam"
    params:
            pipeline_home = config["pipeline_home"],
            batch    = config["cluster"]["job_xenofilter"],
            rulename = "XenofilteR",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            output_dir = lambda wildcards: "mapping/xenofilter." + wildcards.sample + "." + wildcards.genome + "." + wildcards.xeno_genome + "." + wildcards.annotation + "." + wildcards.ref_type
    shell:
            """
            module load R/4.0.0            
            export R_LIBS={params.pipeline_home}/Rlibs/4.0/sample
            rm -rf {params.output_dir}
            mkdir -p {params.output_dir}
            echo -e "{input.ref_genome}\t{input.xeno_genome}" > {params.output_dir}/sample_list.txt
            Rscript {params.pipeline_home}/scripts/runXenofilteR.R -s {params.output_dir}/sample_list.txt \
                    -o {params.output_dir} \
                    -n {wildcards.sample}.{wildcards.genome}.{wildcards.xeno_genome}.{wildcards.annotation}.{wildcards.ref_type}
            """
 
############
#   seq2HLA
############
rule seq2HLA:
    input:  
            lambda wildcards: FASTQS[wildcards.sample]
    output:
            "{sample}/HLA/seq2HLA/{sample}-ClassI.HLAgenotype4digits"
    params:
            rulename= "seq2HLA",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            app_home = config["app_home"],
            python_version  = config["version_common"]['python2'],
            R_version = config["version_common"]['R'],
            bowtie_version = config["version"]['bowtie'],
            batch = config["cluster"]["job_seq2hla"],
            hla=lambda wildcards: config[samples[wildcards.sample]["Genome"]]["hla_ref"],
            fastqs = lambda wildcards: " -1 " + FASTQS[wildcards.sample][0] + " -2 " + FASTQS[wildcards.sample][1] if len(FASTQS[wildcards.sample])==2 else " -1 " + FASTQS[wildcards.sample][0]
    benchmark:
            "{sample}/benchmark/seq2HLA.{sample}.benchmark.txt"
    shell:
            """
            module load bowtie/{params.bowtie_version} python/{params.python_version} R/{params.R_version}
            python {params.app_home}/seq2HLA/seq2HLA.py {params.hla}/seq2HLA/ -1 {input[0]} -2 {input[1]}  -p ${{THREADS}} -r {wildcards.sample}/HLA/seq2HLA/{wildcards.sample}
            """
            
############
#   HLAminer
############
rule HLAminer:
    input:
            lambda wildcards: FASTQS[wildcards.sample]
    output: 
            "{sample}/HLA/HLAminer/HLAminer_HPTASR.csv"
    version:
            config["version"]["HLAminer"]
    params:
            app_home = config["app_home"],
            rulename="HLAminer",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            batch=config["cluster"]["job_hlaminer"],
            hla=lambda wildcards: config[samples[wildcards.sample]["Genome"]]["hla_ref"],
            pipeline_home=pipeline_home
    benchmark:
            "{sample}/benchmark/HLAminer.{sample}.benchmark.txt"
    shell: 
            """
            echo {work_dir}/{input[0]} >{wildcards.sample}/HLA/HLAminer/patient.fof
            echo {work_dir}/{input[1]} >>{wildcards.sample}/HLA/HLAminer/patient.fof
            bash {params.app_home}/HLAminer_v{version}/bin/HPTASRwgs_classI.sh {params.app_home}/HLAminer_v{version}/bin {wildcards.sample}/HLA/HLAminer/
            """
############
##  MergeHLA Calls
#############
rule MergeHLA:
    input:
        A="{sample}/HLA/HLAminer/HLAminer_HPTASR.csv",
        B="{sample}/HLA/seq2HLA/{sample}-ClassI.HLAgenotype4digits",
    output:
        "{sample}/HLA/{sample}.Calls.txt"
    params:
        rulename = "MergeHLA",
        log_dir = lambda wildcards: wildcards.sample + '/log',
        script=pipeline_home + "/scripts/consensusHLA.pl"
    shell: 
        """
        export LC_ALL=C
        perl {params.script} {input.B} {input.A} | sort > {output}
        """

rule StringTie:
    input:
            "{sample}/STAR_{genome}_{annotation}/{sample}.{genome}.{annotation}.star.genome.bam",
    output: 
            gtf="{sample}/StringTie_{genome}_{annotation}/{sample}.{genome}.{annotation}.stringtie.gtf",
            tab="{sample}/StringTie_{genome}_{annotation}/{sample}.{genome}.{annotation}.stringtie.tab",
    version:
            config["version"]["stringtie"]
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_stringtie"],
            rulename = "StringTie",
            strandness = lambda wildcards: "--rf" if samples[wildcards.sample]['Enrichment step'] != "polya" else "",
            ref_gtf = lambda wildcards: config[wildcards.genome]["gtf_" + wildcards.annotation],
            log_dir = lambda wildcards: wildcards.sample + '/log',
    benchmark:
            "{sample}/benchmark/stringtie.{sample}.{genome}.{annotation}.benchmark.txt"
    shell:
            """
            module load stringtie/{version} samtools
            mkdir -p {wildcards.sample}/StringTie_{wildcards.genome}_{wildcards.annotation}
            num_cores=$SLURM_CPUS_ON_NODE
            echo $num_cores
            echo ${{THREADS}}
            stringtie {input} -p ${{THREADS}} -o {output.gtf} -G {params.ref_gtf} -A {output.tab} {params.strandness}
            """

    
rule STAR:
    input:
            lambda wildcards: FASTQS[wildcards.sample]
    output: 
            genome_bam="{sample}/STAR_{genome}_{annotation}/{sample}.{genome}.{annotation}.star.genome.bam",
            genome_bai="{sample}/STAR_{genome}_{annotation}/{sample}.{genome}.{annotation}.star.genome.bam.bai",
            trans_bam="{sample}/STAR_{genome}_{annotation}/{sample}.{genome}.{annotation}.star.trans.bam"
    version:
            config["version"]["star"]
    params:
            star_ref = lambda wildcards: config[wildcards.genome]["star_index_" + wildcards.annotation],
            fastqs = lambda wildcards: " ".join(work_dir + '/' + x for x in FASTQS[wildcards.sample]),
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_star"],
            rulename = "STAR",
            log_dir = lambda wildcards: wildcards.sample + '/log',
    benchmark:
            "{sample}/benchmark/star.{sample}.{genome}.{annotation}.benchmark.txt"
    shell:
            """
            module load STAR/{version} samtools
            cd ${{LOCAL}}/
            num_cores=$SLURM_CPUS_ON_NODE
            echo $num_cores
            echo ${{THREADS}}
            STAR 	--outTmpDir twopass \
                    --genomeDir {params.star_ref} \
                    --readFilesIn {params.fastqs} \
                    --readFilesCommand zcat\
                    --outSAMtype BAM Unsorted\
                    --twopassMode Basic \
                    --outFileNamePrefix {wildcards.sample}.{wildcards.genome}.{wildcards.annotation} \
                    --runThreadN ${{THREADS}} \
                    --outFilterMismatchNmax 2\
                    --outSAMunmapped Within \
                    --quantMode TranscriptomeSAM \
                    --outSAMattrIHstart 0 \
                    --outSAMattributes NM
            echo "Finished STAR twopass mapping"
            mkdir -p {wildcards.sample}/STAR_{wildcards.genome}.{wildcards.annotation}
            #
            mv -f {wildcards.sample}.{wildcards.genome}.{wildcards.annotation}Aligned.toTranscriptome.out.bam {params.work_dir}/{output.trans_bam}
            samtools sort -@ ${{THREADS}} -o {wildcards.sample}.{wildcards.genome}.{wildcards.annotation}Aligned.sortedByCoord.out.bam -O bam {wildcards.sample}.{wildcards.genome}.{wildcards.annotation}Aligned.out.bam
            mv -f {wildcards.sample}.{wildcards.genome}.{wildcards.annotation}Aligned.sortedByCoord.out.bam {params.work_dir}/{output.genome_bam}
            samtools index -@ ${{THREADS}} {params.work_dir}/{output.genome_bam}
            outf={params.work_dir}/{output.genome_bam}
            fn="${{outf%.bam}}"
            samtools flagstat {params.work_dir}/{output.genome_bam} > $fn.flagstat.txt
            samtools flagstat {params.work_dir}/{output.trans_bam} > {params.work_dir}/{output.trans_bam}.flagstat.txt
            #module load deeptools
            #bamCoverage -b {params.work_dir}/{output.genome_bam} -o $fn.bw
            """
            
rule xengsort_pe:
    input:
            R1="{sample}/DATA/{sample}" + suffix_R1,
            R2="{sample}/DATA/{sample}" + suffix_R2,
    output: 
            "{sample}/DATA/{sample}.filtered" + suffix_R1,
            "{sample}/DATA/{sample}.filtered" + suffix_R2,
    version:
            config["version"]["xenome"]
    benchmark:
            "{sample}/benchmark/xengsort.{sample}.benchmark.txt"
    params:
            xengsort_ref = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["xengsort_ref"],            
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_xengsort"],
            rulename = "xengsort_pe",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            pipeline_home=pipeline_home
    shell:
            """
            module load xengsort
            cd ${{LOCAL}}/
            xengsort classify --index {params.xengsort_ref} -T ${{THREADS}} --prefix {wildcards.sample} --classification new \
                --fastq <(zcat {params.work_dir}/{input.R1}) --pairs <(zcat {params.work_dir}/{input.R2})
            perl {params.pipeline_home}/scripts/replaceNbase.pl -i {wildcards.sample}-graft.1.fq > {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R1.fastq
            #perl {params.pipeline_home}/scripts/replaceNbase.pl -i {wildcards.sample}/DATA/{wildcards.sample}-both.1.fq >> {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R1.fastq
            #perl {params.pipeline_home}/scripts/replaceNbase.pl -i {wildcards.sample}/DATA/{wildcards.sample}-neither.1.fq >> {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R1.fastq
            perl {params.pipeline_home}/scripts/replaceNbase.pl -i {wildcards.sample}-graft.2.fq > {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R2.fastq
            #perl {params.pipeline_home}/scripts/replaceNbase.pl -i {wildcards.sample}/DATA/{wildcards.sample}-both.2.fq >> {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R2.fastq
            #perl {params.pipeline_home}/scripts/replaceNbase.pl -i {wildcards.sample}/DATA/{wildcards.sample}-neither.2.fq >> {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R2.fastq            
            gzip {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R1.fastq
            gzip {params.work_dir}/{wildcards.sample}/DATA/{wildcards.sample}.filtered_R2.fastq
            echo -e "sample\thost\tgraft\tambiguous\tboth\tneither" > {params.work_dir}/{wildcards.sample}/DATA/classification.tsv
            cat classification.tsv >> {params.work_dir}/{wildcards.sample}/DATA/classification.tsv
            
            """

rule prepareFASTQ:
    input: 
            lambda wildcards: data_dir + "/" + samples[wildcards.sample]["SampleFiles"] + "/" + samples[wildcards.sample]["SampleFiles"] + "_" + wildcards.suffix
    output: 
            "{sample}/DATA/{sample}_{suffix}",
    shell:
            """
            mkdir -p {wildcards.sample}
            mkdir -p {wildcards.sample}/log
            mkdir -p {wildcards.sample}/DATA
            ln -s {input} {output}
            """
