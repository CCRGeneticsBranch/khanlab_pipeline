try:
    #prepare targets
    FASTQS = {}
    FASTQCS = []
    BAMS = []
    BWS= []
    TDFS = []
    SPP = []
    MACS2 = []
    ANNOTATION = []
    ROSE = []
    MOTIFS = []
    EDENS = []
    RNAseq = []
    COLTRONS = []
    QC = []
    MACS_BAMS = {}
    for sample_id, sample in samples.items():
        ##### Test no input
        sample["Matched normal"] = "."
        #####
        sample["Genome"] = config["genome"]
        SAMPLES.append(sample_id)
        FASTQS[sample_id] = []
        MACS_BAMS[sample_id] = []
        #add default values
        if not "PeakCalling" in sample:
            sample["PeakCalling"] = "narrow"
        if not "LibrarySize" in sample:
            sample["LibrarySize"] = "."
        if not "Matched RNA-seq lib" in sample:
            sample["Matched RNA-seq lib"] = "."
        if not "Matched normal" in sample:
            sample["Matched normal"] = "."
        if not "EnhancePipe" in sample:
            sample["EnhancePipe"] = "no"
        #add FASTQ targets
        if not "SampleFiles" in sample:
            raise Exception('SampleFiles not found in sample sheet')
        sample_file = sample["SampleFiles"]
        # check source FASTQ files exist and determine if this is paired_end
        if not os.path.exists(data_dir + "/" + sample_file + "/" + sample_file + suffix_R1):
            raise Exception(data_dir + "/" + sample_file + "/" + sample_file + suffix_R1 + ' not found')
        samples[sample_id]["PE"] = False
        if os.path.exists(data_dir + "/" + sample_file + "/" + sample_file + suffix_R2):
            samples[sample_id]["PE"] = True
            FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + suffix_R1)
            FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + suffix_R2)
            FASTQCS.append(sample_id + "/qc/" + sample_id + '_R1_fastqc.html')
            FASTQCS.append(sample_id + "/qc/" + sample_id + '_R2_fastqc.html')
        else:
            FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + suffix_SE)
            FASTQCS.append(sample_id + "/qc/" + sample_id + '_R1_fastqc.html')
        if not "Genome" in sample:
            raise Exception('Genome not found in sample sheet')
        if "SpikeIn" in sample and sample["SpikeIn"] == "yes":
            #only add script targets once
            if not "SpikeIn_genome" in sample:
                raise Exception('SpikeIn_genome not found in sample sheet')
            spike_in_genome = sample["SpikeIn_genome"]
            BWS.append(sample_id + "/" + sample_id + "." + config["bin_size"] + ".scaled.bw")
            TDFS.append(sample_id + "/" + sample_id + "." + config["bin_size"] + ".scaled.tdf")
        else:
            sample["SpikeIn"] == "no"
            sample["SpikeIn_genome"] == ""
        if not "LibrarySize" in sample or samples[sample_id]["LibrarySize"] == "." or samples[sample_id]["LibrarySize"] == "":
            samples[sample_id]["LibrarySize"] = config["default_lib_length"]
        BAMS.append(sample_id + "/" + sample_id + ".bam")        
        MACS_BAMS[sample_id].append(sample_id + "/" + sample_id + ".bam")
        BWS.append(sample_id + "/" + sample_id + "." + config["bin_size"] + ".RPM.bw")
        TDFS.append(sample_id + "/" + sample_id + "." + config["bin_size"] + ".RPM.tdf")
        SPP.append(sample_id + "/qc/" + sample_id + ".spp.pdf")
        SPP.append(sample_id + "/qc/" + sample_id + ".spp.txt")
        QC.append(sample_id + "/qc/fingerPrint.pdf")
        QC.append(sample_id + "/qc/fingerPrint.tab")
        has_exp = "Matched RNA-seq lib" in sample and sample["Matched RNA-seq lib"] != "." and sample["Matched RNA-seq lib"] != ""
        if has_exp:
            RNAseq.append(sample_id + "/RNAseq/" + sample["Matched RNA-seq lib"] + "." +  sample["Genome"] + ".ucsc.genes.TPM.txt")
        cutoff_types = ["p","q"]
        for cutoff_type in cutoff_types:
            for cutoff_value in config["macs2"][cutoff_type + "values"]:
                cutoff = "{:.0e}".format(cutoff_value)
                peak_type = sample["PeakCalling"]
                macs_out_prefix = "/MACS_Out_" + cutoff_type + "_"
                MACS2.append(sample_id + macs_out_prefix + cutoff + "/" + sample_id + "_peaks." + peak_type + "Peak.nobl.bed")
                if peak_type == "narrow":
                    MACS2.append(sample_id + macs_out_prefix + cutoff + "/" + sample_id + "." + peak_type + "_summits.bed")
                    MOTIFS.append(sample_id + macs_out_prefix + cutoff + "/motif_" + peak_type + "/knownResults.html")
                if has_exp:
                    EDENS.append(sample_id + macs_out_prefix + cutoff + "/" + sample_id + "_peaks." + peak_type + "Peak.nobl_TPM" + str(config["eden"]["TPM_cutoff"]) + "_multi-genes.txt")
                ANNOTATION.append(sample_id + macs_out_prefix + cutoff + "/" + sample_id + "_peaks." + peak_type + "Peak.nobl.bed.annotation.txt")
                ANNOTATION.append(sample_id + macs_out_prefix + cutoff + "/" + sample_id + "_peaks." +peak_type + "Peak.nobl.bed.annotation.summary")
                
                if sample["EnhancePipe"] == "yes":
                    ROSE.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/" + sample_id + "_peaks_AllEnhancers.table.super.bed")
                    ROSE.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/" + sample_id + "_peaks_AllEnhancers.table.regular.bed")
                    if peak_type == "narrow":
                        ROSE.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/" + sample_id + ".super_summits.bed")
                        ROSE.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/" + sample_id + ".regular_summits.bed")
                        MOTIFS.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/motif_super/knownResults.html")
                        MOTIFS.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/motif_regular/knownResults.html")
                    if has_exp:
                        COLTRONS.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/coltron/" + sample_id + "_CLIQUES_RANKED.txt")
                        EDENS.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/" + sample_id + "_peaks_AllEnhancers.table.super_TPM" + str(config["eden"]["TPM_cutoff"]) + "_multi-genes.txt")
                        EDENS.append(sample_id + macs_out_prefix + cutoff + "/ROSE_out_" + str(config["rose"]["stitch_distance"]) + "/" + sample_id + "_peaks_AllEnhancers.table.regular_TPM" + str(config["eden"]["TPM_cutoff"]) + "_multi-genes.txt")
        
        # if input sample not exists, generate bam file
        input_sample_id = sample["Matched normal"]
        input_sample_file = config["FASTQ_prefix"] + input_sample_id
        if input_sample_id != "." and input_sample_id != "":
            input_sample = {}
            input_sample["Genome"] = sample["Genome"]
            input_sample["LibrarySize"] = sample["LibrarySize"]
            input_sample["SampleFiles"] = input_sample_file
            input_sample["SpikeIn"] = "no"
            input_sample["SpikeIn_genome"] = ""
            samples[input_sample_id] = input_sample
            BAMS.append(input_sample_id + "/" + input_sample_id + ".bam")
            BWS.append(input_sample_id + "/" + input_sample_id + "." + config["bin_size"] + ".RPM.bw")
            MACS_BAMS[sample_id].append(input_sample_id + "/" + input_sample_id + ".bam")
            FASTQS[input_sample_id] = []
            if os.path.exists(data_dir + "/" + input_sample_file + "/" + input_sample_file + suffix_R2):
                FASTQS[input_sample_id].append(input_sample_id + "/DATA/" + input_sample_id + suffix_R1)
                FASTQS[input_sample_id].append(input_sample_id + "/DATA/" + input_sample_id + suffix_R2)
                FASTQCS.append(input_sample_id + "/qc/" + input_sample_id + '_R1_fastqc.html')
                FASTQCS.append(input_sample_id + "/qc/" + input_sample_id + '_R2_fastqc.html')
            else:
                FASTQCS.append(input_sample_id + "/qc/" + input_sample_id + '_R1_fastqc.html')
                FASTQS[input_sample_id].append(input_sample_id + "/DATA/" + input_sample_id + suffix_SE)
        
        
except Exception as err:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    output = io.StringIO()
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=output)
    contents = output.getvalue()
    output.close()
    print(contents)    
    shell("echo 'ChIPseq pipeline has exception: reason " + contents + ". Working Dir:  {work_dir}' |mutt -e 'my_hdr From:chouh@nih.gov' -s 'Khanlab ChIPseq Pipeline Status' `whoami`@mail.nih.gov {emails} ")
    sys.exit()
    
TARGETS = FASTQCS + BAMS + BWS + TDFS + SPP + MACS2 + ANNOTATION + ROSE + MOTIFS + EDENS + RNAseq + COLTRONS + QC

localrules: all, prepareFASTQ, EDEN, prepareSummit, prepareRoseSummit

rule RNAseq_pipeline:
    input:
            config["work_dir"] + "/../../RNAseq"
    output:
            "{sample}/RNAseq/{rnaseq_sample}.{genome}.ucsc.genes.TPM.txt"
    version:
            config["version_common"]["snakemake"]
    params:
            work_dir = config["work_dir"],
            data_dir = config["data_dir"],
            now = config["now"],
            batch    = config["cluster_common"]["subflow"],
            pipeline_home = config["pipeline_home"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            rulename = "RNAseq_pipeline",
            rnaseq_sample = lambda wildcards: samples[wildcards.sample]["Matched RNA-seq lib"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            gene_bed = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["gene_bed"],
            version_python=config["version_common"]["python3"]
    shell:
            """
            mkdir -p {wildcards.sample}/RNAseq
            module load python/{params.version_python}
            python {params.pipeline_home}/scripts/sampleToYaml.py -s {params.rnaseq_sample} -o {wildcards.sample}/RNAseq/{params.rnaseq_sample}.rnaseq.yaml
            module load snakemake/{version}
            {params.pipeline_home}/launch -t rnaseq -w {params.work_dir}/../../RNAseq -s {wildcards.sample}/RNAseq/{params.rnaseq_sample}.rnaseq.yaml -local
            echo -e "Chr\\tStart\\tStop\\tGeneID\\tTPM" > {output}
            join -t $'\\t' -1 4 -2 1 <(sort -k4,4 {params.pipeline_home}/{params.gene_bed} ) <(grep -v gene_id {params.work_dir}/../../RNAseq/{params.rnaseq_sample}/RSEM_{params.genome}_ucsc/{params.rnaseq_sample}.{params.genome}.ucsc.genes.results|sort -k1,1) | awk -F'\\t' 'OFS="\\t"{{print $2,$3,$4,$1,$10}}' | sort -k 1,1 -k2,2n >> {output}
            """

rule coltron:
    input:
            rose_out="{sample}/{macs_dir}/{rose_dir}/{sample}_peaks_AllEnhancers.table.txt",
            bam="{sample}/{sample}.bam",
            exp_file=lambda wildcards: wildcards.sample + "/RNAseq/" +  samples[wildcards.sample]["Matched RNA-seq lib"] + "." + samples[wildcards.sample]["Genome"] + ".ucsc.genes.TPM.txt"
    output:
            "{sample}/{macs_dir}/{rose_dir}/coltron/{sample}_CLIQUES_RANKED.txt"
    benchmark:
            "{sample}/benchmark/coltron.{sample}.{macs_dir}.{rose_dir}.benchmark.txt"
    version:
            config["version"]["bamliquidator"],
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster_common"]["small"],
            pipeline_home = config["pipeline_home"],
            coltron_bin = lambda wildcards: config["coltron"]["path"] if samples[wildcards.sample]["Matched RNA-seq lib"] != "." and samples[wildcards.sample]["Matched RNA-seq lib"] != "" else config["coltron"]["path_noexp"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"].upper(),
            nearest_tf_distance_cutoff = str(config["coltron"]["nearest_tf_distance_cutoff"]),
            exp_cutoff=str(config["coltron"]["nearest_tf_distance_cutoff"]),
            rulename = "coltron",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            R_version=config["version_common"]["R"]
    shell:
            """
            module load meme
            module load bamliquidator/{version}
            module load R/{params.R_version}
            {params.coltron_bin} -e {input.rose_out} -b {input.bam} -g {params.genome} -o {wildcards.sample}/{wildcards.macs_dir}/{wildcards.rose_dir}/coltron/ \
                -n {wildcards.sample} -d {params.nearest_tf_distance_cutoff} -a {input.exp_file} -x {params.exp_cutoff}
            
            """

rule EDEN:
    input:
            bed=lambda wildcards: wildcards.sample + "/" + wildcards.folders + "/" +  wildcards.sample + "_" + wildcards.prefix + ".bed",
            exp_file=lambda wildcards: wildcards.sample + "/RNAseq/" +  samples[wildcards.sample]["Matched RNA-seq lib"] + "." + samples[wildcards.sample]["Genome"] + ".ucsc.genes.TPM.txt"
    output:
            "{sample}/{folders}/{sample}_{prefix}_TPM{TPM_cutoff}_multi-genes.txt"
    wildcard_constraints:
            folders = "MACS.+"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster_common"]["small"],
            pipeline_home = config["pipeline_home"],
            tad_file = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["tad_file"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            super_loci_distance_cutoff = str(config["eden"]["super_loci_distance_cutoff"]),
            nearest_gene_distance_cutoff = str(config["eden"]["nearest_gene_distance_cutoff"]),
            exp_file = lambda wildcards: samples[wildcards.sample]["Matched RNA-seq lib"] + ".gene.TPM.txt",
            gene_bed = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["gene_bed"],
            rulename = "eden",
            log_dir = lambda wildcards: wildcards.sample + '/log',
    shell:
            """
            out_dir=$(dirname "{output}")
            {params.pipeline_home}/scripts/EDEN.pl -e {input.exp_file} -o ${{out_dir}} -x {wildcards.sample}_{wildcards.prefix} -t TPM -c -d {params.pipeline_home}/{params.tad_file} -b {input.bed} -f {wildcards.TPM_cutoff} \
                -s {params.super_loci_distance_cutoff} -n {params.nearest_gene_distance_cutoff} -r {params.pipeline_home}/{params.gene_bed}
            """
            
rule findMotif:
    input:
            "{sample}/{folders}/{sample}.{peak_type}_summits.bed",
    output:
            "{sample}/{folders}/motif_{peak_type}/knownResults.html",
    wildcard_constraints:
            folders = "MACS.+"
    version:
            config["version"]["homer"]
    benchmark:
            "{sample}/benchmark/findMotif.{sample}.{folders}.{peak_type}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_motif"],
            pipeline_home=config["pipeline_home"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            rulename = "findMotif",
            motif_size = str(config["homer"]["motif_size"]),
            log_dir = lambda wildcards: wildcards.sample + '/log',
    shell:
            """  
            module load homer/{version}
            out_dir=$(dirname {output})
            if [ -s {input} ];then
                findMotifsGenome.pl {input} {params.genome} ${{out_dir}} -size {params.motif_size} -p ${{THREADS}} -preparsedDir {params.pipeline_home}/ref/preparsedDir
            else
                touch {output}
            fi
            """

rule prepareRoseSummit:
    input:
            summit=lambda wildcards: wildcards.sample + "/MACS_Out_" + wildcards.cutoff + "/" + wildcards.sample + "." + samples[wildcards.sample]["PeakCalling"] + "_summits.bed",
            se_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.bed",
            re_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.bed",
    output:
            se_summit="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}.super_summits.bed",
            re_summit="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}.regular_summits.bed",
    version:
            config["version_common"]["bedtools"]
    params:
            work_dir = config["work_dir"],
            pipeline_home=config["pipeline_home"],
            rulename = "prepareRoseSummit",
            log_dir = lambda wildcards: wildcards.sample + '/log',
    shell:
            """  
            module load bedtools/{version}
            bedtools intersect -wa -a {input.summit} -b {input.se_table} > {output.se_summit}
            bedtools intersect -wa -a {input.summit} -b {input.re_table} > {output.re_summit}
            """

rule rose:
    input:
            bed=lambda wildcards: wildcards.sample + "/MACS_Out_" + wildcards.cutoff + "/" + wildcards.sample + "_peaks." + samples[wildcards.sample]["PeakCalling"] + "Peak.nobl.bed",
            #bed="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS.bed",
            bam="{sample}/{sample}.bam"
    output:
            rose_out="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.txt",
            se_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.bed",
            re_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.bed",
            se_great="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.GREAT.bed",
            re_great="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.GREAT.bed",
            #no_tss="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS.bed",
            #no_tss_stitched="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS_{stitch_distance}.bed",
    version:
            config["version"]["rose"],            
    benchmark:
            "{sample}/benchmark/rose.{sample}.{cutoff}.{stitch_distance}benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_rose"],
            pipeline_home=config["pipeline_home"],
            peak_type = lambda wildcards: samples[wildcards.sample]["PeakCalling"],
            tss_distance=config["rose"]["tss_distance"],
            tss_bed=lambda wildcards: config["pipeline_home"] + "/" + config[samples[wildcards.sample]["Genome"]]["tss_bed"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"].upper(),
            annotation = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["rose"],
            input_control = lambda wildcards: " -c " + config["work_dir"] + "/" + samples[wildcards.sample]["Matched normal"] + "/" + samples[wildcards.sample]["Matched normal"] + ".bam" if samples[wildcards.sample]["Matched normal"] != "." and samples[wildcards.sample]["Matched normal"] != "" else "",
            rulename = "rose",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            version_R=config["version_common"]["R"],
            version_python=config["version_common"]["python2"],
            version_bedtools=config["version_common"]["bedtools"],
    shell:
            """
            module load bedtools/{params.version_bedtools}
            bedtools intersect -a {input.bed} -b {params.tss_bed} -v > {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS.bed
            bedtools merge -i {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS.bed -d {wildcards.stitch_distance} -c 4,5,6 -o distinct,sum,distinct > {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS_{wildcards.stitch_distance}.bed
            module load rose/{version}
            #module load bamliquidator
            module load python/{params.version_python}
            module load R/{params.version_R}
            #cd /usr/local/apps/bamliquidator/pipeline
            cd $(dirname `which ROSE_main.py`)
            #cd {params.pipeline_home}/apps/rose
            #export PATH=$PATH:{params.pipeline_home}/apps/rose
            ./ROSE_main.py -i {params.work_dir}/{wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS_{wildcards.stitch_distance}.bed -g {params.genome} -r {params.work_dir}/{input.bam} {params.input_control} -t {params.tss_distance} -s {wildcards.stitch_distance} -o {params.work_dir}/{wildcards.sample}/MACS_Out_{wildcards.cutoff}/ROSE_out_{wildcards.stitch_distance}
            {params.pipeline_home}/scripts/roseTable2Bed.sh {params.work_dir}/{output.rose_out} {params.work_dir}/{output.re_table} {params.work_dir}/{output.se_table}
            cd {params.work_dir}
            cut -f1-3 {output.se_table} > {output.se_great}
            cut -f1-3 {output.re_table} > {output.re_great}
            """

rule peakAnnotation:
    input:
            "{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.bed"
    output:
            annotation="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.bed.annotation.txt",
            annotation_summary="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.bed.annotation.summary",
    version:
            config["version"]["homer"]
    benchmark:
            "{sample}/benchmark/peakAnnotation.{sample}.{cutoff}.{peak_type}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster_common"]["small"],
            pipeline_home=config["pipeline_home"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            rulename = "peakAnnotation",
            log_dir = lambda wildcards: wildcards.sample + '/log',
    shell:
            """  
            module load homer/{version}
            annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
            """

rule prepareSummit:
    input:
            nobl="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.bed",
    output:
            summit="{sample}/MACS_Out_{cutoff}/{sample}.{peak_type}_summits.bed",
    wildcard_constraints:
            cutoff = "[p|q]_1e\-\d+"
    version:
            config["version_common"]["bedtools"],
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_macs2"],
            pipeline_home=config["pipeline_home"],
            rulename = "prepareSummit",
            black_list= lambda wildcards: config["pipeline_home"] + "/" + config[samples[wildcards.sample]["Genome"]]["black_list"],
            exclude_chr_list = lambda wildcards: config["pipeline_home"] + "/" + config[samples[wildcards.sample]["Genome"]]["exclude_chr_list"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
    shell:
            """  
            module load bedtools/{version}
            #remove unused chromosomes and blacklist in summit
            {params.pipeline_home}/scripts/removeChr.pl -i {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_summits.bed -o {output.summit}.tmp -l {params.exclude_chr_list}
            sed -i 's/{wildcards.sample}//g' {output.summit}.tmp
            bedtools intersect -a {output.summit}.tmp -b {params.black_list} -v > {output.summit}
            rm {output.summit}.tmp
            """
            
rule macs2:
    input:
            lambda wildcards: MACS_BAMS[wildcards.sample]
    output:
            peak="{sample}/MACS_Out_{cutoff_type}_{cutoff}/{sample}_peaks.{peak_type}Peak",
            renamed_peak="{sample}/MACS_Out_{cutoff_type}_{cutoff}/{sample}_peaks.{peak_type}Peak.bed",
            nobl="{sample}/MACS_Out_{cutoff_type}_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.bed",
            great="{sample}/MACS_Out_{cutoff_type}_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.GREAT.bed",
    wildcard_constraints:
            cutoff_type = "[p|q]",
            cutoff = "1e\-\d+"
    version:
            config["version"]["macs"],
    benchmark:
            "{sample}/benchmark/macs2.{sample}.{cutoff_type}.{cutoff}.{peak_type}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_macs2"],
            pipeline_home=config["pipeline_home"],
            genome = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["macs2_gsize"],
            rulename = "macs2",
            ext_option = lambda wildcards: "" if samples[wildcards.sample]["PE"] else " --extsize " + str(samples[wildcards.sample]["LibrarySize"]),
            bam_format = lambda wildcards: " BAMPE " if samples[wildcards.sample]["PE"] else "BAM",
            peak_type = lambda wildcards: " --broad " if wildcards.peak_type=="broad" else "",
            dup_cutoff = config["dup_cutoff"],
            cutoff_option = lambda wildcards: " -" + wildcards.cutoff_type + " " + wildcards.cutoff + " --broad-cutoff " + wildcards.cutoff if wildcards.peak_type=="broad" else " -" + wildcards.cutoff_type + " " + wildcards.cutoff,
            input_control = lambda wildcards: " -c " + samples[wildcards.sample]["Matched normal"] + "/" + samples[wildcards.sample]["Matched normal"] + ".bam" if samples[wildcards.sample]["Matched normal"] != "." and samples[wildcards.sample]["Matched normal"] != "" else "",
            black_list= lambda wildcards: config["pipeline_home"] + "/" + config[samples[wildcards.sample]["Genome"]]["black_list"],
            exclude_chr_list = lambda wildcards: config["pipeline_home"] + "/" + config[samples[wildcards.sample]["Genome"]]["exclude_chr_list"],
            log_dir = lambda wildcards: wildcards.sample + '/log',
            version_bedtools=config["version_common"]["bedtools"],
    shell:
            """  
            module load macs/{version}
            module load bedtools/{params.version_bedtools}
            dup=`cut -f3 {wildcards.sample}/qc/mapping_summary.txt`
            keep_dup="all"
            if (( $(echo "$dup > {params.dup_cutoff}" |bc -l) )); then
                keep_dup="1"
            fi            
            macs2 callpeak {params.peak_type} -t {wildcards.sample}/{wildcards.sample}.bam {params.input_control} --nomodel  -g {params.genome} -f {params.bam_format}  --keep-dup ${{keep_dup}} --name {wildcards.sample} --outdir {wildcards.sample}/MACS_Out_{wildcards.cutoff_type}_{wildcards.cutoff} {params.ext_option} {params.cutoff_option}
            #remove unused chromosomes and blacklist in peak
            {params.pipeline_home}/scripts/removeChr.pl -i {output.peak} -o {output.renamed_peak} -l {params.exclude_chr_list}
            sed -i 's/{wildcards.sample}//g' {output.renamed_peak}
            bedtools intersect -a {output.renamed_peak} -b {params.black_list} -v > {output.nobl}
            cut -f1-3 {output.nobl} > {output.great}            
            """

rule plotFingerPrint:
    input:
            lambda wildcards: MACS_BAMS[wildcards.sample]
    output:
            pdf="{sample}/qc/fingerPrint.pdf",
            table="{sample}/qc/fingerPrint.tab"
    version:
            config["version"]["deeptools"],
    benchmark:
            "{sample}/benchmark/plotFingerPrint.{sample}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_macs2"],
            pipeline_home=config["pipeline_home"],
            rulename = "plotFingerPrint",
            bam_list = lambda wildcards: ' '.join(MACS_BAMS[wildcards.sample]),
            ext_len = lambda wildcards: str(samples[wildcards.sample]["LibrarySize"]),
            log_dir = lambda wildcards: wildcards.sample + '/log'            
    shell:
            """  
            module load deeptools/{version}
            plotFingerprint \
                -b {params.bam_list} \
                --smartLabels \
                --minMappingQuality 30 --skipZeros \
                -e {params.ext_len} \
                -T "Fingerprints"  \
                -p ${{THREADS}} \
                --plotFile {output.pdf} \
                --outRawCounts {output.table}
            """ 
rule crossCorrelation:
    input:
            "{sample}/{sample}.bam"
    output:
            pdf="{sample}/qc/{sample}.spp.pdf",
            txt="{sample}/qc/{sample}.spp.txt"
    version:
            config["version_common"]["R"]
    benchmark:
            "{sample}/benchmark/crossCorrelation.{sample}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster_common"]["medium"],
            pipeline_home=config["pipeline_home"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            rulename = "crossCorrelation",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            version_samtools=config["version_common"]["samtools"],
    shell:
            """
            module load R/{version}
            module load samtools/{params.version_samtools}
            Rscript {params.pipeline_home}/scripts/run_spp.R -c={input} -odir={wildcards.sample}/spp -savp={output.pdf} -out={output.txt}
            """

rule makeNormTDF:
    input:
            "{sample}/{sample}.{bin_size}.bedgraph"
    output:
            "{sample}/{sample}.{bin_size}.{norm_type}.tdf"
    version:
            config["version"]["igvtools"]
    benchmark:
            "{sample}/benchmark/makeNormTDF.{sample}.{bin_size}.{norm_type}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster_common"]["small"],
            pipeline_home=config["pipeline_home"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            rulename = "makeNormTDF",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            dup_cutoff = config["dup_cutoff"],
    shell:
            """  
            module load igvtools/{version}            
            if [[ "{wildcards.norm_type}" == "scaled" ]]; then
                norm_count=`tail -1 {wildcards.sample}/SpikeIn/spike_map_summary | cut -f2`
            else
                norm_count=`grep mapped {wildcards.sample}/{wildcards.sample}.flagstat.txt | head -1 | cut -d' ' -f1`
            fi
            factor=`echo "1000000/${{norm_count}}" | bc -l`
            {params.pipeline_home}/scripts/scaleTDF.pl -i {input} -o {output} -f ${{factor}} -g {params.genome} -t bedgraph
            """
            
rule makeTDF:
    input:
            "{sample}/{sample}.bam"
    output:
            bg="{sample}/{sample}.{bin_size}.bedgraph"
    version:
            config["version"]["igvtools"]
    benchmark:
            "{sample}/benchmark/makeTDF.{sample}.{bin_size}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster_common"]["small"],
            lib_size = lambda wildcards: samples[wildcards.sample]["LibrarySize"],
            smooth_size = config["smooth_size"],            
            pair = lambda wildcards: " --pairs " if samples[wildcards.sample]["PE"] else "",
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            rulename = "makeTDF",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            dup_cutoff = config["dup_cutoff"],
    shell:
            """  
            module load igvtools/{version}
            dup=`cut -f3 {wildcards.sample}/qc/mapping_summary.txt`
            keep_dup="--includeDuplicates"
            if (( $(echo "$dup > {params.dup_cutoff}" |bc -l) )); then
                keep_dup=""
            fi   
            igvtools count -w {wildcards.bin_size} -e {params.lib_size} ${{keep_dup}} {params.pair}{input} {wildcards.sample}/{wildcards.sample}.{wildcards.bin_size}.tdf {params.genome}
            igvtools tdftobedgraph {wildcards.sample}/{wildcards.sample}.{wildcards.bin_size}.tdf {output.bg}
            """

rule makeBigWig:
    input:
            "{sample}/{sample}.bam"
    output:
            "{sample}/{sample}.{bin_size}.{norm_type}.bw"
    version:
            config["version"]["deeptools"]
    benchmark:
            "{sample}/benchmark/makeBigWig.{sample}.{bin_size}.{norm_type}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster_common"]["medium"],
            lib_size = lambda wildcards: samples[wildcards.sample]["LibrarySize"],
            smooth_size = config["smooth_size"],
            rulename = "makeBigWig",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            dup_cutoff = config["dup_cutoff"],
    shell:
            """  
            module load deeptools/{version}
            dup=`cut -f3 {wildcards.sample}/qc/mapping_summary.txt`
            keep_dup=""
            if (( $(echo "$dup > {params.dup_cutoff}" |bc -l) )); then
                keep_dup="--ignoreDuplicates"
            fi               
            norm_param=" --normalizeUsing CPM "
            if [[ "{wildcards.norm_type}" == "scaled" ]]; then
                spike_count=`tail -1 {wildcards.sample}/SpikeIn/spike_map_summary | cut -f2`
                ref_count=`tail -1 {wildcards.sample}/SpikeIn/spike_map_summary | cut -f1`
                factor=`echo "1000000/${{spike_count}}" | bc -l`
                norm_param=" --scaleFactor ${{factor}} "
            fi
            bamCoverage -b {input} -o {output} -of bigwig -p ${{THREADS}} --smoothLength {params.smooth_size} ${{keep_dup}} --binSize {wildcards.bin_size} -e {params.lib_size} ${{norm_param}}
            
            """
rule BWA:
    input:
            lambda wildcards: FASTQS[wildcards.sample]
    output: 
            bam="{sample}/{sample}.bam",
            bai="{sample}/{sample}.bam.bai",
    version:
            config["version"]["bwa"]
    params:
            bwa_idx = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["bwa_index"] + ".fa" if samples[wildcards.sample]["SpikeIn"] != "yes" else config[samples[wildcards.sample]["Genome"]]["bwa_index"] + "." + samples[wildcards.sample]["SpikeIn_genome"] + ".fa",
            fastqs = lambda wildcards: " ".join(work_dir + '/' + x for x in FASTQS[wildcards.sample]),
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_bwa"],
            rulename = "BWA",
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            spikeIn = lambda wildcards: samples[wildcards.sample]["SpikeIn"],
            SpikeIn_genome = lambda wildcards: samples[wildcards.sample]["SpikeIn_genome"],
            min_mapq = config["min_mapq"],            
            log_dir = lambda wildcards: wildcards.sample + '/log',
            version_samtools=config["version_common"]["samtools"],
            version_picard=config["version_common"]["picard"],
            version_deeptools=config["version"]["deeptools"],
    benchmark:
            "{sample}/benchmark/bwa.{sample}.benchmark.txt"
    shell:
            """  
            module load bwa/{version}
            module load samtools/{params.version_samtools}
            cd ${{LOCAL}}
            bwa mem -t ${{THREADS}} {params.bwa_idx} {params.fastqs} | samtools view -bhq {params.min_mapq} | samtools sort -@ ${{THREADS}} -o {wildcards.sample}.bam -O bam -
            samtools index -@ ${{THREADS}} {wildcards.sample}.bam
            module load picard/{params.version_picard}
            module load deeptools/{params.version_deeptools}
            mkdir -p {wildcards.sample}
            if [[ "{params.spikeIn}" == "yes" ]];then
                for c in `samtools idxstats {wildcards.sample}.bam | cut -f1 | grep -v '*'`;do 
                    samtools view {wildcards.sample}.bam ${{c}} -b > {wildcards.sample}.${{c}}.bam &
                done
                wait
                echo "done split finished"
                mkdir -p {params.work_dir}/{wildcards.sample}/SpikeIn
                samtools merge -@ ${{THREADS}} -f {wildcards.sample}.{params.SpikeIn_genome}.bam {wildcards.sample}.{params.SpikeIn_genome}_*.bam
                samtools merge -@ ${{THREADS}} -f {wildcards.sample}.{params.genome}.bam {wildcards.sample}.chr*.bam
                samtools index -@ ${{THREADS}} {wildcards.sample}.{params.genome}.bam
                samtools index -@ ${{THREADS}} {wildcards.sample}.{params.SpikeIn_genome}.bam
                java -jar $PICARDJARPATH/picard.jar MarkDuplicates AS=true M={wildcards.sample}.{params.SpikeIn_genome}.metrix.txt I={wildcards.sample}.{params.SpikeIn_genome}.bam O={wildcards.sample}.{params.SpikeIn_genome}.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
                java -jar $PICARDJARPATH/picard.jar MarkDuplicates AS=true M={wildcards.sample}.{params.genome}.metrix.txt I={wildcards.sample}.{params.genome}.bam O={wildcards.sample}.{params.genome}.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
                samtools index -@ ${{THREADS}} {wildcards.sample}.{params.SpikeIn_genome}.dd.bam
                samtools index -@ ${{THREADS}} {wildcards.sample}.{params.genome}.dd.bam
                samtools flagstat -@ ${{THREADS}} {wildcards.sample}.{params.SpikeIn_genome}.dd.bam > {wildcards.sample}.{params.SpikeIn_genome}.flagstat.txt
                samtools flagstat -@ ${{THREADS}} {wildcards.sample}.{params.genome}.dd.bam > {wildcards.sample}.{params.genome}.flagstat.txt
                mapped_spikeIn=`grep mapped {wildcards.sample}.{params.SpikeIn_genome}.flagstat.txt | head -1 | cut -d' ' -f1`
                mapped_ref=`grep mapped {wildcards.sample}.{params.genome}.flagstat.txt | head -1 | cut -d' ' -f1`
                echo -e "{params.genome}\t{params.SpikeIn_genome}\n${{mapped_ref}}\t${{mapped_spikeIn}}" > {params.work_dir}/{wildcards.sample}/SpikeIn/spike_map_summary
                #mv {wildcards.sample}.{params.SpikeIn_genome}.* {params.work_dir}/{wildcards.sample}/SpikeIn/
                #mv {wildcards.sample}.{params.genome}.bam {params.work_dir}/{wildcards.sample}/{wildcards.sample}.bam
                #mv {wildcards.sample}.{params.genome}.bam.bai {params.work_dir}/{wildcards.sample}/{wildcards.sample}.bam.bai
                
            else
                java -jar $PICARDJARPATH/picard.jar MarkDuplicates AS=true M={wildcards.sample}.{params.genome}.metrix.txt I={wildcards.sample}.bam O={wildcards.sample}.{params.genome}.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT                                
                samtools index -@ ${{THREADS}} {wildcards.sample}.{params.genome}.dd.bam
                samtools flagstat -@ ${{THREADS}} {wildcards.sample}.{params.genome}.dd.bam > {wildcards.sample}.{params.genome}.flagstat.txt
                #mv {wildcards.sample}.bam {params.work_dir}/{wildcards.sample}
                #mv {wildcards.sample}.bam.bai {params.work_dir}/{wildcards.sample}
            fi
            #move output files
            bam_f={params.work_dir}/{wildcards.sample}/{wildcards.sample}.bam            
            flag_f={params.work_dir}/{wildcards.sample}/{wildcards.sample}.flagstat.txt
            cp {wildcards.sample}.{params.genome}.dd.bam ${{bam_f}}
            cp {wildcards.sample}.{params.genome}.dd.bam.bai ${{bam_f}}.bai
            cp {wildcards.sample}.{params.genome}.metrix.txt {params.work_dir}/{wildcards.sample}/{wildcards.sample}.metrix.txt
            cp {wildcards.sample}.{params.genome}.flagstat.txt ${{flag_f}}
            
            dup=`grep duplicate ${{flag_f}} | cut -d' ' -f1`
            mapped=`grep 'mapped (' ${{flag_f}} | cut -d' ' -f1`
            pect=`echo $dup/$mapped | bc -l`
            echo -e "${{flag_f}}\t${{dup}}\t${{mapped}}\t${{pect}}" > {params.work_dir}/{wildcards.sample}/qc/mapping_summary.txt
            """

rule fastqc:
    input:
            lambda wildcards: wildcards.sample + "/DATA/" + wildcards.sample + "_" + wildcards.suffix + ".fastq.gz"
    output:
            "{sample}/qc/{sample}_{suffix}_fastqc.html"
    params:
            work_dir = config["work_dir"],
            rulename = "fastqc",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            batch    = config["cluster_common"]["medium"]
    benchmark:
            "{sample}/benchmark/fastqc.{suffix}.benchmark.txt"
    version:
            config["version_common"]["fastqc"]
    shell:
            """
            module load fastqc/{version}
            mkdir -p {wildcards.sample}/qc
            fastqc -t ${{THREADS}} -o {wildcards.sample}/qc {input}
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
            if [[ ! -L {output} ]];then
                ln -s {input} {output}
            fi
            """