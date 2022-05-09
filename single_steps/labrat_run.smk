rule labrat_salmon_align:
    input:
        bam = STAR_OUTDIR + "{sample}.Aligned.sorted.out.bam",
        idx = STAR_OUTDIR + "{sample}.Aligned.sorted.out.bam.bai"

    output:
        os.path.join(RSEQC_OUTDIR, "{sample}.geneBodyCoverage.txt")

    params:
        # samples = lambda wildcards, input: ",".join(input.bams),
        librarytype = config[librarytype]

    conda:
        "../env/labratenv.yml"
    
    threads:
        4
    shell:
        """
        LABRAT.py --mode runSalmon\
        --librarytype {params.librarytype} \
        --txfasta <output of makeTFfasta> \
        --reads1 <comma separated list of forward read files> \
        --reads2 <Optional, comma separated list of reverse read files>\
        --samplename <comma separated list of sample names> \
        --threads {threads}
        """
