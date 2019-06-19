workflow myWorkflow {
    
    # java -jar -Dconfig.file=application.conf ./cromwell/cromwell-41.jar run bigdata.wdl --options options.json --inputs inputs.json
    
    File fabz2NR1
    File fabz2NR2
    File fabz2TR1
    File fabz2TR2
    String ReadGroupHeader1
    String ReadGroupHeader2    
    
    call TaskBWA as TaskBWAN { input: fqbz2R1=fabz2NR1, fqbz2R2=fabz2NR2, ReadGroupHeader=ReadGroupHeader1 }
    call TaskBWA as TaskBWAT { input: fqbz2R1=fabz2TR1, fqbz2R2=fabz2TR2, ReadGroupHeader=ReadGroupHeader2 }
    call TaskSamtoolsSordex as TaskSamtoolsSordexN { input: bam_in=TaskBWAN.bam }
    call TaskSamtoolsSordex as TaskSamtoolsSordexT { input: bam_in=TaskBWAT.bam }
    
    call TaskSamtoolStats as TaskSamtoolStatsN { input:
        bam=TaskSamtoolsSordexN.bam,
        bai=TaskSamtoolsSordexN.bai,
    }
    call TaskSamtoolStats as TaskSamtoolStatsT { input:
        bam=TaskSamtoolsSordexT.bam,
        bai=TaskSamtoolsSordexT.bai,
    }
    
    call TaskMantaSomatic { input:
        bamN=TaskSamtoolsSordexN.bam,
        baiN=TaskSamtoolsSordexN.bai,
        bamT=TaskSamtoolsSordexT.bam,
        baiT=TaskSamtoolsSordexT.bai,
    }
    
    call TaskStrelkaSomatic { input:
        bamN=TaskSamtoolsSordexN.bam,
        baiN=TaskSamtoolsSordexN.bai,
        bamT=TaskSamtoolsSordexT.bam,
        baiT=TaskSamtoolsSordexT.bai,
        vcfIndel=TaskMantaSomatic.vcf,
        tbiIndel=TaskMantaSomatic.tbi,
    }
    
    call TaskVarScanCNV { input:
        bamN=TaskSamtoolsSordexN.bam,
        baiN=TaskSamtoolsSordexN.bai,
        bamT=TaskSamtoolsSordexT.bam,
        baiT=TaskSamtoolsSordexT.bai,
    }
    
    
    output {
        File bamN = TaskSamtoolsSordexN.bam
        File baiN = TaskSamtoolsSordexN.bai
        File bamT = TaskSamtoolsSordexT.bam
        File baiT = TaskSamtoolsSordexT.bai
        File bamplotN = TaskSamtoolStatsN.bamplot
        File bamplotT = TaskSamtoolStatsT.bamplot
        File snv = TaskStrelkaSomatic.vcf
        File tbi = TaskStrelkaSomatic.tbi
        File cnv = TaskVarScanCNV.cnv
        File cnv_seg = TaskVarScanCNV.cnv_seg
        File cnv_merged = TaskVarScanCNV.cnv_merged
        File cnv_seg_png_p = TaskVarScanCNV.cnv_seg_png_p
        File cnv_seg_png_s = TaskVarScanCNV.cnv_seg_png_s
        File cnv_seg_png_w = TaskVarScanCNV.cnv_seg_png_w
    }
    
}


task TaskBWA {

    File fqbz2R1
    File fqbz2R2
    String ReadGroupHeader

    command {
        bunzip2 -c ${fqbz2R1} > R1.fq
        bunzip2 -c ${fqbz2R2} > R2.fq
        bwa mem -M -R "${ReadGroupHeader}" -t 8 -K 10000000 /data/Reference/GRCh37-lite.fa R1.fq R2.fq | samblaster -M -e | samtools view -bS -o bwa.bam -
    }
    
    output {
        File bam = "bwa.bam"
    }
}


task TaskSamtoolsSordex {
    
    File bam_in
    
    command {
        samtools sort -@ 8 -o sort.bam ${bam_in}
        samtools index -@ 8 sort.bam
    }
    
    output {
        File bam = "sort.bam"
        File bai = "sort.bam.bai"
    }

}


task TaskSamtoolStats {
    
    File bam
    File bai

    command {
        samtools stats -@ 8 -r /data/Reference/GRCh37-lite.fa ${bam} > bam.stats
        plot-bamstats -p bamplot/ bam.stats
        tar zcf bamplot.tar.gz bamplot/
    }
    
    output {
        File bamplot = "bamplot.tar.gz"
    }

}


task TaskMantaSomatic {
    
    File bamN
    File bamT
    File baiN
    File baiT
    
    command {
        configManta.py --normalBam=${bamN} --tumorBam=${bamT} --referenceFasta=/data/Reference/GRCh37-lite.fa --exome --runDir=MantaSomaticWorkflow
        python MantaSomaticWorkflow/runWorkflow.py --mode=local --jobs=8
    }
    
    output {
        File vcf = "MantaSomaticWorkflow/results/variants/candidateSmallIndels.vcf.gz"
        File tbi = "MantaSomaticWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi"
    }
    
}


task TaskStrelkaSomatic {
    
    File bamN
    File bamT
    File baiN
    File baiT
    File vcfIndel
    File tbiIndel
    
    command {
        configureStrelkaSomaticWorkflow.py --normalBam=${bamN} --tumorBam=${bamT} --referenceFasta=/data/Reference/GRCh37-lite.fa --indelCandidates=${vcfIndel} --exome --runDir=StrelkaSomaticWorkflow
        python StrelkaSomaticWorkflow/runWorkflow.py --mode=local --jobs=8
    }
    
    output {
        File vcf = "StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"
        File tbi = "StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz.tbi"
    }
    
}


task TaskVarScanCNV {
    
    File bamN
    File bamT
    File baiN
    File baiT
    
    command {
        samtools mpileup -q 1 -f /data/Reference/GRCh37-lite.fa ${bamN} ${bamT} | java -jar $VARSCAN_PATH/VarScan.v2.3.9.jar copynumber - varScan --mpileup 1
        java -jar $VARSCAN_PATH/VarScan.v2.3.9.jar copyCaller varScan.copynumber --output-file varScan.copynumber.called
        Rscript /data/Scripts/CBS.R
        perl /data/Scripts/mergeSegments.pl varScan.copynumber.called.seg  --ref-arm-sizes /data/Scripts/armsize.txt --output-basename varScan.copynumber.called.seg.merged
    }
    
    output {
        File cnv = "varScan.copynumber"
        File cnv_seg = "varScan.copynumber.called.seg"
        File cnv_merged = "varScan.copynumber.called.seg.merged.events.tsv"
        File cnv_seg_png_p = "varScan.copynumber.called.seg.p.png"
        File cnv_seg_png_s = "varScan.copynumber.called.seg.s.png"
        File cnv_seg_png_w = "varScan.copynumber.called.seg.w.png"
    }
    
}
