configfile: 'Arch_pipe.config.yaml'
import os

#chromosomes = list(range(1,23))
#chromosomes.append("X")
chromosomes = ["X"]
#chromosomes = [22]

sample_lists = config['sample_lists']

arch_AN_max = os.popen('cat '+config['sample_lists']['arch']+' | wc -l').read().strip() # count max possible arch
arch_AN_max = int(arch_AN_max)*2


def get_bcftools_add_ID_input(wildcards):
    return {'vcf': config['base_vcf'].format(chrom=wildcards.chrom)} # do in this way in learning purposes
    
#def get_seq_sim_input(wildcards):
#    return {'arch_vcf': config['arch_vcf'].format(wildcards)}

def get_seq_sim_input(wildcards):
    return {'arch_vcf': config['arch_vcf'].format(chrom=wildcards.chrom)}

#def get_seq_sim_input(wildcards):
#    return config['arch_vcf'][wildcards]


# for this last rule, I created an env with umap libraray installed, and had not better idea to attach it
# for all R scripts, I installed localy two base libraries (rlang,lifecycle) and attach them from home dir

rule plot_total_umap:
    input:
        "res_concat/seq_sim.collapsed_annotated.tsv.gz"
    output:
        "res_concat/seq_sim.collapsed_annotated.full_plot.umap.png","res_concat/seq_sim.collapsed_annotated.introgres_plot.umap.png"
    envmodules:
        "python/3.9.12","py-mpmath/1.2.1","py-numpy/1.22.4","py-pandas/1.4.2","py-seaborn/0.11.2","py-scikit-learn/1.1.1","py-matplotlib/3.5.2"
    resources:
        mem='20G', time=90
    shell:
        """
        source ~/virr_envs/umap_env/bin/activate; python3 umap_plot.py {input}
	      """
	      
	      
# I just put one file here, while it concated all that matters
rule unite_all_chr:
    input:
        expand("chr_plots/chr{chrom}.chunks.pdf",chrom=chromosomes)
    output:
        "res_concat/seq_sim.collapsed_annotated.tsv.gz"
    params:
      chroms_string = ','.join(map(str,chromosomes))
    envmodules:
        "htslib/1.16"
    shell:
        """
        bash concat_main_results.sh {params.chroms_string}
	      """


rule r2_heatmaps_plots:
    input:
        "res/chr{chrom}.haplotypes.R2_0_8.ILS.seq_sim.tsv","res/chr{chrom}.aSNPs.R2_0_8.tsv"
    output:
        "chr_plots/chr{chrom}.chunks.pdf"
    envmodules:
        "r/4.1.3","r-data-table/1.14.2-r4.1.3","r-dplyr/1.0.7","r-r-utils/2.11.0","r-tidyr/1.1.4","r-gridextra/2.3","r-ggplot2/3.3.5"
    resources:
        mem='50G', time=900
    shell:
        """
        Rscript snippet.plot_chr_R2.R {wildcards.chrom}
	      """


rule seq_sim_wrap_up:
    input:
        file_seq="res/chr{chrom}.seq_sim.tsv",
        arch=sample_lists['arch'],
        target=sample_lists['target']
    output:
        "res/chr{chrom}.seq_sim.collapsed_annotated.tsv","res/chr{chrom}.haplotypes.R2_0_8.ILS.seq_sim.tsv"
    envmodules:
        "r/4.1.3","r-data-table/1.14.2-r4.1.3","r-dplyr/1.0.7","r-r-utils/2.11.0","r-tidyr/1.1.4"
    params:
        denis_sample_name=config['denisova']
    shell:
        """
        Rscript snippet.seq_sim.wrapping_up.R {wildcards.chrom} {input.arch} {input.target} {params.denis_sample_name}
	      """


rule seq_sim:
    input:
        unpack(get_seq_sim_input),
        #get_seq_sim_input,
        #= lambda wildcards: identify_read_groups(f'cram/{wildcards.sample}.bam.cram')
        #arch_vcf = config['arch_vcf'].format(chrom=wildcards.chrom),
        file_ils="res/chr{chrom}.haplotypes.R2_0_8.ILS.tsv",
        arch=sample_lists['arch'],
        target=sample_lists['target'],
        random10_outgr="aux_lists/out_random_10.RS_out.samples.tmp"
    output:
        "res/chr{chrom}.seq_sim.tsv","res/chr{chrom}.seq_sim.wide_debug.tsv"
    log: "logs/seq_sim.{chrom}.log"
    envmodules:
        "bcftools/1.16","r/4.1.3","r-data-table/1.14.2-r4.1.3","r-dplyr/1.0.7","r-r-utils/2.11.0","r-tidyr/1.1.4"
    resources:
        mem='20G', time=2400
    shell:
        """
        (Rscript snippet.seq_sim.R {wildcards.chrom} {input.arch} {input.random10_outgr} {input.target} {input.arch_vcf}) 2>> {log}
	      """


rule ILS_test:
    input:
        file_hap="res/chr{chrom}.haplotypes.R2_0_8.tsv",
        decode=config['rec_maps']['decode'],
        AAmap=config['rec_maps']['AAmap']
    output:
        "res/chr{chrom}.haplotypes.R2_0_8.ILS.tsv","res/chr{chrom}.ILS_raw_values.tsv"
    envmodules:
        "r/4.1.3","r-data-table/1.14.2-r4.1.3","r-dplyr/1.0.7","r-r-utils/2.11.0"
    shell:
        """
        Rscript snippet.ILS_annotating.R {wildcards.chrom} {input.decode} {input.AAmap}
        """


rule collapse_haplo_R2_based:
    input:
        "aux_r2/chr{chrom}.hap.ld.gz"
    output:
        "res/chr{chrom}.haplotypes.R2_0_8.tsv","res/chr{chrom}.aSNPs.R2_0_8.tsv",
        "res/chr{chrom}.linked_haplo_R2_0.5.development.tsv","res/chr{chrom}.shared_positions_R2_0.5.development.tsv"
    envmodules:
        "r/4.1.3","r-data-table/1.14.2-r4.1.3","r-dplyr/1.0.7","r-r-utils/2.11.0"
    resources:
        mem='10G', time=360
    shell:
        """
        Rscript snippet.collapse_R2.R {wildcards.chrom}
        """


rule vcftools_count_r2:
    input:
        file_r2="aux_r2/chr{chrom}.target.for_R2.bcf"        
    output:
        "aux_r2/chr{chrom}.hap.ld.gz"
    log: "logs/vcftools_count_r2.{chrom}.log"
    envmodules:
        "vcftools/0.1.14","htslib/1.16"
    resources:
        mem='20G', time=900
    shell:
        """
        (vcftools --bcf {input.file_r2} --hap-r2 --min-r2 0.05 --out aux_r2/chr{wildcards.chrom} && 
        bgzip -c aux_r2/chr{wildcards.chrom}.hap.ld > {output} && rm aux_r2/chr{wildcards.chrom}.hap.ld) 2>> {log}
	      """
	      

rule bcftools_extract_map:
    input:
        map_id_list="res/chr{chrom}.putative_archaic_variants.list",
        target=sample_lists['target'],
        file="data_id/chr{chrom}.original_vcf_with_IDs.bcf"        
    output:
        "aux_r2/chr{chrom}.target.for_R2.bcf"
    envmodules:
        "bcftools/1.16"
    shell:
        """
      bcftools view --include ID==@{input.map_id_list} \
      -S {input.target} {input.file} -Ob -o {output}
	      """


rule combining_putative_map:
    input:
        "res/chr{chrom}.target.info.tsv","res/chr{chrom}.outgroup.info.tsv","res/chr{chrom}.arch.info.tsv"
    output:
        "res/chr{chrom}.putative_archaic_variants.list","res/chr{chrom}.putative_archaic_variants.tsv"
    envmodules:
        "r/4.1.3","r-data-table/1.14.2-r4.1.3","r-dplyr/1.0.7","r-r-utils/2.11.0","r-tidyr/1.1.4"
    shell:
        """
        Rscript snippet.unit_putative_aSNPs.R {wildcards.chrom}
        """


rule bcftools_putative_map_info:
    input:
        file="data_id/chr{chrom}.original_vcf_with_IDs.bcf",
        arch=sample_lists['arch'],
        outgroup=sample_lists['outgroup'],
        target=sample_lists['target']
    output:
        temp("res/chr{chrom}.target.info.tsv"),temp("res/chr{chrom}.outgroup.info.tsv"),temp("res/chr{chrom}.arch.info.tsv")
    params: max_arch_an=arch_AN_max
    envmodules:
        "bcftools/1.16"
    shell:
        """
        bash bcftools_putative_map.S.sh {input.file} {wildcards.chrom} {input.arch} {input.outgroup} {input.target} {params.max_arch_an}
        """


rule bcftools_add_ID:
    input:
        vcf="data/chr{chrom}.Arch.out_target.bcf"
    output:
        "data_id/chr{chrom}.original_vcf_with_IDs.bcf"
    envmodules:
        "bcftools/1.16"
    shell:
        """
	      bcftools annotate -x ID {input.vcf} | bcftools annotate --set-id +'%CHROM:%POS:%REF:%ALT' \
        -Ob -o {output} && bcftools index {output}
	      """

rule bcftools_preprocess:
    input:
        unpack(get_bcftools_add_ID_input),
        smpl_file="aux_lists/outgroup_target.samples.tmp"
    output:
        "data/chr{chrom}.Arch.out_target.bcf"
    envmodules:
        "bcftools/1.16"
    resources:
        mem='30G', time=1200
    shell:
        """
	      bash preprocess.1kg.sh {wildcards.chrom} {input.smpl_file}
	      """

rule aux_lists:
    input:
        arch=sample_lists['arch'],
        outgroup=sample_lists['outgroup'],
        target=sample_lists['target']
    output:
        out_target_list="aux_lists/outgroup_target.samples.tmp",
        rand_out_list="aux_lists/out_random_10.RS_out.samples.tmp",
        seq_sim_inside_bcftools="aux_lists/arch_10out_target.samples.tmp"
    shell:
        """
        cat {input.arch} {input.outgroup} > aux_lists/arch_outgroup.samples.tmp;
        cat {input.outgroup} {input.target} > {output.out_target_list};
        shuf --random-source {input.outgroup} {input.outgroup} | head -n10 > {output.rand_out_list};
        cat {input.arch} aux_lists/out_random_10.RS_out.samples.tmp {input.target} > {output.seq_sim_inside_bcftools} 
        """
