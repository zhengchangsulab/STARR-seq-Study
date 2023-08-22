#!/bin/bash

crm="ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_36.clean.sort"
gene_tss_info="ensemble_info_2021_5_15.tss.gene_info_clean.bed.sort"
dataset="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge"
non_candidate="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge.non_candidate"


function run_process_data(){
    cp ../count_total_bp.py .
    cp ../filter_chrosome.py .
    folder=$1
    cell_name=${folder/"_STARR"/}
    echo -e ${cell_name}
    
    starr_raw=${cell_name}"_STARR.consolidate.bed"
    starr_sort=${starr_raw}".sort"
    sort -k1,1 -k2,2n ${starr_raw} > ${starr_sort}
    python filter_chrosome.py ${starr_sort}
    
}


function run_process_overlap(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    starr=${cell_name}"_STARR.consolidate.bed.sort"


    bedtools intersect -a ${starr}.filter -b ../${crm} -f 0.5 -F 0.5 -e -wa > ${starr}.filter.fall_in_crm
    bedtools intersect -a ${starr}.filter.fall_in_crm -b ${cell_name}.active_predict.LR.bed -f 0.5 -F 0.5 -e -wa|sort|uniq > ${starr}.filter.fall_in_crm.overlap
    bedtools intersect -a ${cell_name}.active_predict.LR.bed -b ${starr}.filter.fall_in_crm -f 0.5 -F 0.5 -e -wa|sort|uniq > ${cell_name}.active_predict.LR.bed.overlap_with_STARR

    python count_total_bp.py ${starr}.filter.fall_in_crm
    python count_total_bp.py ${starr}.filter.fall_in_crm.overlap
    python count_total_bp.py ${cell_name}.active_predict.LR.bed.overlap_with_STARR
    python count_total_bp.py ${cell_name}.active_predict.LR.bed
}


function run_get_non_activate(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    active_crm=${cell_name}.active_predict.LR.bed
    non_active_crm=${cell_name}.non_active_predict.LR.bed
    bedtools intersect -a ../${crm} -b ${active_crm} -v > ${non_active_crm}
}

function select_enhancers(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    cp ../select_only_enhancers.py .
    python select_only_enhancers.py ${cell_name}
    echo -e ${cell_name}
}

function run_split_starr_abcdef(){
    cp ../count_total_bp.py .
    folder=$1
    cell_name=${folder/"_STARR"/}
    active_crm=${cell_name}.active_predict.LR.bed
    non_active_crm=${cell_name}.non_active_predict.LR.bed
    
    
    starr_peak_merge=${cell_name}"_STARR.consolidate.bed.sort.filter.merge"
    starr_peak=${cell_name}"_STARR.consolidate.bed.sort.filter"
    bedtools merge -i ${starr_peak} > ${starr_peak_merge}
    

    # get non activate CRMs
    #bedtools intersect -a ../${crm} -b ${active_crm} -v > ${non_active_crm}
    ## get a (bp)
    bedtools subtract -a ${non_active_crm} -b ${starr_peak_merge} > ${cell_name}".region.a"
    python count_total_bp.py ${cell_name}".region.a"
    # get non-active crms in a
    bedtools intersect -a ${non_active_crm} -b ${cell_name}".region.a" -f 1 -wa > ${non_active_crm}".fall_in_region_a"
    
    ## get b (bp)
    bedtools subtract -a ${active_crm} -b ${starr_peak_merge} > ${cell_name}".region.b"
    python count_total_bp.py ${cell_name}".region.b"
    # get active crms in b
    bedtools intersect -a ${active_crm} -b ${cell_name}".region.b" -f 1 -wa > ${active_crm}".fall_in_region_b"

    ## get c (bp)
    bedtools intersect -a ${active_crm} -b ${starr_peak_merge} > ${cell_name}".region.c"
    python count_total_bp.py ${cell_name}".region.c"
    # get starr peaks in c
    bedtools intersect -a ${starr_peak} -b ${cell_name}".region.c" -f 1 -wa > ${starr_peak}".fall_in_region_c"


    ## get d (bp)
    bedtools intersect -a ${non_active_crm} -b ${starr_peak_merge} > ${cell_name}".region.d"
    python count_total_bp.py ${cell_name}".region.d"
    # get starr peaks in d 
    bedtools intersect -a ${starr_peak} -b ${cell_name}".region.d" -f 1 -wa > ${starr_peak}".fall_in_region_d"

    ## get e (bp)
    bedtools intersect -a ../${non_candidate} -b ${starr_peak_merge} > ${cell_name}".region.e"
    python count_total_bp.py ${cell_name}".region.e"
    # get starr peaks in e
    bedtools intersect -a ${starr_peak} -b ${cell_name}".region.e" -f 1 -wa> ${starr_peak}".fall_in_region_e"

    starr_peak_number_in_e=(`wc -l ${starr_peak}".fall_in_region_e"`)
    
    ## get f (bp)
    bedtools subtract -a ../${non_candidate} -b ${starr_peak_merge} > ${cell_name}".region.all.f"
    bedtools intersect -a ../${non_candidate} -b ${cell_name}".region.all.f" -f 1 -wa > ${cell_name}".region.non_crmc_fall_in.all.f"
    shuf -n ${starr_peak_number_in_e[0]} ${cell_name}".region.non_crmc_fall_in.all.f" > ${cell_name}".region.f"
    cp ${cell_name}".region.f" ${cell_name}".non_candidates.fall_in_region_f"
    python count_total_bp.py ${cell_name}".region.f"
}

function run_count_bp_f(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    python count_total_bp.py ${cell_name}".region.all.f"

}

function run_reformat_starr(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    active_CRM=${cell_name}.active_predict.LR.bed
    starr=${cell_name}"_STARR.consolidate.bed.sort"
    awk 'BEGIN{print "CRM_id\tgene_id\tgene_symbol\tdistance"} {print $1":"$2"-"$3"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' ${starr}.filter.closest_tss > ${starr}.filter.closest_tss.reformat
}


function run_split_crms(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    active_CRM=${cell_name}.active_predict.LR.bed
    starr=${cell_name}"_STARR.consolidate.bed.sort"
    
    bedtools intersect -a ${active_CRM} -b ${starr}.filter -v > ${active_CRM}".not_ovelrap_with_STARR" #a
    bedtools intersect -a ${starr}.filter -b ${active_CRM} -f 1 -wa > ${starr}".filter.total_fall_in_active_crm" # b                                                                                                                         
    #c
    bedtools intersect -a ${starr}.filter -b ../${crm} -f 1 -wa > ${starr}".filter.total_fall_in_crm"
    bedtools intersect -a ${starr}".filter.total_fall_in_crm" -b ${active_CRM} -f 1 -v > ${starr}".filter.total_fall_in_crm.not_any_fall_in_active_crm"

    #d
    bedtools intersect -a ${starr}.filter -b ../${dataset} -f 1 -wa > ${starr}".filter.total_fall_in_dataset"
    bedtools intersect -a ${starr}".filter.total_fall_in_dataset" -b ../${crm} -f 1 -v > ${starr}".filter.total_fall_in_dataset.not_any_fall_in_crm"


    bedtools closest -a ../${crm} -b ../${gene_tss_info} -d > ${crm}.closest_tss
    bedtools closest -a ${starr}".filter" -b ../${gene_tss_info} -d > ${starr}.filter.closest_tss
    awk 'BEGIN{print "CRM_id\tgene_id\tgene_symbol\tdistance"} {print $1":"$2"-"$3"\t"$7"\t"$8"\t"$9}' ${crm}.closest_tss > ${crm}.closest_tss.reformat
    awk 'BEGIN{print "CRM_id\tgene_id\tgene_symbol\tdistance"} {print $1":"$2"-"$3"\t"$7"\t"$8"\t"$9}' ${starr}.filter.closest_tss > ${starr}.filter.closest_tss.reformat                                                                    


    awk '{print $1":"$2"-"$3}'  ${active_CRM}".not_ovelrap_with_STARR" >  ${active_CRM}".not_ovelrap_with_STARR.a"
    awk '{print $1":"$2"-"$3}' ${starr}".filter.total_fall_in_active_crm" > ${starr}".filter.total_fall_in_active_crm.b"
    awk '{print $1":"$2"-"$3}' ${starr}".filter.total_fall_in_crm.not_any_fall_in_active_crm" > ${starr}".filter.total_fall_in_crm.not_any_fall_in_active_crm.c"
    awk '{print $1":"$2"-"$3}' ${starr}".filter.total_fall_in_dataset.not_any_fall_in_crm" > ${starr}".filter.total_fall_in_dataset.not_any_fall_in_crm.d"
    awk 'BEGIN{print "CRM_id\tgene_id\tgene_symbol\tdistance"} {print $1":"$2"-"$3"\t"$14"\t"$15"\t"$16}' ${starr}.filter.closest_tss > ${starr}.filter.closest_tss.reformat
}




function run_sample_a_region(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    active_crm=${cell_name}.active_predict.LR.bed
    non_active_crm=${cell_name}.non_active_predict.LR.bed
    
    starr_peak_merge=${cell_name}"_STARR.consolidate.bed.sort.filter.merge"
    starr_peak=${cell_name}"_STARR.consolidate.bed.sort.filter"

    shuf -n 10000 ${non_active_crm}".fall_in_region_a" > ${non_active_crm}".fall_in_region_a.sample"
    shuf -n 10000 ${active_crm}".fall_in_region_b" > ${active_crm}".fall_in_region_b.sample"
    
    sort -k1,1 -k2,2n ${non_active_crm}".fall_in_region_a.sample" > ${non_active_crm}".fall_in_region_a.sample.sort"
    sort -k1,1 -k2,2n ${active_crm}".fall_in_region_b.sample" > ${active_crm}".fall_in_region_b.sample.sort"
    sort -k1,1 -k2,2n ${cell_name}".non_candidates.fall_in_region_f" > ${cell_name}".non_candidates.fall_in_region_f.sort"

}

function run_tss_gm12878(){
    starr_peak="GM12878_STARR.consolidate.bed.sort.filter"
    cd GM12878_STARR
    bedtools closest -a ${starr_peak} -b ../${gene_tss_info} -d > ${starr_peak}.closest_tss
    awk 'BEGIN{print "CRM_id\tgene_id\tgene_symbol\tdistance"} {print $1":"$2"-"$3"\t"$14"\t"$15"\t"$16}' ${starr_peak}.closest_tss > ${starr_peak}.closest_tss.reformat
    cd ..
}


function run_tss_non_candidate(){
    non_candidate="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge.non_candidate"    
    bedtools closest -a ${non_candidate} -b ${gene_tss_info} -d > ${non_candidate}.closest_tss
    awk 'BEGIN{print "CRM_id\tgene_id\tgene_symbol\tdistance"} {print $1":"$2"-"$3"\t"$7"\t"$8"\t"$9}' ${non_candidate}.closest_tss > ${non_candidate}.closest_tss.reformat
}


function run_get_conservation(){
    cp ../extract_binding_sites_conservation_0_base_drop_nan.py .
    files=*.region.[abcdef]*
    for f in ${files}
    do
	pbs="run_get_conservation_score_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "python extract_binding_sites_conservation_0_base_drop_nan.py ${f} hg38.phyloP100way.bw phylop100way" >> ${pbs}
	sbatch ${pbs}
    done
}

function run_reformat(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    cp ${cell_name}".region.f" ${cell_name}".non_candidates.fall_in_region_f"
    files=*".fall_in_region_"*
    for f in ${files}
    do
	awk '{print $1":"$2"-"$3}' ${f} > ${f}.reformat
    done

}

#run_tss_non_candidate


function run_plot_abcdef_6_cell_v2(){
    for score in phylop100way #gerp 
    do
	pbs="run_plot_abcdef_score_6_cell."${score}".v2.pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "python plot_conservation_score_density_abcdef_6_cell_v2.py True ${score}" >> ${pbs}
	sbatch ${pbs}
    done

}



function run_count_fall_bp(){
    cp ../run_sort_merge.sh .
    cp ../count_total_bp.py .
    files=*fall_in*[abcdef]
    for f in ${files}
    do
	pbs="run_sort_merge_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "bash run_sort_merge.sh ${f}" >> ${pbs}
	sbatch ${pbs}
    done
}


function run_get_multile_intersect(){
    for region in f #a b c d e f
    do
	pbs="run_get_abcedf_intersection.pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "bash run_get_intersect_abcdef_region.sh ${region}" >> ${pbs}
	sbatch ${pbs}
    done

}

function run_get_starr_fall_in_region(){
    folder=$1


    crm="ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_36.clean.sort"
    gene_tss_info="ensemble_info_2021_5_15.tss.gene_info_clean.bed.sort"
    dataset="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge"
    non_candidate="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge.non_candidate"

    cell_name=${folder/"_STARR"/}
    active_crm=${cell_name}.active_predict.LR.bed
    non_active_crm=${cell_name}.non_active_predict.LR.bed

    starr_peak_merge=${cell_name}"_STARR.consolidate.bed.sort.filter.merge"
    starr_peak=${cell_name}"_STARR.consolidate.bed.sort.filter"

    bedtools intersect -a ${active_crm} -b ${cell_name}".region.c" -f 1  -wa > ${active_crm}".fall_in_region_c"
    bedtools intersect -a ${non_active_crm} -b ${cell_name}".region.d" -f 1 -wa > ${non_active_crm}".fall_in_region_d"
    sort -k1,1 -k2,2n ${active_crm}".fall_in_region_c" > ${active_crm}".fall_in_region_c.sort"
    sort -k1,1 -k2,2n ${non_active_crm}".fall_in_region_d" > ${non_active_crm}".fall_in_region_d.sort"
}





function run_get_length(){
    cp ../extract_crm_length.py .
    files=*fall_in_region_*.sort
    for f in ${files}
    do
	python extract_crm_length.py ${f}
    done

}

function run_plot_enrichheatmap(){
    cp ../plot_enrichHeatmap.R .
    for mark in H3K9me3 H3K27me3 #H3K4me1 H3K27ac ATAC
    do
	regions=(`ls *fall_in*[abcde].sort|grep -v "merge"`)
	for region in ${regions[@]}
	do
	    pbs="run_plot_enrichheatmap_"${mark}"-"${region}".pbs"
	    echo -e ${pbs_header} > ${pbs}
	    echo -e "Rscript plot_enrichHeatmap.R ${region} ${mark}" >> ${pbs}
	    sbatch ${pbs}
	done
    done



}

##################################################################
function run_overlap_crm_starr_in_regions(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    active_crm=${cell_name}.active_predict.LR.bed
    non_active_crm=${cell_name}.non_active_predict.LR.bed
    
    starr_peak_merge=${cell_name}"_STARR.consolidate.bed.sort.filter.merge"
    starr_peak=${cell_name}"_STARR.consolidate.bed.sort.filter"
    bedtools merge -i ${starr_peak} > ${starr_peak_merge}
    

    # get non activate CRMs
    bedtools intersect -a ../${crm} -b ${active_crm} -v > ${non_active_crm}
    ## get a (bp)
    bedtools subtract -a ${non_active_crm} -b ${starr_peak_merge} > ${cell_name}".region.a"
    python count_total_bp.py ${cell_name}".region.a"
    # get non-active crms in a
    bedtools intersect -a ${non_active_crm} -b ${cell_name}".region.a" -f 0.5 -wa > ${non_active_crm}".fall_in_region_a.v0.5"
    
    ## get b (bp)
    bedtools subtract -a ${active_crm} -b ${starr_peak_merge} > ${cell_name}".region.b"
    python count_total_bp.py ${cell_name}".region.b"
    # get active crms in b
    bedtools intersect -a ${active_crm} -b ${cell_name}".region.b" -f 0.5 -wa > ${active_crm}".fall_in_region_b.v0.5"

    ## get c (bp)
    bedtools intersect -a ${active_crm} -b ${starr_peak_merge} > ${cell_name}".region.c"
    python count_total_bp.py ${cell_name}".region.c"
    # get starr peaks in c
    bedtools intersect -a ${starr_peak} -b ${cell_name}".region.c" -f 0.5 -wa > ${starr_peak}".fall_in_region_c.v0.5"
    bedtools intersect -a ${active_crm} -b ${cell_name}".region.c" -f 0.5  -wa > ${active_crm}".fall_in_region_c.v0.5"                                                                                                                              

    ## get d (bp)
    #bedtools intersect -a ${non_active_crm} -b ${starr_peak_merge} > ${cell_name}".region.d"
    #python count_total_bp.py ${cell_name}".region.d"
    # get starr peaks in d 
    bedtools intersect -a ${starr_peak} -b ${cell_name}".region.d" -f 0.5 -wa > ${starr_peak}".fall_in_region_d.v0.5"
    bedtools intersect -a ${non_active_crm} -b ${cell_name}".region.d" -f 0.5 -wa > ${non_active_crm}".fall_in_region_d.v0.5" 

    ## get e (bp)
    bedtools intersect -a ../${non_candidate} -b ${starr_peak_merge} > ${cell_name}".region.e"
    python count_total_bp.py ${cell_name}".region.e"
    # get starr peaks in e
    bedtools intersect -a ${starr_peak} -b ${cell_name}".region.e" -f 0.5 -wa > ${starr_peak}".fall_in_region_e.v0.5"

    #starr_peak_number_in_e=(`wc -l ${starr_peak}".fall_in_region_e"`)
    
    ## get f (bp)
    bedtools subtract -a ../${non_candidate} -b ${starr_peak_merge} > ${cell_name}".region.all.f"
    bedtools intersect -a ../${non_candidate} -b ${cell_name}".region.all.f" -f 1 -wa > ${cell_name}".region.non_crmc_fall_in.all.f"
    shuf -n ${starr_peak_number_in_e[0]} ${cell_name}".region.non_crmc_fall_in.all.f" > ${cell_name}".region.f"
    cp ${cell_name}".region.f" ${cell_name}".non_candidates.fall_in_region_f"
    python count_total_bp.py ${cell_name}".region.f"

}

function run_reformat_v0.5(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    #cp ${cell_name}".region.f" ${cell_name}".non_candidates.fall_in_region_f"
    files=*.fall_in_region*v0.5
    for f in ${files}
    do
	awk '{print $1":"$2"-"$3}' ${f} > ${f}.reformat
    done
}


function run_generate_random_peaks(){
    cp ../generate_random_crm_v3.py .
    files=(`ls *fall_in_region_{a,b,c,d,e,f}.sort`)
    for f in ${files[@]}
    do
	echo -e ${pbs_header} > run_generate_random_peaks.pbs
	echo -e "python generate_random_crm_v3.py ${f}" >> run_generate_random_peaks.pbs
	sbatch run_generate_random_peaks.pbs
    done
}

function run_intersect_TE(){
    files=(`ls *fall_in_region_{a,b,c,d,e,f}.sort`)
    for f in ${files[@]}
    do
	echo -e ${pbs_header} > run_intersect_te_${f}.pbs
	echo -e "bedtools intersect -a ${f} -b ../Repeats/hg38.fa.out.bed -F 1 -f 0.1 -wao > ${f}.overlap_Repeats" >> run_intersect_te_${f}.pbs 
	sbatch run_intersect_te_${f}.pbs 
    done

    files=(`ls *fall_in_region_{a,b,c,d,e,f}.sort.random_peak`)
    for f in ${files[@]}
    do
	echo -e ${pbs_header} > run_intersect_te_${f}.pbs
	echo -e "bedtools intersect -a ${f} -b ../Repeats/hg38.fa.out.bed -F 1 -f 0.1 -wao > ${f}.overlap_Repeats" >> run_intersect_te_${f}.pbs 
	sbatch run_intersect_te_${f}.pbs 
    done
}


function run_intersect_chromHMM(){
    folder=$1
    cell_name=${folder/"_STARR"/}


    files=(`ls *fall_in_region_{a,b,c,d,e,f}.sort`)
    for f in ${files[@]}
    do
	echo -e ${pbs_header} > run_intersect_chromhmm_${f}.pbs 
	echo -e "bedtools intersect -a ${f} -b ../ChromHMM/${cell_name}.bed -wao > ${f}.overlap_chromHMM" >> run_intersect_chromhmm_${f}.pbs 
	sbatch run_intersect_chromhmm_${f}.pbs 
    done

    files=(`ls *fall_in_region_{a,b,c,d,e,f}.sort.random_peak`)
    for f in ${files[@]}
    do
	echo -e ${pbs_header} > run_intersect_chromhmm_${f}.pbs 
	echo -e "bedtools intersect -a ${f} -b ../ChromHMM/${cell_name}.bed -wao > ${f}.overlap_chromHMM" >> run_intersect_chromhmm_${f}.pbs 
	sbatch run_intersect_chromhmm_${f}.pbs 
    done

}

function run_find_active_in_other_cells(){
    folder=$1
    cp ../run_find_activte_starr_in_other_cells.sh .
    for t in 0.5 0.6 0.7 0.8
    do
	bash run_find_activte_starr_in_other_cells.sh ${folder} ${t}
    done
}

function run_compute_gc(){
    cp ../compute_gc.py .
    
    files=(`ls *fall_in_region_{a,b,c,d,e,f}.sort`)
    for f in ${files[@]}
    do
	pbs="run_compute_gc_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "python compute_gc.py ${f}" >> ${pbs}
	sbatch ${pbs}
    done


    files=(`ls *fall_in_region_{a,b,c,d,e,f}.sort.random_peak`)
    for f in ${files[@]}
    do
	pbs="run_compute_gc_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "python compute_gc.py ${f}" >> ${pbs}
	sbatch ${pbs}
    done
    
}


function run_find_full_starr_in_other_cells(){
    folder=$1
    cp ../run_find_full_starr_in_other_cells.sh .
    for t in 0.5 0.6 0.7 0.8
    do
	bash run_find_full_starr_in_other_cells.sh ${folder} ${t}
    done
}


function run_get_tfbs(){
    files=(`ls *fall_in_region_{a,b,c,d,e}.sort`)
    for f in ${files[@]}
    do
	pbs="run_compute_tfbs_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "bedtools intersect -a ${f} -b ../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_0.cre.umotif.sort.uniq -F 1 -wao -sorted > ${f}.umotif.cre" >> ${pbs}
	echo -e "bedtools intersect -a ${f} -b ../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_0.cre.sort.uniq -F 1 -wao -sorted > ${f}.cre" >> ${pbs}
	#sbatch ${pbs}
    done


    files=(`ls *fall_in_region_{a,b,c,d,e}.sort.random_peak`)
    for f in ${files[@]}
    do
	pbs="run_compute_tfbs_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	#echo -e "sort -k1,1 -k2,2n ${f} > ${f}.sort" >> ${pbs}
	echo -e "bedtools intersect -a ${f}.sort -b ../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_0.cre.umotif.sort.uniq -F 1 -wao -sorted > ${f}.sort.umotif.cre" >> ${pbs}
	echo -e "bedtools intersect -a ${f}.sort -b ../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_0.cre.sort.uniq -F 1 -wao -sorted > ${f}.sort.cre" >> ${pbs}
	sbatch ${pbs}
    done

}


function run_get_tfbs_homer(){
    files=(`ls *fall_in_region_{a,b,c,d,e}.sort`)
    for f in ${files[@]}
    do
	pbs="run_compute_tfbs_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "bedtools intersect -a ${f} -b ../homer.KnownMotifs.hg38.191020.bed -F 1 -wao -sorted > ${f}.homer.cre" >> ${pbs}
	sbatch ${pbs}
    done

    files=(`ls *fall_in_region_{a,b,c,d,e}.sort.random_peak.sort`)
    for f in ${files[@]}
    do
	pbs="run_compute_tfbs_"${f}".pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e "bedtools intersect -a ${f} -b ../homer.KnownMotifs.hg38.191020.bed -F 1 -wao -sorted > ${f}.homer.cre" >> ${pbs}
	sbatch ${pbs}
    done

}


function run_create_tfbs_degree(){
    cp ../compute_tfbs_degree.py .

    for f in *sort.umotif.cre #*LR.bed.fall_in_region_b.sort.umotif.cre 
    do
	for ws in 5 10 15
	do
	    pbs="run_compute_tfbs_degree_"${f}"."${ws}".pbs"
	    echo -e ${pbs_header} > ${pbs}
	    echo -e "python compute_tfbs_degree.py ${f} ${ws}" >> ${pbs}
	    sbatch ${pbs}
	done

    done

}


function run_get_ctcf_summit(){

    cell_name=${folder/"_STARR"/}
    echo -e ${cell_name}
    cd CTCF_folder/
    cp ../../get_ctcf_peak_summit.py .
    python get_ctcf_peak_summit.py ${cell_name}
    echo -e ${cell_name}
    pwd
    cd ../
}


function run_overlap_summit(){
    cell_name=${folder/"_STARR"/}
    echo -e ${cell_name}
    cd CTCF_folder/
    cp ../../run_overlap_with_island.sh .

    summit=${cell_name}".ctcf.summit"
    pbs="run_overlap_with_tfbs_island.pbs"
    crm_island="../../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_0.cre.sort.uniq.merge"
    crm="../../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_36.clean.sort"
    echo -e ${pbs_header} > ${pbs}
    echo -e "bedtools intersect -a ${summit} -b ${crm_island} -f 1 -wo > ${summit}.fall_in_island" >> ${pbs}
    echo -e "bedtools intersect -a ${summit} -b ${crm} -f 1 -wo > ${summit}.fall_in_crm" >> ${pbs}
    echo -e "bash run_overlap_with_island.sh ${cell_name}" >> ${pbs}
    sbatch ${pbs}
    cd ../
}



function run_plotenrichheatmap(){
    folder=$1
    cp ../plot_enrichheatmap_v6.R .  
    cp ../plot_enrichheatmap_v6_no_atac.R . 
    cell_name=${folder/"_STARR"/}
    pbs="run_plot_enrichheatmap.pbs"
    echo -e ${pbs_header} > ${pbs}
    #echo -e "Rscript plot_enrichheatmap_v6.R ${cell_name}" >> ${pbs}
    echo -e "Rscript plot_enrichheatmap_v6_no_atac.R ${cell_name}" >> ${pbs} # for Hela-S3 and POU040MAI
    sbatch ${pbs}
}



function run_extract_length(){
    folder=$1
    cell_name=${folder/"_STARR"/}    
    cp ../extract_crm_length.py .
    starr=${cell_name}"_STARR.consolidate.bed.sort"
    python extract_crm_length.py ${starr}
    python extract_crm_length.py ${cell_name}.active_predict.LR.bed
}


function run_sort_uniq_fall_in_region_abcdef(){
    for f in *fall_in_region_*.v0.5
    do
	sort ${f}|uniq > ${f}.sort.uniq
    done
}


function run_get_length_crms_abcdef_v2(){
    cat *_STARR/*non_active_predict.LR.bed.fall_in_region_a.v0.5.sort.uniq > 8_cell_non_active_predict.LR.bed.fall_in_region_a.v0.5.sort
    cat *_STARR/*.active_predict.LR.bed.fall_in_region_b.v0.5.sort.uniq > 8_cell_active_predict.LR.bed.fall_in_region_b.v0.5.sort
    cat *_STARR/*active_predict.LR.bed.fall_in_region_c.v0.5.sort.uniq > 8_cell_active_predict.LR.bed.fall_in_region_c.v0.5.sort
    cat *_STARR/*non_active_predict.LR.bed.fall_in_region_d.v0.5.sort.uniq > 8_cell_non_active_predict.LR.bed.fall_in_region_d.v0.5.sort


    sort 8_cell_non_active_predict.LR.bed.fall_in_region_a.v0.5.sort|uniq > 8_cell_non_active_predict.LR.bed.fall_in_region_a.v0.5.sort.uniq
    sort 8_cell_active_predict.LR.bed.fall_in_region_b.v0.5.sort|uniq > 8_cell_active_predict.LR.bed.fall_in_region_b.v0.5.sort.uniq
    sort 8_cell_active_predict.LR.bed.fall_in_region_c.v0.5.sort|uniq > 8_cell_active_predict.LR.bed.fall_in_region_c.v0.5.sort.uniq
    sort 8_cell_non_active_predict.LR.bed.fall_in_region_d.v0.5.sort|uniq > 8_cell_non_active_predict.LR.bed.fall_in_region_d.v0.5.sort.uniq


    python extract_crm_length.py 8_cell_non_active_predict.LR.bed.fall_in_region_a.v0.5.sort.uniq
    python extract_crm_length.py 8_cell_active_predict.LR.bed.fall_in_region_b.v0.5.sort.uniq
    python extract_crm_length.py 8_cell_active_predict.LR.bed.fall_in_region_c.v0.5.sort.uniq
    python extract_crm_length.py 8_cell_non_active_predict.LR.bed.fall_in_region_d.v0.5.sort.uniq
}

function run_merge_crms_abcdef_regions(){
    cat *_STARR/*non_active_predict.LR.bed.fall_in_region_a|cut -f1,2,3 - > 8_cell_non_active_predict.LR.bed.fall_in_region_a.cat
    cat *_STARR/*.active_predict.LR.bed.fall_in_region_b|cut -f1,2,3 - > 8_cell_active_predict.LR.bed.fall_in_region_b.cat
    cat *_STARR/*active_predict.LR.bed.fall_in_region_c|cut -f1,2,3 - > 8_cell_active_predict.LR.bed.fall_in_region_c.cat
    cat *_STARR/*non_active_predict.LR.bed.fall_in_region_d|cut -f1,2,3 - > 8_cell_non_active_predict.LR.bed.fall_in_region_d.cat
    cat *_STARR/*.non_candidates.fall_in_region_f|cut -f1,2,3 - > 8_cell_non_candidates.fall_in_region_f.cat
    cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_e|cut -f1,2,3 - > 8_cell_STARR.consolidate.bed.sort.filter.fall_in_region_e.cat
    cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_d|cut -f1,2,3 - > 8_cell_STARR.consolidate.bed.sort.filter.fall_in_region_d.cat
    cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_c|cut -f1,2,3 - > 8_cell_STARR.consolidate.bed.sort.filter.fall_in_region_c.cat
}

function run_sort_merge_8_cell_abcedf_regions(){
    for f in 8_CELL*merge
    do
	cmd1="python extract_binding_sites_conservation_0_base_drop_nan.py ${f} hg38.phyloP100way.bw phylop100way"
	pbs="run_process_8_cell_merge.pbs"
	echo -e ${pbs_header} > ${pbs}
	echo -e ${cmd1} >> ${pbs}
	sbatch ${pbs}
    done
}


function run_get_8_cell_coverage(){
    python get_abcd_region_coverage.py
}

function run_starr_seq_number(){
    python run_get_crm_number_vs_starr_peak_number.py
}

function run_merge_gene(){
    folder=$1
    cell_name=${folder/"_STARR"/}
    cd Gene_Expression
    
    cp ../../merge_gene_expression_replicates.py .
    python merge_gene_expression_replicates.py ${cell_name}
    cd ../
}


function run_merge_abcdef_region(){
    for i in f #a b c d e f
    do
	bash run_get_intersect_abcdef_region.sh ${i}
    done
}

function run_plot_conservation_8_cell(){
    pbs="run_plot_abcde_8_cell.pbs"
    echo -e ${pbs_header} > ${pbs}
    echo -e "python plot_conservation_score_density_abcdef_8_cell.py" >> ${pbs}
    sbatch ${pbs}

}

function run_get_full_crm_in_starr_in_cd_regions(){
    folder=$1

    crm="ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_36.clean.sort"
    gene_tss_info="ensemble_info_2021_5_15.tss.gene_info_clean.bed.sort"
    dataset="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge"
    non_candidate="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge.non_candidate"

    cell_name=${folder/"_STARR"/}
    active_crm=${cell_name}.active_predict.LR.bed
    non_active_crm=${cell_name}.non_active_predict.LR.bed
    
    
    starr_peak_merge=${cell_name}"_STARR.consolidate.bed.sort.filter.merge"
    starr_peak=${cell_name}"_STARR.consolidate.bed.sort.filter"
    
    
    bedtools intersect -a ${active_crm}".fall_in_region_c.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_c.v0.5.sort.uniq" -wb > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_partial_crm.v0.5"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_d.v0.5.sort.uniq" -wb > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_partial_crm.v0.5"

    bedtools intersect -a ${active_crm}".fall_in_region_c.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_c.v0.5.sort.uniq" -f 1 -wb > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_full_crm.v0.5"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_d.v0.5.sort.uniq" -f 1 -wb > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_full_crm.v0.5"


    bedtools intersect -a ${active_crm}".fall_in_region_c.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.f1bp.sort.uniq" -wb > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_partial_crm.f1bp"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.f1bp.sort.uniq" -wb > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_partial_crm.f1bp"

    bedtools intersect -a ${active_crm}".fall_in_region_c.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.f1bp.sort.uniq" -f 1 -wb > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_full_crm.f1bp"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.f1bp.sort.uniq" -f 1 -wb > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_full_crm.f1bp"


    bedtools intersect -a ${active_crm}".fall_in_region_c.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.v1bp.sort.uniq" -wb > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_partial_crm.v1bp"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.v1bp.sort.uniq" -wb > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_partial_crm.v1bp"

    bedtools intersect -a ${active_crm}".fall_in_region_c.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.v1bp.sort.uniq" -f 1 -wb > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_full_crm.v1bp"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.v1bp.sort.uniq" -f 1 -wb > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_full_crm.v1bp"
    
}





function run_get_full_crm_in_starr_in_cd_regions_wao(){
    folder=$1

    crm="ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_36.clean.sort"
    gene_tss_info="ensemble_info_2021_5_15.tss.gene_info_clean.bed.sort"
    dataset="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge"
    non_candidate="TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge.non_candidate"

    cell_name=${folder/"_STARR"/}
    active_crm=${cell_name}.active_predict.LR.bed
    non_active_crm=${cell_name}.non_active_predict.LR.bed
    
    
    starr_peak_merge=${cell_name}"_STARR.consolidate.bed.sort.filter.merge"
    starr_peak=${cell_name}"_STARR.consolidate.bed.sort.filter"
    

    bedtools intersect -a ${active_crm} -b ${cell_name}".region.c" -f 1  -wa > ${active_crm}".fall_in_region_c"
    bedtools intersect -a ${non_active_crm} -b ${cell_name}".region.d" -f 1 -wa > ${non_active_crm}".fall_in_region_d"
    sort -k1,1 -k2,2n ${active_crm}".fall_in_region_c" > ${active_crm}".fall_in_region_c.sort"
    sort -k1,1 -k2,2n ${non_active_crm}".fall_in_region_d" > ${non_active_crm}".fall_in_region_d.sort"

    bedtools intersect -a ${active_crm}".fall_in_region_c.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_c.v0.5.sort.uniq" -wao > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_partial_crm.v0.5.wao"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_d.v0.5.sort.uniq" -wao > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_partial_crm.v0.5.wao"

    bedtools intersect -a ${active_crm}".fall_in_region_c.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_c.v0.5.sort.uniq" -F 1 -wao > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_full_crm.v0.5.wao"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v0.5.sort.uniq" -b ${starr_peak}".fall_in_region_d.v0.5.sort.uniq" -F 1 -wao > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_full_crm.v0.5.wao"


    bedtools intersect -a ${active_crm}".fall_in_region_c.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.f1bp.sort.uniq" -wao > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_partial_crm.f1bp.wao"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.f1bp.sort.uniq" -wao > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_partial_crm.f1bp.wao"

    bedtools intersect -a ${active_crm}".fall_in_region_c.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.f1bp.sort.uniq" -F 1 -wao > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_full_crm.f1bp.wao"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.f1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.f1bp.sort.uniq" -F 1 -wao > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_full_crm.f1bp.wao"


    bedtools intersect -a ${active_crm}".fall_in_region_c.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.v1bp.sort.uniq" -wao > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_partial_crm.v1bp.wao"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.v1bp.sort.uniq" -wao > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_partial_crm.v1bp.wao"

    bedtools intersect -a ${active_crm}".fall_in_region_c.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_c.v1bp.sort.uniq" -F 1 -wao > ${active_crm}-${starr_peak}".fall_in_region_c.overlap_full_crm.v1bp.wao"
    bedtools intersect -a ${non_active_crm}".fall_in_region_d.v1bp.sort.uniq" -b ${starr_peak}".fall_in_region_d.v1bp.sort.uniq" -F 1 -wao > ${non_active_crm}-${starr_peak}".fall_in_region_d.overlap_full_crm.v1bp.wao"
}


function run_extract_starr_level(){
    cell_name=${folder/"_STARR"/}
    cp ../extract_starr_level.py .
    cp ../run_extract_starr_level.sh .
    
    for f in *fall_in_region_{a,b,c,d,e,f} *fall_in_region_{a,b,c,d,e,f}.v1bp *fall_in_region_{a,b,c,d,e,f}.f1bp *fall_in_region_{a,b,c,d,e,f}.v0.5
    do
	sbatch --time=10:00:00 --mem=10GB run_extract_starr_level.sh ${cell_name} ${f}
    done

}


function run_extract_histone_level(){
    cell_name=${folder/"_STARR"/}
    cp ../extract_histone_level.py .
    cp ../run_extract_histone_level.sh .
    
    for f in *fall_in_region_{a,b,c,d,e,f} *fall_in_region_{a,b,c,d,e,f}.v1bp *fall_in_region_{a,b,c,d,e,f}.f1bp *fall_in_region_{a,b,c,d,e,f}.v0.5
    do
	for histone in *{H3,ATAC}*.bigWig
	do
	    sbatch --time=10:00:00 --mem=10GB run_extract_histone_level.sh ${cell_name} ${f} ${histone}
	done
    done

}


function run_extract_length(){
    cell_name=${folder/"_STARR"/}
    cp ../run_extract_len.sh .
    for f in *fall_in_region_{a,b,c,d,e,f}.v1bp *fall_in_region_{a,b,c,d,e,f}.f1bp *fall_in_region_{a,b,c,d,e,f}.v0.5
    do
	sbatch --time=10:00:00 --mem=10GB run_extract_len.sh ${f}
    done

}

function run_merge_crms_abcdef_regions_sort_uniq(){

    for ext in v0.5
    do
	cat *_STARR/*non_active_predict.LR.bed.fall_in_region_a.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_non_active_predict.LR.bed.fall_in_region_a.${ext}.sort.uniq
	cat *_STARR/*.active_predict.LR.bed.fall_in_region_b.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_active_predict.LR.bed.fall_in_region_b.${ext}.sort.uniq
	cat *_STARR/*active_predict.LR.bed.fall_in_region_c.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_active_predict.LR.bed.fall_in_region_c.${ext}.sort.uniq
	cat *_STARR/*non_active_predict.LR.bed.fall_in_region_d.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_non_active_predict.LR.bed.fall_in_region_d.${ext}.sort.uniq
	cat *_STARR/*.non_candidates.fall_in_region_f.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_non_candidates.fall_in_region_f.${ext}.sort.uniq
	cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_e.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_STARR.consolidate.bed.sort.filter.fall_in_region_e.${ext}.sort.uniq
	cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_d.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_STARR.consolidate.bed.sort.filter.fall_in_region_d.${ext}.sort.uniq
	cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_c.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell_STARR.consolidate.bed.sort.filter.fall_in_region_c.${ext}.sort.uniq
    done
}


function run_merge_crms_abcdef_regions_sort_uniq_new(){

    for ext in v0.5
    do
	cat *_STARR/*non_active_predict.LR.bed.fall_in_region_a.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.non_active_predict.LR.bed.fall_in_region_a.${ext}.sort.uniq
	cat *_STARR/*.active_predict.LR.bed.fall_in_region_b.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.active_predict.LR.bed.fall_in_region_b.${ext}.sort.uniq
	cat *_STARR/*active_predict.LR.bed.fall_in_region_c.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.active_predict.LR.bed.fall_in_region_c.${ext}.sort.uniq
	cat *_STARR/*non_active_predict.LR.bed.fall_in_region_d.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.non_active_predict.LR.bed.fall_in_region_d.${ext}.sort.uniq
	cat *_STARR/*.non_candidates.fall_in_region_f.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.non_candidates.fall_in_region_f.${ext}.sort.uniq
	cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_e.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.STARR.consolidate.bed.sort.filter.fall_in_region_e.${ext}.sort.uniq
	cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_d.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.STARR.consolidate.bed.sort.filter.fall_in_region_d.${ext}.sort.uniq
	cat *_STARR/*_STARR.consolidate.bed.sort.filter.fall_in_region_c.${ext}.sort.uniq|cut -f1,2,3 -|sort|uniq > 8_cell.STARR.consolidate.bed.sort.filter.fall_in_region_c.${ext}.sort.uniq
    done
}



function run_sort_uniq_fall_full_starr(){
    for f in *crm.{v1bp,v0.5}
    do
	sort ${f}|uniq > ${f}.sort.uniq
    done

}


folders=*_STARR
for folder in ${folders} 
do
    cd ${folder}
    #####################################
    #run_get_non_activate ${folder}
    #select_enhancers ${folder}
    #run_move_bak_crm ${folder}
    #run_process_data ${folder}
    #run_process_overlap ${folder}
    #run_split_starr_abcdef ${folder}
    #run_split_crms ${folder}
    #run_sample_a_region ${folder}
    #run_get_conservation ${folder}
    #run_reformat ${folder}
    #run_plot_abcdef ${folder} #need to run
    #run_count_fall_bp
    #run_get_starr_fall_in_region ${folder}
    #run_get_length
    #run_overlap_crm_starr_in_regions ${folder}
    #run_reformat_v0.5 ${folder}
    #run_reformat_starr ${folder}
    #run_count_bp_f ${folder}
    #run_generate_random_peaks
    #run_intersect_TE
    #run_intersect_chromHMM ${folder}
    #run_find_active_in_other_cells ${folder}
    #run_compute_gc
    #run_find_full_starr_in_other_cells ${folder}
    #run_get_tfbs
    #run_get_tfbs_homer #not run
    #run_create_tfbs_degree
    #run_plotenrichheatmap ${folder}
    #run_extract_length ${folder}

    #run_overlap_crm_starr_in_regions ${folder}
    #run_overlap_crm_starr_in_regions_1bp ${folder}
    #run_overlap_crm_starr_in_regions_f1bp ${folder}
    #run_sort_uniq_fall_in_region_abcdef 
    #run_merge_gene ${folder}
    #run_get_full_crm_in_starr_in_cd_regions ${folder}
    #run_sort_uniq_fall_full_starr
    #run_extract_starr_level ${folder}
    #run_extract_histone_level ${folder}
    #run_extract_length ${folder}
    #run_get_full_crm_in_starr_in_cd_regions_wao ${folder}
    cd ../
done

