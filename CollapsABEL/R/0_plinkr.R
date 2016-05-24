#' A wrapper for plink
#' 
#' @importFrom stringr str_trim
#' 
#' @param	D	Same as plink --D
#' @param	K	Same as plink --K
#' @param	a1_allele	Same as plink --a1-allele
#' @param	a2_allele	Same as plink --a2-allele
#' @param	adjust	Same as plink --adjust
#' @param	all	Same as plink --all
#' @param	all_pheno	Same as plink --all-pheno
#' @param	allele1234	Same as plink --allele1234
#' @param	alleleACGT	Same as plink --alleleACGT
#' @param	allele_count	Same as plink --allele-count
#' @param	allow_extra_chr	Same as plink --allow-extra-chr
#' @param	allow_no_sex	Same as plink --allow-no-sex
#' @param	alt_group	Same as plink --alt-group
#' @param	alt_snp	Same as plink --alt-snp
#' @param	annotate	Same as plink --annotate
#' @param	annotate_snp_field	Same as plink --annotate-snp-field
#' @param	aperm	Same as plink --aperm
#' @param	assoc	Same as plink --assoc
#' @param	attrib	Same as plink --attrib
#' @param	attrib_indiv	Same as plink --attrib-indiv
#' @param	autosome	Same as plink --autosome
#' @param	autosome_num	Same as plink --autosome-num
#' @param	autosome_xy	Same as plink --autosome-xy
#' @param	bcf	Same as plink --bcf
#' @param	bd	Same as plink --bd
#' @param	bed	Same as plink --bed
#' @param	beta	Same as plink --beta
#' @param	bfile	Same as plink --bfile
#' @param	bgen	Same as plink --bgen
#' @param	biallelic_only	Same as plink --biallelic-only
#' @param	bim	Same as plink --bim
#' @param	blocks	Same as plink --blocks
#' @param	blocks_inform_frac	Same as plink --blocks-inform-frac
#' @param	blocks_max_kb	Same as plink --blocks-max-kb
#' @param	blocks_min_maf	Same as plink --blocks-min-maf
#' @param	blocks_recomb_highci	Same as plink --blocks-recomb-highci
#' @param	blocks_strong_highci	Same as plink --blocks-strong-highci
#' @param	blocks_strong_lowci	Same as plink --blocks-strong-lowci
#' @param	bmerge	Same as plink --bmerge
#' @param	border	Same as plink --border
#' @param	bp_space	Same as plink --bp-space
#' @param	case_only	Same as plink --case-only
#' @param	cc	Same as plink --cc
#' @param	cell	Same as plink --cell
#' @param	cfile	Same as plink --cfile
#' @param	chap	Same as plink --chap
#' @param	check_sex	Same as plink --check-sex
#' @param	chr	Same as plink --chr
#' @param	chr_set	Same as plink --chr-set
#' @param	ci	Same as plink --ci
#' @param	clump	Same as plink --clump
#' @param	clump_allow_overlap	Same as plink --clump-allow-overlap
#' @param	clump_annotate	Same as plink --clump-annotate
#' @param	clump_best	Same as plink --clump-best
#' @param	clump_field	Same as plink --clump-field
#' @param	clump_index_first	Same as plink --clump-index-first
#' @param	clump_kb	Same as plink --clump-kb
#' @param	clump_p1	Same as plink --clump-p1
#' @param	clump_p2	Same as plink --clump-p2
#' @param	clump_r2	Same as plink --clump-r2
#' @param	clump_range	Same as plink --clump-range
#' @param	clump_range_border	Same as plink --clump-range-border
#' @param	clump_replicate	Same as plink --clump-replicate
#' @param	clump_snp_field	Same as plink --clump-snp-field
#' @param	clump_verbose	Same as plink --clump-verbose
#' @param	cluster	Same as plink --cluster
#' @param	cluster_missing	Same as plink --cluster-missing
#' @param	cm_map	Same as plink --cm-map
#' @param	cnv_blue	Same as plink --cnv-blue
#' @param	cnv_border	Same as plink --cnv-border
#' @param	cnv_brown	Same as plink --cnv-brown
#' @param	cnv_check_no_overlap	Same as plink --cnv-check-no-overlap
#' @param	cnv_count	Same as plink --cnv-count
#' @param	cnv_del	Same as plink --cnv-del
#' @param	cnv_disrupt	Same as plink --cnv-disrupt
#' @param	cnv_drop_no_segment	Same as plink --cnv-drop-no-segment
#' @param	cnv_dup	Same as plink --cnv-dup
#' @param	cnv_enrichment_test	Same as plink --cnv-enrichment-test
#' @param	cnv_exclude	Same as plink --cnv-exclude
#' @param	cnv_exclude_off_by_1	Same as plink --cnv-exclude-off-by-1
#' @param	cnv_freq_excldue_above	Same as plink --cnv-freq-excldue-above
#' @param	cnv_freq_excldue_below	Same as plink --cnv-freq-excldue-below
#' @param	cnv_freq_excldue_exact	Same as plink --cnv-freq-excldue-exact
#' @param	cnv_freq_exclude_above	Same as plink --cnv-freq-exclude-above
#' @param	cnv_freq_exclude_below	Same as plink --cnv-freq-exclude-below
#' @param	cnv_freq_exclude_exact	Same as plink --cnv-freq-exclude-exact
#' @param	cnv_freq_incldue_exact	Same as plink --cnv-freq-incldue-exact
#' @param	cnv_freq_include_exact	Same as plink --cnv-freq-include-exact
#' @param	cnv_freq_method2	Same as plink --cnv-freq-method2
#' @param	cnv_freq_overlap	Same as plink --cnv-freq-overlap
#' @param	cnv_green	Same as plink --cnv-green
#' @param	cnv_indiv_perm	Same as plink --cnv-indiv-perm
#' @param	cnv_intersect	Same as plink --cnv-intersect
#' @param	cnv_kb	Same as plink --cnv-kb
#' @param	cnv_list	Same as plink --cnv-list
#' @param	cnv_make_map	Same as plink --cnv-make-map
#' @param	cnv_max_kb	Same as plink --cnv-max-kb
#' @param	cnv_max_score	Same as plink --cnv-max-score
#' @param	cnv_max_sites	Same as plink --cnv-max-sites
#' @param	cnv_overlap	Same as plink --cnv-overlap
#' @param	cnv_red	Same as plink --cnv-red
#' @param	cnv_region_overlap	Same as plink --cnv-region-overlap
#' @param	cnv_report_regions	Same as plink --cnv-report-regions
#' @param	cnv_score	Same as plink --cnv-score
#' @param	cnv_seglist	Same as plink --cnv-seglist
#' @param	cnv_sites	Same as plink --cnv-sites
#' @param	cnv_subset	Same as plink --cnv-subset
#' @param	cnv_test	Same as plink --cnv-test
#' @param	cnv_test_1sided	Same as plink --cnv-test-1sided
#' @param	cnv_test_2sided	Same as plink --cnv-test-2sided
#' @param	cnv_test_region	Same as plink --cnv-test-region
#' @param	cnv_test_window	Same as plink --cnv-test-window
#' @param	cnv_track	Same as plink --cnv-track
#' @param	cnv_union_overlap	Same as plink --cnv-union-overlap
#' @param	cnv_unique	Same as plink --cnv-unique
#' @param	cnv_verbose_report_regions	Same as plink --cnv-verbose-report-regions
#' @param	cnv_write	Same as plink --cnv-write
#' @param	cnv_write_freq	Same as plink --cnv-write-freq
#' @param	complement_sets	Same as plink --complement-sets
#' @param	compound_genotypes	Same as plink --compound-genotypes
#' @param	compress	Same as plink --compress
#' @param	condition	Same as plink --condition
#' @param	condition_list	Same as plink --condition-list
#' @param	consensus_match	Same as plink --consensus-match
#' @param	const_fid	Same as plink --const-fid
#' @param	control	Same as plink --control
#' @param	counts	Same as plink --counts
#' @param	covar	Same as plink --covar
#' @param	covar_name	Same as plink --covar-name
#' @param	covar_number	Same as plink --covar-number
#' @param	cow	Same as plink --cow
#' @param	d	Same as plink --d
#' @param	data	Same as plink --data
#' @param	debug	Same as plink --debug
#' @param	decompress	Same as plink --decompress
#' @param	dfam	Same as plink --dfam
#' @param	distance	Same as plink --distance
#' @param	distance_exp	Same as plink --distance-exp
#' @param	distance_matrix	Same as plink --distance-matrix
#' @param	dog	Same as plink --dog
#' @param	dominant	Same as plink --dominant
#' @param	dosage	Same as plink --dosage
#' @param	double_id	Same as plink --double-id
#' @param	dprime	Same as plink --dprime
#' @param	dummy	Same as plink --dummy
#' @param	dummy_coding	Same as plink --dummy-coding
#' @param	each_versus_others	Same as plink --each-versus-others
#' @param	each_vs_others	Same as plink --each-vs-others
#' @param	epistasis	Same as plink --epistasis
#' @param	epistasis_summary_merge	Same as plink --epistasis-summary-merge
#' @param	exclude	Same as plink --exclude
#' @param	exclude_before_extract	Same as plink --exclude-before-extract
#' @param	exclude_snp	Same as plink --exclude-snp
#' @param	exclude_snps	Same as plink --exclude-snps
#' @param	extract	Same as plink --extract
#' @param	fam	Same as plink --fam
#' @param	family	Same as plink --family
#' @param	fast_epistasis	Same as plink --fast-epistasis
#' @param	fid	Same as plink --fid
#' @param	file	Same as plink --file
#' @param	fill_missing_a2	Same as plink --fill-missing-a2
#' @param	filter	Same as plink --filter
#' @param	filter_cases	Same as plink --filter-cases
#' @param	filter_controls	Same as plink --filter-controls
#' @param	filter_females	Same as plink --filter-females
#' @param	filter_founders	Same as plink --filter-founders
#' @param	filter_males	Same as plink --filter-males
#' @param	filter_nonfounders	Same as plink --filter-nonfounders
#' @param	fisher	Same as plink --fisher
#' @param	flip	Same as plink --flip
#' @param	flip_scan	Same as plink --flip-scan
#' @param	flip_scan_threshold	Same as plink --flip-scan-threshold
#' @param	flip_scan_verbose	Same as plink --flip-scan-verbose
#' @param	flip_scan_window	Same as plink --flip-scan-window
#' @param	flip_scan_window_kb	Same as plink --flip-scan-window-kb
#' @param	flip_subset	Same as plink --flip-subset
#' @param	freq	Same as plink --freq
#' @param	freqx	Same as plink --freqx
#' @param	from	Same as plink --from
#' @param	from_bp	Same as plink --from-bp
#' @param	from_kb	Same as plink --from-kb
#' @param	from_mb	Same as plink --from-mb
#' @param	frqx	Same as plink --frqx
#' @param	fst	Same as plink --fst
#' @param	gap	Same as plink --gap
#' @param	gates	Same as plink --gates
#' @param	gc	Same as plink --gc
#' @param	gen	Same as plink --gen
#' @param	gene	Same as plink --gene
#' @param	gene_all	Same as plink --gene-all
#' @param	gene_list	Same as plink --gene-list
#' @param	gene_list_border	Same as plink --gene-list-border
#' @param	gene_report	Same as plink --gene-report
#' @param	gene_report_empty	Same as plink --gene-report-empty
#' @param	gene_report_snp_field	Same as plink --gene-report-snp-field
#' @param	gene_subset	Same as plink --gene-subset
#' @param	genedrop	Same as plink --genedrop
#' @param	genepi	Same as plink --genepi
#' @param	geno	Same as plink --geno
#' @param	genome	Same as plink --genome
#' @param	genome_full	Same as plink --genome-full
#' @param	genome_lists	Same as plink --genome-lists
#' @param	genome_minimal	Same as plink --genome-minimal
#' @param	genotypic	Same as plink --genotypic
#' @param	gfile	Same as plink --gfile
#' @param	gplink	Same as plink --gplink
#' @param	grm	Same as plink --grm
#' @param	grm_bin	Same as plink --grm-bin
#' @param	grm_gz	Same as plink --grm-gz
#' @param	group_avg	Same as plink --group-avg
#' @param	groupdist	Same as plink --groupdist
#' @param	gxe	Same as plink --gxe
#' @param	hap...	Same as plink --hap...
#' @param	hap	Same as plink --hap
#' @param	hap_assoc	Same as plink --hap-assoc
#' @param	hap_freq	Same as plink --hap-freq
#' @param	hap_impute	Same as plink --hap-impute
#' @param	hap_max_phase	Same as plink --hap-max-phase
#' @param	hap_min_phase_prob	Same as plink --hap-min-phase-prob
#' @param	hap_miss	Same as plink --hap-miss
#' @param	hap_phase	Same as plink --hap-phase
#' @param	hap_phase_wide	Same as plink --hap-phase-wide
#' @param	hap_pp	Same as plink --hap-pp
#' @param	hap_snps	Same as plink --hap-snps
#' @param	hap_tdt	Same as plink --hap-tdt
#' @param	hap_window	Same as plink --hap-window
#' @param	hard_call_threshold	Same as plink --hard-call-threshold
#' @param	hardy2	Same as plink --hardy2
#' @param	hardy	Same as plink --hardy
#' @param	help	Same as plink --help
#' @param	het	Same as plink --het
#' @param	hethom	Same as plink --hethom
#' @param	hide_covar	Same as plink --hide-covar
#' @param	homog	Same as plink --homog
#' @param	homozyg	Same as plink --homozyg
#' @param	homozyg_density	Same as plink --homozyg-density
#' @param	homozyg_gap	Same as plink --homozyg-gap
#' @param	homozyg_group	Same as plink --homozyg-group
#' @param	homozyg_het	Same as plink --homozyg-het
#' @param	homozyg_include_missing	Same as plink --homozyg-include-missing
#' @param	homozyg_kb	Same as plink --homozyg-kb
#' @param	homozyg_match	Same as plink --homozyg-match
#' @param	homozyg_snp	Same as plink --homozyg-snp
#' @param	homozyg_verbose	Same as plink --homozyg-verbose
#' @param	homozyg_window_het	Same as plink --homozyg-window-het
#' @param	homozyg_window_kb	Same as plink --homozyg-window-kb
#' @param	homozyg_window_missing	Same as plink --homozyg-window-missing
#' @param	homozyg_window_snp	Same as plink --homozyg-window-snp
#' @param	homozyg_window_threshold	Same as plink --homozyg-window-threshold
#' @param	horse	Same as plink --horse
#' @param	hwe	Same as plink --hwe
#' @param	hwe_all	Same as plink --hwe-all
#' @param	ibc	Same as plink --ibc
#' @param	ibm	Same as plink --ibm
#' @param	ibs_matrix	Same as plink --ibs-matrix
#' @param	ibs_test	Same as plink --ibs-test
#' @param	id_delim	Same as plink --id-delim
#' @param	id_dict	Same as plink --id-dict
#' @param	id_match	Same as plink --id-match
#' @param	iid	Same as plink --iid
#' @param	impossible	Same as plink --impossible
#' @param	impute_sex	Same as plink --impute-sex
#' @param	ind_major	Same as plink --ind-major
#' @param	indep	Same as plink --indep
#' @param	indep_pairphase	Same as plink --indep-pairphase
#' @param	indep_pairwise	Same as plink --indep-pairwise
#' @param	independent_effect	Same as plink --independent-effect
#' @param	indiv_sort	Same as plink --indiv-sort
#' @param	inter_chr	Same as plink --inter-chr
#' @param	interaction	Same as plink --interaction
#' @param	je_cellmin	Same as plink --je-cellmin
#' @param	keep	Same as plink --keep
#' @param	keep_allele_order	Same as plink --keep-allele-order
#' @param	keep_autoconv	Same as plink --keep-autoconv
#' @param	keep_before_remove	Same as plink --keep-before-remove
#' @param	keep_cluster_names	Same as plink --keep-cluster-names
#' @param	keep_clusters	Same as plink --keep-clusters
#' @param	keep_fam	Same as plink --keep-fam
#' @param	lambda	Same as plink --lambda
#' @param	lasso	Same as plink --lasso
#' @param	lasso_select_covars	Same as plink --lasso-select-covars
#' @param	ld	Same as plink --ld
#' @param	ld_snp	Same as plink --ld-snp
#' @param	ld_snp_list	Same as plink --ld-snp-list
#' @param	ld_snps	Same as plink --ld-snps
#' @param	ld_window	Same as plink --ld-window
#' @param	ld_window_kb	Same as plink --ld-window-kb
#' @param	ld_window_r2	Same as plink --ld-window-r2
#' @param	ld_xchr	Same as plink --ld-xchr
#' @param	lfile	Same as plink --lfile
#' @param	liability	Same as plink --liability
#' @param	linear	Same as plink --linear
#' @param	list	Same as plink --list
#' @param	list_23_indels	Same as plink --list-23-indels
#' @param	list_all	Same as plink --list-all
#' @param	logistic	Same as plink --logistic
#' @param	lookup...	Same as plink --lookup...
#' @param	lookup	Same as plink --lookup
#' @param	lookup_gene	Same as plink --lookup-gene
#' @param	lookup_list	Same as plink --lookup-list
#' @param	loop_assoc	Same as plink --loop-assoc
#' @param	maf	Same as plink --maf
#' @param	maf_succ	Same as plink --maf-succ
#' @param	make_bed	Same as plink --make-bed
#' @param	make_founders	Same as plink --make-founders
#' @param	make_grm	Same as plink --make-grm
#' @param	make_grm_bin	Same as plink --make-grm-bin
#' @param	make_grm_gz	Same as plink --make-grm-gz
#' @param	make_just_bim	Same as plink --make-just-bim
#' @param	make_just_fam	Same as plink --make-just-fam
#' @param	make_perm_pheno	Same as plink --make-perm-pheno
#' @param	make_pheno	Same as plink --make-pheno
#' @param	make_rel	Same as plink --make-rel
#' @param	make_set	Same as plink --make-set
#' @param	make_set_border	Same as plink --make-set-border
#' @param	make_set_collapse_group	Same as plink --make-set-collapse-group
#' @param	make_set_complement_all	Same as plink --make-set-complement-all
#' @param	make_set_complement_group	Same as plink --make-set-complement-group
#' @param	map	Same as plink --map
#' @param	mat	Same as plink --mat
#' @param	match	Same as plink --match
#' @param	match_type	Same as plink --match-type
#' @param	matrix	Same as plink --matrix
#' @param	max	Same as plink --max
#' @param	max_maf	Same as plink --max-maf
#' @param	mc	Same as plink --mc
#' @param	mcc	Same as plink --mcc
#' @param	mcovar	Same as plink --mcovar
#' @param	mds_cluster	Same as plink --mds-cluster
#' @param	mds_plot	Same as plink --mds-plot
#' @param	me	Same as plink --me
#' @param	me_exclude_one	Same as plink --me-exclude-one
#' @param	memory	Same as plink --memory
#' @param	mendel	Same as plink --mendel
#' @param	mendel_duos	Same as plink --mendel-duos
#' @param	mendel_multigen	Same as plink --mendel-multigen
#' @param	merge	Same as plink --merge
#' @param	merge_equal_pos	Same as plink --merge-equal-pos
#' @param	merge_list	Same as plink --merge-list
#' @param	merge_mode	Same as plink --merge-mode
#' @param	merge_x	Same as plink --merge-x
#' @param	meta_analysis	Same as plink --meta-analysis
#' @param	meta_analysis_..._field	Same as plink --meta-analysis-...-field
#' @param	mfilter	Same as plink --mfilter
#' @param	mh	Same as plink --mh
#' @param	mhf	Same as plink --mhf
#' @param	min	Same as plink --min
#' @param	mind	Same as plink --mind
#' @param	mishap_window	Same as plink --mishap-window
#' @param	missing	Same as plink --missing
#' @param	missing_code	Same as plink --missing-code
#' @param	missing_genotype	Same as plink --missing-genotype
#' @param	missing_phenotype	Same as plink --missing-phenotype
#' @param	missing_var_code	Same as plink --missing-var-code
#' @param	mlma	Same as plink --mlma
#' @param	mlma_loco	Same as plink --mlma-loco
#' @param	mlma_no_adj_covar	Same as plink --mlma-no-adj-covar
#' @param	model	Same as plink --model
#' @param	model_dom	Same as plink --model-dom
#' @param	model_gen	Same as plink --model-gen
#' @param	model_rec	Same as plink --model-rec
#' @param	model_trend	Same as plink --model-trend
#' @param	mouse	Same as plink --mouse
#' @param	mperm	Same as plink --mperm
#' @param	mperm_save	Same as plink --mperm-save
#' @param	mperm_save_all	Same as plink --mperm-save-all
#' @param	mpheno	Same as plink --mpheno
#' @param	must_have_sex	Same as plink --must-have-sex
#' @param	mwithin	Same as plink --mwithin
#' @param	neighbour	Same as plink --neighbour
#' @param	no_fid	Same as plink --no-fid
#' @param	no_parents	Same as plink --no-parents
#' @param	no_pheno	Same as plink --no-pheno
#' @param	no_sex	Same as plink --no-sex
#' @param	no_snp	Same as plink --no-snp
#' @param	no_x_sex	Same as plink --no-x-sex
#' @param	nonfounders	Same as plink --nonfounders
#' @param	nop	Same as plink --nop
#' @param	not_chr	Same as plink --not-chr
#' @param	nudge	Same as plink --nudge
#' @param	null_group	Same as plink --null-group
#' @param	null_snp	Same as plink --null-snp
#' @param	oblig_cluster	Same as plink --oblig-cluster
#' @param	oblig_clusters	Same as plink --oblig-clusters
#' @param	oblig_missing	Same as plink --oblig-missing
#' @param	out	Same as plink --out
#' @param	output_chr	Same as plink --output-chr
#' @param	output_missing_genotype	Same as plink --output-missing-genotype
#' @param	output_missing_phenotype	Same as plink --output-missing-phenotype
#' @param	oxford_pheno_name	Same as plink --oxford-pheno-name
#' @param	parallel	Same as plink --parallel
#' @param	parameters	Same as plink --parameters
#' @param	parentdt1	Same as plink --parentdt1
#' @param	parentdt2	Same as plink --parentdt2
#' @param	pat	Same as plink --pat
#' @param	pca	Same as plink --pca
#' @param	pca_cluster_names	Same as plink --pca-cluster-names
#' @param	pca_clusters	Same as plink --pca-clusters
#' @param	ped	Same as plink --ped
#' @param	pedigree	Same as plink --pedigree
#' @param	perm	Same as plink --perm
#' @param	perm_batch_size	Same as plink --perm-batch-size
#' @param	perm_count	Same as plink --perm-count
#' @param	pfilter	Same as plink --pfilter
#' @param	pheno	Same as plink --pheno
#' @param	pheno_merge	Same as plink --pheno-merge
#' @param	pheno_name	Same as plink --pheno-name
#' @param	pick1	Same as plink --pick1
#' @param	plist	Same as plink --plist
#' @param	poo	Same as plink --poo
#' @param	pool_size	Same as plink --pool-size
#' @param	ppc	Same as plink --ppc
#' @param	ppc_gap	Same as plink --ppc-gap
#' @param	proxy_...	Same as plink --proxy-...
#' @param	proxy_assoc	Same as plink --proxy-assoc
#' @param	proxy_b_kb	Same as plink --proxy-b-kb
#' @param	proxy_b_maxsnp	Same as plink --proxy-b-maxsnp
#' @param	proxy_b_r2	Same as plink --proxy-b-r2
#' @param	proxy_b_threshold	Same as plink --proxy-b-threshold
#' @param	proxy_b_window	Same as plink --proxy-b-window
#' @param	proxy_dosage	Same as plink --proxy-dosage
#' @param	proxy_drop	Same as plink --proxy-drop
#' @param	proxy_flanking	Same as plink --proxy-flanking
#' @param	proxy_geno	Same as plink --proxy-geno
#' @param	proxy_genotypic_concordance	Same as plink --proxy-genotypic-concordance
#' @param	proxy_glm	Same as plink --proxy-glm
#' @param	proxy_impute	Same as plink --proxy-impute
#' @param	proxy_impute_threshold	Same as plink --proxy-impute-threshold
#' @param	proxy_kb	Same as plink --proxy-kb
#' @param	proxy_list	Same as plink --proxy-list
#' @param	proxy_maf	Same as plink --proxy-maf
#' @param	proxy_maxsnp	Same as plink --proxy-maxsnp
#' @param	proxy_mhf	Same as plink --proxy-mhf
#' @param	proxy_r2	Same as plink --proxy-r2
#' @param	proxy_r2_no_filter	Same as plink --proxy-r2-no-filter
#' @param	proxy_replace	Same as plink --proxy-replace
#' @param	proxy_show_proxies	Same as plink --proxy-show-proxies
#' @param	proxy_sub_maxsnp	Same as plink --proxy-sub-maxsnp
#' @param	proxy_sub_r2	Same as plink --proxy-sub-r2
#' @param	proxy_tdt	Same as plink --proxy-tdt
#' @param	proxy_verbose	Same as plink --proxy-verbose
#' @param	proxy_window	Same as plink --proxy-window
#' @param	prune	Same as plink --prune
#' @param	q_score_file	Same as plink --q-score-file
#' @param	q_score_range	Same as plink --q-score-range
#' @param	qfam...	Same as plink --qfam...
#' @param	qmatch	Same as plink --qmatch
#' @param	qq_plot	Same as plink --qq-plot
#' @param	qt	Same as plink --qt
#' @param	qt_means	Same as plink --qt-means
#' @param	qual_geno_...	Same as plink --qual-geno-...
#' @param	qual_geno_max_threshold	Same as plink --qual-geno-max-threshold
#' @param	qual_geno_scores	Same as plink --qual-geno-scores
#' @param	qual_geno_threshold	Same as plink --qual-geno-threshold
#' @param	qual_max_threshold	Same as plink --qual-max-threshold
#' @param	qual_scores	Same as plink --qual-scores
#' @param	qual_threshold	Same as plink --qual-threshold
#' @param	r2	Same as plink --r2
#' @param	r	Same as plink --r
#' @param	range	Same as plink --range
#' @param	rank	Same as plink --rank
#' @param	read_dists	Same as plink --read-dists
#' @param	read_freq	Same as plink --read-freq
#' @param	read_genome	Same as plink --read-genome
#' @param	read_genome_list	Same as plink --read-genome-list
#' @param	read_genome_minimal	Same as plink --read-genome-minimal
#' @param	recessive	Same as plink --recessive
#' @param	recode12	Same as plink --recode12
#' @param	recode	Same as plink --recode
#' @param	recodeA	Same as plink --recodeA
#' @param	recodeAD	Same as plink --recodeAD
#' @param	recodeHV	Same as plink --recodeHV
#' @param	recode_allele	Same as plink --recode-allele
#' @param	recode_beagle	Same as plink --recode-beagle
#' @param	recode_bimbam	Same as plink --recode-bimbam
#' @param	recode_fastphase	Same as plink --recode-fastphase
#' @param	recode_lgen	Same as plink --recode-lgen
#' @param	recode_rlist	Same as plink --recode-rlist
#' @param	recode_structure	Same as plink --recode-structure
#' @param	recode_vcf	Same as plink --recode-vcf
#' @param	recode_whap	Same as plink --recode-whap
#' @param	reference	Same as plink --reference
#' @param	reference_allele	Same as plink --reference-allele
#' @param	regress_distance	Same as plink --regress-distance
#' @param	regress_pcs	Same as plink --regress-pcs
#' @param	regress_rel	Same as plink --regress-rel
#' @param	rel_check	Same as plink --rel-check
#' @param	rel_cutoff	Same as plink --rel-cutoff
#' @param	remove	Same as plink --remove
#' @param	remove_cluster_names	Same as plink --remove-cluster-names
#' @param	remove_clusters	Same as plink --remove-clusters
#' @param	remove_fam	Same as plink --remove-fam
#' @param	rerun	Same as plink --rerun
#' @param	rice	Same as plink --rice
#' @param	sample	Same as plink --sample
#' @param	score	Same as plink --score
#' @param	score_no_mean_imputation	Same as plink --score-no-mean-imputation
#' @param	script	Same as plink --script
#' @param	seed	Same as plink --seed
#' @param	set	Same as plink --set
#' @param	set_by_all	Same as plink --set-by-all
#' @param	set_collapse_all	Same as plink --set-collapse-all
#' @param	set_hh_missing	Same as plink --set-hh-missing
#' @param	set_max	Same as plink --set-max
#' @param	set_me_missing	Same as plink --set-me-missing
#' @param	set_missing_nonsnp_ids	Same as plink --set-missing-nonsnp-ids
#' @param	set_missing_snp_ids	Same as plink --set-missing-snp-ids
#' @param	set_missing_var_ids	Same as plink --set-missing-var-ids
#' @param	set_names	Same as plink --set-names
#' @param	set_p	Same as plink --set-p
#' @param	set_r2	Same as plink --set-r2
#' @param	set_r2_phase	Same as plink --set-r2-phase
#' @param	set_table	Same as plink --set-table
#' @param	set_test	Same as plink --set-test
#' @param	sex	Same as plink --sex
#' @param	sheep	Same as plink --sheep
#' @param	show_tags	Same as plink --show-tags
#' @param	silent	Same as plink --silent
#' @param	simulate	Same as plink --simulate
#' @param	simulate_haps	Same as plink --simulate-haps
#' @param	simulate_label	Same as plink --simulate-label
#' @param	simulate_missing	Same as plink --simulate-missing
#' @param	simulate_n	Same as plink --simulate-n
#' @param	simulate_ncases	Same as plink --simulate-ncases
#' @param	simulate_ncontrols	Same as plink --simulate-ncontrols
#' @param	simulate_prevalence	Same as plink --simulate-prevalence
#' @param	simulate_qt	Same as plink --simulate-qt
#' @param	simulate_tags	Same as plink --simulate-tags
#' @param	snp	Same as plink --snp
#' @param	snps	Same as plink --snps
#' @param	snps_only	Same as plink --snps-only
#' @param	specific_haplotype	Same as plink --specific-haplotype
#' @param	split_x	Same as plink --split-x
#' @param	standard_beta	Same as plink --standard-beta
#' @param	subset	Same as plink --subset
#' @param	swap_parents	Same as plink --swap-parents
#' @param	swap_sibs	Same as plink --swap-sibs
#' @param	swap_unrel	Same as plink --swap-unrel
#' @param	tab	Same as plink --tab
#' @param	tag_kb	Same as plink --tag-kb
#' @param	tag_mode2	Same as plink --tag-mode2
#' @param	tag_r2	Same as plink --tag-r2
#' @param	tail_pheno	Same as plink --tail-pheno
#' @param	tdt	Same as plink --tdt
#' @param	test_all	Same as plink --test-all
#' @param	test_mishap	Same as plink --test-mishap
#' @param	test_missing	Same as plink --test-missing
#' @param	test_snp	Same as plink --test-snp
#' @param	tests	Same as plink --tests
#' @param	tfam	Same as plink --tfam
#' @param	tfile	Same as plink --tfile
#' @param	thin	Same as plink --thin
#' @param	thin_count	Same as plink --thin-count
#' @param	threads	Same as plink --threads
#' @param	to	Same as plink --to
#' @param	to_bp	Same as plink --to-bp
#' @param	to_kb	Same as plink --to-kb
#' @param	to_mb	Same as plink --to-mb
#' @param	tped	Same as plink --tped
#' @param	transpose	Same as plink --transpose
#' @param	trend	Same as plink --trend
#' @param	tucc	Same as plink --tucc
#' @param	twolocus	Same as plink --twolocus
#' @param	unbounded	Same as plink --unbounded
#' @param	unrelated_heritability	Same as plink --unrelated-heritability
#' @param	update_alleles	Same as plink --update-alleles
#' @param	update_chr	Same as plink --update-chr
#' @param	update_cm	Same as plink --update-cm
#' @param	update_ids	Same as plink --update-ids
#' @param	update_map	Same as plink --update-map
#' @param	update_name	Same as plink --update-name
#' @param	update_parents	Same as plink --update-parents
#' @param	update_sex	Same as plink --update-sex
#' @param	vcf	Same as plink --vcf
#' @param	vcf_filter	Same as plink --vcf-filter
#' @param	vcf_half_call	Same as plink --vcf-half-call
#' @param	vcf_idspace_to	Same as plink --vcf-idspace-to
#' @param	vcf_min_qual	Same as plink --vcf-min-qual
#' @param	vegas	Same as plink --vegas
#' @param	version	Same as plink --version
#' @param	vif	Same as plink --vif
#' @param	whap	Same as plink --whap
#' @param	window	Same as plink --window
#' @param	with_freqs	Same as plink --with-freqs
#' @param	with_phenotype	Same as plink --with-phenotype
#' @param	with_reference	Same as plink --with-reference
#' @param	within	Same as plink --within
#' @param	write_cluster	Same as plink --write-cluster
#' @param	write_covar	Same as plink --write-covar
#' @param	write_dosage	Same as plink --write-dosage
#' @param	write_set	Same as plink --write-set
#' @param	write_set_r2	Same as plink --write-set-r2
#' @param	write_snplist	Same as plink --write-snplist
#' @param	xchr_model	Same as plink --xchr-model
#' @param	zero_cluster	Same as plink --zero-cluster
#' @param	zero_cms	Same as plink --zero-cms
#' @param	one	Same as plink --1
#' @param	twothreefile	Same as plink --23file
#' @param	wait	Logical. If FALSE, the plink process will fork into the background.
#' @param	stdout	Passed to system2, see its documentation.
#' @param	stderr	Passed to system2, see its documentation.
#' @export
plinkr = function(
	D=NULL,
	K=NULL,
	a1_allele=NULL,
	a2_allele=NULL,
	adjust=NULL,
	all=NULL,
	all_pheno=NULL,
	allele1234=NULL,
	alleleACGT=NULL,
	allele_count=NULL,
	allow_extra_chr=NULL,
	allow_no_sex=NULL,
	alt_group=NULL,
	alt_snp=NULL,
	annotate=NULL,
	annotate_snp_field=NULL,
	aperm=NULL,
	assoc=NULL,
	attrib=NULL,
	attrib_indiv=NULL,
	autosome=NULL,
	autosome_num=NULL,
	autosome_xy=NULL,
	bcf=NULL,
	bd=NULL,
	bed=NULL,
	beta=NULL,
	bfile=NULL,
	bgen=NULL,
	biallelic_only=NULL,
	bim=NULL,
	blocks=NULL,
	blocks_inform_frac=NULL,
	blocks_max_kb=NULL,
	blocks_min_maf=NULL,
	blocks_recomb_highci=NULL,
	blocks_strong_highci=NULL,
	blocks_strong_lowci=NULL,
	bmerge=NULL,
	border=NULL,
	bp_space=NULL,
	case_only=NULL,
	cc=NULL,
	cell=NULL,
	cfile=NULL,
	chap=NULL,
	check_sex=NULL,
	chr=NULL,
	chr_set=NULL,
	ci=NULL,
	clump=NULL,
	clump_allow_overlap=NULL,
	clump_annotate=NULL,
	clump_best=NULL,
	clump_field=NULL,
	clump_index_first=NULL,
	clump_kb=NULL,
	clump_p1=NULL,
	clump_p2=NULL,
	clump_r2=NULL,
	clump_range=NULL,
	clump_range_border=NULL,
	clump_replicate=NULL,
	clump_snp_field=NULL,
	clump_verbose=NULL,
	cluster=NULL,
	cluster_missing=NULL,
	cm_map=NULL,
	cnv_blue=NULL,
	cnv_border=NULL,
	cnv_brown=NULL,
	cnv_check_no_overlap=NULL,
	cnv_count=NULL,
	cnv_del=NULL,
	cnv_disrupt=NULL,
	cnv_drop_no_segment=NULL,
	cnv_dup=NULL,
	cnv_enrichment_test=NULL,
	cnv_exclude=NULL,
	cnv_exclude_off_by_1=NULL,
	cnv_freq_excldue_above=NULL,
	cnv_freq_excldue_below=NULL,
	cnv_freq_excldue_exact=NULL,
	cnv_freq_exclude_above=NULL,
	cnv_freq_exclude_below=NULL,
	cnv_freq_exclude_exact=NULL,
	cnv_freq_incldue_exact=NULL,
	cnv_freq_include_exact=NULL,
	cnv_freq_method2=NULL,
	cnv_freq_overlap=NULL,
	cnv_green=NULL,
	cnv_indiv_perm=NULL,
	cnv_intersect=NULL,
	cnv_kb=NULL,
	cnv_list=NULL,
	cnv_make_map=NULL,
	cnv_max_kb=NULL,
	cnv_max_score=NULL,
	cnv_max_sites=NULL,
	cnv_overlap=NULL,
	cnv_red=NULL,
	cnv_region_overlap=NULL,
	cnv_report_regions=NULL,
	cnv_score=NULL,
	cnv_seglist=NULL,
	cnv_sites=NULL,
	cnv_subset=NULL,
	cnv_test=NULL,
	cnv_test_1sided=NULL,
	cnv_test_2sided=NULL,
	cnv_test_region=NULL,
	cnv_test_window=NULL,
	cnv_track=NULL,
	cnv_union_overlap=NULL,
	cnv_unique=NULL,
	cnv_verbose_report_regions=NULL,
	cnv_write=NULL,
	cnv_write_freq=NULL,
	complement_sets=NULL,
	compound_genotypes=NULL,
	compress=NULL,
	condition=NULL,
	condition_list=NULL,
	consensus_match=NULL,
	const_fid=NULL,
	control=NULL,
	counts=NULL,
	covar=NULL,
	covar_name=NULL,
	covar_number=NULL,
	cow=NULL,
	d=NULL,
	data=NULL,
	debug=NULL,
	decompress=NULL,
	dfam=NULL,
	distance=NULL,
	distance_exp=NULL,
	distance_matrix=NULL,
	dog=NULL,
	dominant=NULL,
	dosage=NULL,
	double_id=NULL,
	dprime=NULL,
	dummy=NULL,
	dummy_coding=NULL,
	each_versus_others=NULL,
	each_vs_others=NULL,
	epistasis=NULL,
	epistasis_summary_merge=NULL,
	exclude=NULL,
	exclude_before_extract=NULL,
	exclude_snp=NULL,
	exclude_snps=NULL,
	extract=NULL,
	fam=NULL,
	family=NULL,
	fast_epistasis=NULL,
	fid=NULL,
	file=NULL,
	fill_missing_a2=NULL,
	filter=NULL,
	filter_cases=NULL,
	filter_controls=NULL,
	filter_females=NULL,
	filter_founders=NULL,
	filter_males=NULL,
	filter_nonfounders=NULL,
	fisher=NULL,
	flip=NULL,
	flip_scan=NULL,
	flip_scan_threshold=NULL,
	flip_scan_verbose=NULL,
	flip_scan_window=NULL,
	flip_scan_window_kb=NULL,
	flip_subset=NULL,
	freq=NULL,
	freqx=NULL,
	from=NULL,
	from_bp=NULL,
	from_kb=NULL,
	from_mb=NULL,
	frqx=NULL,
	fst=NULL,
	gap=NULL,
	gates=NULL,
	gc=NULL,
	gen=NULL,
	gene=NULL,
	gene_all=NULL,
	gene_list=NULL,
	gene_list_border=NULL,
	gene_report=NULL,
	gene_report_empty=NULL,
	gene_report_snp_field=NULL,
	gene_subset=NULL,
	genedrop=NULL,
	genepi=NULL,
	geno=NULL,
	genome=NULL,
	genome_full=NULL,
	genome_lists=NULL,
	genome_minimal=NULL,
	genotypic=NULL,
	gfile=NULL,
	gplink=NULL,
	grm=NULL,
	grm_bin=NULL,
	grm_gz=NULL,
	group_avg=NULL,
	groupdist=NULL,
	gxe=NULL,
	hap...=NULL,
	hap=NULL,
	hap_assoc=NULL,
	hap_freq=NULL,
	hap_impute=NULL,
	hap_max_phase=NULL,
	hap_min_phase_prob=NULL,
	hap_miss=NULL,
	hap_phase=NULL,
	hap_phase_wide=NULL,
	hap_pp=NULL,
	hap_snps=NULL,
	hap_tdt=NULL,
	hap_window=NULL,
	hard_call_threshold=NULL,
	hardy2=NULL,
	hardy=NULL,
	help=NULL,
	het=NULL,
	hethom=NULL,
	hide_covar=NULL,
	homog=NULL,
	homozyg=NULL,
	homozyg_density=NULL,
	homozyg_gap=NULL,
	homozyg_group=NULL,
	homozyg_het=NULL,
	homozyg_include_missing=NULL,
	homozyg_kb=NULL,
	homozyg_match=NULL,
	homozyg_snp=NULL,
	homozyg_verbose=NULL,
	homozyg_window_het=NULL,
	homozyg_window_kb=NULL,
	homozyg_window_missing=NULL,
	homozyg_window_snp=NULL,
	homozyg_window_threshold=NULL,
	horse=NULL,
	hwe=NULL,
	hwe_all=NULL,
	ibc=NULL,
	ibm=NULL,
	ibs_matrix=NULL,
	ibs_test=NULL,
	id_delim=NULL,
	id_dict=NULL,
	id_match=NULL,
	iid=NULL,
	impossible=NULL,
	impute_sex=NULL,
	ind_major=NULL,
	indep=NULL,
	indep_pairphase=NULL,
	indep_pairwise=NULL,
	independent_effect=NULL,
	indiv_sort=NULL,
	inter_chr=NULL,
	interaction=NULL,
	je_cellmin=NULL,
	keep=NULL,
	keep_allele_order=NULL,
	keep_autoconv=NULL,
	keep_before_remove=NULL,
	keep_cluster_names=NULL,
	keep_clusters=NULL,
	keep_fam=NULL,
	lambda=NULL,
	lasso=NULL,
	lasso_select_covars=NULL,
	ld=NULL,
	ld_snp=NULL,
	ld_snp_list=NULL,
	ld_snps=NULL,
	ld_window=NULL,
	ld_window_kb=NULL,
	ld_window_r2=NULL,
	ld_xchr=NULL,
	lfile=NULL,
	liability=NULL,
	linear=NULL,
	list=NULL,
	list_23_indels=NULL,
	list_all=NULL,
	logistic=NULL,
	lookup...=NULL,
	lookup=NULL,
	lookup_gene=NULL,
	lookup_list=NULL,
	loop_assoc=NULL,
	maf=NULL,
	maf_succ=NULL,
	make_bed=NULL,
	make_founders=NULL,
	make_grm=NULL,
	make_grm_bin=NULL,
	make_grm_gz=NULL,
	make_just_bim=NULL,
	make_just_fam=NULL,
	make_perm_pheno=NULL,
	make_pheno=NULL,
	make_rel=NULL,
	make_set=NULL,
	make_set_border=NULL,
	make_set_collapse_group=NULL,
	make_set_complement_all=NULL,
	make_set_complement_group=NULL,
	map=NULL,
	mat=NULL,
	match=NULL,
	match_type=NULL,
	matrix=NULL,
	max=NULL,
	max_maf=NULL,
	mc=NULL,
	mcc=NULL,
	mcovar=NULL,
	mds_cluster=NULL,
	mds_plot=NULL,
	me=NULL,
	me_exclude_one=NULL,
	memory=NULL,
	mendel=NULL,
	mendel_duos=NULL,
	mendel_multigen=NULL,
	merge=NULL,
	merge_equal_pos=NULL,
	merge_list=NULL,
	merge_mode=NULL,
	merge_x=NULL,
	meta_analysis=NULL,
	meta_analysis_..._field=NULL,
	mfilter=NULL,
	mh=NULL,
	mhf=NULL,
	min=NULL,
	mind=NULL,
	mishap_window=NULL,
	missing=NULL,
	missing_code=NULL,
	missing_genotype=NULL,
	missing_phenotype=NULL,
	missing_var_code=NULL,
	mlma=NULL,
	mlma_loco=NULL,
	mlma_no_adj_covar=NULL,
	model=NULL,
	model_dom=NULL,
	model_gen=NULL,
	model_rec=NULL,
	model_trend=NULL,
	mouse=NULL,
	mperm=NULL,
	mperm_save=NULL,
	mperm_save_all=NULL,
	mpheno=NULL,
	must_have_sex=NULL,
	mwithin=NULL,
	neighbour=NULL,
	no_fid=NULL,
	no_parents=NULL,
	no_pheno=NULL,
	no_sex=NULL,
	no_snp=NULL,
	no_x_sex=NULL,
	nonfounders=NULL,
	nop=NULL,
	not_chr=NULL,
	nudge=NULL,
	null_group=NULL,
	null_snp=NULL,
	oblig_cluster=NULL,
	oblig_clusters=NULL,
	oblig_missing=NULL,
	out=NULL,
	output_chr=NULL,
	output_missing_genotype=NULL,
	output_missing_phenotype=NULL,
	oxford_pheno_name=NULL,
	parallel=NULL,
	parameters=NULL,
	parentdt1=NULL,
	parentdt2=NULL,
	pat=NULL,
	pca=NULL,
	pca_cluster_names=NULL,
	pca_clusters=NULL,
	ped=NULL,
	pedigree=NULL,
	perm=NULL,
	perm_batch_size=NULL,
	perm_count=NULL,
	pfilter=NULL,
	pheno=NULL,
	pheno_merge=NULL,
	pheno_name=NULL,
	pick1=NULL,
	plist=NULL,
	poo=NULL,
	pool_size=NULL,
	ppc=NULL,
	ppc_gap=NULL,
	proxy_...=NULL,
	proxy_assoc=NULL,
	proxy_b_kb=NULL,
	proxy_b_maxsnp=NULL,
	proxy_b_r2=NULL,
	proxy_b_threshold=NULL,
	proxy_b_window=NULL,
	proxy_dosage=NULL,
	proxy_drop=NULL,
	proxy_flanking=NULL,
	proxy_geno=NULL,
	proxy_genotypic_concordance=NULL,
	proxy_glm=NULL,
	proxy_impute=NULL,
	proxy_impute_threshold=NULL,
	proxy_kb=NULL,
	proxy_list=NULL,
	proxy_maf=NULL,
	proxy_maxsnp=NULL,
	proxy_mhf=NULL,
	proxy_r2=NULL,
	proxy_r2_no_filter=NULL,
	proxy_replace=NULL,
	proxy_show_proxies=NULL,
	proxy_sub_maxsnp=NULL,
	proxy_sub_r2=NULL,
	proxy_tdt=NULL,
	proxy_verbose=NULL,
	proxy_window=NULL,
	prune=NULL,
	q_score_file=NULL,
	q_score_range=NULL,
	qfam...=NULL,
	qmatch=NULL,
	qq_plot=NULL,
	qt=NULL,
	qt_means=NULL,
	qual_geno_...=NULL,
	qual_geno_max_threshold=NULL,
	qual_geno_scores=NULL,
	qual_geno_threshold=NULL,
	qual_max_threshold=NULL,
	qual_scores=NULL,
	qual_threshold=NULL,
	r2=NULL,
	r=NULL,
	range=NULL,
	rank=NULL,
	read_dists=NULL,
	read_freq=NULL,
	read_genome=NULL,
	read_genome_list=NULL,
	read_genome_minimal=NULL,
	recessive=NULL,
	recode12=NULL,
	recode=NULL,
	recodeA=NULL,
	recodeAD=NULL,
	recodeHV=NULL,
	recode_allele=NULL,
	recode_beagle=NULL,
	recode_bimbam=NULL,
	recode_fastphase=NULL,
	recode_lgen=NULL,
	recode_rlist=NULL,
	recode_structure=NULL,
	recode_vcf=NULL,
	recode_whap=NULL,
	reference=NULL,
	reference_allele=NULL,
	regress_distance=NULL,
	regress_pcs=NULL,
	regress_rel=NULL,
	rel_check=NULL,
	rel_cutoff=NULL,
	remove=NULL,
	remove_cluster_names=NULL,
	remove_clusters=NULL,
	remove_fam=NULL,
	rerun=NULL,
	rice=NULL,
	sample=NULL,
	score=NULL,
	score_no_mean_imputation=NULL,
	script=NULL,
	seed=NULL,
	set=NULL,
	set_by_all=NULL,
	set_collapse_all=NULL,
	set_hh_missing=NULL,
	set_max=NULL,
	set_me_missing=NULL,
	set_missing_nonsnp_ids=NULL,
	set_missing_snp_ids=NULL,
	set_missing_var_ids=NULL,
	set_names=NULL,
	set_p=NULL,
	set_r2=NULL,
	set_r2_phase=NULL,
	set_table=NULL,
	set_test=NULL,
	sex=NULL,
	sheep=NULL,
	show_tags=NULL,
	silent=NULL,
	simulate=NULL,
	simulate_haps=NULL,
	simulate_label=NULL,
	simulate_missing=NULL,
	simulate_n=NULL,
	simulate_ncases=NULL,
	simulate_ncontrols=NULL,
	simulate_prevalence=NULL,
	simulate_qt=NULL,
	simulate_tags=NULL,
	snp=NULL,
	snps=NULL,
	snps_only=NULL,
	specific_haplotype=NULL,
	split_x=NULL,
	standard_beta=NULL,
	subset=NULL,
	swap_parents=NULL,
	swap_sibs=NULL,
	swap_unrel=NULL,
	tab=NULL,
	tag_kb=NULL,
	tag_mode2=NULL,
	tag_r2=NULL,
	tail_pheno=NULL,
	tdt=NULL,
	test_all=NULL,
	test_mishap=NULL,
	test_missing=NULL,
	test_snp=NULL,
	tests=NULL,
	tfam=NULL,
	tfile=NULL,
	thin=NULL,
	thin_count=NULL,
	threads=NULL,
	to=NULL,
	to_bp=NULL,
	to_kb=NULL,
	to_mb=NULL,
	tped=NULL,
	transpose=NULL,
	trend=NULL,
	tucc=NULL,
	twolocus=NULL,
	unbounded=NULL,
	unrelated_heritability=NULL,
	update_alleles=NULL,
	update_chr=NULL,
	update_cm=NULL,
	update_ids=NULL,
	update_map=NULL,
	update_name=NULL,
	update_parents=NULL,
	update_sex=NULL,
	vcf=NULL,
	vcf_filter=NULL,
	vcf_half_call=NULL,
	vcf_idspace_to=NULL,
	vcf_min_qual=NULL,
	vegas=NULL,
	version=NULL,
	vif=NULL,
	whap=NULL,
	window=NULL,
	with_freqs=NULL,
	with_phenotype=NULL,
	with_reference=NULL,
	within=NULL,
	write_cluster=NULL,
	write_covar=NULL,
	write_dosage=NULL,
	write_set=NULL,
	write_set_r2=NULL,
	write_snplist=NULL,
	xchr_model=NULL,
	zero_cluster=NULL,
	zero_cms=NULL,
	one = NULL,
	twothreefile = NULL,
	stdout=collenv$.plink_stdout,
	stderr=collenv$.plink_stderr,
	wait=TRUE
) {
	paramList = mget(names(formals()),sys.frame(sys.nframe()))
	# Should I wait for the process to finish?
	wait = paramList$wait
	paramList$wait = NULL

	# stdout and stderr settings, default is to the R console
	stdout = paramList$stdout
	paramList$stdout = NULL
	stderr = paramList$stderr
	paramList$stderr = NULL
	
	paramVector = unlist(paramList)
	paramVector = paramVector[!is.null(paramVector)]
#	paramVector = stringr::str_trim(paramVector)
	
	
	paramName = names(paramVector)
	names(paramVector) = NULL
	paramName = gsub("_", "-", paramName)
	paramName = paste("--", paramName, sep="")
	
	
	if("--one" %in% paramName) {
		idx = which(paramName == "--one")
		paramName[idx] = "--1"
	}
	if("--twothreefile" %in% paramName) {
		idx = which(paramName == "--twothreefile")
		paramName[idx] = "--23file"
	}


	nParam = length(paramName)
	idxOdd = seq(1, nParam * 2, 2)
	idxEven = seq(2, nParam * 2, 2)
	paramNameWithValue = character(nParam * 2)
	paramNameWithValue[idxOdd] = paramName
	paramNameWithValue[idxEven] = paramVector
#	paramNameWithValue = ifelse(paramNameWithValue == "", "", paste0("'", paramNameWithValue, "'"))

	
	ret = system2('plink', paramNameWithValue, wait=wait, stdout=stdout, stderr=stderr)
    if(ret != 0) {
		warning("plink failed.")
		invisible(FALSE)
    } else {
		invisible(TRUE)
	}
}

