linkR_examples <- function(name='salmon', fdir=NULL){

	# RETURN LIST
	rlist <- list()

	# GET EXTDATA DIRECTORY
	if(is.null(fdir)){
		fdir <- paste0(path.package("linkR"), "/extdata/")
	}

	# CONVERT NAME TO LOWERCASE
	name <- tolower(name)
	
	# GET SPECIES NAME AND FILE NAME
	if(name == 'salmon'){
		fname <- 'salmo_salar_1'
		name <- 'salmon'
	}
	if(grepl('salmo_salar', name, ignore.case=TRUE)){
		fname <- name
		name <- 'salmon'
	}
	if(name == 'owl'){
		fname <- 'bubo_virginianus_FMNH488595'
		name <- 'owl'
	}
	if(grepl('FMNH', name, ignore.case=TRUE)){
		fname <- gsub('fmnh', 'FMNH', name)
		name <- 'owl'
	}

	if(name == 'salmon'){

		# GET LANDMARKS
		rlist$landmarks <- as.matrix(read.table(file=paste0(fdir, fname, '.txt'), row.names=1))

		# GET LANDMARK-LINK ASSOCIATIONS
		rlist$lm.assoc.ref <- as.matrix(read.table(file=paste0(fdir, 'salmon_link_associations.txt')))

		# SET PATH CONNECTIONS
		rlist$path.connect <- list(
			c('SUSP_NC_ANT_R', 'SUSP_NC_POS_R', 'HYD_SUSP_R', 'LJ_QUAD_R', 'SUSP_NC_ANT_R'),
			c('HYD_SUSP_R', 'URO_HYP_DOR_R', 'HYD_BHYL_R', 'HYD_SUSP_R'),
			c('LJ_QUAD_R', 'LJ_SYMPH_VEN', 'LJ_SYMPH_DOR', 'LJ_QUAD_R'),
			c('URO_HYP_DOR_R', 'HYD_MID_DOR', 'HYD_MID_VEN', 'HYD_BHYL_R', 'URO_HYP_DOR_R'),
			c("NC_SAG_CURVE_POS0030", "NC_VT", "SUSP_NC_POS_R", paste0("ORB", 8:2, "_R"), "ORB_ANT_R", "PMX_TIP", "NC_SAG_ANT[0-9]+", "NC_SAG_CURVE_POS[0-9]+"),
			c("NC_SAG_CURVE_POS0030", "NC_VT", "SUSP_NC_POS_L", paste0("ORB", 8:2, "_L"), "ORB_ANT_L", "PMX_TIP", "NC_SAG_ANT[0-9]+", "NC_SAG_CURVE_POS[0-9]+"),
			c('SUSP_NC_ANT_L', 'SUSP_NC_POS_L', 'HYD_SUSP_L', 'LJ_QUAD_L', 'SUSP_NC_ANT_L'),
			c('HYD_SUSP_L', 'URO_HYP_DOR_L', 'HYD_BHYL_L', 'HYD_SUSP_L'),
			c('URO_HYP_DOR_L', 'HYD_MID_DOR', 'HYD_MID_VEN', 'HYD_BHYL_L', 'URO_HYP_DOR_L'),
			c('BHYL', 'LJ_SYMPH_VEN'),
			c('LJ_QUAD_L', 'LJ_SYMPH_VEN', 'LJ_SYMPH_DOR', 'LJ_QUAD_L'),
			c("ORB_ANT_L", paste0("ORB", 2:8, "_L"))
		)
	}

	if(name == 'owl'){

		# GET LANDMARKS
		rlist$landmarks <- as.matrix(read.table(file=paste0(fdir, fname, '.txt'), row.names=1))

		# GET LANDMARK-LINK ASSOCIATIONS
		rlist$lm.assoc.ref <- as.matrix(read.table(file=paste0(fdir, 'owl_link_associations.txt')))

		rlist$path.connect <- list(
			c("upper_bill_culmen[0-9]+", "nasalfrontalhinge_cranium_L", "jugal_upperbeak_L", "upper_bill_tomium_L[0-9]+", "upperbeak_tip", "upper_bill_culmen[0]+[1]"),
			c("upper_bill_culmen[0-9]+", "nasalfrontalhinge_cranium_R", "jugal_upperbeak_R", "upper_bill_tomium_R[0-9]+", "upperbeak_tip", "upper_bill_culmen[0]+[1]"),
			c("quadrate_cranium_mc_R", "quadrate_cranium_lc_R", "quadrate_jugal_R", "mand_condyle_quadrate_uni_R", "mand_condyle_quadrate_lat_R", "mand_condyle_quadrate_med_R", "quadrate_pterygoid_R", "quadrate_cranium_mc_R"),
			c("quadrate_cranium_mc_R", "orbital_proc_quadrate_sup_base_R", "orbital_proc_quadrate_distal_R", "orbital_proc_quadrate_inf_base_R", "quadrate_pterygoid_R"),
			c("quadrate_cranium_lc_R", "orbital_proc_quadrate_sup_base_R", "quadrate_cranium_mc_R"),
			c("palatine_lat_crest_ant_R", "palatine_lat_crest_R[0-9]+", "pterygoid_palatine_R", paste0("palatine_med_crest_R00", 30:10), paste0("palatine_med_crest_R000", 9:1), "palatine_med_crest_ant_R"),
			#c("palatine_med_crest_ant_R", "palatine_med_crest_R[0-9]+", "pterygoid_palatine_R"),
			c("palatine_lat_crest_ant_L", "palatine_lat_crest_L[0-9]+", "pterygoid_palatine_L"),
			c("palatine_med_crest_ant_L", "palatine_med_crest_L[0-9]+", "pterygoid_palatine_L"),
			c("quadrate_cranium_mc_L", "orbital_proc_quadrate_sup_base_L", "orbital_proc_quadrate_distal_L", "orbital_proc_quadrate_inf_base_L", "quadrate_pterygoid_L"),
			c("quadrate_cranium_lc_L", "orbital_proc_quadrate_sup_base_L", "quadrate_cranium_mc_L"),
			c("quadrate_cranium_mc_L", "quadrate_cranium_lc_L", "quadrate_jugal_L", "mand_condyle_quadrate_uni_L", "mand_condyle_quadrate_lat_L", "mand_condyle_quadrate_med_L", "quadrate_pterygoid_L", "quadrate_cranium_mc_L"),
			c("orbital_proc_quadrate_inf_base_L", "quadrate_jugal_L"),
			c("orbital_proc_quadrate_inf_base_R", "quadrate_jugal_R", "mand_condyle_quadrate_uni_R", "mand_condyle_quadrate_lat_R", "mand_condyle_quadrate_med_R", "quadrate_pterygoid_R"),
			c("orbital_proc_quadrate_inf_base_L", "quadrate_jugal_L", "mand_condyle_quadrate_uni_L", "mand_condyle_quadrate_lat_L", "mand_condyle_quadrate_med_L", "quadrate_pterygoid_L"),
			c("nasalfrontalhinge_cranium_R", "jugal_upperbeak_R"),
			c("nasalfrontalhinge_cranium_L", "jugal_upperbeak_L"),
			c("nasalfrontalhinge_cranium_R", "nasalfrontalhinge_cranium_L"),
			c("cranium_sagittal[0-9]+", "cranium_occipital"),
			c("quadrate_pterygoid_R", "pterygoid_palatine_R", "upperbeak_palatine_R"),
			c("quadrate_pterygoid_L", "pterygoid_palatine_L", "upperbeak_palatine_L"),
			c("quadrate_jugal_R", "jugal_upperbeak_R"),
			c("quadrate_jugal_L", "jugal_upperbeak_L")
		)
	}

	# MATCH LANDMARK AND LINK ASSOCIATIONS
	rlist$lm.assoc <- rep(NA, nrow(rlist$landmarks))
	for(i in 1:nrow(rlist$landmarks)){
		for(j in 1:nrow(rlist$lm.assoc.ref)){
			if(grepl(rlist$lm.assoc.ref[j, 1], rownames(rlist$landmarks)[i])){
				rlist$lm.assoc[i] <- rlist$lm.assoc.ref[j, 2]
				break
			}
		}
	}

	rlist
}
