library(testthat)
data(UPS2MS)

# AbsoluteQuantification.default
test_that("AbsoluteQuantification.default", {
	expect_that(AbsoluteQuantification.default(ProteinInference.default(UPS2_SRM, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE))[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=1.670342,"r.squared"=0.9512209,"calibration_covar"=16.9878), tolerance = .001))
	expect_that(AbsoluteQuantification.default(ProteinInference.default(UPS2_LFQ, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE))[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=2.152926,"r.squared"=0.8009153,"calibration_covar"=11.64113), tolerance = .001))
	expect_that(AbsoluteQuantification.default(ProteinInference.default(UPS2_SC, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE))[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=5.405285,"r.squared"=0.5646759,"calibration_covar"=86.19907), tolerance = .001))
})

# predict.AbsoluteQuantification
test_that("predict.AbsoluteQuantification", {
	expect_that(predict.AbsoluteQuantification(AbsoluteQuantification.default(ProteinInference.default(UPS2_SRM, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE)))$prediction$concentration[1:5],equals(c(0.29041894,4.77261553,1.57177383,5.84919886,5.46618250), tolerance = .001))
	expect_that(predict.AbsoluteQuantification(AbsoluteQuantification.default(ProteinInference.default(UPS2_LFQ, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE)))$prediction$concentration[1:5],equals(c(4.382068,5.895762,5.882757,5.963327,3.818318), tolerance = .001))
	expect_that(predict.AbsoluteQuantification(AbsoluteQuantification.default(ProteinInference.default(UPS2_SC, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE)))$prediction$concentration[1:5],equals(c(3.295982,3.295982,7.622184,4.665851,6.148056), tolerance = .001))
})

# folderror.AbsoluteQuantification
test_that("folderror.AbsoluteQuantification", {
	expect_that(folderror.AbsoluteQuantification(c(10,9,8),c(5,4,2)),equals(c(2.00,2.25,4.00)))
})

# cval.AbsoluteQuantification
test_that("cval.AbsoluteQuantification: Monte Carlo", {
	set.seed(131)
	expect_that(cval.AbsoluteQuantification(AbsoluteQuantification.default(ProteinInference.default(UPS2_SRM, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE)),cval_method="mc",mcx=2)$cval[c("r.squared","mfe")],equals(list("r.squared"=0.9292566,"mfe"=1.761073), tolerance = .001))

	set.seed(131)
	expect_that(cval.AbsoluteQuantification(AbsoluteQuantification.default(ProteinInference.default(UPS2_LFQ, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE)),cval_method="mc",mcx=2)$cval[c("r.squared","mfe")],equals(list("r.squared"=0.5639523,"mfe"=2.804344), tolerance = .001))

	set.seed(131)
	expect_that(cval.AbsoluteQuantification(AbsoluteQuantification.default(ProteinInference.default(UPS2_SC, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE)),cval_method="mc",mcx=2)$cval[c("r.squared","mfe")],equals(list("r.squared"=0.4131045,"mfe"=6.591413), tolerance = .001))
})

# AbsoluteQuantification.default
test_that("AbsoluteQuantification.default", {
	set.seed(131)
	REP1<-UPS2_SRM[sample(854,800),]
	REP1[which(REP1$protein_id %in% unique(REP1$protein_id)[1:10]),]$concentration<-"?"
	REP1$run_id<-"UPS2_SRM_REP1"

	REP2<-UPS2_SRM[sample(854,800),]
	REP2[which(REP2$protein_id %in% unique(REP2$protein_id)[1:10]),]$concentration<-"?"
	REP2$run_id<-"UPS2_SRM_REP2"

	REP3<-UPS2_SRM[sample(854,800),]
	REP3[which(REP3$protein_id %in% unique(REP3$protein_id)[1:10]),]$concentration<-"?"
	REP3$run_id<-"UPS2_SRM_REP3"

	UPS2_SRM_REP<-rbind(REP1,REP2,REP3)

	data_abs<-AbsoluteQuantification.default(ProteinInference.default(UPS2_SRM_REP, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE))
	data_abs_predict<-predict(data_abs)
	data_abs_cval<-cval(data_abs,mcx=2)
	data_abs_cval_predict<-predict(data_abs_cval)

	expect_that(data_abs[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=1.820787,"r.squared"=0.9425091,"calibration_covar"=17.84446), tolerance = .001))
	expect_that(data_abs_predict[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=1.820787,"r.squared"=0.9425091,"calibration_covar"=17.84446), tolerance = .001))
	expect_that(data_abs_cval$cval[c("mfe","r.squared")],equals(list("mfe"=1.887385,"r.squared"=0.9339684), tolerance = .001))
	expect_that(data_abs_cval_predict$cval[c("mfe","r.squared")],equals(list("mfe"=1.887385,"r.squared"=0.9339684), tolerance = .001))
	expect_that(pivot(data_abs),throws_error("Apply predict before pivot to AbsoluteQuantification object"))
	expect_that(pivot(data_abs_predict),equals(structure(c(1.53455646337375, 145.93611087053, 5.64284491012796, 435.816627098295, 295.299516247426, 0.683691282425213, 178.206176091835, 0.419191108154238, 39.1400080492642, 0.832309209124821, 99.7714602473884, 6.43660876450224, 866.239803375142, 36.9583880103523, 5.44454279866606, 5.01041177221375, 24.1166828284089, 5.85570926647292, 486.44630632458, 0.167820084579975, 2.47506259687646, 0.248808640247866, 128.165736354375, 297.633563973043, 45.106249617401, 274.135089492113, 191.93098648792, 40.9880277905149, 1.53455646337375, 145.93611087053, 5.64284491012796, 435.816627098295, 295.299516247426, 0.683691282425213, 178.206176091835, 0.419191108154238, 39.1400080492642, 0.832309209124821, 99.7714602473884, 6.43660876450224, 866.239803375142, 36.9583880103523, 5.44454279866606, 5.01041177221375, 24.1166828284089, 5.85570926647292, 486.44630632458, 0.167820084579975, 2.47506259687646, 0.248808640247866, 128.165736354375, 297.633563973043, 45.106249617401, 274.135089492113, 191.93098648792, 40.9880277905149, 1.53455646337375, 145.93611087053, 5.64284491012796, 435.816627098295, 295.299516247426, 0.683691282425213, 178.206176091835, 0.419191108154238, 39.1400080492642, 0.832309209124821, 99.7714602473884, 6.43660876450224, 866.239803375142, 36.9583880103523, 5.44454279866606, 5.01041177221375, 24.1166828284089, 5.85570926647292, 486.44630632458, 0.167820084579975, 2.47506259687646, 0.248808640247866, 128.165736354375, 297.633563973043, 45.106249617401, 274.135089492113, 191.93098648792, 40.9880277905149), .Dim = c(28L, 3L), .Dimnames = list(c("O76070ups|SYUG_HUMAN_UPS", "P00167ups|CYB5_HUMAN_UPS", "P00709ups|LALBA_HUMAN_UPS", "P00915ups|CAH1_HUMAN_UPS", "P00918ups|CAH2_HUMAN_UPS", "P01008ups|ANT3_HUMAN_UPS", "P01031ups|CO5_HUMAN_UPS", "P01127ups|PDGFB_HUMAN_UPS", "P01133ups|EGF_HUMAN_UPS", "P01344ups|IGF2_HUMAN_UPS", "P02144ups|MYG_HUMAN_UPS", "P02753ups|RETBP_HUMAN_UPS", "P02768ups|ALBU_HUMAN_UPS", "P04040ups|CATA_HUMAN_UPS", "P06732ups|KCRM_HUMAN_UPS", "P12081ups|SYHC_HUMAN_UPS", "P15559ups|NQO1_HUMAN_UPS", "P16083ups|NQO2_HUMAN_UPS", "P41159ups|LEP_HUMAN_UPS", "P55957ups|BID_HUMAN_UPS", "P61626ups|LYSC_HUMAN_UPS", "P61769ups|B2MG_HUMAN_UPS", "P62937ups|PPIA_HUMAN_UPS", "P62988ups|UBIQ_HUMAN_UPS", "P63165ups|SUMO1_HUMAN_UPS", "P68871ups|HBB_HUMAN_UPS", "P69905ups|HBA_HUMAN_UPS", "Q06830ups|PRDX1_HUMAN_UPS"), c("UPS2_SRM_REP1", "UPS2_SRM_REP2", "UPS2_SRM_REP3"))), tolerance = .001))

	data_absnc<-AbsoluteQuantification.default(ProteinInference.default(UPS2_SRM_REP, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose",peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose",transition_summary = "sum",combine_precursors = FALSE, consensus_proteins = FALSE, consensus_peptides = FALSE, consensus_transitions = FALSE))
	data_absnc_predict<-predict(data_absnc)
	data_absnc_cval<-cval(data_absnc,mcx=2)
	data_absnc_cval_predict<-predict(data_absnc_cval)

	expect_that(data_absnc[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=1.765394,"r.squared"=0.9472684,"calibration_covar"=17.88808), tolerance = .001))
	expect_that(data_absnc_predict[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=1.765394,"r.squared"=0.9472684,"calibration_covar"=17.88808), tolerance = .001))
	expect_that(data_absnc_cval$cval[c("mfe","r.squared")],equals(list("mfe"=1.812872,"r.squared"=0.9388508), tolerance = .001))
	expect_that(data_absnc_cval_predict$cval[c("mfe","r.squared")],equals(list("mfe"=1.812872,"r.squared"=0.9388508), tolerance = .001))
	expect_that(pivot(data_absnc),throws_error("Apply predict before pivot to AbsoluteQuantification object"))
	expect_that(pivot(data_absnc_predict),equals(structure(c(1.38223739187041, 126.564657393243, 5.02799902824017, 374.543541167615, 254.605017913418, 0.619982830852972, 611.234999275816, 0.381679074417682, 34.3181129754078, 0.895186461344148, 93.7527986592859, 5.95305205773649, 740.211406818313, 46.9742161587062, 4.85274741132051, 6.18408673559932, 21.2309239420806, 5.21606323620805, 417.673100504483, 0.153970070428787, 2.2205449586935, 0.410444430505147, 111.383768391298, 256.600616226537, 43.1914863336051, 236.503425764807, 221.213305211335, 51.2702026363797, 1.38223739187041, 126.564657393243, 5.02799902824017, 374.543541167615, 254.605017913418, 0.619982830852972, 154.294681019423, 0.381679074417682, 34.3181129754078, 0.936237707155779, 86.8019358149232, 5.95305205773649, 740.211406818313, 38.6132089077754, 4.85274741132051, 4.46889127816852, 21.2309239420806, 5.21606323620805, 432.018798210824, 0.153970070428787, 2.2205449586935, 0.227528422249866, 111.273192807584, 256.600616226537, 69.9631961910401, 241.832257481942, 166.075414709717, 35.9246789226138, 1.38223739187041, 126.564657393243, 5.02799902824017, 374.543541167615, 254.605017913418, 0.619982830852972, 528.141323004369, 0.381679074417682, 34.3181129754078, 0.936237707155779, 93.7527986592859, 5.72900052763951, 740.211406818313, 46.9742161587062, 4.85274741132051, 6.18408673559932, 21.2309239420806, 5.21606323620805, 432.018798210824, 0.153970070428787, 2.2205449586935, 0.410444430505147, 111.383768391298, 256.600616226537, 69.9631961910401, 241.832257481942, 187.545323236375, 42.5152748447174), .Dim = c(28L, 3L), .Dimnames = list(c("O76070ups|SYUG_HUMAN_UPS", "P00167ups|CYB5_HUMAN_UPS", "P00709ups|LALBA_HUMAN_UPS", "P00915ups|CAH1_HUMAN_UPS", "P00918ups|CAH2_HUMAN_UPS", "P01008ups|ANT3_HUMAN_UPS", "P01031ups|CO5_HUMAN_UPS", "P01127ups|PDGFB_HUMAN_UPS", "P01133ups|EGF_HUMAN_UPS", "P01344ups|IGF2_HUMAN_UPS", "P02144ups|MYG_HUMAN_UPS", "P02753ups|RETBP_HUMAN_UPS", "P02768ups|ALBU_HUMAN_UPS", "P04040ups|CATA_HUMAN_UPS", "P06732ups|KCRM_HUMAN_UPS", "P12081ups|SYHC_HUMAN_UPS", "P15559ups|NQO1_HUMAN_UPS", "P16083ups|NQO2_HUMAN_UPS", "P41159ups|LEP_HUMAN_UPS", "P55957ups|BID_HUMAN_UPS", "P61626ups|LYSC_HUMAN_UPS", "P61769ups|B2MG_HUMAN_UPS", "P62937ups|PPIA_HUMAN_UPS", "P62988ups|UBIQ_HUMAN_UPS", "P63165ups|SUMO1_HUMAN_UPS", "P68871ups|HBB_HUMAN_UPS", "P69905ups|HBA_HUMAN_UPS", "Q06830ups|PRDX1_HUMAN_UPS"), c("UPS2_SRM_REP1", "UPS2_SRM_REP2", "UPS2_SRM_REP3"))), tolerance = .001))
})
