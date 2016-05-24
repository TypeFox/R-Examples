#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% test file 1.2 : pipeline of the pre processing
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# setwd("V:/etude25/2891")
# cd ../../data1/etude/etude25/2891/

require(MRIaggr)

#### 0 loading (cf demo : ExtractNifti) ####
path <- system.file(file.path("nifti"), package = "MRIaggr")

param_nifti <- unlist(strsplit(list.files(path), split=c(".nii"), fixed = TRUE))
n.files <- length(param_nifti)

ls.array <- list()
for(iter_file in 1:n.files){
  ls.array[[iter_file]] <- readMRI(paste(path, param_nifti[iter_file], sep = "/"),
                                   format = "nifti")  
}
MRIaggr.Pat1 <- constMRIaggr(ls.array = ls.array,
                             identifier = "Pat1", param = param_nifti,
                             default_value = NULL,pos_default_value = c(1,1,1),
                             rm.ls.array = TRUE)


#### Step A pre-rocessing : work on the raw contrast parameters ####

#### 1- compression ####
MRIaggr.Pat1red <- constCompressMRIaggr(MRIaggr.Pat1, compression.factor = 2, verbose = TRUE)

par(mfrow=c(2,4),mar=rep(2.5,4))
multiplot(MRIaggr.Pat1, param = "DWI_t0", window = NULL, breaks = seq(0,300,1), legend = TRUE)
multiplot(MRIaggr.Pat1red, param = "DWI_t0", window = NULL, breaks = seq(0,300,1), legend = TRUE)

#### 2- allocation of the clinical data ####
df.clinical <- data.frame(Age = 68,
                          Sex = "male",
                          NIHSS_H0 = 8,
                          NIHSS_t2 = 7,
                          NIHSS_D2 = 4,
                          NIHSS_W1 = NA,
                          NIHSS_M1 = 2,
                          M1_Recanalization = 0,
                          M1_Hemorrhage = 0,
                          Symptoms_admission = 330,
                          Treatment = 0,
                          FinalStroke_volume = 29.526,
                          AcuteStroke_volume = 25.2598)

allocClinic(MRIaggr.Pat1) <- df.clinical
selectClinic(MRIaggr.Pat1)
selectClinic(MRIaggr.Pat1, param = "Age")

#### 3- computation of the area of interest before brain extraction ####

#### conversion of the mask to logical
allocContrast(MRIaggr.Pat1, param = "MASK_DWI_t0", overwrite = TRUE) <- as.logical(selectContrast(MRIaggr.Pat1, param = "MASK_DWI_t0") > 0)
allocContrast(MRIaggr.Pat1, param = "MASK_T2_FLAIR_t2", overwrite = TRUE) <- as.logical(selectContrast(MRIaggr.Pat1, param = "MASK_T2_FLAIR_t2") > 0)

#### necrosis
res <- calcTableLesion(MRIaggr.Pat1, maskN = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"))
allocTable(MRIaggr.Pat1, type = "lesion", overwrite = TRUE) <- res

calcTableLesion(MRIaggr.Pat1, maskN = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
         numeric2logical = TRUE, update.object = TRUE, overwrite = TRUE, verbose = TRUE)

selectTable(MRIaggr.Pat1, type = "lesion")



#### reperfusion-hypoperfusion
res <- calcTableHypoReperf(MRIaggr.Pat1, param = c("TTP","MTT"), timepoint = "t0",
                           overwrite = TRUE, update.object = TRUE)

selectTable(MRIaggr.Pat1, type = "hypoperfusion")

res <- calcTableHypoReperf(MRIaggr.Pat1, param = c("TTP","MTT"), timepoint = "t0",
                           mask = "MASK_DWI_t0",
                           overwrite = TRUE, update.object = TRUE)

selectTable(MRIaggr.Pat1, type = "hypoperfusion")

res <- calcTableHypoReperf(MRIaggr.Pat1, param = c("TTP","MTT"), timepoint = c("t0","t1"),
                           mask = "MASK_DWI_t0",
                           overwrite = TRUE, update.object = TRUE)


selectTable(MRIaggr.Pat1, type = "reperfusion")

#### threshold of the clinician mask
res <- calcROCthreshold(MRIaggr.Pat1, mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), param = c("DWI_t0","T2_FLAIR_t2"),
         numeric2logical = TRUE, plot = "boxplot_Youden")

calcROCthreshold(MRIaggr.Pat1, mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), param = c("DWI_t0","T2_FLAIR_t2"),
         numeric2logical = TRUE, plot = "ROC_Youden", update.object = TRUE, overwrite = TRUE)

selectDescStats(MRIaggr.Pat1, name = "Mask_threshold")


#### 4- filtering a contrast parameter ####
res <- calcFilter(MRIaggr.Pat1, param = c("T2_FLAIR_t2","DWI_t0","TTP_t0"),
                 filter="2D_M3", bilateral=FALSE, na.rm = FALSE, update.object = TRUE, overwrite = TRUE)

calcFilter(MRIaggr.Pat1, param = c("T2_FLAIR_t2","DWI_t0","TTP_t0"),
                 filter="2D_G3", bilateral=FALSE, na.rm = FALSE, update.object = TRUE, overwrite = TRUE)

par(mfrow=c(2,2),mar=rep(2,4))
multiplot(MRIaggr.Pat1, param = "T2_FLAIR_t2", main = "raw", num.main = FALSE,
              num = 1, window = NULL, breaks = c(-100,seq(0, 600),900))
multiplot(MRIaggr.Pat1, param = "T2_FLAIR_t2_2D_G3", main = "2D_G3", num.main = FALSE,
              num = 1, window = NULL, legend = FALSE, breaks = c(-100,seq(0, 600),900))
multiplot(MRIaggr.Pat1, param = "T2_FLAIR_t2_2D_M3", main = "2D_M3", num.main = FALSE,
              num = 1, window = NULL, legend = FALSE, breaks = c(-100,seq(0, 600),900))

#### 5- computation of the lesion spatial groups ####
res <- calcGroupsMask(MRIaggr.Pat1, mask = "MASK_DWI_t0", W.spatial_res = c(1.875,1.875,6),
         numeric2logical = TRUE, W.range = 6)

res <- calcGroupsMask(MRIaggr.Pat1, mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
         W.spatial_res = c(1.875,1.875,6), numeric2logical = TRUE,
         W.range = 6, update.object = TRUE, overwrite = TRUE)

selectDescStats(MRIaggr.Pat1, name = "GroupsLesion")

#### 6- brain extraction ####

#### calc.mask
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "threshold",
               th.select_optima = 2,
               overwrite = TRUE, update.object = TRUE)

multiplot(MRIaggr.Pat1, param = "mask")
multiplot(MRIaggr.Pat1, param = "T2_GRE_t0")

res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "kmeans",
                      kmeans.n_groups = 2:4,
                      update.object = TRUE, overwrite = TRUE)

multiplot(MRIaggr.Pat1, param = "mask")
multiplot(MRIaggr.Pat1, param = "mask", num = 1, legend = FALSE)

multiplot(MRIaggr.Pat1, param = "T2_GRE_t0",index1 = "mask")

#### liss.mask
res <- calcSmoothMask(MRIaggr.Pat1, update.object = TRUE, overwrite = TRUE)
            
multiplot(MRIaggr.Pat1, param = "mask", legend = FALSE)

table(selectContrast(MRIaggr.Pat1, param = "mask"))


#### update TableLesion
res <- calcTableLesion(MRIaggr.Pat1, maskN = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), mask = "mask",
     numeric2logical = TRUE, update.object=FALSE, overwrite = TRUE)
selectTable(MRIaggr.Pat1, "lesion")

#### Step B pre-processing  : work on the reduced contrast parameters ####

#### 0- contruct.red
MRIaggr.Pat1_red <- constReduceMRIaggr(MRIaggr.Pat1, mask = "mask")
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2")

#### 1- update of the area of interest ####
res <- calcTableHypoReperf(MRIaggr.Pat1_red,
                            param = c("TTP","MTT"),
                            timepoint = c("t0","t1"),
                            update.object = TRUE,
                            overwrite = TRUE
)

#### 2- probabilistic segmentation in tissue class ####
res <- calcTissueType(MRIaggr.Pat1_red, param = "T1_t0", update.object = TRUE,
                niter = 100)

multiplot(MRIaggr.Pat1_red, num = 1,
              param = c("CSF","WM","GM"), legend = FALSE,
              palette = "rgb")

multiplot(MRIaggr.Pat1_red,
              param = c("CSF","WM","GM"), legend = FALSE,
              palette = "hsv", cex = 1.25)



#### 3- identification of the hemispheres ####
res <- calcHemisphere(MRIaggr.Pat1_red, param = "T2_GRE_t0",
                mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
                verbose = TRUE, update.object = TRUE, overwrite = TRUE)

multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",
              midplane = TRUE)


#### 4- calculation of the contralateral values ####
res <- calcContralateral(MRIaggr.Pat1_red, param = c("DWI_t0","MTT_t0","TTP_t0","T2_FLAIR_t2"),
                         num = NULL, type = "mean", param.ref = "T1_t0",
                         distband = 1, lambda = 1,
                         verbose = TRUE, update.object = TRUE, overwrite = TRUE)

carto_test <- selectContrast(MRIaggr.Pat1_red, param = c("DWI_t0","DWI_t0_contro","T2_FLAIR_t2","T2_FLAIR_t2_contro","MASK_DWI_t0","MASK_T2_FLAIR_t2"))
carto_test$MASK_DWI_t0 <- as.logical(carto_test$MASK_DWI_t0)
carto_test$MASK_T2_FLAIR_t2 <- as.logical(carto_test$MASK_T2_FLAIR_t2)


par(mfrow = c(2,4), mar = rep(1.5,4), mgp  = c(2,0.5,0))
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", num = 1:3,
             window = NULL, hemisphere = "left", main = "raw - slice ")
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2_contro", num = 1:3,
             window = NULL, hemisphere = "lesion", main = "normalized - slice ")

res <- calcContralateral(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"),
                         num = NULL, type = "mean", param.ref = "T1_t0", distband = 1, lambda = 1,
                         verbose = TRUE)

multiplot(res$data[,c("i","j","k")],
             contrast = res$data[,"DWI_t0_contro", drop = FALSE],
             index1 = res$data[res$index_plot$index.plot_lesionR,c("i","j","k"), drop = FALSE]
)

multiplot(res$data[,c("i","j","k")],
             contrast = res$data[,"DWI_t0_contro", drop = FALSE],
              index1 = res$data[res$index_plot$index.plot_lesionR,c("i","j","k"), drop = FALSE],
              index2 = res$data[res$index_plot$index.plot_lesionL,c("i","j","k"), drop = FALSE]
)




#### 5- computation of the neighborhood matrix and the regional values ####
res <- calcW(MRIaggr.Pat1_red, range = 10, spatial_res = c(1.875,1.875,6),
              upper = TRUE,
              update.object = TRUE, overwrite = TRUE)

res  <- calcRegionalContrast(MRIaggr.Pat1_red, param = c("DWI_t0","MTT_t0","TTP_t0","T2_FLAIR_t2"),
							 bandwidth = 1.875, update.object = TRUE, overwrite = TRUE)

par(mfrow = c(2,4), mar = rep(1.5,4), mgp  = c(2,0.5,0))
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", num = 1:3,
             window = NULL, main = "raw - slice ")
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2_regional", num = 1:3,
             window = NULL, main = "regional - slice ")


#### 5- normalisation and tissue distribution #### 
res <- calcNormalization(MRIaggr.Pat1_red, param = c("DWI_t0","MTT_t0","TTP_t0","T2_FLAIR_t2","DWI_t0_regional","MTT_t0_regional","TTP_t0_regional"),
                         update.object = TRUE, overwrite = TRUE)

par(mfrow = c(2,4), mar = rep(1.5,4),mgp  = c(2,0.5,0))
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", num = 1:3,
             window = NULL, main = "raw - slice ")
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", num = 1:3,
             window = NULL, main = "normalized norm",
             norm_mu = "contralateral", norm_sigma = "contralateral")

names(selectNormalization(MRIaggr.Pat1_red))

selectNormalization(MRIaggr.Pat1_red, type = "global", mu = TRUE, sigma = FALSE)

res <- calcDistTissues(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"),
                       param.membership  = c("CSF","WM","GM","MASK_DWI_t0")
)

