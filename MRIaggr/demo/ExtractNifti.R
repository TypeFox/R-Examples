#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Fichier test 1.1 : loading NIFTI files and convert them to carto3D and MRIaggr object
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# setwd("V:/etude25/2891")

require(MRIaggr)

path <- system.file(file.path("nifti"), package = "MRIaggr")

#### A- conversion to Carto3D then to IRMaggr ####
#### 1- Import ####

#### manual ####

nifti.Pat1_MTT_t0 <- readMRI(file.path(path, "MTT_t0"), format = "nifti")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path, "TTP_t0"), format = "nifti")

nifti.Pat1_MTT_t1 <- readMRI(paste(path, "MTT_t1", sep = "/"), format = "nifti")
nifti.Pat1_TTP_t1 <- readMRI(paste(path, "TTP_t1", sep = "/"), format = "nifti")

nifti.Pat1_DWI_t0 <- readMRI(paste(path, "DWI_t0", sep = "/"), format = "nifti")
nifti.Pat1_T1_t0 <- readMRI(paste(path, "T1_t0", sep = "/"), format = "nifti")
nifti.Pat1_T2_GRE_t0 <- readMRI(paste(path, "T2_GRE_t0", sep = "/"), format = "nifti")
nifti.Pat1_T2_FLAIR_t2 <- readMRI(paste(path, "T2_FLAIR_t2", sep = "/"), format = "nifti")

nifti.Pat1_MASK_DWI_t0 <- readMRI(paste(path, "MASK_DWI_t0", sep = "/"), format = "nifti")
nifti.Pat1_MASK_T2_FLAIR_t2 <- readMRI(paste(path, "MASK_T2_FLAIR_t2", sep = "/"), format = "nifti")


#### automatic ####

param_nifti <- unlist(strsplit(list.files(path), split = ".nii" , fixed = TRUE))
n.files <- length(param_nifti)

for(iter_file in 1:n.files){
  assign(x = paste("nifti.Pat1_", param_nifti[iter_file], sep = ""),
         value = readMRI(paste(path, param_nifti[iter_file], sep = "/"),
                         format = "nifti")
  )  
}

#### display  ####
fields::image.plot(nifti.Pat1_DWI_t0[,,1,1])


#### 2- conversion from nifti to Carto3D ####

#### manual ####
ls.Carto3D.Pat1 <- list()

ls.Carto3D.Pat1$MTT_t0 <- constCarto3D(nifti.Pat1_MTT_t0,
                                    identifier = "Pat1", param = "MTT_t0")
ls.Carto3D.Pat1$TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,
                                     identifier = "Pat1", param = "TTP_t0")

ls.Carto3D.Pat1$MTT_t1 <- constCarto3D(nifti.Pat1_MTT_t1,
                                     identifier = "Pat1", param = "MTT_t1")
ls.Carto3D.Pat1$TTP_t1 <- constCarto3D(nifti.Pat1_TTP_t1,
                                     identifier = "Pat1", param = "TTP_t1")

ls.Carto3D.Pat1$MASK_DWI_t0 <- constCarto3D(nifti.Pat1_MASK_DWI_t0,
                                     identifier = "Pat1", param = "MASK_DWI_t0")
ls.Carto3D.Pat1$MASK_T2_FLAIR_t2 <- constCarto3D(nifti.Pat1_MASK_T2_FLAIR_t2,
                                     identifier = "Pat1", param = "MASK_T2_FLAIR_t2")

ls.Carto3D.Pat1$T1.H0 <- constCarto3D(nifti.Pat1_T1_t0,
                                     identifier = "Pat1", param = "T1_t0")
ls.Carto3D.Pat1$T2_GRE.H0 <- constCarto3D(nifti.Pat1_T2_GRE_t0,
                                      identifier = "Pat1", param = "T2_GRE_t0")

multiplot(ls.Carto3D.Pat1$T1.H0)

#### automatic ####

ls.Carto3D.Pat1 <- list()
for(iter_file in 1:n.files){
  eval(parse(text=paste(
    "ls.Carto3D.Pat1[[iter_file]] <- constCarto3D(nifti.Pat1_", param_nifti[iter_file],",
                        identifier = \"Pat1\", param = param_nifti[iter_file])",
    sep = "")))
}
names(ls.Carto3D.Pat1) <- param_nifti

#### display 
multiplot(ls.Carto3D.Pat1[[1]], num = 1,
                window = FALSE)
  
multiplot(ls.Carto3D.Pat1[[1]], num = 1:3,
                window = FALSE)
  
multiplot(ls.Carto3D.Pat1[[1]], num = 1:3, axes = FALSE, legend = FALSE,
                window = FALSE)
  


#### 3- conversion from Carto3D to MRIaggr ####

MRIaggr.Pat1 <- Carto3D2MRIaggr(ls.Carto3D = ls.Carto3D.Pat1,
                                rm.Carto3D = TRUE)

#### display
multiplot(MRIaggr.Pat1, num = 1:3, param = "DWI_t0",
            window = FALSE)

multiplot(MRIaggr.Pat1, num = 1:3, param = "DWI_t0",
            window = FALSE, axes = FALSE)



#### B- Direct conversion to MRIaggr ####


#### B1- import ####
param_nifti <- unlist(strsplit(list.files(path), split = ".nii", fixed = TRUE))
n.files <- length(param_nifti)

ls.array <- list()
for(iter_file in 1:n.files){
  ls.array[[iter_file]] <- readMRI(paste(path, param_nifti[iter_file], sep = "/"),
                                   format = "nifti")  
}

#### display 
fields::image.plot(ls.array[[1]][,,1,1])


#### B2- conversion to MRIaggr ####
MRIaggr.Pat1 <- constMRIaggr(ls.array = ls.array,
                              identifier = "Pat1", param = param_nifti,
                              default_value = NULL, pos_default_value = c(1,1,1),
                              rm.ls.array = TRUE)

#### display
multiplot(MRIaggr.Pat1, num = 1:3, param = "DWI_t0",
               window = FALSE, legend = FALSE)



