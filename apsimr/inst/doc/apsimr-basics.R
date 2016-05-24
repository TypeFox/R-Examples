## ----run,eval=FALSE------------------------------------------------------
#  apsimExe <-"C:/Program Files (x86)/Apsim75-r3008/Model/Apsim.exe"
#  apsimWd <- "~/APSIM"
#  toRun <- c("Canopy.apsim", "Continuous Wheat.apsim")
#  results <- apsim(exe = apsimExe, wd = apsimWd, files = toRun)

## ----p1,eval=FALSE-------------------------------------------------------
#  plot(results$"Continuous Wheat", geom = 'line')

## ----p2,eval=FALSE-------------------------------------------------------
#  plot(results$"Continuous Wheat", one_plot = TRUE, geom = 'line') + theme_bw()

## ----p3,eval=FALSE-------------------------------------------------------
#  plot(results$"Continuous Wheat", y = 'yield') + geom_line(colour = 'red') + theme_bw()

## ----file,eval=FALSE-----------------------------------------------------
#  apsimFile <- "Canopy.apsim"
#  apsimWd <- "~/APSIM"

## ----vars,eval=FALSE-----------------------------------------------------
#  apsimVar <- c("SoilWater/Thickness", "SoilOrganicMatter/SoilCN")
#  apsimValue <- list(c(rep(200, 2), rep(300, 9)), 10)

## ----edit,eval=FALSE-----------------------------------------------------
#  edit_apsim(file = apsimFile, wd = apsimWd, var = apsimVar,
#             value = apsimValue, overwrite = FALSE)

## ----simfileedit,eval=FALSE----------------------------------------------
#  simFile <- "Soil.xml"
#  simVar <- c("nitrification_pot", "dnit_nitrf_loss","wfnit_values")
#  simValue <- list(abs(rnorm(1)), abs(rnorm(1)), c(0,2,2,1))
#  edit_sim_file(file = simFile, wd = apsimWd, var = simVar,
#                value = simValue, overwrite = FALSE)

## ----apsimSA,eval=FALSE,size="footnotesize"------------------------------
#  library(sensitivity)
#  
#  meanYield<-function(x){
#    return(mean(x$lai_cowpea))
#  }
#  
#  vars <- c("SoilOrganicMatter/SoilCN", "SoilWater/DiffusConst", "SoilWater/CNCov")
#  
#  n <- 20
#  X1 <- data.frame(SoilCN = runif(n, 5, 25),
#                   DiffusConst = runif(n, 20, 50), CNCov = runif(n, 0, 1))
#  X2 <- data.frame(SoilCN = runif(n, 5, 25),
#                   DiffusConst = runif(n, 20, 50), CNCov = runif(n, 0, 1))
#  
#  sobolResults <- soboljansen(model = apsim_vector, X1, X2, exe = apsimExe, wd = apsimWd,
#                          vars = vars, to.run = apsimFile, g = meanYield, overwrite = TRUE)
#  plot(sobolResults)

## ----apsimSAEmulator,eval=FALSE,size="footnotesize"----------------------
#  n <- 61
#  emulX <- data.frame(SoilCN = runif(n, 5, 25),
#                   DiffusConst = runif(n, 20, 50), CNCov = runif(n, 0, 1))
#  emulRes <- apsim_emul_sa(model = apsim_vector, X = emulX, method = "singleGAM",
#                           exe = apsimExe, wd = apsimWd, vars = vars, to.run = apsimFile,
#                           g = meanYield, overwrite = TRUE)
#  plot(emulRes)

## ----apsimMGSA,eval=FALSE,size="footnotesize"----------------------------
#  
#  rawYield <- function(x){
#    return(x$lai_cowpea)
#  }
#  
#  cowpeaY <- apsim_vector(X = emulX, exe = apsimExe, wd = apsimWd,
#                          vars = vars, to.run = apsimFile, g = rawYield,
#                          overwrite = TRUE)
#  
#  vCowpea<-var(cowpeaY)
#  pcsY<-svd(vCowpea)$u
#  N<-nrow(cowpeaY)
#  ones<-matrix(rep(1,N),ncol=1)
#  yBar<-(1/N)*ones%*%t(ones)%*%cowpeaY
#  cowpeaPCS<-(cowpeaY-yBar)%*%pcsY
#  
#  mgsaRes <- vector("list",5)
#  
#  for(i in 1:5){
#    mgsaRes[[i]] <- apsim_emul_sa(y = cowpeaPCS[,i], X = emulX, method = "separateGAM")
#  }
#  

