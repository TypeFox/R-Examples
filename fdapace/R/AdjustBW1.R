# This function adjusts the bandwidth for smoothing of mean function in PCA
##########################################################################
# Input:  - kernel: 'gauss', 'epan' or other kernel type
#         - bopt: working bandwidth for mean function smoothing
#         - npoly: (Not used)
#         - nder: degree of derivative need to be estimated
#         - dataType: whether the functional data is dense and dataType ("Dense")
#         - verbose: 
##########################################################################
# Output: - bopt: bandwidth after adjustment
##########################################################################

AdjustBW1 <- function(kernel, bopt, npoly, nder, dataType, verbose){
  # for Gaussian kernel
  if(kernel == 'gauss'){
    if(dataType == "Dense"){
      bwmu_fac = c(1.1, 0.8, 0.8)
    }
    else bwmu_fac = c(1.1, 1.2, 2)
    
    if(nder > 2){ facID = 3 }
    else if(nder >= 0 && nder <= 2){facID = nder + 1 }
    else facID = 1
    
    bopt = bopt * bwmu_fac[facID]
  }
  # for Epanechnikov kernel
  if(kernel == 'epan'){
    if(dataType == "Dense"){
      bwmu_fac = c(1.1, 1.0, 1.1)
    }
    else bwmu_fac = c(1.1, 1.2, 1.5)
    
    if(nder > 2){ facID = 3 }
    else if(nder >= 0 && nder <= 2){ facID = nder + 1 }
    else facID = 1
    
    bopt = bopt * bwmu_fac[facID]
  }
  bopt
}
