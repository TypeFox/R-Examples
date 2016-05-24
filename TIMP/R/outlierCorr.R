"outlierCorr" <- function(oldRes, fence=3, saturCorr=FALSE, saturThresh=.05,
                          saturMin=NA, saturDivMax=3, outlierCorr=TRUE,
                          newM = TRUE) {
  ## this function takes the result returned by fitModel
  ## and returns the list of datasets with saturation correction applied 
  m <- oldRes$currModel@modellist
  res <- oldRes$currModel@fit@resultlist
  dt <- oldRes$currModel@data
  weightList <- vector("list", length(m))
  if(outlierCorr) {
    cntall <- 0
    for(i in 1:length(m)){
      if(identical(dt[[i]]@outMat,matrix() ) || newM ) 
        dt[[i]]@outMat <- matrix(FALSE, nrow=nrow(dt[[i]]@psi.df),
                                ncol=ncol(dt[[i]]@psi.df))
      cnt <- 0
      x2 <- m[[i]]@x2
      x <- m[[i]]@x
      weightList[[i]] <- matrix(1, nrow=m[[i]]@nt, ncol=m[[i]]@nl)
      fitted <- resid <- matrix(nrow = m[[i]]@nt, ncol = m[[i]]@nl)
      for(j in 1:length(x2)){
        if(m[[i]]@clpType == "x2"){
          fitted[,j] <- res[[i]]@fitted[[j]]
          resid[,j] <- res[[i]]@resid[[j]]  
        }
        if(m[[i]]@clpType == "x"){
          fitted[j,] <- res[[i]]@fitted[[j]]
          resid[j,] <- res[[i]]@resid[[j]]  
        }
        
      }
      for(j in 1:length(x2)) {
        rr <- resid[,j]
        iq <- IQR(rr)
        fH <- quantile(rr,3/4)
        fL <- quantile(rr,1/4)
        outlow <-  which(rr > (fH+fence*iq))
        outhigh <- which(rr < (fL-fence*iq))
        outind <- append(outlow, outhigh)
        if(length(outind)>0) {
          ## dt[[i]]@psi.df[outind,j] <- fitted[outind,j]
          ## above line was for plotting but now want to plot outliers w/sym.
          dt[[i]]@outMat[outind,j] <- TRUE
          weightList[[i]][outind,j] <- 1e-10
          cnt <- cnt + length(outind)
        }
      }
      cntall <- cntall + cnt 
      cat("Removed", cnt, "outliers in dataset", i,"\n") 
    }
    cat("Removed", cntall, "outliers in total\n")
  }
  if(saturCorr) {
    cntall <- 0
    for(i in 1:length(m)){
      if(identical(dt[[i]]@outMat,matrix() ) || newM )
        dt[[i]]@satMat <- matrix(FALSE, nrow=nrow(dt[[i]]@psi.df),
                                 ncol=ncol(dt[[i]]@psi.df))
      if(is.na(saturMin)) 
        saturMin <- max(dt[[i]]@psi.df)/saturDivMax
      cnt <- 0
      x2 <- m[[i]]@x2
      x <- m[[i]]@x
      weightList[[i]] <- matrix(1, nrow=m[[i]]@nt, ncol=m[[i]]@nl)
      fitted <- resid <- matrix(nrow = m[[i]]@nt, ncol = m[[i]]@nl)
      for(j in 1:length(x2)){
        fitted[,j] <- res[[i]]@fitted[[j]]
        resid[,j] <- res[[i]]@resid[[j]]  
      }
      for(j in 1:length(x2)) {
        rr <- resid[,j]
        satind <-  intersect(which(- rr/fitted[,j] > saturThresh),
                             which(dt[[i]]@psi.df[,j] > saturMin))
        if(length(satind)>0) {
          ##dt[[i]]@psi.df[satind,j] <- fitted[satind,j]
          ## above line was for plotting but now want to plot outliers w/sym.
          dt[[i]]@satMat[satind,j] <- TRUE
          weightList[[i]][satind,j] <- 1e-10
          cnt <- cnt + length(satind)
        }
      }
      cntall <- cntall + cnt 
      cat("Removed", cnt, "saturated observations in dataset", i,"\n") 
    }
    cat("Removed", cntall, "saturated observations in total\n")
  }  
  list(dt=dt, weightList = weightList)
}
