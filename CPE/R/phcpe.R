#Qianxing Mo, moq@mskcc.org
#Department of Epidemiology and Biostatistics
#Memorial Sloan-Kettering Cancer Center, NY 10021

#The input for phcpe fuction must be a 'coxph' or 'cph' object
#Note, the default setting for model.matrix has been changed, now it doesn't
#need to exclude the first column of the return matrix.
# old: design = model.matrix(coxfit)[,-1]; updated: design = model.matrix(coxfit)

phcpe <- function(coxfit, CPE.SE=FALSE,out.ties=FALSE) {
  if(class(coxfit)[1] != "coxph" && class(coxfit)[1] != "cph"){
    stop("Error! Input must be a coxph or cph object")
  }

  row <- as.integer(sum(coxfit$n))
  col <- as.integer(length(coxfit$coefficients))
  design <- model.matrix(coxfit)
  design <- as.double(as.vector(t(design)))  
  xbeta <- as.double(as.vector(coxfit$linear.predictors))
  varbeta <- as.double(as.vector(t(coxfit$var)))
  bandwidth <- as.double(0.5*sd(coxfit$linear.predictors)*(row^(-1/3)))

  if(CPE.SE==TRUE){
    if(row >= 3000) {
      cat(c("It may take about n*n minutes to calculate 10000*n rows of data.\n"))
    }
    if(out.ties == FALSE){
      res <- .C("coxcpe",row,col,bandwidth,xbeta,design,varbeta,out=as.double(rep(0, 3)),PACKAGE="CPE")
    }else{
      res <- .C("cpeNoTies",row,col,bandwidth,xbeta,design,varbeta,out=as.double(rep(0, 3)),PACKAGE="CPE")
    }
    return(list(CPE = res$out[1], CPE.SE = res$out[3]))
  }else {
    if(out.ties == FALSE){
      res <- .C("coxcpeOnly",row,xbeta,out=as.double(0), PACKAGE="CPE")
    }else {
      res <- .C("cpeOnlyNoTies",row,xbeta,out=as.double(0), PACKAGE="CPE")
    } 
    return(list(CPE=res$out))
  }
}

