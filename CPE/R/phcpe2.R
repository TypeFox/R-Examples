#Qianxing Mo, moq@mskcc.org
#Department of Epidemiology and Biostatistics
#Memorial Sloan-Kettering Cancer Center, NY 10021

#The input for phcpe fuction must be a 'coxph' or 'cph' object

phcpe2 <- function(coef,coef.var,design,CPE.SE=FALSE,out.ties=FALSE){

  covar = as.matrix(coef.var)
  design = as.matrix(design)
  row <- as.integer(nrow(design))
  col <- as.integer(ncol(design))

  if(dim(covar)[1] != dim(covar)[2] || dim(covar)[1] != length(coef) || length(coef) != col){
    cat("Error: the dimensions of coef, coef.var, or design do not match!\n")
    stop("length(coef) == ncol(design) == dim(coef.var)[1] == dim(coef.var)[2]\n")
  }
  
  xbeta <- as.double(as.vector(design%*%coef))
  design <- as.double(as.vector(t(design)))
  varbeta <- as.double(as.vector(t(covar)))
  bandwidth <- as.double(0.5*sd(xbeta)*(row^(-1/3)))

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

