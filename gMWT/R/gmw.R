# Version: 06-07-2013, Df

# Changes:
# 28-06-2013: Added the keepPM option, DF
# 28-06-2013: Changed the layout and added comments, DF
# 06-07-2013: Added the mwAkw option

gmw <- function(X, g, goi=NULL, test="mw", type="permutation", prob="pair", nper=2000, alternative="greater", mc=1, output="min", cluster=NULL, order=TRUE, keepPM= FALSE, mwAkw=FALSE, alg=NULL){

# Input checks
 if(is.null(alg)) alg <- "Cnaive"
 type <- match.arg(type,c("permutation","external","asymptotic"))
 test <- match.arg(test,c("uit","triple","jt","jt*","mw","kw"))
 prob <- match.arg(prob,c("single","pair","triple"))
 alg <- match.arg(alg,c("Cnaive","Rnaive","Csubmat","Rsubmat"))
 alternative <- match.arg(alternative,c("smaller","greater","two.sided"))
 output <- match.arg(output,c("min","full"))
 if(mwAkw==TRUE) mwAkw <- 1
 if(keepPM & type!="permutation") warning("Input mismatch! KeepPM=TRUE requires also type='permutation'. There won't be now a permutation matrix available for W&Y.")


# Adjust the group labels
  g <- relabelGroups(g)

# Get flags, if the input X will be analyzed as vector or matrix and calculate certain constants
  dimX <- dim(X)
  XisVector <- is.null(dimX)
  if(is.null(goi)) goi <- g
  goi <- unique(goi)

# Initialze the class variables for the output
  DNAME <- paste(deparse(substitute(X)),", grouping vector:",deparse(substitute(group)))
  TEST<- switch(test,"uit"="Union Intersection Test",
                     "triple"="Triple Based Test",
 		     "jt"="Jonckheere-Terpstra",
		     "jt*"="Jonckheere-Terpstra *",
		     "mw"="Mann-Whitney Test",
		     "kw"="Kruskal-Wallis Test")
  TYPE <- switch(type,"permutation"="Permutation Test",
                      "asymptotic"="Asymptotic Test",
		      "external"="Included in base/other package Test")
  ALTERNATIVE=""
  STATISTIC=""
  PVAL=""

  PARAMETERS <- list(DNAME,TEST,TYPE,ALTERNATIVE,STATISTIC,PVAL,dimX,XisVector)

  if(test=="uit"){

    res <- uit.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output, keepPM=keepPM)

  } else if(test=="triple"){

    res <- triple.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output,alg, keepPM=keepPM, order=order)

} else if(test=="mw"){

    if(prob=="pair")
    {
      res <- mw.gmw(X=X,g=g,goi=goi,type=type,nper=nper,alternative=alternative,mc=mc,PARAMETERS,output=output,order=order,keepPM=keepPM)
    } else if(prob=="single"){
      res <- mw.single.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output, keepPM=keepPM)
    }

} else if(test=="jt"){

    res <- jt.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output, keepPM=keepPM)

} else if(test=="jt*"){

    res <- jt.star.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output, keepPM=keepPM)

} else if(test=="kw"){
    if(alternative!="two.sided"){
      alternative <- "two.sided"
     # warning("I switched the alternative to 'two.sided' for the Kruskal-Wallis test.")
    }
    if(is.numeric(mwAkw)){
      if(mwAkw<0 || mwAkw >1) stop("The significance level alpha has to be between 0 and 1.")
      resKW <- kw.gmw(X,g,cluster,goi,type,nper,mc,PARAMETERS,output, keepPM=keepPM)
      takeMW <- which(resKW$p.values<mwAkw)
      PARAMETERS[[2]] <- "Mann-Whitney Test"
      res <- mw.gmw(X=X,g=g,goi=goi,type=type,nper=nper,alternative="two.sided",mc=mc,PARAMETERS,output=output,order=TRUE, keepPM=keepPM)
      res$p.values <- rbind(resKW$p.values,res$p.values)
      temp <- rownames(res$p.values)
      res$p.values <- res$p.values[,which((res$p.values[1,]<=mwAkw)==TRUE)]
      rownames(res$p.values) <- temp
    } else if(mwAkw==FALSE){
    res <- kw.gmw(X,g,cluster,goi,type,nper,mc,PARAMETERS,output, keepPM=keepPM)
    } else {
      stop("You picked a non-valid option for 'mwAkf'. Please pick a sig.-level for the pre-KW test, or choose the default 'FALSE'.")
    }
 } else { 

    stop("There is not such test, sorry!")
    res <- NULL

 }
 class(res) <- "gmw"
 attr(res,"keepPM") <- keepPM
 attr(res,"test") <- test
 attr(res,"alternative") <- alternative
 return(res)
}