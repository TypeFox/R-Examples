#<<BEGIN>>
mcratio <- function(x, pcentral = .5, pvar = .975, punc = .975, na.rm=FALSE)
#TITLE Ratio of uncertainty and the variability
#DESCRIPTION
# Provides measures of variability, uncertainty, and both combined for an 
#\samp{mc} or an \samp{mcnode} object.
#KEYWORDS distribution
#INPUTS
#{x}<<an \samp{mc} or an \samp{mcnode} object>>
#[INPUTS]
#{pcentral}<<the quantile for the central tendency.>>         .
#{pvar}<<the quantile for the measure of variability.>>
#{punc}<<the quantile for the measure of uncertainty.>>
#{na.rm}<<a logical value indicating whether NA values should be stripped before the computation proceeds.>>
#VALUE
#A matrix. 
#DETAILS
#The function evaluates three ratios for each \samp{mcnode}. Given:
#{A}<<the \samp{(100 * pcentral)}th percentile of uncertainty (by default the median) for the \samp{(100 * pcentral)}th percentile of variability>>
#{B}<<the \samp{(100 * pcentral)}th percentile of uncertainty for the \samp{(100 * pvar)}th percentile of variability>>
#{C}<<the \samp{(100 * punc)}th percentile of uncertainty for the \samp{(100 * pcentral)}th percentile of variability>>
#{D}<<the \samp{(100 * punc)}th percentile of uncertainty for the \samp{(100 * pvar)}th percentile of variability>>

#The following ratio are estimated
#{*}<<Variability Ratio: B / A>>
#{*}<<Uncertainty Ratio: C / A>>
#{*}<<Overall Uncertainty Ratio: D / A>>

#For multivariate nodes, the statistics are evaluate on each dimension or on statistics according to the corresponding \samp{outm} value. 

#EXAMPLE
#data(total)
#mcratio(total, na.rm=TRUE)

#REFERENCE
#Ozkaynak, H., Frey, H.C., Burke, J., Pinder, R.W. (2009) "Analysis of coupled model uncertainties in source-to-dose modeling
#of human exposures to ambient air pollution: A PM2.5 case study", Atmospheric environment, Volume 43, Issue 9, March 2009, Pages 1641-1649.
#CREATED 10-02-10
#--------------------------------------------

{
  if(inherits(x,"mcnode")) {
    nom <- deparse(substitute(x))
    x <- list(x)
    names(x) <- nom}
  
  outm <- lapply(x,attr,which="outm")

  x <- x[outm!="none"]
  outm <- outm[outm!="none"]
  typen <- sapply(x,attr,which="type")
  dimm <- sapply(x,dim)
  nvariates <- dimm[3,]
  nom <- names(x)
  
  lesnom <- function(outm,nvariates,nom){
    if(outm[1] == "each"){
      if(nvariates==1) nomsortie <- nom
      else nomsortie <- paste(nom,1:nvariates,sep=".")}
    else nomsortie <- paste(nom,": ",outm," of variates",sep="")
    return(nom=nomsortie)}

  nomsortie <- mapply(lesnom,outm,nvariates,nom,SIMPLIFY=FALSE)

  gerout <- function(node,outm,nvariates){
    if(outm[1]=="each"){
      if(nvariates==1) return(node)
      return(lapply(1:nvariates,function(x) node[,,x,drop=FALSE]))
      }
    res <- mapply(function(outm){
                  func <- get(outm,mode="function")
                  node <- apply(node,c(1,2),func)
                  dim(node) <- c(dim(node),1)
                  node},outm,SIMPLIFY=FALSE)
    if(length(outm)==1) return(res[[1]])
    return(res)
  }

  x <- mapply(gerout,x,outm,nvariates,SIMPLIFY=FALSE)

  LESSTAT <- function(x,prob,wdim){
    if(is.list(x)) return(lapply(x,LESSTAT,prob=prob,wdim=wdim))
    return(apply(x,wdim,quantile,prob = prob,na.rm=na.rm))    
        }
        
  res   <- mapply(LESSTAT, x,  MoreArgs = list(prob=c(pcentral,pvar), wdim=2),SIMPLIFY=FALSE)
  res   <- mapply(LESSTAT, res,MoreArgs = list(prob=c(pcentral,punc), wdim=1),SIMPLIFY=FALSE)

  LESRATIOS <- function(x){
    if(is.list(x)) return(lapply(x,LESRATIOS))
    A <- x[1,1]
    B <- x[1,2]
    C <- x[2,1]
    D <- x[2,2]
    VariabilityR <- B/A
    UncertaintyR <- C/A
    OverallUR <- D/A
    return(c(A,B,C,D,VariabilityR,UncertaintyR,OverallUR))
    }

    res   <- lapply(res, LESRATIOS)
    res <- matrix(unlist(res),byrow=TRUE,ncol=7,
                  dimnames=list(unlist(nomsortie),c("A","B","C","D","VariabilityR","UncertaintyR","Over.Unc.R")))

  return(res)}

