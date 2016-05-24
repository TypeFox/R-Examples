#<<BEGIN>>
quantile.mc <- function(x, probs = seq(0, 1, 0.01), lim=c(0.025,0.975), na.rm=TRUE,  ...)
#TITLE Quantiles of a mc Object
#DESCRIPTION
# Evaluates quantiles of a \samp{mc} object. This function is used by \samp{plot.mc}
#KEYWORDS univar
#INPUTS
# {x}<<a \samp{mc} objects>>
#[INPUTS]
#{probs}<<the quantiles to be calculated>>
#{na.rm}<<TRUE or FALSE>>
#{lim}<<a vector of numbers (between 0 and 1) indicating the enveloppe. Maybe \samp{NULL} or empty.>>
#{\dots}<<For generic method consistancy.>>
#DETAILS
#The quantiles are evaluated in the variability dimension.
#Then, the median, the mean and the \samp{lim} quantiles are evaluated for each of these quantiles.
#VALUE
#A list of quantiles.
#SEE ALSO
#\code{\link{plot.mc}}, \code{\link{quantile}}.
#EXAMPLE
#data(total)
#quantile(total$xVUM3)
#quantile(total)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#

{
  lprobs <- length(probs)

  nomprob <- paste(probs*100,"%",sep="")
  nomlim <-  paste(lim*100,"%",sep="")
  nomlim <- c("median","mean",nomlim)
  lim    <- c(0.5,lim)

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

  LESSTAT <- function(x,typen){
    if(is.list(x)) return(lapply(x,LESSTAT,typen=typen))
        if(typen=="0") {
            restmp <- matrix(x[,,1,drop=FALSE],dimnames=list("NoInc","NoVar"))
            }

        else if(typen=="V"){
            restmp <- quantile(x[,,1], probs = probs, na.rm = na.rm, names=FALSE)
            restmp <- matrix(restmp,nrow=1,dimnames=list("NoInc",nomprob))
            }

        else if(typen=="U") {
            restmp <- quantile(x[,,1],probs = lim, na.rm = na.rm, names=FALSE)
            restmp <- c(mean(x[,,1],na.rm=na.rm),restmp)
            restmp[1:2] <- restmp[2:1]
            restmp <- matrix(restmp,ncol=1,dimnames=list(nomlim,"NoVar"))
            }

        else if(typen=="VU") {
            prem <- apply(x[,,1,drop=FALSE],2,quantile,probs = probs, na.rm = na.rm, names=FALSE)
            restmp <-  apply(prem,1,quantile,probs = lim, na.rm = na.rm, names=FALSE)
            restmp <-  rbind(rowMeans(prem,na.rm=na.rm), restmp)
            restmp[1:2,] <- restmp[2:1,]
            dimnames(restmp) <- list(nomlim,nomprob)
         }
#         attr(restmp,which="type") <- typen
         return(restmp)
      }

  res <- mapply(LESSTAT,x,typen,SIMPLIFY=FALSE)
  res <- mapply("attr<-",res,"type",typen,SIMPLIFY=FALSE)

  names(res) <- nom
  class(x) <- "plotmc"
  return(res)}

#<<BEGIN>>
quantile.mcnode <- function(x, ...)
#ISALIAS quantile.mc
#--------------------------------------------
{ nom <- deparse(substitute(x))
  x <- list(x)
  names(x) <- nom
  quantile.mc(x, ...)}