#<<BEGIN>>
summary.mc <- function(object,probs = c(0,0.025,0.25,0.5,0.75,0.975,1),lim=c(0.025,0.975),...)
#TITLE Summary of mcnode and mc Object
#DESCRIPTION
# Provides a summary of a \samp{mcnode}, a \samp{mc} or a \samp{mccut} object.
#KEYWORDS univar
#INPUTS
#{object}<<a \samp{mcnode} or a \samp{mc} object or a \samp{mccut} object.>>
#{x}<<A \samp{summary.mc} object as provided by the \samp{summary.mc} function.>>
#[INPUTS]
#{probs}<<A vector of values used for the quantile function (variability dimension).>>
#{digits}<<Number of digits in the print.>>
#{lim}<<A vector of values used for the quantile function (uncertainty dimension).>>
#{\dots}<<For generic functions consistancy.>>
#VALUE
#a list.
#DETAILS
#The mean, the standard deviation and the \samp{probs} quantiles will be evaluated in the variability dimension.
#The median, the mean and the \samp{lim} quantiles will then be evaluated on these statistics in the uncertainty dimension.
#
#Multivariate nodes:
#
#If the \samp{"outm"} attributes of the mcnode is "none", the node is not evaluated, if it is "each"
#the variates are evaluated one by one, if it is a function (e.g. "mean"), the function is applied on the
#\samp{nvariates} dimension before providing a classical output.
#
#SEE ALSO
# \code{\link{mcnode}} for mcnode objects, \code{\link{mc}} for mc objects, \code{\link{mccut}} for mccut objects, \code{\link{quantile}}
#EXAMPLE
#data(total)
#summary(xVUM3)
#summary(total)

#CREATED 08-01-25
#REVISED 08-01-25
#--------------------------------------------
{
  yaprob <- length(probs) > 0
  if(yaprob) {nom1 <- paste(probs*100,"%",sep="")
  	nom1[probs==0] <- "Min"
  	nom1[probs==1] <- "Max"} else nom1 <- NULL
	nom1 <- c("mean","sd",nom1,"nsv","Na's")

  yalim <- length(lim) > 0
  if(yalim) nom2 <- paste(lim*100,"%",sep="") else nom2 <- NULL
  nom2 <- c("median","mean",nom2)

  warnna <- FALSE
  stat1 <- function(x){	# x is a array(nvar, nunc, 1)
    nsimna <- apply(x,2,function(x) sum(is.na(x)))
    if(yaprob) quant <- t(apply(x,2,quantile,na.rm=TRUE,probs=probs,names = FALSE))
      else quant <- NULL
    if(warnna==FALSE && any(nsimna!=0)) {
        warning("NA's value in mc object",call.=FALSE)
        eval(expression(warnna <- TRUE),sys.parent())}

        cbind(    colMeans(x,na.rm=TRUE),
                  apply(x,c(2,3),sd,na.rm=TRUE),	#compliance with deprecated sd(<matrix>)
                  quant,
                  apply(x,2,length),
                  nsimna)
     }

    stat2 <- function(output){
        output <- rbind( rowMeans(output,na.rm=TRUE),
                         apply(output,1,quantile,na.rm=TRUE,probs=c(0.5,lim),names = FALSE))
   	    output[c(1,2),] <- output[c(2,1),]
   	    return(output)
     }

  outm <- lapply(object,attr,which="outm")
  object <- object[outm != "none",drop=FALSE]
  if(length(object)==0) return(NULL)
  outm <- lapply(object,attr,which="outm")
  nom <- names(object)
  typen <- sapply(object,attr,which="type")
  dimm <- sapply(object,dim)
  nvariates <- dimm[3,]

    # Deal with multivariate nodes

  gerout <- function(node,name,outm,nvariates){
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

  object <- mapply(gerout,object,nom,outm,nvariates,SIMPLIFY=FALSE)

  # Function for stat

  lesstat <-function(object,typen){

        if(is.list(object)) return(mapply(lesstat,object,typen,SIMPLIFY=FALSE))

        if(typen=="VU") {
          out <- stat2(t(stat1(object[,,1,drop=FALSE])))
          nomtmp <- list(nom2,nom1)}

        else if(typen=="V") {
          out <- stat1(object[,,1,drop=FALSE])
          nomtmp <- list("NoUnc",nom1)}

        else if(typen=="U") {
          out <- stat2(object[,,1,drop=FALSE])
          nomtmp <- list(nom2,"NoVar")}

        else if(typen=="0") {
          out <- matrix(object[,,1])
          nomtmp <- list("NoVar","NoUnc")
        }
      dimnames(out) <- nomtmp
      return(out)
      }

  output <- mapply(lesstat,object,typen,SIMPLIFY=FALSE)
  output <- mapply("attr<-",output,"type",typen,SIMPLIFY=FALSE)

	class(output) <- c("summary.mc","listof")

  return(output)
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




#<<BEGIN>>
summary.mcnode <- function(object,probs = c(0,0.025,0.25,0.5,0.75,0.975,1),lim=c(0.025,0.975),digits = 3,...)
#ISALIAS summary.mc
#--------------------------------------------
{
  summary.mc(list(node=object), probs = probs, lim=lim, digits = 3,... )
}

#<<BEGIN>>
print.summary.mc <- function(x,digits=3,...)
#ISALIAS summary.mc
#--------------------------------------------
{
  x <- lapply(x,function(y) if(is.list(y)) lapply(y,"attr<-","type",NULL) else {attr(y,"type") <- NULL;y})
  NextMethod(x,digits=digits,...)
  }
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

