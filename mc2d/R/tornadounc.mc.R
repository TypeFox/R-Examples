#<<BEGIN>>
tornadounc <- function(mc,...) 
#ISALIAS tornadounc.mc
#--------------------------------------------
 UseMethod("tornadounc")
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
tornadounc.default <- function(mc,...) 
#ISALIAS tornadounc.mc
#--------------------------------------------
 tornadounc.mc(mc,...)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
tornadounc.mc <- function(mc,output = length(mc), quant=c(0.5,0.75,0.975),use = "all.obs",	method=c("spearman","kendall","pearson"),...)
#TITLE Computes Correlation between Inputs and Output in a mc Object (tornado) in the Uncertainty Dimension
#KEYWORDS univar
#NAME tornadounc
#DESCRIPTION Provides statistics for a tornado chart. Evaluates correlations between output and inputs of a \samp{mc} object in the uncertainty dimension.
#INPUTS
#{mc}<<a \samp{mc} object.>>
#{x}<<a \samp{tornadounc} object.>>
#[INPUTS]
#{output}<<The rank or the name of the output to be considered. Should be a \samp{"VU"} or a \samp{"U" type mcnode}. By default: the last element of \samp{mc}.>>
#{quant}<<The vector of quantiles used in the variability dimension.>>
#{use}<<An optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "all.obs", "complete.obs" or "pairwise.complete.obs" (see \code{\link{cor}}).>>
#{method}<<A character string indicating which correlation coefficient (or covariance) is to be computed. One of "spearman" (default), "kendall" or "pearson", can be abbreviated (see \code{\link{cor}}). Warning : "pearson" is the default for \code{\link{cor}}).>>
#{\dots}<<Further arguments to be passed to the final print function.>>
#DETAILS
# The \samp{tornadounc.mc} function computes the spearman's rho statistic between 
#{*}<<values (\samp{"U" type mcnode}) or statistics calculated in the variability dimension (\samp{"VU" type mcnode}) of inputs and>>
#{*}<<values (\samp{"U" type mcnode}) or statistics calculated in the variability dimension (\samp{"VU" type mcnode}) of one output.>>
#The statistics are the mean, the median and the quantiles specified by \samp{quant}.
#
#It is useful to estimate a rank-based measure of association between one set of
# random variable of a \samp{mc} object (the output) and the others in the uncertainty dimension.</>
#\samp{tornadounc.mccut} may be applied on a \code{\link{mccut}} object if a \samp{summary.mc} function was used in the third block of the
#\code{\link{evalmccut}} call.
#
#If output refers to a \samp{"0"} or \samp{"V" mcnode}, it is an error.
#
#If use is "all.obs", then the presence of missing observations will produce an error.
#If use is "complete.obs" then missing values are handled by casewise deletion.
#Finally, if use has the value "pairwise.complete.obs" then the correlation between each pair of variables
#is computed using all complete pairs of observations on those variables.
#VALUE
#An invisible object of class \samp{tornadounc}.
#A \samp{tornadounc} object is a list of objects containing the following objects:
#{value}<<a matrix of values of correlation coefficients. Each row are the value
#or the statistics of inputs, each columns the value or the statistics of outputs.>>
#{output}<<the name of the output>>
#{method}<<the method used>>
#{use}<<the \samp{use} parameter>>
#SEE ALSO
# \code{\link{cor}}.</>
# \code{\link{tornado}} for tornado in the variability dimension.</>
# \code{\link{plot.tornadounc}} to draw the results.</>
#EXAMPLE
#data(total)
#tornadounc(total,3)
#tornadounc(total,4,use="complete")
#tornadounc(total,7,use="complete.obs")
#tornadounc(total,8,use="complete.obs")
#(y <- tornadounc(total,10,use="complete.obs"))
#plot(y,1,1)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
	method <- match.arg(method)
 	na.method <- pmatch(use, lesmet <- c("all.obs", "complete.obs", "pairwise.complete.obs"))

  if(length(output) > 1) stop("Only one output permitted")
  if(is.numeric(output)) output <- names(mc)[output]
	if(!(output %in% names(mc)))	stop("Output is not a valid value")

  typen <- sapply(mc,attr,which="type")
  outm <-  lapply(mc,attr,which="outm")
  typeout <- typen[output]
	if(typeout!="U" && typeout!="VU")	stop("Output is not a 'U' or a 'VU' node")
  if(outm[output]=="none") stop("Output has a 'none' outm attribute")

  #Select nodes according to the type of the output
  quel <- outm!="none" & (typen == "U" | (typeout=="VU" & typen=="VU"))
  if(sum(quel) < 2) stop("No valid pairs of mcnode")

  # Select the nodes
  mc <- mc[quel]
  outm <- outm[quel]
  typen <- typen[quel]
  dimm <- sapply(mc,dim)
  nvariates <- dimm[3,]
  nom <- names(mc)
  nomi <- nom[nom != output]

   # Which col are complete
  quelk <- !apply(sapply(mc,function(x) apply(x,2,function(x) any(is.na(x)))),1,any)
  nco <- sum(quelk)              # remaining dimension of uncertainty

  if(!all(quelk)){
    if(na.method==1) stop("NA values. Change 'use' option")
    if(na.method==2){
      mc <- lapply(mc,"[",,quelk,,drop=FALSE)    # ! perte structure
      use <-  "all.obs"}
  }

  # Real name of input / output
  
  yaprob <- length(quant) > 0
  if(yaprob) {nom1 <- paste(quant*100,"%",sep="")
  	nom1[quant==0] <- "Min"
  	nom1[quant==1] <- "Max"} else nom1 <- NULL
  nom1 <- c("mean","sd",nom1)

  lesnom <- function(nom,outm,nvariates,typen){
    if(outm[1] == "each"){
      if(nvariates==1) nomsortie <- nom
      else nomsortie <- paste(nom,1:nvariates,sep=".")}
    else nomsortie <- paste(nom,": ",outm," of variates",sep="")

    if(typen == "VU")
      nomsortie <- lapply(nomsortie, function(x) paste(nom1,x))

    return(nom=nomsortie)}

  nomsortie <- mapply(lesnom,nom,outm,nvariates,typen,SIMPLIFY=FALSE)
  nomlistfin <- lesnom(output,outm[[output]],nvariates[[output]],"U")           # Name final list

  # Deal with multivariate nodes

  gerout <- function(node,name,outm,nvariates){
    if(outm[1]=="each")
      res <- lapply(1:nvariates,function(x) node[,,x,drop=FALSE])
    else
        res <- mapply(function(outm){
                  func <- get(outm,mode="function")
                  node <- apply(node,c(1,2),func)
                  dim(node) <- c(dim(node),1)
                  node},outm,SIMPLIFY=FALSE)

  return(res)
  }

  mc <- mapply(gerout,mc,nom,outm,nvariates,SIMPLIFY=FALSE)

  # Build Statistics on inputs and output

  stat <- function(node,outm,typen,nomsortie){
    if(typen=="U")
      tmp <- lapply(node,function(x) x[1,,])
    else {
      tmp <- lapply(node,function(x)
                              apply(x,c(2,3),function(x) c(mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE), # Compliance with the deprecation of sd(<matrix>)
                                              quantile(x,na.rm=TRUE,prob=quant))))
                                              }
      node <- lapply(tmp,matrix,nrow=nco,byrow=TRUE)

  return(node)
  }

  mc <- mapply(stat,mc,outm,typen,nomsortie,SIMPLIFY=FALSE)

  nomout <- nomsortie[[output]]
  out <- mc[[output]]

  nomin <- unlist(nomsortie[nomi])
  mc <- matrix(unlist(mc[nomi]),nrow=nco)
  
  lescorr <- mapply(function(x) as.matrix(cor(x,mc,method=method,use=use)),out,SIMPLIFY=FALSE)
  lescorr <- lapply(lescorr,"colnames<-",value=nomin)
  lescorr <- mapply("rownames<-",lescorr,value=nomout,SIMPLIFY=FALSE,USE.NAMES=TRUE)
  names(lescorr) <- nomlistfin

  tc <- list(value = lescorr, output = output, method = method, use = use)
	class(tc) <- "tornadounc"
  return(tc)
	}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
print.tornadounc <- function(x, ...)
#ISALIAS tornadounc.mc
#--------------------------------------------
{
tmethod <- c("Spearman's rho statistic","Kendall's tau statistic","Pearson correlation")
	tmethod <- tmethod[x$method==c("spearman","kendall","pearson")]
  cat("Tornado on uncertainty\n")
  cat(tmethod,"\n")
	cat("Output: ",x$output,"\n")
	print(x$value,...)
 }
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

