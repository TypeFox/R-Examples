#<<BEGIN>>
tornado <- function(mc,output = length(mc),use = "all.obs",	method=c("spearman","kendall","pearson"),lim=c(0.025,0.975))
#TITLE Computes Correlation between Inputs and Output in a mc Object (tornado) in the Variability Dimension;
#DESCRIPTION Provides statistics for a tornado chart. Evaluates correlations between output and inputs of a \samp{mc} object.
#KEYWORDS univar
#INPUTS
#{mc}<<a \code{\link{mc}} object or a \code{\link{mccut}} object.>>
#{x}<<A \samp{tornado} object as provided by the \samp{tornado} function.>>
#[INPUTS]
#{output}<<(for \samp{mc} objects only). The rank or the name of the output to be considered. By default: the last element of the \samp{mc}.>>
#{use}<<(for \samp{mc} objects only). An optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "all.obs", "complete.obs" or "pairwise.complete.obs" (see \code{\link{cor}}).>>
#{method}<<(for \samp{mc} objects only). A character string indicating which correlation coefficient (or covariance) is to be computed. One of "spearman" (default), "kendall" or "pearson", can be abbreviated (see \code{\link{cor}}). Warning : the default is not the same in \code{\link{cor}}.>>
#{lim}<<A vector of quantiles used to compute the credible interval in two-dimensional models.>>
#{\dots}<<Further arguments to be passed to the final print function.>>
#DETAILS
# The tornado function computes the spearman's rho statistic. It is used to estimate a rank-based measure of association between one set of 
#random variable of a \samp{mc} object (the output) and the others (the inputs).</>
#\samp{tornado} may be applied on a \samp{mccut} object if a \samp{tornado} function was used in the third block of the
#\code{\link{evalmccut}} call.
#VALUE
#An invisible object of class tornado.
#A tornado object is a list of objects containing the following objects:
#{value}<<the value of correlation coefficients>>
#{output}<<the name of the output>>
#{method}<<the method used>>
#{use}<<the use parameter>>
#DETAILS
#If "output" refers to a \samp{"0" mcnode}, it is an error.
#If "output" refers to a \samp{"V" mcnode}, correlations are only provided for other \samp{"V" mcnode}s.
#If "output" refers to a \samp{"U" mcnode}, correlations are only provided for other \samp{"U" mcnode}s.
#If "output" refers to a \samp{"VU" mcnode}, correlations are only provided for other \samp{"VU" mcnode}s and \samp{"V" mcnode}s.
#
#If use is "all.obs", then the presence of missing observations will produce an error.
#If use is "complete.obs" then missing values are handled by casewise deletion.
#Finally, if use has the value "pairwise.complete.obs" then the correlation between each pair of variables
#is computed using all complete pairs of observations on those variables.
#SEE ALSO
# \code{\link{cor}}.</>
# \code{\link{plot.tornado}} to draw the results.</>
#EXAMPLE
#data(total)
#tornado(total,2,"complete.obs","spearman",c(0.025,0.975))
#tornado(total,4,"pairwise.complete.obs","spearman",c(0.025,0.975))
#tornado(total,6,"complete.obs","kendall",c(0.025,0.975))
#tornado(total,8,"complete.obs","spearman",c(0.025,0.975))
#(y <- tornado(total,10,"complete.obs","spearman",c(0.025,0.975)))
#plot(y)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
  if(inherits(mc,"mccut")){
    summ <- function(x) {
      y <-c(mean=mean(x,na.rm=TRUE),quantile(x,probs=c(0.5,lim),na.rm=TRUE))
      y[1:2] <- y[2:1]
      return(y)}
    mc <- mc[[which(sapply(mc,inherits,what="tornado.mccut"))[1]]]
    if(length(mc)==0) stop("tornado was not evaluated in evalmccut : impossible to evaluate")
    mc$value <- lapply(mc$value,function(x) drop(apply(x,c(1,3),summ)))
    class(mc) <- "tornado"
    return(mc)
  }
  
	method <- match.arg(method)
 	na.method <- pmatch(use, lesmet <- c("all.obs", "complete.obs", "pairwise.complete.obs"))

  if(!inherits(mc,"mc")) stop("tornado is not implemented for this object")
  if(is.numeric(output)) output <- names(mc)[output]
	if(!(output %in% names(mc)))	stop("Output is not a valid value")

  typen <- sapply(mc,attr,which="type")
  outm <-  lapply(mc,attr,which="outm")
  typeout <- typen[output]
  nameout <- names(mc[[output]])
	if(typeout!="V" && typeout!="VU")	stop("Output is not a 'V' or a 'VU' node")

  #Select nodes according to the type of the output
  quel <- outm!="none" & (typen == "V" | (typeout=="VU" & typen=="VU"))
  if(sum(quel) < 2) stop("No valid pairs of mcnode")      # Verify latter if one is the output

  # Select the nodes
  mc <- mc[quel]
  outm <- outm[quel]
  typen <- typen[quel]
  nom <- names(mc)
  dimm <- sapply(mc,dim)
  nvariates <- dimm[3,]
  nom <- names(mc)
  nomi <- nom[nom != output]

  
  
  # Which row are complete
  quelk <- !apply(sapply(mc,function(x) apply(x,1,function(x) any(is.na(x)))),1,any)

  if(!all(quelk)){
    if(na.method==1) stop("NA values. Change 'use' option")
    if(na.method==2){
      mc <- lapply(mc,"[",quelk,,,drop=FALSE)    # ! perte structure
      use <-  "all.obs"}
  }

  lesnom <- function(outm,nvariates,nom){
    if(outm[1] == "each"){
      if(nvariates==1) nomsortie <- nom
      else nomsortie <- paste(nom,1:nvariates,sep=".")}
    else nomsortie <- paste(nom,": ",outm," of variates",sep="")
    return(nom=nomsortie)}

  nomsortie <- mapply(lesnom,outm,nvariates,nom,SIMPLIFY=FALSE)
  rnbinput <- length(nomsortie)

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

  mc <- mapply(gerout,mc,outm,nvariates,SIMPLIFY=FALSE)

  # Deal with outm for the output
  nomout <- nomsortie[[output]]
  out <- mc[[output]]
  nco <-  dimm[2,output]

  nomin <- unlist(nomsortie[nomi])
  mc <- mc[nomi]


  # Build correlations
  calcorr <- function(inp,out){
    nunci <- dim(inp)[2]
   sapply(1:nco,function(y) cor(out[,y,],inp[,ifelse(nunci==1,1,y),],method=method,use=use))
  }

# Alternative bien plus lente  
#  calcorr <- function(inp,out){
#    nunci <- dim(inp)[2]
#    res <- cor(out[,,1],inp[,,1],method=method,use=use)
#    if(nunci==1) return(as.vector(res))
#    return(diag(res))
#    }

  if(!is.list(out)) out <- list(out)
  res <- rapply(out,function(x) matrix(rapply(mc,calcorr,how="unlist",out=x),nrow=nco,dimnames=list(NULL,nomin)),how="replace")

  if(typeout=="VU") {
    res <- lapply(res,function(x) {
                  tmp <- rbind(mean=colMeans(x,na.rm=TRUE),
                                  apply(x,2,quantile,probs=c(0.5,lim),na.rm=TRUE))
                  tmp[1:2,] <- tmp[2:1,]
                  rownames(tmp)[1:2] <- c("median","mean")
                  return(tmp)})
    }
  else   res <- lapply(res,"rownames<-",value="corr")

  names(res) <- nomout

  
  tc <- list(value = res, output = output, method = method, use = use)
	class(tc) <- "tornado"
  return(tc)
	}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#<<BEGIN>>
print.tornado <- function(x, ...)
#ISALIAS tornado
#--------------------------------------------
{	tmethod <- c("Spearman's rho statistic","Kendall's tau statistic","Pearson correlation")
	tmethod <- tmethod[x$method==c("spearman","kendall","pearson")]
  cat(tmethod,"\n")
	cat("Output: ",x$output,"\n")
	print(x$value,...)
 }
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


