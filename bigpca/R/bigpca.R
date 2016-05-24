###NAMESPACE ADDITIONS###
# Depends: R (>= 3.0), grDevices, graphics, stats, utils, reader (>= 1.0.1), NCmisc (>= 1.1), bigmemory (>= 4.0.0), biganalytics
# Imports: parallel, methods, bigmemory.sri, irlba
# Suggests:
# importFrom(parallel, mclapply)
# importFrom(irlba, irlba)
# import(methods, bigmemory.sri, grDevices, graphics, stats, utils, reader, NCmisc, bigmemory, biganalytics)
###END NAMESPACE###


.onLoad <- function(libname, pkgname) {
  options(deleteFileBacked=TRUE) # generally will be nice to set this to TRUE to avoid be prevented from overwriting files
}

# may want to use updated irlba2() function in order to allow bigger Matrices to be pca-ed #

#setMethod("print", signature=signature(x="big.matrix"), function(x,...) prv.big.matrix(x,...))

# internal function for testing: 
# use: keeper <- manage.test.files() at start of example
# then: unlink(bigpca:::manage.test.files(FALSE,keepers)) at end of examples
# then inspect the example output for warnings with unlink()
manage.test.files <- function(start=TRUE,keepers=NULL) {
  test.files <- c("testMyBig.bck","t.bigMat.bck","splitmatR1.bck","split1.bck","sel.bck","functestdn.bck","functest.bck","funclongcol.bck","t.bigMat.RData",
    "sel.RData","functestdn.RData","fn.RData","testMyBig.dsc","t.bigMat.dsc","splitmatR1.dsc","split1.dsc","sel.dsc","functestdn_file_rowname_list_check_this.txt",
    "functestdn.dsc","functest.dsc","funclongcol.dsc","bigrowstemp.txt","bigcolstemp.txt","test.dsc","test.bck")
  if(start) { 
    ii <- which(list.files() %in% test.files) 
    if(length(ii)>0) { keepers <- list.files()[ii] }
    return(keepers)
  } else {
    ii <- which(list.files() %in% test.files[!test.files %in% keepers]) 
    if(length(ii)>0) { 
      killers <- list.files()[ii] 
      warning("unlink(",paste(killers,collapse=","),")")
    } else { killers <- NULL }
  }
  return(killers)
}




#' Tidier display function for big matrix objects
#' 
#' This function prints the first and last columns and rows of a big matrix, and
#' a few more than this if desired. Allows previewing of a big.matrix without 
#' overloading the console. 
#'
#' @param bigMat the description file, big.matrix object, or big.matrix.descriptor object, 
#'  anything that can be read by get.big.matrix()
#' @param dir the directory containing the big.matrix backing/description files
#' @param name logical, whether to print a name for the matrix
#' @param dat logical, whether to print any of the matrix contents (overrides row/col)
#' @param descr character, optional name of the description file, which if not null will be displayed
#' @param bck character, optional name of the backing file, which if not null will be displayed
#' @param mem logical, whether to display the amount of memory used by the object
#' @param rows integer, number of rows to display
#' @param cols integer, number of columns to display
#' @param rcap character, caption to display for the rows
#' @param ccap character, caption to display for the columns
#' @param ... additional arguments to prv.large (from NCmisc) which displays the end result
#' @return Prints to console a compact representation of the bigMat matrix, 
#'  with the first few rows and columns, and the last row and column. Note that sometimes
#'  the initial printing of a big.matrix can take a little while. But subsequently the printout
#'  should be almost instantaneous.
#' @seealso \code{\link{get.big.matrix}}
#' @export
#' @examples 
#' bM <- filebacked.big.matrix(20, 50,
#'        dimnames = list(paste("r",1:20,sep=""), paste("c",1:50,sep="")),
#'        backingfile = "test.bck",  backingpath = getwd(), descriptorfile = "test.dsc")
#' bM[1:20,] <- replicate(50,rnorm(20))
#' prv.big.matrix(bM)
#' prv.big.matrix(bM,rows=10,cols=4)
#' unlink(c("test.dsc","test.bck"))  # clean up files
prv.big.matrix <- function(bigMat,dir="",rows=3,cols=2,name=NULL,dat=TRUE,
                           descr=NULL,bck=NULL,mem=FALSE,rcap="",ccap="",...) {
  # print summary of big matrix contents
  #must.use.package("bigmemory")
  if(is.null(name) & is.character(bigMat)) { name <- basename(bigMat[1]) } # if it's a file name
  bigMat <- get.big.matrix(bigMat,dir=dir)
  if(!is.null(name)) { 
    cat("\nBig matrix; '",name,"', with: ",sep="") 
  } else { cat("Big matrix with: ") }
  nC <- ncol(bigMat); nR <- nrow(bigMat)
  bcap <- paste(" (",c(rcap[1],ccap[1]),")",sep="") 
  if(rcap!="") { rcap <- bcap[1] }; if(ccap!="") { ccap <- bcap[2] }
  cat(nR," rows",rcap,", ",nC," columns",ccap,"\n",sep="") 
  if(is.sub.big.matrix(bigMat)) { cat("[a sub.big.matrix object]\n")}
  cat(" - data type:",is(bigMat[1,1])[1],"\n")
  if(!is.null(descr)) { cat(" - descriptor file;",descr,"\n") }
  if(is.filebacked(bigMat)) { 
    if(!is.null(bck)) {
      cat(" - backing file;",bck,"\n") }
  } else {
    cat(" - not filebacked")
  }
  if(dat) {
    prv.large(bigMat,rows=rows,cols=cols,...)
  } else {
    if(!is.null(colnames(bigMat))) {
      cat(" - columns:",paste(colnames(bigMat)[1:max(rows,cols)],collapse=", "),
          "...",colnames(bigMat)[nC],"\n")
    } else { cat(" - no column names\n") }
    if(!is.null(rownames(bigMat))) {
      cat(" -    rows:",paste(rownames(bigMat)[1:max(rows,cols)],collapse=", "),
          "...",rownames(bigMat)[nR],"\n")  
    } else { cat(" - no row names\n") }
  }
  if(mem) {
    total.datapoints <- nR*nC
    disk.est <- round(estimate.memory(bigMat))
    cat("Total of",total.datapoints,"data-points, using",disk.est,"GB estimated disk space\n")
  }
  cat("\n")
}




#' Estimate the variance percentages for uncalculated eigenvalues
#'
#' If using a function like irlba() to calculate PCA, then you can choose (for speed) 
#' to only calculate a subset of the eigenvalues. So there is no exact percentage of variance explained 
#' by the PCA, or by each component as you will get as output from other routines.
#' This code uses a linear, or b*1/x model, to estimate the AUC for the unknown eigenvalues, providing
#' a reasonable estimate of the variances accounted for by each unknown eigenvalue, and
#' the predicted eigenvalue sum of the unknown eigenvalues.
#'
#' @param eigenv the vector of eigenvalues actually calculated
#' @param min.dim the size of the smaller dimension of the matrix submitted to singular
#'  value decomposition, e.g, number of samples - i.e, the max number of possible eigenvalues,
#'  alternatively use 'M'.
#' @param M optional enter the original dataset 'M'; simply used to derive the dimensions,
#'  alternatively use 'min.dim'.
#' @param elbow the number of components which you think explain the important portion
#'  of the variance of the dataset, so further components are assumed to be reflecting
#'  noise or very subtle effects, e.g, often the number of components used is decided
#'  by the 'elbow' in  a scree plot (see 'pca.scree.plot')
#' @param linear whether to use a linear model to model the 'noise' eigenvalues; alternative
#'  is a 1/x model with no intercept.
#' @param estimated logical, whether to return the estimated variance percentages for unobserved eigenvalues
#'  along with the real data; will also generate a factor describing which values in the returned
#'  vector are observed versus estimated.
#' @param print.est whether to output the estimate result to the console
#' @param print.coef whether to output the estimate regression coefficients to the console
#' @param add.fit.line logical, if there is an existing scree plot, adds the fit line from this estimate
#'  to the plot ('pca.scree.plot' can use this option using the parameter of the same name)
#' @param col colour for the fit line
#' @param ignore.warn ignore warnings when an estimate is not required (i.e, all eigenvalues present)
#' @return By default returns a list where the first element ''variance.pcs' are the known variance
#'  percentages for each eigenvalue based on the estimated divisor, the second element
#'  'tail.auc' is the area under the curve for the estimated eigenvalues. If estimate
#'  =TRUE then a third element is return with separate variance percentages for
#'  each of the estimated eigenvalues.
#' @seealso \code{\link{pca.scree.plot}}
#' @export
#' @examples
#' nsamp <- 100; nvar <- 300; subset.size <- 25; elbow <- 6
#' mat <- matrix(rnorm(nsamp*nvar),ncol=nsamp) 
#' # or use: # mat <- crimtab-rowMeans(crimtab) ; subset.size <- 10 # crimtab centred
#' prv.large(mat)
#' pca <- svd(mat,nv=subset.size,nu=0) # calculates subset of V, but all D
#' require(irlba)
#' pca2 <- irlba(mat,nv=subset.size,nu=0) # calculates subset of V & D
#' pca3 <- princomp(mat,cor=TRUE) # calculates all
#' # number of eigenvalues for svd is the smaller dimension of the matrix
#' eig.varpc <- estimate.eig.vpcs(pca$d^2,M=mat)$variance.pcs
#' cat("sum of all eigenvalue-variances=",sum(eig.varpc),"\n")
#' print(eig.varpc[1:elbow])
#' # number of eigenvalues for irlba is the size of the subset if < min(dim(M))
#' eig.varpc <- estimate.eig.vpcs((pca2$d^2)[1:subset.size],M=mat)$variance.pcs
#' print(eig.varpc[1:elbow])  ## using 1/x model, underestimates total variance
#' eig.varpc <- estimate.eig.vpcs((pca2$d^2)[1:subset.size],M=mat,linear=TRUE)$variance.pcs
#' print(eig.varpc[1:elbow])  ## using linear model, closer to exact answer
#' eig.varpc <- estimate.eig.vpcs((pca3$sdev^2),M=mat)$variance.pcs
#' print(eig.varpc[1:elbow])  ## different analysis, but fairly similar var.pcs
estimate.eig.vpcs <- function(eigenv=NULL,min.dim=length(eigenv),M=NULL,elbow=NA,linear=TRUE,
                              estimated=FALSE,print.est=TRUE,print.coef=FALSE,
                              add.fit.line=FALSE,col="blue",ignore.warn=FALSE) {
  ## if matrix is optionally inputted, calculate the minimum dim automatically
  if(length(eigenv)<2) { warning("was passed insufficient eigenvalues (",length(eigenv),")") ; return(eigenv) }
  if(!is.null(M)) { if(!is.null(dim(M))) { min.dim <- min(dim(M),na.rm=T) } }
  if(all(is.na(min.dim))) { min.dim <- length(eigenv) }
  n.comp <- length(eigenv) # max(c(min.dim,length(eigenv)),na.rm=T)
  if(all(is.na(elbow))) { 
    if(n.comp==min.dim) {
      elbow <- quick.elbow(eigenv)
    } else {
      elbow <- 3 
    }
  }
  elbow <- round(min(n.comp,elbow,na.rm=T)) # make sure not > n.comp
  if(!is.numeric(eigenv)) { warning("eigenv not numeric"); return(NULL) }
  #preview(c("elbow","min.dim","n.comp","eigenv"))
  if(is.na(min.dim) | ((min.dim-n.comp)<2) | ((n.comp-elbow)<(min(20,min.dim/20,na.rm=T))) ) {
    # if most/all eigenvalues already present this is not needed, or if parameters insufficient
    # then don't try to calculate the AUC of the remaining eigenvalues
    if(n.comp==min.dim) {
      if(!ignore.warn) { cat("All eigenvalues present, estimate not required\n") }
      estimated <- F
    } else {
      warning("didn't attempt to estimate eigenvalues as there were",
        " very few unknowns compared to the number of samples,",
        " or not enough eigenvalues between the elbow and 'min.dim'")
    }
    var.pcs <- eigenv[1:n.comp]/(sum(eigenv)); tail.var <- 0
  } else {
    # estimate combined variance of eigenvalues not calculated by irlba using 1/x model
    if(!linear) {
      xx <- 1/(1:length(eigenv))
      ab <- lm(eigenv[elbow:n.comp]~0+xx[elbow:n.comp])$coef
      tail.var <- ((log(min.dim)-log(n.comp))*ab[1]) # integral evaluated
      mod.txt <- "[b/x, no intercept]"
      predy <- ab[1]*(1/c((elbow+1):min.dim))
    } else {
      xx <- 1:length(eigenv)
      ab <- lm(eigenv[elbow:n.comp]~xx[elbow:n.comp])$coef
      zeropoint <- round(ab[1]/abs(ab[2])); zeropoint <- min(c(min.dim,zeropoint),na.rm=T)
      tail.var <- (zeropoint-n.comp)*(ab[1]+((n.comp)*ab[2]))*.5
     # tail.var <- ((ab[1]*(zeropoint-n.comp))+((((n.comp-elbow)^2)-((zeropoint-elbow)^2))*ab[2]*.5)) # integral evaluated
      mod.txt <- "[a + bx]"
      predy <- ab[1]+(ab[2]*c((elbow+1):min.dim))
      predy[(zeropoint-elbow):(min.dim-elbow)] <- 0
    }
    # intercept ignored as the asymptote should theoretically be zero so assumption
    # values > this reflect noise variance that might dissipate as x--> n.samp
    if(print.est) {
      not.calc <- min.dim-length(eigenv)
      cat(" estimate of eigenvalue sum of",not.calc,"uncalculated eigenvalues:",(as.numeric(tail.var)),"\n")
    }
    if(print.coef) {
      cat(" slope",mod.txt,":",as.numeric(tail(ab,1)),"\n")
    }
    if(add.fit.line) {
      # add fitted line to scree plot if one exists
      predx <- c((elbow+1):min.dim); print(predy)
      #prv(c("predx","predy"))
      print(length(predx)); print(length(predy))
      try(lines(predx,predy,col=col),T)
    }
    var.pcs <- eigenv[1:n.comp]/(sum(eigenv)+tail.var)
  }
  if(estimated) {
    var.pcs.comb <- numeric(min.dim); realornot <- rep("estimated",times=min.dim)
    var.pcs.comb[(1+min.dim-length(predy)):min.dim] <- predy/(sum(eigenv)+tail.var)
    var.pcs.comb[1:length(var.pcs)] <- var.pcs
    realornot[1:length(var.pcs)] <- "observed"
    out <- list(var.pcs.comb,tail.var,as.factor(realornot))
    names(out) <- c("variance.pcs","tail.auc","estimated")
  } else {
    out <- list(var.pcs,tail.var)
    names(out) <- c("variance.pcs","tail.auc")
  }
  return(out)
}


#' Make scree plots for any PCA
#'
#' Make a scree plot using eigenvalues from princomp(), prcomp(), svd(), irlba(), big.PCA(), etc.
#' Note that most these return values which need to be squared to be proper eigenvalues.
#' There is also an option to use the estimate.eig.vpcs() function to estimate any missing
#' eigenvalues (e.g, if using a function like irlba' to calculate PCA) and then to visualise
#' the fitline of the estimate on the scree plot.
#'
#' @param eigenv the vector of eigenvalues actually calculated
#' @param elbow the number of components which you think explain the important chunk
#'  of the variance of the dataset, so further components are modelled as reflecting
#'  noise or very subtle effects, e.g, often the number of components used is decided
#'  by the 'elbow' in  a scree plot (see 'pca.scree.plot')
#' @param min.dim the size of the smaller dimension of the matrix submitted to singular
#'  value decomposition, e.g, number of samples - i.e, the max number of possible eigenvalues,
#'  alternatively use 'M'.
#' @param M optional enter the original dataset 'M'; simply used to derive the dimensions,
#'  alternatively use 'min.dim'.
#' @param linear whether to use a linear model to model the 'noise' eigenvalues; alternative
#'  is a 1/x model with no intercept.
#' @param printvar logical, whether to print summary of variance calculations
#' @param add.fit.line logical, if there is an existing scree plot, adds the fit line from this estimate
#'  to the plot ('pca.scree.plot' can use this option using the parameter of the same name)
#' @param n.xax number of components to include on the x-axis
#' @param verbose logical, whether to display additional output
#' @param return.data logical, whether to return the percentages of variance explained for each component, or nothing (just plot)
#' @param ... further arguments to the plot function
#' @return Either a vector of variance percentages explained, or nothing (just a plot), depending on value of 'return.data'
#' @seealso \code{\link{pca.scree.plot}}
#' @export
#' @examples
#' require(irlba)
#' nsamp <- 100; nvar <- 300; subset.size <- 25; elbow <- 6
#' mat <- matrix(rnorm(nsamp*nvar),ncol=nsamp) 
#' #this gives the full solution
#' pca <- svd(mat,nv=subset.size,nu=0)
#  # test with larger and smaller subset, larger gives 1/x better fit, smaller, x
#' pca2 <- irlba(mat,nv=subset.size,nu=0)
#' # show alternate fits for linear versus 1/x fit
#' pca.scree.plot((pca2$d^2)[1:subset.size],n.xax=100,add.fit.line=TRUE,
#'                min.dim=min(dim(mat)),linear=TRUE, elbow=6, ylim=c(0,1400))
#' pca.scree.plot((pca2$d^2)[1:subset.size],n.xax=100,add.fit.line=TRUE,
#'               min.dim=min(dim(mat)),linear=FALSE, elbow=40, ylim=c(0,1400))
#' subset.size <- 75
#' pca2 <- irlba(mat,nv=subset.size,nu=0)
#' pca.scree.plot((pca2$d^2)[1:subset.size],n.xax=100,add.fit.line=TRUE,
#'               min.dim=min(dim(mat)),linear=TRUE, elbow=6, ylim=c(0,1400))
#' pca.scree.plot((pca2$d^2)[1:subset.size],n.xax=100,add.fit.line=TRUE,
#'               min.dim=min(dim(mat)),linear=FALSE, elbow=40, ylim=c(0,1400))
pca.scree.plot <- function(eigenv,elbow=NA,printvar=TRUE,min.dim=NA,M=NULL,add.fit.line=FALSE,
                           n.xax=max(30,length(eigenv)),linear=TRUE,verbose=FALSE,return.data=FALSE,...) 
{
  # do SCREE PLOTS AND calculate EIGENVALUE VARIANCE after a PCA
  if(!is.null(M)) { if(!is.null(dim(M))) { min.dim <- min(dim(M),na.rm=T) } }
  n.comp <- length(eigenv)
  if(all(is.na(elbow))) { 
    if(n.comp==min.dim) {
      elbow <- quick.elbow(eigenv)
    } else {
      elbow <- 3 
    }
  }
  elbow <- round(min(n.comp,elbow,na.rm=T))
  n.xax <- min(n.xax,n.comp) # constrain n.xax to not be greater than number of components
  plot(eigenv[1:n.xax],bty="l",xlab="number of principle components",
       ylab="eigenvalues",bg="green",pch=21,...)
  abline(v=(elbow+.5),lty="dashed")
  legend("topright",legend=c("Principle components","scree plot 'elbow' cutoff"),
         pt.bg=c("green",NA),pch=c(21,NA),lty=c(NA,"dashed"),bty="n")
  scree.calc <- estimate.eig.vpcs(eigenv=eigenv,min.dim=min.dim,elbow=elbow,ignore.warn=!verbose,
                  print.est=T,print.coef=T,add.fit.line=add.fit.line,col="blue",linear=linear)
  if(printvar & verbose) {
    cat(" sum of eigen-variance:",round(sum(eigenv)+scree.calc$tail.auc,2),"\n")
    cat(" variance % estimates: \n ",round(scree.calc$variance.pcs,2),"\n")
  }
  if(return.data) {
    return(scree.calc$variance.pcs)
  } else {
    return(invisible())
  }
}


# quickly choose an elbow for a PC. 
# at variance below 5% per component, choose the largest % drop
# designed for variance percentages, but will also work given a full set of Evalues
#' Quickly estimate the 'elbow' of a scree plot (PCA)
#' 
#' This function uses a rough algorithm to estimate a sensible 'elbow' to
#' choose for a PCA scree plot of eigenvalues. The function looks at an initial arbitrarily 'low'
#' level of variance and looks for the first eigenvalue lower than this. If the very first eigenvalue
#' is actually lower than this (i.e, when the PCs are not very explanatory) then this 'low' value is
#' iteratively halved until this is no longer the case. After starting below this arbitrary threshold
#' the drop in variance explained by each pair of consecutive PCs is standardized by dividing over the 
#' larger of the pair. The largest percentage drop in the series below 'low' % is selected as the 'elbow'.
#' @param varpc numeric, vector of eigenvalues, or 'percentage of variance' explained datapoints for
#'  each principle component. If only using a partial set of components, should first pass to 
#'  estimate.eig.vpcs() to estimate any missing eigenvalues.
#' @param low numeric, between zero and one, the threshold to define that a principle component
#'  does not explain much 'of the variance'.
#' @param max.pc maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n')
#' @return The number of last principle component to keep, prior to the determined elbow cutoff
#' @export
#' @seealso \code{\link{estimate.eig.vpcs}}
#' @author Nicholas Cooper 
#' @examples
#' # correlated data
#' mat <- sim.cor(100,50)
#' result <- princomp(mat)
#' eig <- result$sdev^2
#' elb.a <- quick.elbow(eig)
#' pca.scree.plot(eig,elbow=elb.a,M=mat) 
#' elb.b <- quick.elbow(eig,low=.05) # decrease 'low' to select more components
#' pca.scree.plot(eig,elbow=elb.b,M=mat) 
#' # random (largely independent) data, usually higher elbow #
#' mat2 <- generate.test.matrix(5,3)
#' result2 <- princomp(mat2)
#' eig2 <- result2$sdev^2
#' elb2 <- quick.elbow(result2$sdev^2)
#' pca.scree.plot(eig2,elbow=elb2,M=mat2)
quick.elbow <- function(varpc,low=.08,max.pc=.9) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones 
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      #print(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else { 
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("no eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee) 
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
  	elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}




#' A multicore 'apply' function for big.matrix objects
#'
#' # to put into NCmisc
#' Multicore method to run a function for a big.matrix that could be run using 'apply'
#' on a regular matrix (when parameter use.apply=T [default]). Otherwise for a
#' function that might be more efficient done in done in chunks (e.g, utilising vectorised 
#'  functions) use.apply=F can be set so that processing is done on larger submatrices, rather
#' than 1 row/column at a time. Input to specify whether to perform the function
#' row or columnwise is equivalent to 'apply' syntax, 1=by-rows, 2=by-columns.
#' This function is useful for big.matrix processing even without multiple cores, particulary
#' when MARGIN=1 (row-wise). While native colmean, colmin and colsd functions for 
#' big.matrix objects are very fast (and will probably outperform bmcapply even with 1
#'  core versus many), these are only natively implemented for column-wise operations and 
#' the equivalent operations if needing to be row-wise should be faster with bmcapply for
#' matrices larger than available RAM.
#' Can also be used for regular matrices although there is unlikely to be a speed advantage.
#' @seealso \code{\link{get.big.matrix}}
#' @param bigMat the big.matrix object to apply the function upon, can enter as a filename,
#'  description object or any other valid parameter to get.big.matrix(). Can also use with a standard matrix
#' @param MARGIN 1=row-wise, 2=column-wise, see same argument for base:::apply()
#' @param FUN the function to apply, should return a result with 1 dimension that has the
#'  same length as dim(bigMat)[MARGIN]=L; i.e, a vector length L, matrix (L,x) or (x,L) or list[[L]].
#'  Note that using a custom 'combine.fn' parameter might allow exceptions to this.
#' @param dir directory argument for get.big.matrix(), ie. the location of the bigMat backing file if 
#'  not in the current working directory.
#' @param by integer, the number of rows/columns to process at once. The default should work in most
#'  situations however, if the dimension not specified by MARGIN is very large, this might need
#'  to be smaller, or if the function being applied is much more efficient performed 
#'  on a large matrix than several smaller ones then this 'by' parameter should be increased
#'  within memory contraints. You should make sure 'estimate.memory(c(by,dim(bigMat)[-MARGIN]))'
#'  doesn't exceed available RAM.
#' @param n.cores integer, the number of parallel cores to utilise; note that sometimes if a machine
#'  has only a few cores this can result in slower performance by tying up resources
#'  which should be available to perform background and system operations.
#' @param use.apply logical, if TRUE then use the 'apply' function to apply FUN to
#'  each submatrix, or if FALSE, then directly apply FUN to submatrices, which
#'  means that FUN must return results with at least 1 dimension the same as the input, 
#'  or you can use a custom 'combine.fn' parameter to recombine results from submatrices.
#' @param convert logical, only need to change this parameter when use.apply=FALSE. If use are using a 
#'  function that can natively run on big.matrix objectsthen you can increase speed 
#'  by setting convert=FALSE. Most functions will expect a regular matrix
#'   and may fail with a big.matrix, so default convert=TRUE behaviour
#'  will convert submatrices to a regular matrix just before processing.
#' @param combine.fn a custom function to recombine input from sub.matrix processing. Default
#'  combine functions are list(), cbind() and rbind(); so a custom function should
#'  expect the same input as these; ie., a list of unspecified length, which will be
#'  the list of results from parallel calls on submatrices of bigMat, usually of size by*X.
#' @param ... if use.apply=TRUE, then additional arguments for apply(); else additional arguments
#'  for FUN.
#' @return Result depends on the function 'FUN' called, and the parameter 'combine.fn', but if MARGIN=1 usually is a vector 
#'  of length nrow(bigMat), or if MARGIN=2 a vector of length ncol(bigMat).
#' @export
#' @examples
#' # set up a toy example of a big.matrix (functions most relevant when matrix is huge)
#' bM <- filebacked.big.matrix(20, 50,
#'        dimnames = list(paste("r",1:20,sep=""), paste("c",1:50,sep="")),
#'        backingfile = "test.bck",  backingpath = getwd(), descriptorfile = "test.dsc")
#' bM[1:20,] <- replicate(50,rnorm(20))
#' prv.big.matrix(bM)
#' # compare native bigmemory column-wise function to multicore [native probably faster]
#' v1 <- colsd(bM) # native bigmemory function
#' v2 <- bmcapply(bM,2,sd,n.cores=2) # use up to 2 cores if available
#' print(all.equal(v1,v2))
#' # compare row-means approaches
#' v1 <- rowMeans(as.matrix(bM))
#' v2 <- bmcapply(bM,1,mean,n.cores=2) # use up to 2 cores if available
#' v3 <- bmcapply(bM,1,rowMeans,use.apply=FALSE)
#' print(all.equal(v1,v2)); print(all.equal(v2,v3))
#' # example using a custom combine function; taking the mean of column means
#' weight.means.to.scalar <- function(...) { X <- list(...); mean(unlist(X)) }
#' v1 <- bmcapply(bM, 2, sd, combine.fn=weight.means.to.scalar)
#' v2 <- mean(colsd(bM))
#' print(all.equal(v1,v2))
#' ## note that this function works with normal matrices, however, multicore
#' # operation is only likely to benefit speed when operations take more than 10 seconds
#' # so this function will mainly help using large matrices or intensive functions
#' test.size <- 5 # try increasing this number, or use more intensive function than sd()
#' # to test relative speed for larger matrices
#' M <- matrix(runif(10^test.size),ncol=10^(test.size-2)) # normal matrix
#' system.time(bmcapply(M,2,sd,n.cores=2)) # use up to 2 cores if available
#' system.time(apply(M,2,sd)) # 
#' unlink(c("test.bck","test.dsc"))
bmcapply <- function(bigMat,MARGIN,FUN,dir=NULL,by=200,n.cores=1,
                     use.apply=TRUE,convert=!use.apply,combine.fn=NULL,...) {
  # multicore way of calculating a function (e.g, dlrs) for a big.matrix,
  # when use.apply=T, a function that could be done with apply(); when use.apply=F, one that
  # is best done in chunks rather than row by row or col by col
  # can do it column wise (bycol=T, eg for all samples), or row-wise, eg for all snps
  # 'by' is number of rows/cols to process in each chunk
  if(is.null(dir)) { dir <- getwd() } else { if(!file.exists(dir)) { dir <- getwd() } }
  if(!as.numeric(MARGIN)[1] %in% c(1,2)) { stop("MARGIN must be 1=rows, or 2=columns") }
  if(!is.function(FUN)) { stop("FUN must be a function") }
  if(!is.numeric(by)) { by <- 200 } else { by <- round(by) }
  bycol <- as.logical(as.numeric(MARGIN[1])-1)
  #must.use.package("multicore")
  if(is.null(dim(bigMat))) { stop("this function only works on matrix objects") }
  if(length(dim(bigMat))!=2) { stop("this function only works on matrix objects") }
  if(bycol) { tot.main <- ncol(bigMat) } else { tot.main <- nrow(bigMat)}
  stepz <- round(seq(from=1,to=(tot.main+1),by=by))
  # if the last step is too small, merge with previous
  sc.lst <- head(tail(stepz,2),1); lst <- tail(stepz,1)
  if((tail(stepz,1)) != (tot.main+1)) { stepz <- c(stepz,(tot.main+1)) }
  if(lst-sc.lst <=2) { ll <- length(stepz); if(ll>2) { stepz <- stepz[-(ll-1)] } else { warning("number of rows/columns quite small, may cause issues") } }
  #
  split.to <- length(stepz)-1
  result <- numeric(tot.main)
  ## define the function
  big.fnc <- function(dd,func,stepz,bigMat,dir,bycol=T,...)
  {
    big.type <- is.big.matrix(bigMat)
    x1 <- stepz[dd]; x2 <- stepz[dd+1]-1 #subset row selection
    if(bycol) {
      if(big.type) {
        next.block <- sub.big.matrix(bigMat, firstCol=x1, lastCol=x2, backingpath=dir )
      } else {
        next.block <- bigMat[,x1:x2]
      }
      dm <- 2
    } else {
      if(big.type) {
        next.block <- sub.big.matrix(bigMat, firstRow=x1, lastRow=x2, backingpath=dir )
      } else {
        next.block <- bigMat[x1:x2,]
      }
      dm <- 1
    }
    if(convert | is.null(dim(next.block))) { 
      next.block <- bigmemory::as.matrix(next.block)
#      next.block <- next.block[1:nrow(next.block),1:ncol(next.block)]  # conv to standard matrix if req'd
    } 
    if(!use.apply) { out <- func(next.block,...)  } else {  out <- biganalytics::apply(next.block,dm,func,...) }
    rm(next.block) ; gc() # remove the sub-matrix pointer each iteration  
    return(out)
  }
  ## run function as mclapply()
  if(split.to>=1) {
    #print("z1")
    result.list <- parallel::mclapply(1:split.to, FUN=big.fnc, func=FUN, stepz=stepz, 
                                       bigMat=bigMat, dir=dir, bycol=bycol, mc.cores=n.cores,...)
    #print("z2")
    if(is.function(combine.fn)) {
      result <- do.call(combine.fn,args=result.list)
    } else {
      comb.fn <- choose.comb.fn(result.list,stepz) # automatically decide best combining function based on largest and most common result
      result <- do.call(comb.fn,args=result.list)
      #result <- unlist(result.list,recursive=unlist.results) 
    }
  } else {
    result <- NULL; warning("matrix had insufficient columns returning NULL")
  }
  return(result)  
}


## internal function
choose.comb.fn <- function(result.list,stepz) {
  ### NB: when NCmisc is a package this should require(NCmisc)
  # choose the best function to combine multicore results from bigmcapply
  ## DEFINE function to choose best based on 1 example result:
  sub.comb.fn <- function(dm,ls,ll) {
    # evaluate output based on dimensions of an example result
    if(ls==ll) { comb.fn <- "c" } else {
      if(!is.null(dm)){
        if(dm[1]==dm[2] & dm[1]==ls) {
          warning("looks like results have dim[1]==dim[2], suggest using 'combine.fn' parameter to obtain desired result")
          comb.fn <- c("list","cbind","rbind")
        } else {
          if(dm[1]==ls | dm[2]==ls) {
            if(dm[2]==ls) {
              comb.fn <- "cbind"
            } else {
              comb.fn <- "rbind"
            }
          } else {
            warning("mapping from function results was not obvious, suggest using 'combine.fn' parameter")
            comb.fn <- "list"
          }
        }
      } else {
        comb.fn <- "list"
      }
    }
    return(comb.fn)
  }
  ## process max and mode results using function defined above
  rawlens <- sapply(result.list,function(x) { length(unlist(x)) })
  #print(is(rawlens)); print(dim(rawlens)); print(head(rawlens))
  lL <- length(rawlens); if(length(rawlens)>1) { rawlens <- rawlens[-lL] }
  max.res <- which(rawlens==max(rawlens,na.rm=T))[1] # choose the largest result to define structure
  mode.res <- which(rawlens==Mode(rawlens,multi=T))[1] # choose the largest result to define structure
  dm <- dim(result.list[[max.res]]); ll <- length(result.list[[max.res]]); ls <- length(stepz[max.res]:stepz[max.res+1])-1
  fn.max <- sub.comb.fn(dm,ls,ll)
  dm <- dim(result.list[[mode.res]]); ll <- length(result.list[[mode.res]]); ls <- length(stepz[mode.res]:stepz[max.res+1])-1
 # print(dm); print(ls); print(ll); print(mode.res); print(max.res)
  fn.mode <- sub.comb.fn(dm,ls,ll)
  if(fn.max[1]!=fn.mode[1]) {
    # resolve any conflict between results based on max or mode
    allf <- c(fn.max,fn.mode); 
    if(length(allf)==4) { 
      fn.max <- Mode(allf) # by chance one had equal dimensions, but this should pick the best by convergence
    } else {
      warning("differential testing of result output from multiple cores gave different structures, suggest using 'combine.fn' to ensure correct output")
    }
  }
  return(fn.max[1])
}


#' Attempt to install the bigalgebra package using SVN
#'
#' The bigalgebra package for efficient algebraic operations on big.matrix objects
#' has now been submitted to CRAN, so this function is now mostly redundant.
#' It used to require installation from SVN and some tinkering, such as changing the 
#' description file to add the dependency, and linking 'BH' to allow the package to work.
#' This may still be required on older versions of R that do not support the bigalgebra
#' package uploaded to CRAN, but I cannot confirm this.
#' This function automatically performs these corrections. First, it attempts to check-out
#'  the latest version of bigalgebra from SVN version management system and then corrects 
#' the description file, then tries to install the package.
#' Note you must also have 'BLAS' installed on your system to utilise this package
#' effectively. PCA functions in the present package are better with bigalgebra installed,
#' but will still run without it. For more information on installation alternatives, 
#' type big.algebra.install.help().
#' Returns TRUE if bigalgebra is already installed.
#' @seealso \code{\link{big.algebra.install.help}}
#' @param verbose whether to report on installation progress/steps
#' @return If SVN is installed on your system, along with BLAS, this function should install the bigalgebra package,
#'  else it will return instructions on what to do to fix the issue 
#' @export
#' @examples
#' # not run # svn.bigalgebra.install(TRUE)
svn.bigalgebra.install <- function(verbose=FALSE) {
  # this is a major hack to install bigalgebra from SVN,
  # manually modifying the DESCRIPTION file to depend and link to 'BH'
  warning("bigalgebra is now on CRAN and so this function should be unnecessary")
  cur.dir <- getwd()
  cat("\nAttempting to install the bigalgebra package using SVN")
  if(verbose) { cat("\n") } else { cat(" .") }
  #my.fn <- file("bigalgebra.install.log",open="w")
  #sink(file=my.fn,type="message")
  if(!check.linux.install("svn")) { return(F) }
  nons <- cat.path(getwd(),"tempfdsg345t")
  dir.create(nons)
  setwd(nons)
  system("svn checkout svn://scm.r-forge.r-project.org/svnroot/bigmemory",
         intern=!verbose, ignore.stderr=!verbose)
  cat(".")
  setwd("./bigmemory/pkg")
  a.mod <- F
  des.fn <- "./bigalgebra/DESCRIPTION"
  if(is.file(des.fn,dir=getwd())) {
    DES <- readLines(des.fn); cat(".")
    l1 <- grep("Depends: bigmemory",DES)
    l2 <- grep("LinkingTo: bigmemory",DES)
    if(length(l1)==1 & length(l2)==1) {
      if(length(grep("BH",l1))==0) {
        if(verbose) { cat("modifying bigalgebra DESCRIPTION file to depend on BH\n") }
        DES[l1] <- gsub("Depends: bigmemory","Depends: BH, bigmemory",DES[l1])
        a.mod <- T; cat(".")
      }
      if(length(grep("BH",l2))==0) {
        if(verbose) { cat("modifying bigalgebra DESCRIPTION file to link to BH\n") }
        DES[l2] <- gsub("LinkingTo: bigmemory","LinkingTo: BH, bigmemory",DES[l2])
        a.mod <- T  ; cat(".")
      }
    }
    if(a.mod) {
      writeLines(DES,con=des.fn); cat(".")
    }
    system("REFBLAS=1 R CMD INSTALL bigalgebra",intern=!verbose, ignore.stderr=!verbose)
    cat(". done\n")
    suc <- T
  } else {
    warning("bigalgebra DESCRIPTION file not found, installation failed")
    suc <- F
  }
  setwd(nons)
  system("rm -rf bigmemory", intern=!verbose, ignore.stderr=!verbose)
  setwd(cur.dir)
  unlink(nons)
  return(suc)
}


#' Attempt to install the bigalgebra package
#'
#' The bigalgebra package  has now been submitted to CRAN, so this function is now 
#' mostly redundant. It may still be useful for some, and it will still work,
#' as the first step to check CRAN, so at the risk of affecting existing code
#' I will leave the function here for now.
#' This function attempts to see whether bigalgebra is installed, then checks CRAN in case it 
#' has been updated, then check RForge. Failing that, it will attempt to install
#' using svn.bigalgebra.install(). Returns TRUE if already installed.
#' The bigalgebra package for efficient algebraic operations on big.matrix objects
#' was not currently on CRAN, and used to fail a check on dependencies. Changing the 
#' description file was needed to add the dependency, and linking 'BH' allow3e the package to work.
#' This function attempts to check-out the latest version of bigalgebra from SVN
#' version management system and corrects the description file then installs.
#' Note you must also have 'BLAS' installed on your system to utilise this package
#' effectively. PCA functions in the present package are better with bigalgebra installed,
#' but will still run without it. For more information on installation alternatives, 
#' type big.algebra.install.help().
#' @seealso \code{\link{svn.bigalgebra.install}}
#' @param verbose whether to report on installation progress/steps
#' @return If bigalgebra is already installed, or can be installed from RForge or SVN,
#'  this should load or install the bigalgebra package,
#'  else will return instructions on what to do next to fix the issue 
#' @export
#' @examples 
#' # not run # big.algebra.install.help(TRUE)
big.algebra.install.help <- function(verbose=FALSE) {
  ## bigalgebra package doesn't install easily using the regular R way of installing packages
  # here try a simple way that might work, and if not, provide links and instructions to 
  # guide a manual installation
  try({ if(do.call("require",args=list("bigalgebra"))) { return(T) } })
  if("bigalgebra" %in% search.cran("big")[[1]]) { must.use.package("bigalgebra",T); return() }
  warning("bigalgebra is now on CRAN and so this function should be unnecessary")
  cat("\nbigalgebra installation not found, will attempt to install now, but it can be tricky\n")
  do.call("install.packages",args=list("bigalgebra", repos="http://R-Forge.R-project.org"))
  if(do.call("require",args=list("bigalgebra"))) {
    cat("bigalgebra seems to have installed successfully\n")
    return(T)
  } else {
    tt <- svn.bigalgebra.install(verbose=verbose)
    if(!tt) {
      cat("standard bigalgebra installation has failed\n")
      cat("go to: http://connectir.projects.nitrc.org/download/\n")
      cat("for installation tips\n")
      cat("can use a command line like: REFBLAS=1 R CMD INSTALL bigalgebra\n")
      cat("where bigalgebra is the source package (for example, bigalgebra_0.8.1.tar.gz)\n")
      cat("downloaded from https://r-forge.r-project.org/R/?group_id=556\n")
      cat("Your system may also be missing a BLAS installation, in which case you might try\n")
      cat("installing OpenBLAS; see instructions at http://xianyi.github.com/OpenBLAS/\n")
      return(F)
    } else {
      return(T)
    }
  }
}


#' Retrieve a big.matrix object
#'
#' This function can load a big.matrix object using a big.matrix.descriptor object, the
#' name of a description file, the name of a binary file containing a big.matrix.descriptor
#' or if passed a big.matrix object, it will just return that object. Only the object or
#' file name plus the directory containing the backing file are required.
#' @param fn the name of a description file, the name of a binary file containing a
#'  big.matrix.descriptor, a big.matrix object or a big.matrix.descriptor object.
#' @param dir directory containing the backing file (if not the working directory)
#' @param verbose whether to display information on method being used, or minor warnings
#' @return Returns a big.matrix object, regardless of what method was used as reference/input
#' @export
#' @examples 
#' # set up a toy example of a big.matrix 
#' bM <- filebacked.big.matrix(20, 50,
#'        dimnames = list(paste("r",1:20,sep=""), paste("c",1:50,sep="")),
#'        backingfile = "test.bck",  backingpath = getwd(), descriptorfile = "test.dsc")
#' bM[1:20,] <- replicate(50,rnorm(20))
#' # Now have a big matrix which can be retrieved using this function in 4 ways:
#' d.bM <- describe(bM)
#' save(d.bM,file="fn.RData")
#' bM1 <- get.big.matrix("test.dsc")
#' bM2 <- get.big.matrix(d.bM)
#' bM3 <- get.big.matrix("fn.RData")
#' bM4 <- get.big.matrix(bM)
#' prv.big.matrix(bM)
#' prv.big.matrix(bM1)
#' prv.big.matrix(bM2)
#' prv.big.matrix(bM3)
#' prv.big.matrix(bM4)
#' unlink(c("fn.RData","test.bck","test.dsc"))
get.big.matrix <- function(fn,dir="",verbose=FALSE)
{
  # loads a big.matrix either using an big.matrix description object
  # , or this object in a binary file or text file, or points to a bigmatrix or matrix
  if(all(dir=="")) { dir <- getwd() }
  if(exists("validate.dir.for",mode="function")) {
    ## plumbCNV specific code ##
    dir <- do.call("validate.dir.for",list(dir=dir,elements="big",warn=F))  
    dir.big <- dir$big
  } else {
    # otherwise
    dir.big <- dir
    if(is.list(dir)) { if(!is.null(dir[["big"]])) { dir.big <- dir$big } }
  }
  if(is(fn)[1]=="big.matrix.descriptor")
  {
    bigMat2 <- attach.big.matrix(fn,path=dir.big)
  } else {
    if(is(fn)[1]=="big.matrix" | is(fn)[1]=="matrix")
    {
      if(is(fn)[1]=="matrix") {
        bigMat2 <- as.big.matrix(fn,descriptorfile="TEMPBIG",backingpath=dir.big)
      } else { bigMat2 <- fn }
    } else {
      lastchar <- substr(dir.big,nchar(dir.big),nchar(dir.big))
      if (length(grep(".RData",fn))==0) {
        fn <- basename(fn)
        if(!fn %in% list.files(dir.big)) { 
          stop(paste("Error: big.matrix file '",fn,"' not in 'dir.big'",sep=""))
        }
        if(verbose) { cat(" loading big matrix using text description\n") }
        if(lastchar=="/") { dir.big <- substr(dir.big,1,nchar(dir.big)-1) }
        bigMat2 <- attach.big.matrix(fn,path=dir.big)
      } else {
        if(verbose) { cat(" loading big matrix using RData description\n") }
        if(lastchar!="/") { dir.big <- paste(dir.big,"/",sep="") }
        filenm <- cat.path(dir.big,fn,must.exist=T)
        dscnm <- paste(load(filenm))
        big.fn <- NULL
        for (ii in 1:length(dscnm)) {
          if("big.matrix.descriptor" %in% is(get(dscnm[ii])))
          { big.fn <- dscnm[ii] } 
        }
        if(!is.null(big.fn)) {
          descr <- get(big.fn) 
        } else {
          stop(paste("Error: didn't find bigmatrix descriptor in file",fn))
        }
        bigMat2 <- attach.big.matrix(descr,path=dir.big) 
      }
    }
  }
  return(bigMat2)
}





# internal? function to quickly extract whether a delimited matrix
# has row or column names, and what the column names are
quick.mat.format <- function(fn) {
  temp.fn <- "twmpwerw123.txt"
  txtx <- readLines(fn,n=11) # change this to like 5 when reader is updated to save memory
  writeLines(txtx,con=temp.fn)
  first2 <- reader(temp.fn)
  unlink(temp.fn) # remove this temporary file straight away
  if(is.null(dim(first2))) {
    ## assuming no header row #
    return(list(rownames=F,colnames=F,ncol=1,cnames=NULL))
  }
  rn <- rownames(first2); if(rn[1]=="1") { rn <- NULL }
  if(!is.null(rn)) { ern <- T } else { ern <- F }
  cn <- colnames(first2); if(cn[1]=="V1") { cn <- NULL }
  if(!is.null(cn)) { ecn <- T } else { ecn <- F }
  ncl <- ncol(first2)
  return(list(rownames=ern,colnames=ecn,ncol=ncl,cnames=cn))
}


## internal function for import.big.matrix
check.text.matrix.format <- function(fn,ncol=NA,header=NULL,row.names=NULL,sep="\t") 
{
  # read the first two lines of a text matrix file and determine
  # whether it has row and or column names, return T/F indices for this
  # plus a further 
  # assumes tab as separator
  fn <- fn[1]
  if(!is.character(fn)) { warning("file name 'fn' should be a character()") }
  if(!file.exists(fn)) { stop(paste("Error: file",fn,"not found")) } 
  dat.file <- file(fn); headrow <- NULL
  rnames.in.file <- F; first.is.head <- F; name.lookup <- F
  # read first two lines of matrix format datafile
  open(con=dat.file,open="r")
  next.line <- readLines(dat.file,n=1)
  line1 <- strsplit(next.line,sep,fixed=T)[[1]]
  next.line <- readLines(dat.file,n=1)
  line2 <- strsplit(next.line,sep,fixed=T)[[1]]
  close(con=dat.file)
  ## check whether first row is header or just plain data
  frst <- length(line1)
  scnd <- length(line2)
  if(is.na(ncol)) {
    if((scnd-frst)==1) { rnames.in.file <- T; first.is.head <- T; name.lookup <- F;
                         headrow <- line1 }
  } else {
   # prv(c("frst","ncol"))
    if (frst!=ncol & frst!=ncol+1) {
      #preview(c("frst","ncol"))
      stop("dimensions of import file do not match id lists specified, exiting")
      break; break; 
    } else {
      if(length(which(paste(line1[1:10]) %in% paste(header)))>8) {
        # first line seems to be the header
        if(!all(header==paste(line1[(frst-ncol+1):frst]))) {
          stop("Error: ID list did not match file header. Please check the files (bash$ head -1 <filename>) and fix this")
        } else {
          # will need to go to next line to avoid reading header as data
          first.is.head <- T
        }
      } else {
        first.is.head <- F
      }
    }
    ## check whether second row starts with a rowname, or just plain data
    if(length(line2)==ncol+1) {
      rnames.in.file <- T
    } else {
      # seem to be no rownames in file
      rnames.in.file <- F
    }
  }
  # if colnames row found + header specified, check they are the same
  if(first.is.head & !is.null(header)) {
    if(!all(paste(line1)==paste(header))) {
      if(!paste(line2[1]) %in% paste(header)) {
        warning("there seems to be a header of column labels in the raw file but not the header specified\n",
                "proceeding with these, but please check this mismatch is expected (e.g, might be col numbers)")
        name.lookup <- F
      } else {
        warning("order of header (column labels) in data file does not match order of inputted 'header' list\n",
                "proceeding, but using name lookup method will make import very slow and possibly unstable\n",
                "If your source file has an inconsistent order this is the only supported import method.\n",
                "Otherwise, recommend to cancel (ctrl-C), amend this discrepancy and run again")
        name.lookup <- T
      }
    }
  }
  # if rownames row found + rownames specified, check they are the same
  if(rnames.in.file & !is.null(row.names)) {
   # prv(c("line2","row.names"),counts=list(1,1))
    if(!paste(line2[1]) %in% paste(row.names[1:2])) {
      if(!paste(line2[1]) %in% paste(row.names)) {
        warning("there seem to be row labels in the raw file but not the rownames specified\n",
                "proceeding with these, but please check this mismatch is expected (e.g, might be row numbers)")
        name.lookup <- F
      } else {
        warning("order of row labels in data file does not match order of inputted rowname list\n",
                "proceeding, but using name lookup method will make import very slow and possibly unstable\n",
                "If your source file has an inconsistent order this is the only supported import method.\n",
                "Otherwise, recommend to cancel (ctrl-C), amend this discrepancy and run again")
        name.lookup <- T
      }
    }
  }
  out <- list(first.is.head,rnames.in.file,name.lookup,headrow)
  names(out) <- c("header","rnames","match","colnames")
  return(out)
}



# # patch some reader functions that cause trouble for long format files :(
# # fix problem in reader:::get.delim :(
# get.delim <- function(...,delims=c("\t"," ","\t| +",";",",")) {
#   return(reader:::get.delim(...,delims=delims))
# }
# # fix problem in reader:::file.ncol :(
# file.ncol <- function(fn,...) { reader:::file.ncol(fn,del=get.delim(fn),...) }
# 
# # add some more known extensions to rmv.ext
# rmv.ext <- function(...) {
#   return(reader:::rmv.ext(...,more.known=c("RDA","DSC","BCK")))
# }


#' Generate a test matrix of random data
#' 
#' Generates a test matrix of easily specified size and type. Options allow
#' automated row and column names (which might resemble labels for a SNP analysis)
#' and return of several different formats, matrix, data.frame or big.matrix.
#' You can specify the randomisation function (e.g, rnorm, runif, etc), as well
#' as parameters determining the matrix size. Can also generate big.matrix objects,
#' and an important feature is that the method to generate big.matrix objects is
#' scalable so that very large matrices for simulation can be generated only limited
#' by disk space and not by RAM.
#' @param size 10^size is the total number of datapoints simulated. 6 or less are fairly quick to generate,
#'  while 7 takes a few seconds. 8 will take under a minute, 9 around ten minutes, 10, perhaps over an hour.
#'  Values are coerced to the range of integers c(2:10).
#' @param row.exp similar to 'nrow' when creating a matrix, except this is exponential, giving 10^row.exp rows.
#' @param rand a function, must return 'n' values, when rand(n) is called, eg., rnorm(), runif(), numeric()
#' @param dimnames logical, whether to generate some row and column names
#' @param data.frame logical, whether to return as a data.frame (FALSE means return a matrix)
#' @param big.matrix logical, whether to return as a big.matrix (overrides data.frame). If a file.name
#'  is used then the big.matrix will be filebacked and this function returns a list with a 
#'  a big.matrix, and the description and backing filenames.
#' @param file.name if a character, then will write the result to tab file instead of returning
#'  the object, will return the filename; overrides data.frame. Alternatively, if big.matrix=TRUE,
#'  then this provides the basename for a filebacked big.matrix.
#' @param tracker logical, whether to display a progress bar for large matrices (size>7) where progress will be slow
#' @return Returns a random matrix of data for testing/simulation, can be a data.frame or big.matrix if those options are selected
#' @export
#' @author Nicholas Cooper 
#' @examples 
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' mat <- (generate.test.matrix(5)); prv(mat)
#' lst <- (generate.test.matrix(5,3,big.matrix=TRUE,file.name="bigtest"))
#' mat <- lst[[1]]; prv(mat); headl(lst[2:3]); 
#' unlink(unlist(lst[2:3]))
#' setwd(orig.dir) # reset working dir to original
generate.test.matrix <- function(size=5,row.exp=2,rand=rnorm,dimnames=TRUE,
                                 data.frame=FALSE,big.matrix=FALSE,file.name=NULL, 
                                 tracker=TRUE) {
  test.size <- min(max(2,round(size)),10) # try increasing this number for larger matrices
  dim.dif <- max(min(row.exp,size-2),1)
  big.d <- 10^(test.size-dim.dif); bigger.d <- (10^(test.size-(dim.dif-1)))-1
  other.d <- 10^(dim.dif) ; another.d <- (10^(dim.dif+1))-1
  if(!is.function(rand)) { stop("rand must be a function that generates 'n' simulated values") }
  nr <- (10^test.size)/big.d; nc <- big.d
  if(dimnames) {
    rown <- paste("rs",sample(10:99,nr,replace=T),sample(other.d:another.d,nr),sep="")
    coln <- paste("ID",sample(1:9,nc,replace=T),sample(big.d:bigger.d,nc),sep="")
    dnn <- list(rown,coln)
  } else { dnn <- NULL }
  M <- NULL
  if(big.matrix) {
    if(test.size>7) { if(is.null(file.name)) { file.name="big_test" } }
    # automatically extract the path if one is included in file.name
    if(!is.null(file.name)) {
      if(dirname(file.name)!=".") { pth <- dirname(file.name); file.name <- basename(file.name) } else { pth <- NULL}
      bck <- cat.path(fn=file.name,ext="bck")
      descr <- cat.path(fn=file.name,ext="dsc")
    } else { bck <- descr <- NULL }
    if(test.size>7) {
     M <- big.matrix(nrow=nr,ncol=nc,backingfile=bck,descriptorfile=descr,dimnames=dnn,backingpath=pth)
     if(nr>nc) {
       for (cc in 1:nr) { M[cc,] <- rand(nc); if(tracker) { loop.tracker(cc,nr) }  }
     } else {
       for (cc in 1:nc) { M[,cc] <- rand(nr);  if(tracker) { loop.tracker(cc,nc) } }
     }
    }
  }
  if(is.null(M)) {
    M <- matrix(rand(10^test.size),ncol=nc) # normal matrix
    colnames(M) <- coln; rownames(M) <- rown
  }
  #prv(big.d,bigger.d,other.d,another.d,M)
  if(is.character(file.name) & !big.matrix) {
    write.table(M,sep="\t",col.names=dimnames,row.names=dimnames,file=file.name,quote=F) # no dimnames
    return(file.name)
  } 
  if(data.frame & !big.matrix) {
    return(as.data.frame(M))
  }
  if(big.matrix) {
    if(is.character(file.name)) {
      if(!is.big.matrix(M)) { M <- as.big.matrix(M,descriptorfile=descr,backingfile=bck) }
      return(list(M=M,descr=descr,bck=bck))
    } else {
      return(as.big.matrix(M))
    }
  }
  return(M)
}



#' Internal function used by import.big.data
dir.force.slash <- function(dir) {
  # make sure 'dir' directory specification ends in a / character
  if(!is.null(dim(dir))) { stop("dir should be a vector") }
  dir <- paste(dir)
  dir.ch <- .Platform$file.sep
  the.test <- (dir!="" & substr(dir,nchar(dir),nchar(dir))!=dir.ch)
  dir[the.test] <- paste(dir[the.test],dir.ch,sep="")
  return(dir)
}


#' Load a text file into a big.matrix object
#'
#' This provides a faster way to import text data
#' into a big.matrix object than bigmemory::read.big.matrix(). The
#' method allows import of a data matrix with size exceeding RAM limits.
#' Can import from a matrix delimited file with or without row/column names,
#' or from a long format dataset with no row/columns names (these should be
#' specified as separate lists).
#' @param input.fn character, or list, either a single file name of the data, or
#'  a list of multiple file name if the data is stored as multiple files. If multiple,
#'  then the corresponding list of row or column names that is unique between files
#'  should be a list of the same length.
#' @param dir character, the directory containing all files. Or, if files are split between
#'  directories, then either include the directories explicitly in the filenames, or
#'  multiple directories can be entered as a list, with names 'big', 'ano' and 'col', where
#'  big is the location for big.matrix objects to file-back to, 'ano' is the location
#'  of row and column names, and 'col' is the location of the raw text datafiles.
#' @param long logical, if TRUE, then the data is assumed to be in long format, where
#'  each datapoint is on a new line, and the file is structured so that the data for
#'  each case/sample/id is consecutive and ordered consistently between samples. If using
#'  long format the file should contain no row or column names, these should be specified
#'  in either rows.fn/cols.fn file name arguments, or row.names/col.names vector arguments.
#'  If long=FALSE, then the dimensions of the file will be automatically detected; including
#'  if the file is in long format, however, if you know the data is in long format, specifying
#'  this explicitly will be quicker and guarantees the correct import method.
#' @param rows.fn character, with the name of a text file containing the list of row labels
#'  for the dataset. Unnecessary if importing from a matrix with row/column names in the file,
#'  or if using the row.names parameter. Must be a list of filenames if row names are split
#'  across multiple input.fn files.
#' @param cols.fn character, with the name of a text file containing the list of column labels
#'  for the dataset. Unnecessary if importing from a matrix with row/column names in the file,
#'  or if using the col.names parameter. Must be a list of filenames if column names are split
#'  across multiple input.fn files.
#' @param pref character, optional prefix to use in naming the big.matrix files (description/backing files)
#' @param delete.existing logical, if a big.matrix already exists with the same name as implied
#'  by the current 'pref' and 'dir' arguments, then default behaviour (FALSE) is to return an error.
#'  to overwrite any existing big.matrix file(s) of the same name(s), set this parameter to TRUE.
#' @param ret.obj logical, whether to return a big.matrix.descriptor object (TRUE), or
#'  just the file name of the big.matrix description file of the imported dataset.
#' @param verbose logical, whether to display extra information about import progress and
#'  notifications.
#' @param row.names character vector, optional alternative to specifying rows.fn file name(s),
#'  directly specify row names as a single vector, or a list of vectors if multiple input files
#'  with differing row names are being imported.
#' @param col.names character vector, optional alternative to specifying cols.fn file name(s),
#'  directly specify oclumn names as a single vector, or a list of vectors if multiple input files
#'  with differing column names are being imported.
#' @param dat.type character, data type being imported, default is "double", but can specify any type
#'  supported by a filebacked.big.matrix(), namely, "integer","char","short"; note these
#'  are C-style data types; double=numeric, char=character, integer=integer, short=numeric (although
#'  will be stored with less precision in the C-based big.matrix object).
#' @param ram.gb numeric, the number of gigabytes of free RAM that it is ok for the import
#'  to use. The higher this amount, the quicker the import will be, as flushing RAM contents
#'  to the hard drive more regularly slows down the process. Setting this lower
#'  will reduce the RAM footprint of the import. Note that if you set it too high, it can't
#'  be guaranteed, but usually R and bigmemory will do a reasonable job of managing the memory,
#'  and it shouldn't crash your computer.
#' @param hd.gb numeric, the amount of free space on your hard disk; if you set this
#'  parameter accurately the function will stop if it believes there is insufficient
#'  disk space to import the object you have specified. By default this is set to 
#'  1 terabyte, so if importing an object larger than that, you will have to increase
#'  this parameter to make it work.
#' @param tracker logical, whether to display a progress bar for the importing process
#' @return Returns a big.matrix containing the data imported (single big.matrix even
#'  when text input is split across multiple files)
#' @export
#' @examples 
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' # Collate all file names to use in this example #
#' all.fn <- c("rownames.txt","colnames.txt","functestdn.txt","funclongcol.txt","functest.txt",
#'  paste("rn",1:3,".txt",sep=""),paste("cn",1:3,".txt",sep=""),
#'  paste("split",1:3,".txt",sep=""),
#'  paste("splitmatCd",1:3,".txt",sep=""),paste("splitmatRd",1:3,".txt",sep=""),
#'  paste("splitmatC",1:3,".txt",sep=""), paste("splitmatR",1:3,".txt",sep=""))
#' any.already <- file.exists(all.fn)
#' if(any(any.already)) { 
#'  warning("files already exist in the working directory with the same names as some example files") }
#' # SETUP a test matrix and reference files # 
#' test.size <- 4 # try increasing this number for larger matrices
#' M <- matrix(runif(10^test.size),ncol=10^(test.size-2)) # normal matrix
#' write.table(M,sep="\t",col.names=FALSE,row.names=FALSE,
#'  file="functest.txt",quote=FALSE) # no dimnames
#' rown <- paste("rs",sample(10:99,nrow(M),replace=TRUE),sample(10000:99999,nrow(M)),sep="")
#' coln <- paste("ID",sample(1:9,ncol(M),replace=TRUE),sample(10000:99999,ncol(M)),sep="")
#' r.fn <- "rownames.txt"; c.fn <- "colnames.txt"
#' Mdn <- M; colnames(Mdn) <- coln; rownames(Mdn) <- rown
#' # with dimnames
#' write.table(Mdn,sep="\t",col.names=TRUE,row.names=TRUE,file="functestdn.txt",quote=FALSE) 
#' prv.large(Mdn)
#' writeLines(paste(as.vector(M)),con="funclongcol.txt")
#' in.fn <- "functest.txt"
#' 
#' ### IMPORTING SIMPLE 1 FILE MATRIX ##
#' writeLines(rown,r.fn); writeLines(coln,c.fn)
#' #1. import without specifying row/column names
#' ii <- import.big.data(in.fn); prv.big.matrix(ii) # SLOWER without dimnames!
#' #2. import using row/col names from file
#' ii <- import.big.data(in.fn,cols.fn="colnames.txt",rows.fn="rownames.txt")
#' prv.big.matrix(ii)
#' #3. import by passing colnames/rownames as objects
#' ii <- import.big.data(in.fn, col.names=coln,row.names=rown)
#' prv.big.matrix(ii)
#' 
#' ### IMPORTING SIMPLE 1 FILE MATRIX WITH DIMNAMES ##
#' #1. import without specifying row/column names, but they ARE in the file
#' in.fn <- "functestdn.txt"
#' ii <- import.big.data(in.fn); prv.big.matrix(ii)
#' 
#' ### IMPORTING SIMPLE 1 FILE MATRIX WITH MISORDERED rownames ##
#' rown2 <- rown; rown <- sample(rown);
#' # re-run test3 using in.fn with dimnames
#' ii <- import.big.data(in.fn, col.names=coln,row.names=rown)
#' prv.big.matrix(ii)
#' # restore rownames: 
#' rown <- rown2
#' 
#' ### IMPORTING SIMPLE 1 FILE LONG FORMAT by columns ##
#' in.fn <- "funclongcol.txt"; #rerun test 2 # 
#' ii <- import.big.data(in.fn,cols.fn="colnames.txt",rows.fn="rownames.txt")
#' prv.big.matrix(ii)
#' 
#' ### IMPORTING multifile LONG by cols ##
#' # create the dataset and references
#' splF <- factor(rep(c(1:3),ncol(M)*c(.1,.5,.4)))
#' colnL <- split(coln,splF); MM <- as.data.frame(t(M))
#' Ms2 <- split(MM,splF)
#' Ms2 <- lapply(Ms2,
#'    function(X) { X <- t(X); dim(X) <- c(nrow(M),length(X)/nrow(M)); X } )
#' # preview Ms2 - not run # lapply(Ms2,prv.large)
#' colfs <- paste("cn",1:length(colnL),".txt",sep="")
#' infs <- paste("split",1:length(colnL),".txt",sep="")
#' # create multiple column name files and input files
#' for(cc in 1:length(colnL)) { writeLines(colnL[[cc]],con=colfs[cc]) }
#' for(cc in 1:length(infs)) { 
#'   writeLines(paste(as.vector((Ms2[[cc]]))),con=infs[cc]) }
#'   
#' # Now test the import using colnames and rownames lists
#' ii <- import.big.data(infs, col.names=colnL,row.names=rown)
#' prv.big.matrix(ii)
#' 
#' ### IMPORTING multifile MATRIX by rows ##
#' # create the dataset and references
#' splF <- factor(rep(c(1,2,3),nrow(M)*c(.1,.5,.4)))
#' rownL <- split(rown,splF)
#' Ms <- split(M,splF)
#' Ms <- lapply(Ms,function(X) { dim(X) <- c(length(X)/ncol(M),ncol(M)); X } )
#' # preview Ms - not run # lapply(Ms,prv.large)
#' # create multiple row name files and input files
#' rowfs <- paste("rn",1:length(rownL),".txt",sep="")
#' for(cc in 1:length(rownL)) { writeLines(rownL[[cc]],con=rowfs[cc]) }
#' infs <- paste("splitmatR",1:length(colnL),".txt",sep="")
#' for(cc in 1:length(infs)) { 
#'  write.table(Ms[[cc]],sep="\t",col.names=FALSE,row.names=FALSE,file=infs[cc],quote=FALSE) }
#'  
#' # Now test the import using colnames and rownames files
#' ii <- import.big.data(infs, col.names="colnames.txt",rows.fn=rowfs)
#' prv.big.matrix(ii)
#' 
#' # DELETE ALL FILES ##
#' unlink(all.fn[!any.already]) # prevent deleting user's files
#' ## many files to clean up! ##
#' unlink(c("funclongcol.bck","funclongcol.dsc","functest.bck","functest.dsc",
#'  "functestdn.RData","functestdn.bck","functestdn.dsc","functestdn_file_rowname_list_check_this.txt",
#'  "split1.bck","split1.dsc","splitmatR1.bck","splitmatR1.dsc"))
#' setwd(orig.dir) # reset working dir to original
import.big.data <- function(input.fn=NULL, dir=getwd(), long=FALSE, rows.fn=NULL, cols.fn=NULL, 
                              pref="", delete.existing=TRUE, ret.obj=FALSE, verbose=TRUE, row.names=NULL, col.names=NULL,
                              dat.type="double", ram.gb=2, hd.gb=1000, tracker=TRUE)
{
  # import from a text (hopefully long format) datafile to a big.matrix
  # return bigmatrix description 'object' or name of description file according to 'ret.obj'
  ### CALIBRATE OPTIONS AND PERFORM CHECKS ###
  ## dir.force.slash <- reader:::dir.force.slash # dont use internal function from 'reader'
  if(all(dir=="")) { dir <- getwd() }
  ## for compatibility with plumbCNV directory object
  if(exists("validate.dir.for",mode="function")) {
    ## plumbCNV specific code ##
    dir <- do.call("validate.dir.for",list(dir=dir,elements=c("big","ano","col"),warn=FALSE))  
  } else {
    # otherwise
    dir <- list(big=dir,ano=dir,col=dir)
    if(is.list(dir)) { if(!is.null(dir[["big"]])) { dir.big <- dir$big } }
  }
  dat.type <- tolower(dat.type)
  if(!dat.type %in% c("double","short","character","integer","char","numeric")) {
    dat.type <- options()$bigmemory.default.type; as.type <- dat.type
  } else {
    as.type <- dat.type
    if(dat.type=="short") { as.type <- "numeric" }
    if(dat.type=="char") { as.type <- "character" }
    if(dat.type=="character") { dat.type <- "char" }
    if(dat.type=="numeric") { dat.type <- "double" }
  }
  ## Define data types for big.matrix
    # label
  input.fn <- cat.path(dir$col,input.fn)
  file.rn <- character()
  spec.rn <- spec.cn <- T; miswarn <- F
  #### GET THE SHAPED INPUT FILE NAME(S) ####
  # print(input.fn)
  if(length(input.fn)>1 & ((is.null(rows.fn)&is.null(row.names)) | (is.null(cols.fn)&is.null(col.names)))) {
    stop(paste("When using multiple input files, you must specify row and column names using",
               "row.names + col.names, or rows.fn + cols.fn.",
               "If the input files split the dataset by columns, then the column file names need",
               "to be entered as a list corresponding to unique column names in each input.fn entry,",
               "OR equivalently, a equal length list of row names if the data is split by rows"))
  }
  #### GET THE CORRESPONDING column LIST(S) ####
  if((length(col.names)>length(input.fn)) | is.list(col.names)) {
    if(is.list(col.names))  {
      if(all(sapply(col.names,is)[1,]=="character")) {
        ID.list <- col.names
      }
    } else {
      ID.list <- list(col.names)
    }
  } else {
    if(all(!is.null(cols.fn))) {
      cols.fn <- find.file(cols.fn,dir$ids)  # column id file name
      cat("Reading column and row names...\n")
      cat(paste(" reading column names from",cols.fn,"\n"),sep="")
      ID.list <- lapply(cols.fn,readLines)  #readLines(cols.fn)
    } else {
      cat("no column names specified\n"); spec.cn <- F
      ll <- max(1,length(cols.fn))
      ID.list <- vector("list",ll)
      qmf <- quick.mat.format(input.fn[1])
      if(qmf$ncol==1 | long)  { stop("If reading a long format file, must provide row and column names directly, or via filenames") }
      ecn <- qmf$colnames
      ern <- qmf$rownames; #prv("ern")
      if(ecn) {
        for (cc in 1:ll) { ID.list[[cc]] <- quick.mat.format(input.fn[cc])$cnames }
      } else {
        for (cc in 1:ll) { ID.list[[cc]] <- paste("col",1:file.ncol(input.fn[cc],excl.rn=ern),sep="") }
      }
    }
  }
  ##print(headl(ID.list))
  cmb.ID.list <- paste(unlist(do.call("c",ID.list)))
  if(anyDuplicated(cmb.ID.list)) { stop("Cannot have duplicated column names") }
  ##print(length(ID.list[[1]]))
  num.c.fls <- length(ID.list)
  ### GET THE ORDERED row LIST ###
  if((length(row.names)>length(input.fn)) | is.list(row.names)) {
    if(is.list(row.names))  {
      if(all(sapply(row.names,is)[1,]=="character")) {
        rows.list <- row.names
      }
    } else {
      rows.list <- list(row.names)
    }
  } else {
    if(all(!is.null(rows.fn))) {
      rows.fn <- find.file(rows.fn,dir$ano,dir)  # row annotation file name
      cat(paste(" reading row names from",rows.fn,"\n"),sep="")
      rows.list <- lapply(rows.fn,readLines)  #readLines(cols.fn)
    } else {
      cat("no row names specified\n"); spec.rn <- F
      ll <- max(1,length(cols.fn))
      rows.list <- vector("list",ll)
      ecn <- as.numeric(quick.mat.format(input.fn[1])$colnames)
      for (cc in 1:ll) { 
        nr <- (file.nrow(input.fn[cc])-ecn)
        rows.list[[cc]] <- paste("row",1:nr,sep="") 
      }
    }
  }
  cmb.row.list <- paste(unlist(do.call("c",rows.list)))
  if(anyDuplicated(cmb.row.list)) { stop("Cannot have duplicated row names") }
  if(length(rows.list)>1 & length(ID.list)>1) {
    warning("cannot enter both multiple row and column file names")
    return(NULL)
  }
  num.r.fls <- length(rows.list)
  numfls <- num.r.fls*num.c.fls # one of these should be 1
  if(length(rows.list)>1) {
    multi.mode <- T; col.mode <- F
  } else {
    if(length(ID.list)>1) {
      multi.mode <- T; col.mode <- T
    } else {
      multi.mode <- F; col.mode <- T
    }
  }
  # determine what set of input files have been specified
  if(length(input.fn)>1) {
    if(length(input.fn)==numfls) {
      if(verbose) {
        if(col.mode) {
          cat(paste("reading a single cohort from",numfls,"source files.\n"))
        } else {
          cat(paste("reading a single varset from",numfls,"source files.\n"))
        }
      }
    } else {
      stop("Error: when reading from multiple source files, need an equal number of row or col id files")
    }
  } else { if(numfls!=1) { warning(paste("length of ID list was",numfls,"but only 1 input file")) } } 
  #### DETERMINE FILE DIMENSIONS ####
  num.sub <- length(cmb.ID.list) #ID.list)
  smp.szs <- sapply(ID.list,length)
  num.row <- length(cmb.row.list)
  row.szs <- sapply(rows.list,length)
  #prv(c("num.sub","num.row","cmb.ID.list","cmb.row.list"))
  if(col.mode) {
    fil.ofs <- c(0,cumsum(smp.szs)) #note last element is the end of the last file
  } else {
    fil.ofs <- c(0,cumsum(row.szs))
  }
  cat(paste(" found",num.sub,"column names and",num.row,"marker names\n"))
  cls <- num.sub; rws <- num.row
  #cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
  #memory.estimate <- as.double((as.double(rws)*as.double(cls))/cells.per.gb)
  em <- estimate.memory(c(rws,cls))
  if (em > hd.gb) {
    stop(paste("Insufficient disk space availability (hd.gb) expected for this import method\n",
     "Please free up some disk space and try again",sep=""))
  } else {
    divs <- em/(ram.gb/2) # aim to flush once memory reaches half of RAM
    if(divs<1) { do.fl <- F } else { do.fl <- T }
    divs <- round(divs)
  }
  # use 'pref' as the name of the big.matrix backing files for this cohort
  if(pref=="") { pref <- rmv.ext(basename(input.fn))[1] }
  bck.fn <- paste(pref,"bck",sep=".")
  des.fn <- paste(pref,"dsc",sep=".")
  rda.fn <- paste(pref,"RData",sep=".")
  ### DELETE EXISTING FILE IF HAS SAME NAME ###
  if ((!des.fn %in% list.files(dir$big)) | delete.existing )
  {
    if(delete.existing & (des.fn %in% list.files(dir$big)))
    {
      dfn <- cat.path(dir$big,des.fn)
      cat("\n deleting",dfn,"\n")
      unlink(dfn)
    } else {
      #all clear, no files already exist with same name
    }
  } else {
    cat(paste("\nWarning: Big matrix description file",des.fn,"already exists in",dir$big,"\n"))
    cat("You may wish to delete, rename or move this file, or use option 'delete.existing'=T, before re-running this script\n")
    #stop()
  }
  ### MAIN IMPORT LOOP ###
  cat("\nCreating big matrix object to store group data")
  cat("\n predicted disk use: ",round(em,1),"GB\n")
  if(is.character(dir.force.slash(dir$big))) { if(dir$big=="") { dir$big <- getwd() } }
  if(delete.existing & file.exists(cat.path(dir$big,bck.fn))) { unlink(cat.path(dir$big,bck.fn)) }
  bigVar <- big.matrix(nrow=rws,ncol=cls, backingfile=bck.fn, dimnames=list(cmb.row.list,cmb.ID.list),
                       type=dat.type, backingpath=dir.force.slash(dir$big),
                       descriptorfile=des.fn)
  for(ff in 1:numfls) {
    ifn <- cat.path(dir$col,input.fn[ff],must.exist=T)
    if(col.mode) { ffc <- ff; ffr <- 1 } else { ffc <- 1; ffr <- ff }
    # previously: create file name depending on whether in baf or lrr directory
    #test file type if matrix
    #preview(c("ff","numfls","ifn","col.mode","ffr","ffc"))
    if(long) { input.is.vec <- T } else {
      ncifn <- file.ncol(ifn); if(is.na(ncifn)) { stop("input file ",ifn," was empty") }
      if(ncifn>1) { input.is.vec <- F } else { input.is.vec <- T }
    }
    nxt.rng <- (fil.ofs[ff]+1):(fil.ofs[ff+1])
    #preview("nxt.rng")
    if(col.mode) { cls1 <- length(nxt.rng); rws1 <- rws } else { cls1 <- cls; rws1 <- length(nxt.rng) }
    #preview(c("cls1","rws1","cls","rws"))
    if(!input.is.vec) {
      if(spec.rn & spec.cn) {
        frm <- check.text.matrix.format(fn=ifn,ncol=cls1,header=ID.list[[ffc]],row.names=rows.list[[ffr]])
        #prv(c("cls1","rws1","cls","frm"))
      } else {
        frm <- check.text.matrix.format(fn=ifn)
        if(!is.null(frm$colnames)) { ID.list[[ff]][1:length(frm$colnames)] <- frm$colnames }
      }
      if(frm$rname & !frm$match) { file.rn <- character(rws) } # recording rownames as we go
      if(frm$match & !col.mode) { stop("cannot use matching method with separate files by rows") }
    }
    dat.file <- file(ifn)
    open(con=dat.file,open="r")
    cat(paste(" opening connection to ",c("matrix","long")[1+input.is.vec],
              " format datafile (",ff,"/",numfls,"): ",basename(ifn),"\n",sep=""))
    cat("\nLoading text data into big matrix object:\n")
    if(!input.is.vec)
    {
      ## read from matrix format tab file
      twty.pc <- round(rws1/divs) # flush data every 'n' rows
      for (cc in 1:rws1) {
        if (do.fl & (cc %% twty.pc)==0)  { fl.suc <- bigmemory::flush(bigVar) ; if(!fl.suc) { cat("flush failed\n") } }
        if(tracker) { loop.tracker(cc,rws1) }
        next.line <- readLines(dat.file,n=1)
        next.row <- strsplit(next.line,"\t",fixed=T)[[1]]
        if (cc==1 & frm$header) { 
          # need to go to next line to avoid reading header as data
          next.line <- readLines(dat.file,n=1); next.row <- strsplit(next.line,"\t",fixed=T)[[1]]
        }
        if (frm$rnames) {
          if(frm$match) {
            lbl <- next.row[1]; 
            row.sel <- match(lbl,cmb.row.list)
            bigVar[row.sel,nxt.rng] <- next.row[-1]
            file.rn[row.sel] <- lbl
          } else {
            if(col.mode) {
              #prv("bigVar",counts=list(cc=cc,nxt.rng=nxt.rng[1]))
              #prv(c("next.row","nxt.rng"))
              file.rn[cc] <- next.row[1]; bigVar[cc,nxt.rng] <- next.row[-1]
            } else {
              selc <- 1:(length(next.row)-1)
              file.rn[nxt.rng[cc]] <- next.row[1]; bigVar[nxt.rng[cc],selc] <- next.row[-1]
            }
          }
        } else {
          if(col.mode) {
            bigVar[cc,nxt.rng] <- next.row
          } else {
            bigVar[nxt.rng[cc],] <- next.row
          }
        }
      }
    } else {
      ## read from (long) vector format tab file
      if(col.mode) {
        ## col-wise file splits
        twty.pc <- round(cls1/divs) # flush data every 'n' cols
        for (cc in 1:cls1) {
          if(tracker) { loop.tracker(cc,cls1) }
          if (do.fl & (cc %% twty.pc)==0)  { fl.suc <- bigmemory::flush(bigVar) ; if(!fl.suc) { cat("flush failed\n") } }
          bigVar[,(cc+fil.ofs[ff])] <- as(readLines(dat.file,n=rws),as.type)
        }
      } else {
        ## row-wise file splits
        twty.pc <- round(rws1/divs)
        for (cc in 1:cls) {
          if(tracker) { loop.tracker(cc,cls) }
          if (do.fl &(cc %% twty.pc)==0)  { fl.suc <- bigmemory::flush(bigVar) ; if(!fl.suc) { cat("flush failed\n") } }
          bigVar[nxt.rng,cc] <- as(readLines(dat.file,n=rws1),as.type)
        }
      }
    }
    close(dat.file)
  }
  ### FINISH UP, RETURN BIGMATRIX DESCRIPTION ###
  # in the case of different rownames found in matrix, then show following warning text:
  if(!input.is.vec) {
    if(frm$rname & !frm$match) {
      if(!spec.rn) {
        if(nrow(bigVar)==length(file.rn)) {
          options(bigmemory.allow.dimnames=TRUE)
          rownames(bigVar) <- file.rn; cat("updated big.matrix rownames from names in file(s)\n")
          bigmemory::flush(bigVar)
          big.des <- bigmemory::describe(bigVar)
          des.fn <- cat.path(dir=dir$big,fn=des.fn,ext="RData")
          save(big.des,file=des.fn)
          warning("Had to change description file to a binary file to update rownames. This can be read in with get.big.matrix() [and should be faster to load]")
        } else {
          #prv(c("file.rn"))
          miswarn <- T
        }
      } else {
        #print("gere")
      }
    } else {
      if(frm$match) {
        miswarn <- T
      } else {
        #print("hereg")
      }
    }
  }
  if(miswarn) {
    ofn <- cat.path(dir$ano,pref=pref,"_file_rowname_list_check_this.txt")
    if(exists("file.rn")) { 
      warning("rownames didn't match what was in filecheck the list in the file at:\n ",ofn) 
      writeLines(paste(file.rn),con=ofn) ; cat("\n file preview:\n"); print(head(file.rn,10)); cat("\n")
    } else {
      warning(paste("rownames didn't match what was in file"))
    }
  }
  cat("\n")
  cat(paste(" created big.matrix description file:",des.fn,"\n"))
  cat(paste(" created big.matrix backing file:",bck.fn,"\n"))
  if(ret.obj) {
    return(bigmemory::describe(bigVar))
  } else {
    # convert to RData file regardless as the old descriptor files are now very slow
    big.des <- bigmemory::describe(bigVar)
    des.fn <- cat.path(dir=dir$big,fn=des.fn,ext="RData")
    save(big.des,file=des.fn)
    return(des.fn)
  }
  cat("...complete!\n")
}


# mean replacement code internal
row.rep <- function(X) { X[is.na(X)] <- mean(X,na.rm=T); X }


#internal, analog of 'select.samp.snp.custom' from older versions of plumbCNV
select.col.row.custom <- function(bigMat,row,col,verbose=T)
{
  # based on files/vectors of row-ids and column-ids create selection
  # vectors to select only the ids in these lists for a matrix
  if(verbose) { cat(" calculating selections for rows\n") }
  # try to detect whether a vector of IDs, or file names
  row.ref <- rownames(bigMat)  ; col.ref <- colnames(bigMat) 
  row.sel <- col.sel <- NULL
  byname <- T
  if (length(row)==1 & length(col)==1 & is.character(row) & is.character(col))
  {
    if(verbose) { cat(" [assuming 'col' and 'row' are file names containing column and row ids]") }
    if(file.exists(row)) {
      row.sel <- readLines(row)
    } else {
      if(row=="") {
        if(verbose) { cat(c(" row subset file was empty, selecting all\n")) }
        row.sel <- row.ref
      } else {
        stop("Error: argument 'row' should be a vector of rows names length>1 or a filename with a list of rows (no header)")
      }
    }
    if(file.exists(col)) {
      col.sel <- readLines(col)
    } else {
      if(col=="") {
        if(verbose) { cat(c(" column subset file was empty, selecting all\n")) }
        col.sel <- col.ref
      } else {
        stop("Error: argument 'col' should be a vector of column names length>1 or a filename with a list of rows (no header)")
      }
    }
  } else { 
    if(is.character(row) | is.character(col)) {
      #cat("[assuming 'col' and 'row' are vectors of column and row ids]")
      # if blank then assign all ids
      if(all(row=="")) {
        row.sel <- row.ref
      } else {
        row.sel <- row
      }
      if(all(col=="")) {
        col.sel <- col.ref
      } else { 
        col.sel <- col  
      }
    } else {
      byname <- F
      # assume numeric or logical selection
      if(is.logical(row)) { row <- which(row) }
      if(is.logical(col)) { col <- which(col) }
      to.order.r <- row[row<=nrow(bigMat)]
      to.order.c <- col[col<=ncol(bigMat)]
    }
  }
  # use sort/exclusion lists to get reordering vectors
  #row.sel <- row.sel ; col.sel <- col.sel
  
  #print(head(row.sel));print(head(col.sel))
  if(byname) {
    # selection is based on row/col names
    to.order.r <- narm(match(row.sel,rownames(bigMat)))
    to.order.c <- narm(match(col.sel,colnames(bigMat)))
  } else {
    # selection is based on row/col numbers or logical
  }
  if (!(length(to.order.r[!is.na(to.order.r)])>0 & length(to.order.c[!is.na(to.order.c)])>0))
  { warning("selection of rows and/or columns has resulted in an empty dataset",
            "\ncheck rownames, column names and selection lists for errors") }
  
  out.list <- list(to.order.r,to.order.c,row.sel,col.sel)
  names(out.list) <- c("to.order.r","to.order.c","row.list","column.list")
  return(out.list)
}


#' Reduce one dimension of a large matrix in a strategic way
#' 
#' Thin the rows (or columns) of a large matrix or big.matrix in order to reduce the size of the
#' dataset while retaining important information. Percentage of the original size or a new number 
#' of rows/columns is selectable, and then there are four methods to choose the data subset.
#' Simple uniform and random selection can be specified. Other methods look at the correlation
#' structure of a subset of the data to derive non-arbitrary selections, using correlation, PCA,
#' or association with a phenotype or some other categorical variable. Each of the four methods
#' has a separate function in this package, which you can see for more information, this function
#' is merely a wrapper to select one of the four.
#' @param bigMat a big.matrix object, or any argument accepted by get.big.matrix(), which includes
#'  paths to description files or even a standard matrix object.
#' @param keep numeric, by default a proportion (decimal) of the original number of rows/columns to choose
#'  for the subset. Otherwise if an integer>2 then will assume this is the size of the desired subset,
#'  e.g, for a dataset with 10,000 rows where you want a subset size of 1,000 you could set 'keep' as
#'  either 0.1 or 1000.
#' @param how character, only the first two characters are required and they are not case sensitive,
#'  select what method to use to perform subset selection, options are:
#' 'uniform': evenly spaced selection when random=FALSE, or random selection otherwise;
#'  see uniform.select().
#' 'correlation': most correlated subset when hi.cor=TRUE, least correlated otherwise;
#'  see subcor.select().
#' 'pca': most representative variables of the principle components of a subset;
#'  see subpc.select().
#' 'association': most correlated subset with phenotype if least=FALSE, or least correlated otherwise;
#'  see select.least.assoc().
#' @param dir directory containing the filebacked.big.matrix, same as 'dir' for get.big.matrix.
#' @param rows logical, whether to choose a subset of rows (TRUE), or columns (FALSE). rows is always
#'  TRUE when using 'association' methods.
#' @param random logical, whether to use random selections and subsets (TRUE), or whether to use uniform
#'  selections that should give the same result each time for the same dataset (FALSE)
#' @param hi.cor logical, if using 'correlation' methods, then whether to choose the most correlated (TRUE)
#'  or least correlated (FALSE).
#' @param least logical, if using 'association' methods, whether to choose the least associated (TRUE) or 
#'  most associated variables with phenotype
#' @param pref character, a prefix for big.matrix backing files generated by this selection
#' @param verbose logical, whether to display more information about processing
#' @param ret.obj logical, whether to return the result as a big.matrix object (TRUE), or as a reference
#'  to the binary file containing the big.matrix.descriptor object [either can be read with get.big.matrix() or
#'  prv.big.matrix()]
#' @param ... other arguments to be passed to uniform.select, subpc.select, subcor.select, or select.least.assoc
#' @return A smaller big.matrix with fewer rows and/or columns than the original matrix
#' @export
#' @seealso \code{\link{uniform.select}}, \code{\link{subpc.select}}, \code{\link{subcor.select}}, 
#' \code{\link{select.least.assoc}}, \code{\link{big.select}}, \code{\link{get.big.matrix}}
#' @author Nicholas Cooper 
#' @examples 
#' bmat <- generate.test.matrix(5,big.matrix=TRUE)
#' prv.big.matrix(bmat)
#' # make 5% random selection:
#' lmat <- thin(bmat)
#' prv.big.matrix(lmat)
#' # make 10% most orthogonal selection (lowest correlations):
#' lmat <- thin(bmat,.10,"cor",hi.cor=FALSE)
#' prv.big.matrix(lmat)
#' # make 10% most representative selection:
#' lmat <- thin(bmat,.10,"PCA",ret.obj=FALSE) # return file name instead of object
#' print(lmat)
#' prv.big.matrix(lmat)
#' # make 25% selection most correlated to phenotype
#' # create random phenotype variable
#' pheno <- rep(1,ncol(bmat)); pheno[which(runif(ncol(bmat))<.5)] <- 2
#' lmat <- thin(bmat,.25,"assoc",phenotype=pheno,least=FALSE,verbose=TRUE)
#' prv.big.matrix(lmat)
#' # tidy up temporary files:
#' unlink(c("thin.bck","thin.dsc","thin.RData"))
thin <- function(bigMat,keep=0.05,how=c("uniform","correlation","pca","association"),
                 dir="",rows=TRUE,random=TRUE,hi.cor=TRUE,least=TRUE,
                 pref="thin",verbose=FALSE,ret.obj=TRUE,...) {
    how <- toupper(substr(how[1],1,2)); meth <- "none"
    if(is.character(dir)) { if(dir=="") { dir <- getwd() } } else { dir <- getwd() }
    if(how=="UN") {
      meth <- "uniform"
      rc <- uniform.select(bigMat,keep=keep,dir=dir,rows=rows,random=random,...)
      if(rows) { subrc <- rc[[1]] } else { subrc <- rc[[2]] }
    }
    if(how=="CO") {
      meth <- "correlation"
      #hi.cor is the extra parameter for this
      subrc <- subcor.select(bigMat,keep=keep,dir=dir,rows=rows,random=random,hi.cor=hi.cor,...)
    }
    if(how=="PC") {
      meth <- "pca"
      subrc <- subpc.select(bigMat,keep=keep,dir=dir,rows=rows,random=random,...)
    }
    if(how=="AS") {
      meth <- "association"
      rows <- T
      #least is the extra parameter for this
      test.args <- list(...)
      if(!"phenotype" %in% names(test.args)) { stop("must include argument 'phenotype' when how='association'") }
      subrc <- select.least.assoc(bigMat,keep=keep,dir=dir,verbose=verbose,least=least,...)
    }
    #prv(subrc)
    if(meth=="none") { stop("invalid subsetting method (",how,") specified") }
    if(rows) { sr <- subrc; sc <- 1:ncol(bigMat) } else { sr <- 1:nrow(bigMat); sc <- subrc }
    bigSubMat <- big.select(bigMat, select.rows=sr, select.cols=sc, dir=dir, 
                           deepC=T, pref=pref, verbose=verbose , delete.existing=T)
    if(ret.obj) {
      return(get.big.matrix(bigSubMat))
    } else {
      return(bigSubMat)
    }
}




#' Select a subset of a big.matrix
#' 
#' Select a subset of big.matrix using indexes for a subset of rows and columns.
#' Essentially a wrapper for bigmemory::deepcopy, but with slightly more flexible
#' parameters. bigMat can be entered in any form accepted by get.big.matrix(), row and
#' column selections can be vectors of indexes, names or file.names containing indexes.
#' Default is to process using deepcopy, but processing without using bigmemory native
#' methods is a faster option when matrices are small versus available RAM. File names
#' for backing files are managed only requiring you to enter a prefix, or optionally
#' use the default and gain filebacked functionality without having to bother choosing
#' filename parameters.
#' @param bigMat a big.matrix, matrix or any object accepted by get.big.matrix()
#' @param select.rows selection of rows of bigMat, can be numbers, logical, rownames, or a file with names. 
#'  If using a filename argument, must also use a filename argument for select.cols (cannot mix)
#' @param select.cols selection of columns of bigMat, can be numbers, logical, colnames, or a file with names
#' @param dir the directory containing the bigMat backing file (e.g, parameter for get.big.matrix()).
#' @param deepC logical, whether to use bigmemory::deepcopy, which is slowish, but scalable, or 
#'  alternatively to use standard indexing which converts the result to a regular matrix object,
#'  and is fast, but only feasible for matrices small enough to fit in memory.
#' @param pref character, prefix for the big.matrix backingfile and descriptorfile, and optionally
#'  an R binary file containing a big.matrix.descriptor object pointing to the big.matrix result.
#' @param verbose whether to display extra information about processing and progress
#' @param delete.existing logical, if a big.matrix already exists with the same name as implied
#'  by the current 'pref' and 'dir' arguments, then default behaviour (FALSE) is to return an error.
#'  to overwrite any existing big.matrix file(s) of the same name(s), set this parameter to TRUE.
#' @return A big.matrix with the selected (in order) rows and columns specified
#' @export
#' @author Nicholas Cooper 
#' @examples 
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' bmat <- generate.test.matrix(5,big.matrix=TRUE)
#' # take a subset of the big.matrix without using deepcopy
#' sel <- big.select(bmat,c(1,2,8),c(2:10),deepC=FALSE,verbose=TRUE)
#' prv.big.matrix(sel)
#' # now select the same subset using row/column names from text files
#' writeLines(rownames(bmat)[c(1,2,8)],con="bigrowstemp.txt")
#' writeLines(colnames(bmat)[c(2:10)],con="bigcolstemp.txt")
#' sel <- big.select(bmat, "bigrowstemp.txt","bigcolstemp.txt", delete.existing=TRUE)
#' prv.big.matrix(sel)
#' unlink(c("bigcolstemp.txt","bigrowstemp.txt","sel.RData","sel.bck","sel.dsc"))
#' setwd(orig.dir) # reset working dir to original
big.select <- function(bigMat, select.rows=NULL, select.cols=NULL, dir=getwd(), 
                       deepC=TRUE, pref="sel", delete.existing=FALSE, verbose=FALSE)
{
  # sort and exclude snps/samples from a big.matrix
  if(exists("validate.dir.for",mode="function")) {
    ## plumbCNV specific code ##
    dir <- do.call("validate.dir.for",list(dir=dir,elements=c("big","pc"),warn=F))  
  } else {
    # otherwise
    dir <- list(big=dir,pc=dir)
  }
  # bigmatrix file names for re-ordered filtered matrix (which is the final output of this script)
  if(is.character(bigMat) & length(bigMat)==1) { Fn <- gsub("descrFile","",bigMat) } else { Fn <- "" }
  bck.fn.o <- cat.path(fn=Fn,pref=pref,ext="bck")
  des.fn.o <- cat.path(fn=Fn,pref=pref,ext="dsc")
  R.descr <- cat.path(dir$big,rmv.ext(des.fn.o),ext=".RData"); print(rmv.ext(des.fn.o))
  
  bigMat <- get.big.matrix(bigMat,dir)
  if(verbose) { cat(paste(" attached matrix with dims:",paste(dim(bigMat),collapse=","),"\n")) }
  # get list of deleting/reordering vectors using annotation files
  trans.list <- select.col.row.custom(bigMat,row=select.rows,col=select.cols,verbose=verbose)
  if(verbose) { cat(paste(" selected",length(trans.list[[2]]),"listed samples and",length(trans.list[[1]]),"variables\n")) }
  wrn <- "trans.list was already attached, detaching now..\n"
  #while("trans.list" %in% search()) { detach(trans.list); warning(wrn) }
  #attach(trans.list)
  to.order.r <- trans.list$to.order.r; to.order.c <- trans.list$to.order.c
  if(verbose) {
    cat("\nReordering Variables and Samples...\n")
    cat("\nINDEXES SUMMARY\n")
    cat(paste(length(to.order.r),"row indexes range is from",min(to.order.r),"to",max(to.order.r),"\n"))
    cat("-->",head(to.order.r),sep=", "); cat ("\n")
    cat(paste(length(to.order.c),"col indexes range is from",min(to.order.c),"to",max(to.order.c),"\n"))
    cat("-->",head(to.order.c),sep=", ")
    cat("\n\n raw big.matrix summary before selection/ordering:\n\n")
    prv.big.matrix(bigMat,"bigMat")
  }
  if(!deepC)
  {
    # this is fast with available RAM (like 20 secs)
    if(verbose) { cat(" running reorder in system memory\n") }
    bigMat1 <- bigMat[to.order.r,to.order.c]
    if(verbose) {
      cat(" adding colnames\n") ; colnames(bigMat1) <- colnames(bigMat)[to.order.c] 
      cat(" adding rownames\n") ; rownames(bigMat1) <- rownames(bigMat)[to.order.r] 
      cat(" converting matrix to big.matrix\n") 
    }
    if(delete.existing & file.exists(cat.path(dir$big,basename(bck.fn.o)))) { unlink(cat.path(dir$big,basename(bck.fn.o))) }
    bigMat2 <- as.big.matrix(bigMat1, backingfile=basename(bck.fn.o),
                             descriptorfile=basename(des.fn.o),backingpath=dir$big)
    if(verbose) { cat(paste(" matrix descr saved as standard description file:",des.fn.o,"\n")) }
    descr <- bigmemory::describe(bigMat2)
  } else {
    #this is slow but creates backing file and will speed up ops later
      ### DELETE EXISTING FILE IF HAS SAME NAME ###
      if ((!basename(des.fn.o) %in% list.files(dir$big)) | delete.existing )
      {
        if(delete.existing & (basename(des.fn.o) %in% list.files(dir$big)))
        {
          dfn <- cat.path(dir$big,basename(des.fn.o))
          cat("\n deleting",dfn,"\n")
          unlink(dfn)
        } else {
          #all clear, no files already exist with same name
        }
      } else {
        cat(paste("\nWarning: Big matrix description file",basename(des.fn.o),"already exists in",dir$big,"\n"))
        cat("You may wish to delete, rename or move this file, or use option 'delete.existing'=T, before re-running this script\n")
        #stop()
      }
    if(verbose) { cat(" starting deep copy...") }
    if(delete.existing & file.exists(cat.path(dir$big,basename(bck.fn.o)))) { unlink(cat.path(dir$big,basename(bck.fn.o))) }
    bigMat2 <- deepcopy(bigMat, cols = to.order.c, rows = to.order.r,
                        backingfile=basename(bck.fn.o),backingpath=dir$big,
                        descriptorfile=basename(des.fn.o) )
    if(verbose) { cat("done\n") }
    if(verbose) { cat("\nAdding names\n") }
    options(bigmemory.allow.dimnames=TRUE)
    colnames(bigMat2) <- colnames(bigMat)[to.order.c]
    if(verbose) { cat(" added colnames\n") }
    rownames(bigMat2) <- rownames(bigMat)[to.order.r]  
    if(verbose) { cat(" added rownames\n") }
    descr <- bigmemory::describe(bigMat2)
    bigmemory::flush(bigMat2) # hopefully this will ensure the row/colnames are added to the file backing
    if(verbose) { cat(paste(" due to use of deep copy option, recommend only to use descr saved as rbinary description file\n")) }
  }
  save(descr,file=R.descr)
  if(verbose) {
    cat(paste(" created big.matrix description file:",des.fn.o,"\n"))
    cat(paste(" created big.matrix backing file:",bck.fn.o,"\n"))
    cat(paste(" created big.matrix binary description file:",basename(R.descr),"\n"))
  } 
  #while("trans.list" %in% search()) { detach(trans.list) }
  #return(descr) 
  return(R.descr)
}



#' Selection of a representative variable subset
#' 
#' Returns a subset (size='keep') of row or column numbers that are most representative of a dataset.
#' This function performs PCA on a small subset of columns and all rows (when rows=TRUE, or vice
#'  -versa when rows=FALSE), and selects rows (rows=TRUE) most correlated to the first 'n' principle
#' components, where 'n' is chosen by the function quick.elbow(). The number of variables selected
#' corresponding to each component is weighted according to how much of the variance is explained
#' by each component.
#' @param bigMat a big.matrix, matrix or any object accepted by get.big.matrix()
#' @param keep numeric, by default a proportion (decimal) of the original number of rows/columns to choose
#'  for the subset. Otherwise if an integer>2 then will assume this is the size of the desired subset,
#'  e.g, for a dataset with 10,000 rows where you want a subset size of 1,000 you could set 'keep' as
#'  either 0.1 or 1000.
#' @param rows logical, whether the subset should be of the rows of bigMat. If rows=FALSE, then 
#'  the subset is chosen from columns, would be equivalent to calling subpc.select(t(bigMat)),
#'  but avoids actually performing the transpose which can save time for large matrices.
#' @param dir the directory containing the bigMat backing file (e.g, parameter for get.big.matrix()).
#' @param random logical, passed to uniform.select(), whether to take a random or uniform selection
#'  of columns (or rows if rows=F) to run the subset PCA.
#' @param ram.gb maximum size of the matrix in gigabytes for the subset PCA, 0.1GB is the default 
#'  which should result in minimal processing time on a typical system. Increasing this
#'  increases the processing time, but also the representativeness of the subset chosen. Note
#'  that some very large matrices will not be able to be processed by this function unless 
#'  this parameter is increased; basically if the dimension being thinned is more than 5% of
#'  this memory limit (see estimate.memory() from NCmisc).
#' @param ... further parameters to pass to big.PCA() which performs the subset PCA used to 
#'  determine the most representative rows (or columns).
#' @return A set of row or column indexes (depents on 'rows' parameter) of the most representative
#'  variables in the matrix, as defined by most correlated to principle components
#' @export
#' @seealso \code{\link{thin}}, \code{\link{uniform.select}}, \code{\link{big.PCA}}, \code{\link{get.big.matrix}}
#' @author Nicholas Cooper 
#' @examples 
#' mat <- matrix(rnorm(200*2000),ncol=200) # normal matrix
#' bmat <- as.big.matrix(mat)              # big matrix
#' ii <- subpc.select(bmat,.05,rows=TRUE) # thin down to 5% of the rows
#' ii <- subpc.select(bmat,45,rows=FALSE) # thin down to 45 columns
#' # show that rows=T is equivalent to rows=F of the transpose (random must be FALSE)
#' ii1 <- subpc.select(mat,.4,rows=TRUE,random=FALSE)
#' ii2 <- subpc.select(t(mat),.4,rows=FALSE,random=FALSE)
#' print(all.equal(ii1,ii2))
subpc.select <- function(bigMat,keep=.05,rows=TRUE,dir=getwd(),random=TRUE,ram.gb=0.1,...) {  
  # select a subset of variables based on the most representative variables in the PCs of a subset
  bigMat <- get.big.matrix(bigMat,dir=dir)
  if(rows) { N <- nrow(bigMat) } else { N <- ncol(bigMat) }
  if(keep>2) {
    new.n <- round(keep)
  } else {
    new.n <- round(max(0,min(1,keep))*N)
  }
  min.other <- min(c(dim(bigMat)/2,20)) # e.g, minimum of 20 sample subset to test with very large variable set
  max.other <- min(c(dim(bigMat)/2,200)) # a good amount to test on
  if(estimate.memory(c(N,min.other))>ram.gb) {
    ## too big to do this
    warning("The matrix selected has too many variables to do a mini-PCA on a meaningful subset. Suggest choosing an alternative subset method")
    return(NULL)
  } 
  if(estimate.memory(c(N,max.other))<ram.gb) {
    # go with max
    keeper <- max.other
  } else {
    # go with best we can up to that memory point
    mult <- estimate.memory(c(N,max.other))/ram.gb
    keeper <- (max.other/mult)
  }
  #print(keeper)
  rc <- uniform.select(bigMat,keep=keeper,rows=!rows,dir=dir,random=random)
  ## do it here
  sub.mat <- bigMat[rc[[1]],rc[[2]]]
  #prv.large(sub.mat)
  pt <- "package:"; pkgset <- gsub(pt,"",search()[grep(pt,search())])
  uba <- (all(c("irlba","bigalgebra") %in% pkgset))
  if(rows) {
    quick.pc <- big.PCA(sub.mat,pcs.to.keep=NA,use.bigalgebra=uba,...)
  } else {
    quick.pc <- big.PCA(t(sub.mat),pcs.to.keep=NA,use.bigalgebra=uba,...)
  }
  el <- quick.elbow(quick.pc$Evalues)
  varpcs <- estimate.eig.vpcs(eigenv=quick.pc$Evalues,
        M=sub.mat,elbow=el,print.est=F,estimated=T,ignore.warn=T)$variance.pcs
  el <- quick.elbow(varpcs)
  pc <- quick.pc$PCs[,1:el]; 
  #print(dim(pc))
  if(ncol(sub.mat)==nrow(pc)) {
    cormat <- abs(cor(pc,t(sub.mat)))
  } else {
    cormat <- abs(cor(pc,sub.mat))
  }
  #return(cormat)
  #print(dim(cormat))
  selected <- logical(ncol(cormat)) # number of variables long logical
  new.set <- NULL
  vpc <- varpcs[1:el]/sum(varpcs[1:el]) # to prevent vagaries of roundoff error
  # make sure total number to choose will be exactly right
  if(random) {
    per.pc <- round(new.n*vpc)
    while(sum(per.pc)>new.n) { rr <- sample(1:length(per.pc),1); per.pc[rr] <- max(0,per.pc[rr]-1) }
    while(sum(per.pc)<new.n) { rr <- sample(1:length(per.pc),1); per.pc[rr] <- max(0,per.pc[rr]+1) }
  } else {
    brks <- quantile(1:new.n,cumsum(rev(sort(vpc))))
    per.pc <- table(findInterval(1:new.n,brks,rightmost.closed=T)+1)
  }
  for (cc in 1:nrow(cormat)) {
    nxt.row <- cormat[cc,]; indc <- 1:ncol(cormat)
    nxt.row <- nxt.row[!selected]; indc <- indc[!selected]
    if(per.pc[cc]>0) {
      chc <- (indc[rev(order(nxt.row))])[1:per.pc[cc]]
      #print((nxt.row[rev(order(nxt.row))])[1:per.pc[cc]])
      new.set <- c(new.set,chc)
      selected[new.set] <- T
    }
  }  
  return(sort(new.set))
}  


#' Selection of the most correlated variable subset
#' 
#' Returns a subset (size='keep') of row or column numbers that are most correlated to other
#' variables in the dataset (or if hi.cor=F), then those that are least correlated.
#' This function performs cor() on a small subset of columns and all rows (when rows=TRUE, or vice
#'  -versa when rows=FALSE), and selects rows (rows=TRUE) with greatest/least absolute sum of correlations.
#' @param bigMat a big.matrix, matrix or any object accepted by get.big.matrix()
#' @param keep numeric, by default a proportion (decimal) of the original number of rows/columns to choose
#'  for the subset. Otherwise if an integer>2 then will assume this is the size of the desired subset,
#'  e.g, for a dataset with 10,000 rows where you want a subset size of 1,000 you could set 'keep' as
#'  either 0.1 or 1000.
#' @param rows logical, whether the subset should be of the rows of bigMat. If rows=FALSE, then 
#'  the subset is chosen from columns, would be equivalent to calling subpc.select(t(bigMat)),
#'  but avoids actually performing the transpose which can save time for large matrices.
#' @param hi.cor logical, whether to choose the most correlated (TRUE) or least correlated subset (FALSE).
#' @param dir the directory containing the bigMat backing file (e.g, parameter for get.big.matrix()).
#' @param random logical, passed to uniform.select(), whether to take a random or uniform selection
#'  of columns (or rows if rows=FALSE) to run the subset PCA.
#' @param ram.gb maximum size of the matrix in gigabytes for the subset PCA, 0.1GB is the default 
#'  which should result in minimal processing time on a typical system. Increasing this
#'  increases the processing time, but also the representativeness of the subset chosen. Note
#'  that some very large matrices will not be able to be processed by this function unless 
#'  this parameter is increased; basically if the dimension being thinned is more than 5% of
#'  this memory limit (see estimate.memory() from NCmisc).
#' @return A set of row or column indexes (depents on 'rows' parameter) of the most inter-correlated
#'  (or least) variables in the matrix.
#' @export
#' @seealso \code{\link{thin}}, \code{\link{uniform.select}}, \code{\link{get.big.matrix}}
#' @author Nicholas Cooper 
#' @examples 
#' mat <- matrix(rnorm(200*2000),ncol=200)
#' bmat <- as.big.matrix(mat)
#' ii1 <- subcor.select(bmat,.05,rows=TRUE) # thin down to 5% of the rows
#' ii2 <- subcor.select(bmat,45,rows=FALSE) # thin down to 45 columns
#' prv(ii1,ii2)
#' # show that rows=T is equivalent to rows=F of the transpose (random must be FALSE)
#' ii1 <- subcor.select(mat,.4,rows=TRUE,random=FALSE)
#' ii2 <- subcor.select(t(mat),.4,rows=FALSE,random=FALSE)
#' print(all.equal(ii1,ii2))
subcor.select <- function(bigMat,keep=.05,rows=TRUE,hi.cor=TRUE,dir=getwd(),random=TRUE,ram.gb=0.1) {  
  # select a subset of variables based on the most representative variables in the PCs of a subset
  if(rows) { N <- nrow(bigMat) } else { N <- ncol(bigMat) }
  if(keep>2) {
    new.n <- round(keep)
  } else {
    new.n <- round(max(0,min(1,keep))*N)
  }
  min.other <- min(c(dim(bigMat)/2,20)) # e.g, minimum of 20 sample subset to test with very large variable set
  max.other <- min(c(dim(bigMat)/2,200)) # a good amount to test on
  if(estimate.memory(c(N,min.other))>ram.gb) {
    ## too big to do this
    warning("The matrix selected has too many variables to do a mini-PCA on a meaningful subset. Suggest choosing an alternative subset method")
    return(NULL)
  } 
  if(estimate.memory(c(N,max.other))<ram.gb) {
    # go with max
    keeper <- max.other
  } else {
    # go with best we can up to that memory point
    mult <- estimate.memory(c(N,max.other))/ram.gb
    keeper <- (max.other/mult)
  }
  #print(keeper)
  rc <- uniform.select(bigMat,keep=keeper,rows=!rows,dir=dir,random=random)
  ## do it here
  sub.mat <- bigMat[rc[[1]],rc[[2]]]
  #prv.large(sub.mat)
  #print(dim(pc))
  if(!rows) {
    cormat <- abs(cor(sub.mat))
  } else {
    cormat <- abs(cor(t(sub.mat)))
  }
  sorter <- order(rowSums(cormat))
  if(hi.cor) {
    rank <- rev(sorter)[1:new.n] # select top new.n correlations
  } else {
    rank <- sorter[1:new.n] # select new.n least correlated
  }
  return(sort(rank))
}


#' Derive a subset of a large dataset
#' 
#' Either randomly or uniformly select rows or columns from a large dataset
#' to form a new smaller dataset.
#' @param bigMat a big.matrix object, or any argument accepted by get.big.matrix(), which includes
#'  paths to description files or even a standard matrix object.
#' @param keep numeric, by default a proportion (decimal) of the original number of rows/columns to choose
#'  for the subset. Otherwise if an integer>2 then will assume this is the size of the desired subset,
#'  e.g, for a dataset with 10,000 rows where you want a subset size of 1,000 you could set 'keep' as
#'  either 0.1 or 1000.
#' @param dir directory containing the filebacked.big.matrix, same as dir for get.big.matrix.
#' @param rows logical, whether the subset should be of the rows of bigMat. If rows=FALSE, then 
#'  the subset is chosen from columns, would be equivalent to calling subpc.select(t(bigMat)),
#'  but avoids actually performing the transpose which can save time for large matrices.
#' @param random logical, passed to uniform.select(), whether to take a random or uniform selection
#'  of columns (or rows if rows=FALSE) to run the subset PCA.
#' @param ram.gb maximum size of the matrix in gigabytes for the subset PCA, 0.1GB is the default 
#'  which should result in minimal processing time on a typical system. Increasing this
#'  increases the processing time, but also the representativeness of the subset chosen. Note
#'  that some very large matrices will not be able to be processed by this function unless 
#'  this parameter is increased; basically if the dimension being thinned is more than 5% of
#'  this memory limit (see estimate.memory() from NCmisc).
#' @return A set of row or column indexes (depents on 'rows' parameter) of uniformly distributed
#'  (optionally reproduceable) or randomly selected variables in the matrix.
#' @export
#' @seealso \code{\link{subpc.select}}
#' @author Nicholas Cooper 
#' @examples 
#' mat <- matrix(rnorm(200*100),ncol=200)  # standard matrix
#' bmat <- as.big.matrix(mat)              # big.matrix
#' ii1 <- uniform.select(bmat,.05,rows=TRUE) # thin down to 5% of the rows
#' ii2 <- uniform.select(bmat,45,rows=FALSE,random=TRUE) # thin down to 45 columns
#' prv(ii1,ii2)
## tailor this to make a random uniform selection
uniform.select <- function(bigMat,keep=.05,rows=TRUE,dir="",random=TRUE,ram.gb=0.1) {  
  # select an exactly evenly spaced subset (reproduceable)
  if(rows) { N <- nrow(bigMat) } else { N <- ncol(bigMat) }
  if(keep>2) {
    new.n <- round(keep)
  } else {
    new.n <- round(max(0,min(1,keep))*N)
  }
  if(!random) {
    X <- cut.fac(N,new.n); 
    indx <- tapply(1:N,X,function(x) { round(median(x)) })
    indx <- as.integer(indx)
    # select a randomly uniform subset
  } else {
    indx <- sort(sample(N,new.n))
  }
  if(rows) { 
    outlist <- list(order.r=indx, order.c=1:ncol(bigMat)) 
  } else {
    outlist <- list(order.r=1:nrow(bigMat), order.c=indx)
  }
  return(outlist)
}


#' Quick association tests for phenotype
#' 
#' Simplistic association tests, only meant for purposes of preliminary variable
#' selection or creation of priors, etc. Quickly obtain association p-values for
#' a big.matrix against a list of phenotypes for each row, where columns are 
#' samples and column labels correspond to the rownames of the sample.info dataframe
#' which contains the phenotype information, in a column labelled 'use.col'.
#' @param bigMat a big.matrix object, or any argument accepted by get.big.matrix(), which includes
#'  paths to description files or even a standard matrix object.
#' @param dir directory containing the filebacked.big.matrix, same as dir for get.big.matrix.
#' @param sample.info a data.frame with rownames corresponding to colnames of the bigMat. Must
#'  also contain a column named 'use.col' (default 'phenotype') which contains the categorical
#'  variable to perform the association test for phenotype, etc. This file may contain extra
#'  ids not in colnames(bigMat), although if any column names of bigMat are missing from
#'  sample.info a warning will be given, and the call is likely to give incorrect results.
#' @param use.col the name of the phenotype column in the data.frame 'sample.info'
#' @param p.values logical, whether to return p.values from the associations 
#' @param F.values logical, whether to return F.values from the associations
#' @param n.cores integer, if wanting to process the analysis using multiple cores, specify the number
#' @return Depending on options selected returns either a list of F values and p values, or just F, or just p-values
#'  for association with each variable in the big.matrix.
#' @param verbose logical, whether to display additional output on progress
#' @return If both F.values and p.values are TRUE, returns dataframe of both statistics for each variable, else a vector.
#' If the phenotype has 20 more or more unique categories, it will be assumed to be continuous and the association
#' test applied will be correlation. If there are two categories a t-test will be used, and 3 to 19 categories, an ANOVA#
#' will be used. Regardless of the analysis function, output will be converted to an F statistic and/or associated p-values.
#' Except if p.values and F.values are both set to FALSE and the phenotype is continuous, then pearsons correlation values
#' will be returned
#' @export
#' @seealso \code{\link{get.big.matrix}}
#' @author Nicholas Cooper 
#' @examples 
#' bmat <- generate.test.matrix(5,big.matrix=TRUE)
#' pheno <- rep(1,ncol(bmat)); pheno[which(runif(ncol(bmat))<.5)] <- 2
#' ids <- colnames(bmat); samp.inf <- data.frame(phenotype=pheno); rownames(samp.inf) <- ids
#' both <- quick.pheno.assocs(bmat,samp.inf); prv(both)
#' Fs <- quick.pheno.assocs(bmat,samp.inf,verbose=TRUE,p.values=FALSE); prv(Fs)
#' Ps <- quick.pheno.assocs(bmat,samp.inf,F.values=FALSE); prv(Ps)
quick.pheno.assocs <- function(bigMat,sample.info=NULL,use.col="phenotype",dir="",
                               p.values=TRUE,F.values=TRUE,n.cores=1,verbose=FALSE)
{
  ## use sample.info
  # create list of bigMatrix locations
  # go 'N' snps at time, concatenate into 1 file, run regression assocs
  # for 2 phenotypes/grps gives t values - ordinally equivalent to logistic regression with 2 groups
  if(!all(c(use.col) %in% colnames(sample.info))) { stop(paste("sample.info was invalid for association tests, need ",use.col," column")) }
  if(verbose) {
    cat(" running row-wise tests against",use.col,"\n")
  }
  bigMat <- get.big.matrix(bigMat,dir)
  samp.list <- colnames(bigMat); 
  tot.samps <- length(samp.list)
  # check that sample.info is valid, if not attempt to fix it
  if(!is.data.frame(sample.info)) { stop("sample.info must be a data.frame (containing phenotype information for each row id)") }
  ## ENSURE ONLY USING SAMPLES IN THE BIGMATRIX ##
  cutt <- which(!rownames(sample.info) %in% samp.list)
  if(length(cutt)>0) {  sample.info <- sample.info[-cutt,,drop=F]  }
  if(ncol(bigMat)!=nrow(sample.info)) { 
    n.mis <- length(which(colnames(bigMat) %in% rownames(sample.info)))
    warning(n.mis," samples in BigMat were found in sample.info [failure likely]")
  }
  ## determine test to use based on number of phenotypes ##
  #print(head(sample.info)); iii34 <- sample.info[[paste(use.col)]]; prv(iii34)
  #print(table(sample.info[[paste(use.col)]],useNA=NULL))
  n.phenos <- length(table(sample.info[[paste(use.col)]],useNA=NULL))
  #si <- table(sample.info[[paste(use.col)]],useNA=NULL)
  #prv(si,use.col,sample.info)
  t.type <- "single"
  if(n.phenos==2) { t.type <- "t.test"}
  if(n.phenos>2) { t.type <- "anova"}
  if(n.phenos>19) { t.type <- "cor"; if(!p.values & !F.values) { t.type <- "cor2"; F.values <- T } }
  if(verbose) {
    cat(" found ",n.phenos," ",use.col,"s, ",t.type," will be used to summarise rows most associated with ",use.col,"\n\n",sep="")
  }
  three.test <- function(col,pheno) { return(summary(aov(col~pheno))[[1]][["F value"]][1]) }
  two.test <- function(col,pheno) { return((cor.test(col,pheno)$statistic)^2)  }
  cnt.test <- function(col,pheno) { return(cor(col,pheno,use="pairwise.complete"))  }
  ph.test <- switch(t.type,anova=three.test,t.test=two.test,cor=two.test,cor2=cnt.test,single=NULL)
  deg.fr <- switch(t.type,anova=n.phenos-1,t.test=1,cor=1,single=NULL)
  if(is.null(ph.test)) { stop("Error: used option for association test by ",use.col," but there is only 1 phenotype in file")}
  row.labels <- rownames(bigMat)
  full.size <- length(row.labels)
  good.mem.lim <- 10^7 # allows fast processing at this limit
  if(ncol(bigMat)>10^5) { good.mem.lim <- 10*(ncol(bigMat)) }
  #print(estimate.memory(bigMat))
  if(estimate.memory(bigMat) < .3) { tracker <- F } else { tracker <- T }
  opt.size <- round(good.mem.lim/tot.samps)
  n.segs <- ceiling(full.size/opt.size)
  #prv(n.segs,bigMat)
  seg.starts <- 1+c(0:(n.segs+2))*opt.size
  seg.starts[seg.starts>full.size] <- full.size+1 # +1 ensures final snp included
  seg.starts <- seg.starts[!duplicated(seg.starts)]
  results <- vector("list", n.segs)
  pheno <- sample.info[[use.col]]
  if(!is.numeric(pheno)) { 
    pp <- pheno; pheno <- as.numeric(factor(paste(pp)))
    warning("pheno was not numeric, so recoded:\n")
    ii <- table(pp,pheno); cat(paste(rownames(ii),colnames(ii),sep=" ==> "),"\n")  
  }
  test.seg <- matrix(numeric(),nrow=opt.size,ncol=tot.samps)
  # break analysis into chunks that fit in memory
  # NB: avoided parallel computing - it does not seem to add any speed here?
  kk <- proc.time()[3]
  if(n.cores>1) {
    Fvalues <- bmcapply(bigMat,1,FUN=ph.test,dir=dir,by=200,n.cores=n.cores,pheno=pheno)
    #prv(Fvalues)
  } else {  
    for (dd in 1:n.segs)
    {
      row.subset <- (seg.starts[dd]):(seg.starts[dd+1]-1)
      nr <- length(row.subset)
      # test subset of snps for each segment in turn
      test.seg[1:nr,] <- bigMat[row.subset,] 
      results[[dd]] <- apply(test.seg[1:nr,],1,ph.test,pheno=pheno)
      if(tracker) { loop.tracker(dd,n.segs) }
    }
    Fvalues <- (do.call(c,results))
  }
  #Fvalues <- round(Fvalues,4) # remember t^2 = F
  if(verbose) {
    cat(" took",round((proc.time()[3]-kk)/60,1),"minutes\n")
  }
  ### option to return p values? could also have option for CHI squared.
  ## also should be able to the bigmatrix apply for a speedup
  # p.from.t <- function(t,df) { pt(t, df, lower=FALSE) }
  p.from.f <- function(FF,k,n) {  pf(FF, k, n,lower.tail = FALSE) }
  # pvalues <- sapply(Fvalues^.5,p.from.t,df=tot.samps-1)
  if(is.numeric(Fvalues)) {
    if(p.values) { 
      pvalues <- sapply(Fvalues,p.from.f,k=deg.fr,n=(tot.samps-1))
      out <- pvalues 
    } else { out <- Fvalues }
    if(verbose) {
      cat("\nSummary of",if(p.values & !F.values) { "p-value" } else { "F" },"statistics returned:\n"); print(summary(out)); cat("\n")
    }
    if(length(row.labels)==length(out))
    { names(out) <- row.labels } else {
      stop("Error: analysis failure as list of row labels did not match length of test statistics")
    }
  } else { out <- rep(1,length(Fvalues)) }
  if(!(p.values & F.values)) {
    return(out)
  } else {
    out <- data.frame(F=Fvalues,p=pvalues)
    rownames(out) <- row.labels
    return(out)
  }
}


#' Select subset of rows least associated with a categorical variable
#' 
#' Runs a quick association analysis on the dataset against a phenotype/categorical
#' variable stored in a dataframe, and uses the results as a way to select
#' a subset of the original matrix, so you may wish to select the 'N' least associated
#' variables, or the 'N' most associated.
#' @param bigMat a big.matrix object, or any argument accepted by get.big.matrix(), which includes
#'  paths to description files or even a standard matrix object.
#' @param keep numeric, by default a proportion (decimal) of the original number of rows/columns to choose
#'  for the subset. Otherwise if an integer>2 then will assume this is the size of the desired subset,
#'  e.g, for a dataset with 10,000 rows where you want a subset size of 1,000 you could set 'keep' as
#'  either 0.1 or 1000.
#' @param dir directory containing the filebacked.big.matrix, same as dir for get.big.matrix.
#' @param phenotype a vector which contains the categorical variable to perform an association 
#'  test for phenotype, etc. This should be the same length as the number of columns (e.g, samples)
#'  in bigMat.
#' @param least logical, whether to select TRUE, the top least associated variables, or FALSE, the most 
#'  associated.
#' @param n.cores integer, if wanting to process the analysis using multiple cores, specify the number
#' @param verbose logical, whether to display additional output
#' @return A set of row or column indexes (depents on 'rows' parameter) of the variables most dependent
#'  (or indepent) variables measured by association with a [continuous/categorical] phenotype.
#' @export
#' @seealso \code{\link{quick.pheno.assocs}}
#' @author Nicholas Cooper 
#' @examples 
#' bmat <- generate.test.matrix(5,big.matrix=TRUE)
#' pheno <- rep(1,ncol(bmat)); pheno[which(runif(ncol(bmat))<.5)] <- 2
#' most.correl <- select.least.assoc(bmat,phenotype=pheno,least=FALSE)
#' least.correl <- select.least.assoc(bmat,phenotype=pheno,least=TRUE)
#' cor(bmat[least.correl,][1,],pheno)  # least correlated
#' cor(bmat[most.correl,][1,],pheno)  # most correlated
select.least.assoc <- function(bigMat,keep=.05,phenotype=NULL,least=TRUE,dir="",n.cores=1,verbose=TRUE)
{
  # select a subset of variables based on the most representative variables in the PCs of a subset
  bigMat <- get.big.matrix(bigMat,dir=dir)
  N <- nrow(bigMat) 
  if(keep>2) {
    new.n <- round(keep)
  } else {
    new.n <- round(max(0,min(1,keep))*N)
  }
  if(is.null(colnames(bigMat)) | is.null(rownames(bigMat)) ) {
    options(bigmemory.allow.dimnames=TRUE)
    if(is.null(colnames(bigMat))) { colnames(bigMat) <- paste("c",1:ncol(bigMat),sep="") }
    if(is.null(rownames(bigMat))) { rownames(bigMat) <- paste("r",1:nrow(bigMat),sep="") }
  }
  if(length(phenotype)!=ncol(bigMat)) { stop("phenotype must be the same length as ncol(bigMat)") }
  sample.info <- data.frame(phenotype=phenotype); rownames(sample.info) <- colnames(bigMat)
  ##association version of selecting subset of snps for PCA
  # takes the 'pc.to.keep'% least associated with #sample.info$phenotype'
  myPs <- quick.pheno.assocs(bigMat=bigMat,sample.info=sample.info,
                             dir=dir,F.values=F,n.cores=n.cores,verbose=verbose)
  rank <- order(myPs)
  if(least) { rank <- rev(rank) }
  if(verbose) {
    #cat("\nsummary of p-values for association tests\n"); print(summary(myPs))
  }
  kept.snps <- rank[1:new.n]
  return(sort(kept.snps))
}


# internal function
cut.fac <- function(N,n.grps,start.zero=FALSE,factor=TRUE) {
  X <- findInterval(((1:N)/(N/n.grps)), 1:n.grps, all.inside=F, rightmost.closed=T)+as.numeric(!start.zero)
  #X <- cut(1:N,n.grps,labels=F)
  if(factor) { X <- as.factor(X) }
  return(X)
}





#' PCA/Singular Value Decomposition for big.matrix
#' 
#' At the time of submission there was no existing native method to conduct principle components
#' analysis (PCA) on big.matrix objects. This function allows singular value decomposition (SVD) of
#' very large matrices, very efficiently, versus the default method. The major speed advantages
#' occur when the 'bigalgebra' package is installed, and when the argument for this function
#' 'SVD'=TRUE. Regular PCA can be conducted using SVD=FALSE but it will be slower and the maximum
#' matrix size able to produce a result, given memory limits, will be smaller. SVD is not exactly
#' the same as PCA, but from my testing the components produced will correlate R>.9 with components
#' of PCA on the same matrix. This function is not completely native to big.matrix objects so 
#' there is one step where the matrix submitted needs to be loaded into memory, so if your 
#' big.matrix object is larger than the allowed size of a standard R-matrix [[which is roughly 3GB; 
#' you can check using NCmisc::estimate.memory()], then this function will fail unless you set the 
#' option 'thin' to a percentage that, multiplied by the original matrix memory-load, is under 3GB.
#' For large matrices in my applications, components produced with thinning are still highly 
#' correlated with components produced using the full dataset. For a breakdown of thinning methods,
#' see the description for the function thin() for more  information. Even with medium sized
#' matrices, for instance 15,000 x 50,000 in size, this function is orders of magnitude faster
#' than the standard R PCA functions, usually running in a matter of minutes, rather than hours
#' or days in examples that I have tested, due to much better handling of memory for internal
#' transpose and eigen operations by using the 'bigmemory' architecture.
#'
#' @param bigMat a big.matrix object, or any argument accepted by get.big.matrix(), which includes
#'  paths to description files or even a standard matrix object.
#' @param dir directory containing the filebacked.big.matrix, and also where the output
#'  file will be saved by default if the save.pcs option is set TRUE. 
#' @param pcs.to.keep integer, number of principle components to keep. Singular Value Decomposition
#'  methods are able to run faster if they don't need to calculate every single PC for a large
#'  matrix. Default is to calculate only the first 50; in practice even fewer than this are generally
#'  used directly. Apart from reducing processing time, this can also reduce storage/RAM burden for 
#'  the resulting matrix. Set to NA, or a number >= min(dim(bigMat)) in order to keep all PCs.
#' @param thin decimal, percentage of the original number of rows you want to thin the matrix to.
#'  see function thin() for details of how this can be done, pass arguments to thin() using ...
#'  Even though this PCA function uses mainly 'big.matrix' native methods, there is a step where the
#'  matrix must be stored fully in memory, so this limits the size of what matrix can be processed,
#'  depending on RAM limits. If you want to conduct PCA/SVD on a matrix larger than RAM you can thin
#'  the matrix to a percentage of the original size. Usually such a large matrix will contain correlated
#'  measures and so the exclusion of some data-rows (variables) will have only a small impact on the
#'  resulting principle components. In some applications tested using this function, using only 5% 
#'  of 200,000 variables a PCA gave extremely similar results to using the full dataset.
#' @param SVD logical, whether to use a Singular Value Decomposition method or a PCA method. The 
#'  eigenvalues and eigenvectors of each alternative will be highly correlated so for most applications,
#'  such as PC-correction, this shouldn't make much difference to the result. However, using SVD=TRUE
#'  can provide significant advantages in speed, or even allow a solution on a matrix that would be
#'  to large to practically compute full PCA. Note that only in SVD mode, and with the bigalgebra
#'  package installed will the full speed advantage of this function be utilised.
#' @param LAP logical, whether to use La.svd() instead of svd() when SVD=TRUE, see base:svd for more info.
#' @param center whether to 'centre' the matrix rows by subtracting rowMeans() before conducting the PCA. This is usually
#'  advisable, although you may wish to skip this if the matrix is already centred to save extra processing.
#'  unlike prcomp there is no option to standardize or use the correlation matrix, if this is desired, please
#'  standardize the bigMat object before running this function. Alternatively, 'center' can be a vector of length
#'  nrow(bigMat) which gives an offset to correct each row by.
#' @param save.pcs whether to save the principle component matrix and eigenvalues to a binary file with name pcs.fn
#' @param pcs.fn name of the binary when save.pcs=TRUE
#' @param return.loadings logical, whether to return the 'loadings' (or the other singular vector when SVD=T); could result in a speed decrease
#' @param verbose whether to display detailed progress of the PCA
#' @param use.bigalgebra logical, whether to use the bigalgebra package for algebraic operations. For large
#'  datasets bigalgebra should provide a substantial speedup, and also facilitates use of larger matrices. This 
#'  relies on having bigalgebra installed and loaded, which requires some manual configuration as bigalgebra
#'  is not currently available on CRAN, but only SVN and RForge. See svn.bigalgebra.install() or big.algebra.install.help()
#'  Default is to use bigalgebra if available (TRUE), but setting this FALSE prevents the check for bigalgebra which would be
#'  cleaner if you know that you don't have it installed.
#' @param ... if thin is TRUE, then these should be any additional arguments for thin(), e.g, 'keep', 'how', etc.i
#' @param delete.existing logical, whether to automatically delete filebacked matrices (if they exist) 
#' before rewriting. This is because of an update since 20th October 2015 where bigmemory won't allow
#' overwrite of an existing filebacked matrix. If you wish to set this always TRUE or FALSE, use
#'  options(deleteFileBacked)
#' @return A list of principle components/singular vectors (may be incomplete depending on options selected), and of
#'  the eigenvalues/singular values.
#' @export
#' @seealso \code{\link{get.big.matrix}}, \code{\link{PC.correct}}
#' @author Nicholas Cooper
#' @examples 
## PRELIMINARY EXAMPLES: demonstration of PCA versus SVD ##
#' # create an example matrix and its transpose
#' min.dim <- 200; nvar <- 500; subset.size <- 50
#' mat <- matrix(rnorm(min.dim*nvar),ncol=min.dim) 
#' prv.large(mat)
#' t.mat <- t(mat)
#' # create two alternative covariance matrices
#' MMs <- t.mat %*% mat
#' MsM <- mat %*% t.mat
#' # run singular value decomposition
#' pca <- svd(mat)   
#' D <- pca$d # singular values (=sqrt(eigenvalues))
#' V <- pca$v # right singular vector
#' U <- pca$u # left singular vector
#' sig <- mat-mat; diag(sig) <- D; 
#' MMs2 <- V %*% (t(sig) %*% sig) %*% t(V)
#' sig <- t.mat-t.mat; diag(sig) <- D; 
#' MsM2 <- U %*% (sig %*% t(sig)) %*% t(U)
#' # show that the covariance matrices are equal to the functions of 
#' # the left and right singular vectors
#' prv(MMs,MsM); prv(MMs2,MsM2)
#' pr <- princomp(mat) # PCA using eigendecomposition of cov matrix
#' L <- matrix(rep(0,40000),ncol=200); diag(L) <- pr[[1]]^2 # eigenvalues as diag
#' mat2 <- (pr[[2]]) %*% L %*%  solve(pr[[2]]) # = eigenvectors * eigenvalues * inv(eigenvectors)
#' prv.large(cov(mat)); prv.large(mat2) #  == COVmat (may be slight tolerance differences)
#' ## Now demonstrate the correlation between SVD and PCA ##
#' # the right singular vector is highly correlated with the pca loadings:
#' median(abs(diag(cor(V,pr[["loadings"]]))))
#' # the left singular vector is highly correlated with the pca scores (eigenvectors):
#' median(abs(diag(cor(U,pr[["scores"]]))))
#' cor(pr$sdev,D) # the singular values are equivalent to the eigenvalues
#' 
#' ## MAIN EXAMPLES ##
#' bmat <- as.big.matrix(mat,backingfile="testMyBig.bck",descriptorfile="testMyBig.dsc")
#' result <- big.PCA(bmat) #,verbose=TRUE)
#' headl(result)
#' # plot the eigenvalues with a linear fit line and elbow placed at 13
#' Eigv <- pca.scree.plot(result$Evalues,M=bmat,elbow=6,printvar=FALSE)
#' unlink(c("testMyBig.bck","testMyBig.dsc"))
#' ##  generate some data with reasonable intercorrelations ##
#' mat2 <- sim.cor(500,200,genr=function(n){ (runif(n)/2+.5) })
#' bmat2 <- as.big.matrix(mat2,backingfile="testMyBig.bck",descriptorfile="testMyBig.dsc")
#' # calculate PCA on decreasing subset size 
#' result2 <- big.PCA(bmat2,thin=FALSE)
#' result3 <- big.PCA(bmat2,thin=TRUE,keep=.5)
#' result4 <- big.PCA(bmat2,thin=TRUE,keep=.5,how="cor")
#' result5 <- big.PCA(bmat2,thin=TRUE,keep=.5,how="pca")
#' result6 <- big.PCA(bmat2,thin=TRUE,keep=.2)
#' normal <- result2$PCs
#' thinned <- result3$PCs
#' corred <- result4$PCs
#' pced <- result5$PCs
#' thinner <- result6$PCs
#' ## correlate the resulting PCs with the un-thinned PCs
#' cors.thin.with.orig <- apply(cor(normal,thinned),1,max)
#' cors.corred.with.orig <- apply(cor(normal,corred),1,max)
#' cors.pced.with.orig <- apply(cor(normal,pced),1,max)
#' cors.thinner.with.orig <-apply(cor(normal,thinner),1,max)
#' plot(cors.thin.with.orig,type="l",col="red",ylim=c(0,1))
#' lines(cors.thinner.with.orig,col="orange")
#' lines(cors.corred.with.orig,col="lightblue")
#' lines(cors.pced.with.orig,col="lightgreen")
#' # can see that the first component is highly preserved,
#' # and next components, somewhat preserved; try using different thinning methods
#' unlink(c("testMyBig.bck","testMyBig.dsc"))
big.PCA <- function(bigMat,dir=getwd(),pcs.to.keep=50,thin=FALSE,SVD=TRUE,LAP=FALSE,center=TRUE,
                    save.pcs=FALSE,use.bigalgebra=TRUE,pcs.fn="PCsEVsFromPCA.RData",
                    return.loadings=FALSE,verbose=FALSE,delete.existing=getOption("deleteFileBacked"),...) 
{
  # run principle components analysis on the SNP subset of the LRR snp x sample matrix
  # various methods to choose from with pro/cons of speed/memory, etc.
  #  must use SNP-subset to avoid LD, destroying main effects, +avoid huge memory requirements
  #dir <- validate.dir.for(dir,c("big","pc"))
  if(exists("validate.dir.for",mode="function")) {
    ## plumbCNV specific code ##
    dir <- do.call("validate.dir.for",list(dir=dir,elements=c("big","pc"),warn=F))  
  } else {
    # otherwise
    dir <- list(big=dir,pc=dir)
  }
  #must.use.package(c("irlba"),T)
  if(thin) {
    if(verbose) {  prv.big.matrix(bigMat) }
    bigMat <- thin(bigMat,dir=dir,...)
    if(verbose) {  prv.big.matrix(bigMat) }
  } 
  pcaMat <- get.big.matrix(bigMat,dir)
  #print(dim(pcaMat))
  if(verbose & !thin) { prv.big.matrix(pcaMat,name="pcaMat") }
  est.mem <- estimate.memory(pcaMat)
  if(est.mem>1) {
    cat(" estimated memory required for",nrow(pcaMat),"x",ncol(pcaMat),"matrix:",round(est.mem,2),
      "GB. If this exceeds available,\n  then expect PCA to take a long time or fail!\n")
  }
  #print(packages.loaded())
  subMat <- pcaMat[1:nrow(pcaMat),1:ncol(pcaMat)] # must convert bigmatrix to plain matrix here, no pca yet takes a bigmatrix
  rm(pcaMat)
  # center using row means
  if(length(center)>1) {
    if(length(center)==nrow(subMat)) {
      CM <- center
      center <- TRUE
      rm.sub <- function(X) { 
        mmm <-  matrix(rep(CM,times=ncol(subMat)),ncol=ncol(subMat),byrow=F)
        prv(mmm)
        return(mmm)
      }
    } else {
      rm.sub <- rowMeans ; warning("center vector didn't match number of rows of 'bigMat', data left uncentered")
      center <- FALSE
    }
  } else { 
      rm.sub <- rowMeans # a bit hacky?
  }
  if(center) {
    if(verbose) { cat(" centering data by row means...") }
    subMat <- subMat - rm.sub(subMat)  #matrix(rep(rowMeans(subMat),times=ncol(subMat)),ncol=ncol(subMat))
    subMat[is.na(subMat)] <- 0 # replace missing with the mean
    cat(" means for first 10 snps:\n")
    print(round(head(rowMeans(subMat),10))) # show that centering has worked
  } else {
    subMat <- apply(subMat,1,row.rep)
  }
  if(verbose) { cat(" replaced missing data with mean (PCA cannot handle missing data)\n") }
  #subMat <- t(subMat) # transpose
  dimz <- dim(subMat)
  if(!is.numeric(pcs.to.keep) | is.integer(pcs.to.keep)) { pcs.to.keep <- NA }
  if(is.na(pcs.to.keep)) { pcs.to.keep <- min(dimz) }
  if(pcs.to.keep > min(dimz)) { 
    # ensure not trying to extract too many pcs
    warning(paste("selected too many PCs to keep [",pcs.to.keep,"], changing to ",min(dimz),"\n",sep="")) 
    pcs.to.keep <- min(dimz)
  } 
  if(!SVD & (dimz[2]>dimz[1])) {
    if(verbose) { cat(" PCA using 'princomp' (only for datasets with more samples than markers)\n") }
    print(system.time(result <- princomp(t(subMat))))
    PCs <- result$scores[,1:pcs.to.keep]
    loadings <- result$loadings[,1:pcs.to.keep]
    Evalues <- result$sdev^2 # sds are sqrt of eigenvalues
  } else {
    if(!SVD) {
      if(verbose) {
        cat(" PCA by crossproduct and solving eigenvectors\n")
        cat(" obtaining crossproduct of the matrix and transpose XtX...")
      }
      uu <-(system.time(xtx <- crossprod(subMat)))
      if(verbose) { 
        cat("took",round(uu[3]/60,1),"minutes\n")
        cat(" obtaining eigen vectors of the crossproduct XtX...")
      }
      uu <-(system.time(result <- eigen((xtx/nrow(subMat)),symmetric=T)))
      if(verbose) {  cat("took",round(uu[3]/60,1),"minutes\n") }
      PCs <- result$vectors[,1:pcs.to.keep]
      Evalues <- result$values
      loadings <- NULL
    } else {
      pt <- "package:"; pkgset <- gsub(pt,"",search()[grep(pt,search())])
      do.fast <- (!LAP & (all(c("irlba","bigalgebra") %in% pkgset)))
      if(!use.bigalgebra) { do.fast <- F }
      if(verbose) {
        cat(" PCA by singular value decomposition...") # La.svd gives result with reversed dims. (faster?)
      } 
      if(return.loadings)  { nU <- pcs.to.keep } else { nU <- 0 }
      if(!LAP) {
        if(do.fast) {
          uu <-(system.time(result <- irlba(subMat,nv=pcs.to.keep,nu=nU,mult=matmul))) 
        } else {
          if(use.bigalgebra & verbose) { warning("[without 'bigalgebra' package, PCA runs slowly for large datasets,", 
              "see 'big.algebra.install.help()']\n") }
          uu <-(system.time(result <- svd(subMat,nv=pcs.to.keep,nu=nU)))
        }
        if(verbose) { cat("took",round(uu[3]/60,1),"minutes\n") }
        PCs <- result$v[,1:pcs.to.keep]
        #print("thus ones"); prv(result,return.loadings,nU)
        if(return.loadings) { loadings <- result$u[,1:pcs.to.keep] } else { loadings <- NULL }
        Evalues <- result$d^2 # singular values are the sqrt of eigenvalues 
      } else {
        if(verbose) { cat("\n [using LAPACK alternative with La.svd]") }
        uu <- (system.time(result <- La.svd(subMat,nv=pcs.to.keep,nu=nU)))
        if(verbose) { cat("took",round(uu[3]/60,1),"minutes\n") }
        PCs <- t(result$vt)[,1:pcs.to.keep]  ##?
       # print("thOs ones")
        if(return.loadings) { loadings <- result$u[,1:pcs.to.keep] } else { loadings <- NULL }
        Evalues <- result$d^2 # singular values are the sqrt of eigenvalues
      }
    }
  }
  rownames(PCs) <- colnames(subMat)
  colnames(PCs) <- paste("PC",1:pcs.to.keep,sep="")
  if(save.pcs) {
    ofn <- cat.path(dir$pc,pcs.fn)
    cat(paste("~wrote PC data to file:",ofn,"\n"))
    save(PCs,Evalues,loadings,file=ofn) }
  if(return.loadings & exists("loadings")) {
    colnames(loadings) <- paste("PC",1:pcs.to.keep,sep="")
    rownames(loadings) <- rownames(subMat)
    out.dat <- list(PCs,Evalues,loadings)
    names(out.dat) <- c("PCs","Evalues","loadings")
  } else {
    out.dat <- list(PCs,Evalues)
    names(out.dat) <- c("PCs","Evalues")
  }
  return(out.dat)
}






#' Correct a big.matrix by principle components
#' 
#' Principle components (PC) can be used as a way of capturing bias (when common variance represents bias)
#' and so PC correction is a way to remove such bias from a dataset. Using the first 'n' PCs from an 
#' an analysis performed using big.PCA(), this function will transform the original matrix by regressing
#' onto the 'n' principle components (and optionally gender) and returing the residuals. The result
#' is returned as a big.matrix object, so that objects larger than available RAM can be processed, and
#' multiple processors can be utilised for greater speed for large datasets.
#' 
#' @param pca.result result returned by 'big.PCA()', or a list with 2 elements containing
#'  the principle components and the eigenvalues respectively (or SVD equivalents). Alternatively,
#'  can be the name of an R binary file containing such an object.
#' @param bigMat a big.matrix with exactly corresponding samples (columns) to those submitted to PCA prior to correction
#' @param dir directory containing the big.matrix backing file
#' @param num.pcs number of principle components (or SVD components) to correct for
#' @param n.cores number of cores to use in parallel for processing
#' @param pref prefix to add to the file name of the resulting corrected matrix backing file
#' @param big.cor.fn instead of using 'pref' directly specify the desired file name
#' @param write whether to write the result to a file.backed big.matrix or to simply
#'  return a pointer to the resulting corrected big.matrix
#' @param sample.info if using 'correct.sex=TRUE' then this object should be
#'  a dataframe containing the sex of each sample, with sample names as rownames
#' @param correct.sex if sample.info is a dataframe containing a column named 'gender' or 'sex'
#'  (case insensitive), then add a sex covariate to the PC correction linear model
#' @param add.int logical, whether to maintain the pre-corrected means of each variable, i.e, post-correction
#'  add the mean back onto the residuals which will otherwise have mean zero for each variable.
#' @param preserve.median logical, if add.int=TRUE, then setting this parameter to TRUE will preserve
#'  the median of the original data, instead of the mean. This is because after PC-correction the 
#'  skew may change.
#' @param tracker logical, whether to display a progress bar
#' @param verbose logical, whether to display preview of pre- and post- corrected matrix
#' @param delete.existing logical, whether to automatically delete filebacked matrices (if they exist) 
#' before rewriting. This is because of an update since 20th October 2015 where bigmemory won't allow
#' overwrite of an existing filebacked matrix. If you wish to set this always TRUE or FALSE, use
#'  options(deleteFileBacked)
#' @return A big.matrix of the same dimensions as original, corrected for n PCs and an optional covariate (sex)
#' @export
#' @seealso \code{\link{big.PCA}}
#' @author Nicholas Cooper 
#' @examples 
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' mat2 <- sim.cor(500,200,genr=function(n){ (runif(n)/2+.5) })
#' bmat2 <- as.big.matrix(mat2,backingfile="testMyBig.bck",descriptorfile="testMyBig.dsc")
#' ## calculate PCA ##
#' # result2 <- big.PCA(bmat2,thin=FALSE)
#' # corrected <- PC.correct(result2,bmat2)
#' # corrected2 <- PC.correct(result2,bmat2,n.cores=2)
#' # c1 <- get.big.matrix(corrected) ; c2 <- get.big.matrix(corrected2)
#' # all.equal(as.matrix(c1),as.matrix(c2))
#' unlink(c("testMyBig.bck","testMyBig.dsc"))
#' setwd(orig.dir) # reset working dir to original
PC.correct <- function(pca.result,bigMat,dir=getwd(),num.pcs=9,n.cores=1,pref="corrected",
                            big.cor.fn=NULL,write=FALSE,sample.info=NULL,correct.sex=FALSE,
                            add.int=FALSE,preserve.median=FALSE, tracker=TRUE,verbose=TRUE,
                            delete.existing=getOption("deleteFileBacked"))
{
  ## using results of a PCA analysis, run correction for 'num.pcs' PCs on a dataset
  # uncorrected matrix
  #dir <- validate.dir.for(dir,c("big","pc"))
  if(exists("validate.dir.for",mode="function")) {
    ## plumbCNV specific code ##
    dir <- do.call("validate.dir.for",list(dir=dir,elements=c("big","pc"),warn=F))  
  } else {
    # otherwise
    dir <- list(big=dir,pc=dir)
  }
  origMat <- get.big.matrix(bigMat,dir)
  if(verbose) {
    cat("\nRunning Principle Components correction (PC-correction), using LRR-dataset:\n")
    prv.big.matrix(origMat,name="origMat")
  }
  if(n.cores>1) { multi <- T } else { multi <- F }
  # get filenames now to add to result later
  rN <- rownames(origMat); cN <- colnames(origMat)
  # run pca.correction using vectors (PCs) and values from PCA
  if(!is.list(pca.result)) {
    if(is.character(pca.result)) {
      ofn <- cat.path(dir$pc,pca.result) 
      if(file.exists(ofn))
      {
        pca.file <- get(load(ofn))
        cat(" loaded PCA values and vectors\n")
        nn <- which(toupper(names(pca.result))=="PCS")
        if(length(nn)==0) { nn <- 1 }
        PCs <- pca.result[[nn]]
      } else {
        stop("Error: file",ofn,"does not exist\n")
      }
    } else {
      if(ncol(origMat) %in% dim(pca.result))
      {
        #given dimensions = number of samples, assume PCs entered as matrix
        PCs <- pca.result
      } else {
        stop("Error: expecting file name or PC matrix: pca.result\n")
      }
    } 
  } else {
    nn <- which(toupper(names(pca.result))=="PCS")
    if(length(nn)==0) { nn <- 1 }
    PCs <- pca.result[[nn]]
  }
  # create new matrix same size, ready for corrected values
  nR <- nrow(origMat); nC <- ncol(origMat)
  if(nR<2) { stop("too few rows to run PC correction") }
  if(nC<2) { stop("too few columns to run PC correction") }
  if(verbose) { cat(" creating new file backed big.matrix to store corrected data...") }
  fnm <- paste(pref,"bck",sep=".")
  if(delete.existing & file.exists(cat.path(dir$big,fnm))) { unlink(cat.path(dir$big,fnm))  }
  pcCorMat <- filebacked.big.matrix(nR,nC, backingfile=paste(pref,"bck",sep="."),
                                    backingpath=dir$big, descriptorfile=paste(pref,"dsc",sep="."))
  cat("done\n")
  if(!is.filebacked(pcCorMat) | !is.filebacked(origMat)) {
    warning("at least one of the big.matrices is not filebacked, memory problems may be encountered")
  }
  # in result$vectors, PCs are the columns, only need first 10 or so
  # rows are subjects / samples
  col.sel <- 1:ncol(origMat)
  nPCs <- PCs[,1:num.pcs]; sex.txt <- ""
  if(correct.sex) { 
    ## add a column to the PC matrix which will allow to covary for sex too
    coln <- which(tolower(colnames(sample.info)) %in% c("gender","sex")) 
    if(length(coln)>0) { 
      indx <- match(colnames(origMat),rownames(sample.info))
      if(all(!is.na(indx))) {
        sex.vec <- sample.info[indx,coln[1]] 
        sex.vec[is.na(sex.vec)] <- mean(sex.vec,na.rm=T) # replace missing sex with mean 
        nPCs <- cbind(sex.vec,nPCs) ; sex.txt <- " (covarying for sex)"
      } else { warning("could not correct for sex as sample.info was missing some IDs") }
    }
  }
  nPCs <- cbind(rep(1,nrow(nPCs)),nPCs) # this adds intercept term for lm.fit() [remove if using lm() ]
  cat(" correcting by principle components",sex.txt,", taking the regression-residual for each variable\n",sep="")
  jj <- proc.time()
  nsamples <- ncol(origMat)
  num.snps <- nrow(origMat); sampz <- 1:nsamples
  snps.per.proc <- 400;   flush.freq <- 20
  if(nsamples>10000) { snps.per.proc <- 300 }; if(nsamples>20000) { snps.per.proc <- 150 }
  if(nsamples>40000) { snps.per.proc <- 50 }; if(nsamples>80000) { snps.per.proc <- 10 }
  ## assume if we have lots of cores, we'd also have lots of RAM too
  if(n.cores>5) { snps.per.proc <- snps.per.proc*2; flush.freq <- flush.freq*2 } 
  if(n.cores>15) { snps.per.proc <- snps.per.proc*2 }
  #snps.per.proc <- max(snps.per.proc,n.cores) # at least 1 snp per core as a minimum
  if(snps.per.proc<n.cores) { n.cores <- snps.per.proc %/% 5 } # if really low num, use less cores
  stepz <- round(seq(from=1,to=num.snps+1,by=snps.per.proc))
  if((tail(stepz,1)) != num.snps+1) { stepz <- c(stepz,num.snps+1) }
  # if the last step is too small, merge with previous
  sc.lst <- head(tail(stepz,2),1); lst <- tail(stepz,1)
  if(lst-sc.lst <=4) { ll <- length(stepz); if(ll>2) { stepz <- stepz[-(ll-1)] } else { warning("number of SNPs quite small, may cause issues") } }
  split.to <- length(stepz)-1
  big.extras <- T # flush memory every 'n' iterations.
  # this simple way works (instead of big for-loop) but hogs memory and is no faster
  # [NB: requires transpose of target corrected big matrix dimensions]
  ### pcCorMat <- apply(origMat,1,PC.fn,nPCs=nPCs,col.sel=sampz)
  for (dd in 1:split.to)
  {
    x1 <- stepz[dd]; x2 <- stepz[dd+1]-1 #subset row selection
    # use of this 'sub.big.matrix' structure, stops the memory leak behaviour which spirals
    # the memory relating to 'origMat' out of control. 
    next.rows <- sub.big.matrix(origMat, firstRow=x1, lastRow=x2, backingpath=dir$big )
    if(length(Dim(next.rows))!=2) { stop("expected a matrix, got a vector")}
   # prv(next.rows,x1,x2,nPCs,multi,split.to,stepz)
    # next.rows is now a pointer to a matrix subset, must use 'as.matrix' to coerce to a regular R object 
    if(multi) {
      #pcCorMat[x1:x2,] <- PC.fn.mat.multi(next.rows[1:nrow(next.rows),1:ncol(next.rows)],nPCs,mc.cores=n.cores,add.int=add.int)
      pcCorMat[x1:x2,] <- PC.fn.mat.multi(bigmemory::as.matrix(next.rows),nPCs,mc.cores=n.cores,add.int=add.int, pm=preserve.median)
    } else {
      #pcCorMat[x1:x2,] <- PC.fn.mat.apply(next.rows[1:nrow(next.rows),1:ncol(next.rows)],nPCs,add.int=add.int)
      pcCorMat[x1:x2,] <- PC.fn.mat.apply(bigmemory::as.matrix(next.rows),nPCs,add.int=add.int, pm=preserve.median)
    }
    if(tracker) { loop.tracker(dd,split.to) }
    ## Every 'flush.freq' iterations clean up the memory, remove the 
    ##  big.matrix object 'pcCorMat' and re-attach it 
    if(dd %% flush.freq == 0) {    
      fl.suc <- bigmemory::flush(pcCorMat) & bigmemory::flush(next.rows)
      if(!fl.suc) { cat("flush failed\n") } 
      gc()  # garbage collection
      if(big.extras) {
        RR <- bigmemory::describe(pcCorMat)
        rm(pcCorMat)
        pcCorMat <- attach.big.matrix(RR,path=dir$big)
      }
    }
    rm(next.rows) # remove the sub-matrix pointer each iteration or this memory builds up 
  }
  
  options(bigmemory.allow.dimnames=TRUE)
  rownames(pcCorMat) <- rN;  colnames(pcCorMat) <- cN 
  ll <- proc.time()
  time.taken <- round((ll-jj)[3]/3600,3)
  if(time.taken>1/180) {  cat(paste(" PC-Correction took",time.taken,"hours\n")) }
  bigmemory::flush(pcCorMat) # should allow names to take  
  if(verbose) {
    cat("\nPC-corrected dataset produced:\n")
    prv.big.matrix(pcCorMat,name="pcCorMat")
  }  
  mat.ref <- bigmemory::describe(pcCorMat)
  if(write) {
    if(is.null(big.cor.fn) | !is.character(big.cor.fn)) {
      big.fn <- paste("describePCcorrect",num.pcs,".RData",sep="")
    } else {
      big.fn <- big.cor.fn[1]
    }
    ofn <- cat.path(dir$big,big.fn)
    save(mat.ref,file=ofn)
    cat(paste("~wrote PC-corrected data description file to:\n ",ofn,"\n"))
    return(big.fn)
  } else {
    return(mat.ref)
  }
}



#' Transpose function for big.matrix objects
#'
#' At the time of writing, there is no transpose method for big.matrix()
#' This function returns a new filebacked big.matrix which is the transpose of the input
#' big.matrix. max.gb allows periodic manual flushing of the memory to be conducted in case
#' the built-in memory management of R/bigmemory is not working as desired.
#' This method is a non-native (not using the raw C objects from the package but merely
#' standard R accessors and operations) algorithm to transpose a big matrix efficiently
#' for memory usage and speed. A blank matrix is created on disk and the data is
#' block-wise transposed and buffered into the new matrix.
#'
#' @param bigMat default, a big.matrix(), although if 'file.ok' is set TRUE, then
#'  this can be a big.matrix descriptor, or a file location
#' @param dir the directory for the matrix backing file (preferably for both the original
#'  and the proposed transposed matrix). If this is left NULL and bigMat contains a path,
#'  this path (via dirname(bigMat)) will be used; if it doesn't contain a path the current
#'  working directory will be used
#' @param name the basename of the new transposed matrix
#' @param R.descr the name of a binary file that will store the big.matrix.descriptor
#'  for the transposed matrix. If "" then the descriptor won't be saved. If NULL, then
#'  it will be <name>.RData
#' @param max.gb the maximum number of GB of data to process before flushing the big.matrix
#' @param verbose whether to print messages about each stage of the process
#' @param tracker whether to use a progress bar. NA means it will only be used if the matrix
#'  in question is larger than 1GB.
#' @param file.ok whether to accept big.matrix.descriptors or filenames as input for 
#'  'bigMat'; if T, then anything that works with get.big.matrix(bigMat,dir) is acceptable
#' @param delete.existing logical, whether to automatically delete filebacked matrices (if they exist) 
#' before rewriting. This is because of an update since 20th October 2015 where bigmemory won't allow
#' overwrite of an existing filebacked matrix. If you wish to set this always TRUE or FALSE, use
#'  options(deleteFileBacked)
#' @return A big.matrix that is the transpose (rows and columns switched) of the original matrix
#' @export
#' @examples 
#' bM <- filebacked.big.matrix(200, 500,
#'        dimnames = list(paste("r",1:200,sep=""), paste("c",1:500,sep="")),
#'        backingfile = "test.bck",  backingpath = getwd(), descriptorfile = "test.dsc")
#' bM[1:200,] <- replicate(500,rnorm(200))
#' prv.big.matrix(bM)
#' tbM <- big.t(bM,verbose=TRUE)
#' prv.big.matrix(tbM)
#' unlink(c("t.bigMat.RData","t.bigMat.bck","t.bigMat.dsc","test.bck","test.dsc"))
big.t <- function(bigMat,dir=NULL,name="t.bigMat",R.descr=NULL,max.gb=NA,
                  verbose=F,tracker=NA,file.ok=T,delete.existing=getOption("deleteFileBacked")) {
  #this can be slow!
  if(is.null(R.descr)) { R.descr <- cat.path(dirname(name),name,ext="RData") }
  if(!is.big.matrix(bigMat)) {
    if(is.matrix(bigMat) | is.data.frame(bigMat)) {
      warning("just a regular matrix, used t()") ; return(t(bigMat)) 
    } else {
      if(file.ok) { 
        # bigMat can be a file path
        if(is.null(dir) & file.exists(bigMat)) {
          if(dirname(bigMat)!=".") { dir <- dirname(bigMat) }
        }
        try(bigMat <- get.big.matrix(bigMat,dir)) 
      }
      if(!is.big.matrix(bigMat)) {
        stop("invalid object for big.matrix transposition")    
      }
    }
  } else {
    if(is.null(dir)) { dir <- getwd() }
  }
  slow.is <- 1 #GB  [ above this if tracker = NA, then do a progress bar ]
  nR <- nrow(bigMat); nC <- ncol(bigMat)
  if(is.na(tracker)) { if(estimate.memory(c(nR,nC))>slow.is) { tracker <- T } else { tracker <- F } }
  cN <- colnames(bigMat); rN <- rownames(bigMat)
  if(verbose) { cat(" creating",nC,"x",nR,"target matrix,",name,"...") }
  des <- paste(name,"dsc",sep=".")
  bck <- paste(name,"bck",sep=".")
  if(file.exists(cat.path(dir,bck)) & delete.existing) { unlink(cat.path(dir,bck)) }
  bigTrans <- big.matrix(nrow=nC,ncol=nR, backingfile=bck,
                         backingpath=dir, descriptorfile=des)
  if(verbose) { cat("done\n"); cat("\nAdding names\n") }
  options(bigmemory.allow.dimnames=TRUE)
  colnames(bigTrans) <- rN
  if(verbose) { cat(" added colnames\n") }
  rownames(bigTrans) <- cN
  if(verbose) { cat(" added rownames\n") }
  d2 <- d1 <- 0
  #try({
  split.to <- 10*round(estimate.memory(bigMat)) # split into .1GB chunks, save RAM without creating groups too small to process
  #if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
  stepz <- round(seq(from=1,to=nC+1,length.out=round((split.to+1))))
  # if the last step is too small, merge with previous
  sc.lst <- head(tail(stepz,2),1); lst <- tail(stepz,1)
  if((tail(stepz,1)) != nC+1) { stepz <- c(stepz,nC+1) }
  if(lst-sc.lst <=2) { ll <- length(stepz); if(ll>2) { stepz <- stepz[-(ll-1)] } else { warning("number of columns quite small, may cause issues") } }
  #
  split.to <- length(stepz)-1
  if(verbose) { cat(" transposing 'bigMat' into new big.matrix object:\n") }
  if(is.na(max.gb)) { max.gb <- split.to + 10 }
  for (cc in 1:split.to)
  {
    # within submatrix cols
    c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
    # do the copying
    lilColRange <- c(c1:c2)
    if(tracker) {      loop.tracker(cc,split.to) }
    if(is.finite(sum(as.numeric(lilColRange)))) {
      #cat(range(lilColRange)); cat(dim(bigTrans)); cat(dim(bigMat))
      bigTrans[lilColRange,1:nR] <- t(bigMat[1:nR,lilColRange])
    } else {
      cat(" Warning: empty interval ignored\n")
    }
    if(cc %% (max.gb*10) == 0) {
      # reset memory after every 'max.gb' 1GB chunks to prevent skyrocketing RAM use #
      fl.suc <- bigmemory::flush(bigTrans) ;  if(!fl.suc) { cat("flush failed\n") } ; gc()  
      if(T) {
        RR <- bigmemory::describe(bigTrans); rm(bigTrans); bigTrans <- attach.big.matrix(RR,path=dir)
      }
    }
  }
  #})
  if(verbose) { cat(" combining complete, converting result to big matrix\n") }
  descr <- bigmemory::describe(bigTrans)
  bigmemory::flush(bigTrans) # hopefully this will ensure the row/colnames are added to the file backing
  
  if(verbose) {
    cat(paste(" created big.matrix description file:",des,"\n"))
    cat(paste(" created big.matrix backing file:",bck,"\n"))
    if(length(R.descr)>0 & all(R.descr!="")) { 
      save(descr,file=cat.path(dir,R.descr)) 
      cat(paste(" created big.matrix binary description file:",basename(R.descr),"\n"))
    }
  }
  return(bigTrans)
}


### INTERNAL FUNCTIONS ###
#' Internal
#' pm is 'preserve.median'
PC.fn.mat <- function(next.rows, nPCs, add.int=FALSE, pm=FALSE)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time)
  col.sel <- 1:ncol(next.rows)
  for (dd in 1:nrow(next.rows)) {
    # compiled PC.fn should speed up these ops a little
    next.rows[dd,] <- PC.fn(next.rows[dd,],nPCs,col.sel,add.int=add.int,pm=pm) 
  }  
  return(next.rows)
}


#' Internal
#' pm is 'preserve.median'
PC.fn.mat.apply <- function(nextrows, nPCs, add.int=F, pm=FALSE)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time), vectorized version
  # testing shows the for-loop (non-vectorized) to be slightly faster, maybe because of t()
  # when using PC.fn.2 must pass in vec of 1's if you want the intecept
  nc <- ncol(nextrows); if(nc<1) { warning("nextrows had less than 1 column"); return(rep(0,length(nextrows))) }
  col.sel <- 1:nc
  nextrows <- t(apply(nextrows,1,PC.fn.2,nPCs=nPCs,col.sel=col.sel,add.int=add.int,pm=pm))
  return(nextrows)
}


#' Internal
PC.fn.mat.multi <- function(nextrows, nPCs, mc.cores=1, add.int=F, pm=FALSE)
{
  # matrix version of PC.fn (used to PC-correct one SNP at a time), vectorized version
  # testing shows the for-loop (non-vectorized) to be slightly faster, maybe because of t()
  # when using PC.fn.2 must pass in vec of 1's if you want the intecept
  #mc.cores <- 2
  col.sel <- 1:ncol(nextrows)
  nextrows <- lapply(seq_len(nrow(nextrows)), function(i) nextrows[i,]) # multi slows this down
  #nextrows <- multicore::mclapply(nextrows,PC.fn,nPCs=nPCs,col.sel=col.sel,mc.cores=mc.cores)
 # cat("~")
 # save(nextrows,nPCs,col.sel,mc.cores,add.int,pm,file="TESTOPOO.RData"); stop()
  nextrows <- parallel::mclapply(nextrows,PC.fn.2,nPCs=nPCs,col.sel=col.sel,mc.cores=mc.cores, add.int=add.int, pm=pm)
 # cat("+")
  nextrows <- do.call("rbind",nextrows)
 # cat("`")
  return(nextrows)
}


#' Internal
PC.fn.previous <- function(next.row, nPCs, col.sel, add.int=F)
{
  # apply PC correction for a single SNP, allowing for missing data.
  bad1 <- which(is.na(next.row))
  if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
  if(add.int) {  int <- mean(next.row[sel]) } else { int <- 0 }
  next.row[sel] <- lm(next.row ~ nPCs,na.action="na.exclude")$residuals + int
  return(next.row)
}

PC.fn <- function(next.row, nPCs, col.sel, add.int=F, pm=FALSE)
{
  # apply PC correction for a single SNP, allowing for missing data.
  bad1 <- which(is.na(next.row))
  if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
  if(add.int) { 
    if(pm) { 
      int <- median(next.row[sel])
      next.row[sel] <- lm(next.row ~ nPCs,na.action="na.exclude")$residuals
      next.row[sel] <- next.row[sel] - median(next.row[sel]) + int
    } else { 
      int <- mean(next.row[sel]) 
      next.row[sel] <- lm(next.row ~ nPCs,na.action="na.exclude")$residuals + int
    }
  } else { 
    next.row[sel] <- lm(next.row ~ nPCs,na.action="na.exclude")$residuals
  }
  return(next.row)
}


#' Internal
PC.fn.2 <- function(next.row, nPCs, col.sel, add.int=F, pm=FALSE)
{
  # apply PC correction for a single SNP, allowing for missing data.
  # when using PC.fn.2 must pass in vec of 1's if you want the intecept
  bad1 <- which(is.na(next.row))
  if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
  if(add.int) { 
    if(pm) { 
      int <- median(next.row[sel])
      next.row[sel] <- lm.fit(x=nPCs[sel,],y=next.row[sel])$residuals
      next.row[sel] <- next.row[sel] - median(next.row[sel]) + int
    } else { 
      int <- mean(next.row[sel]) 
      next.row[sel] <- lm.fit(x=nPCs[sel,],y=next.row[sel])$residuals + int
    }
  } else { 
    #prv(nPCs[sel,],next.row[sel])
    next.row[sel] <- stats::lm.fit(x=nPCs[sel,],y=next.row[sel])$residuals
    #next.row[sel] <- stats::lm.fit(x=sel*2,y=sel+rnorm(length(sel)))$residuals
    #next.row[sel] <- rep(1,length(sel))
  }
  return(next.row)
}


#' Internal
PC.fn.2.previous <- function(next.row, nPCs, col.sel, add.int=F, pm=FALSE)
{
  # apply PC correction for a single SNP, allowing for missing data.
  # when using PC.fn.2 must pass in vec of 1's if you want the intecept
  bad1 <- which(is.na(next.row))
  if(length(bad1)>0) { sel <- -bad1 } else { sel <- col.sel }
  if(add.int) { int <- mean(next.row[sel]) } else { int <- 0 }
  next.row[sel] <- lm.fit(x=nPCs[sel,],y=next.row[sel])$residuals + int
  return(next.row)
}


#' Internal
matmul <- function(A, x, transpose=FALSE)
{
  # bigalgebra friendly version of matrix multiplier function
  if(transpose) {
    return(t( t(x) %*% A)) 
  } else {
    return (A %*% x) 
  }
}


