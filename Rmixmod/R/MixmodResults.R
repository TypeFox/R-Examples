###################################################################################
##                                MixmodResults.R                                ##
###################################################################################

###################################################################################
##' @include global.R
##' @include Parameter.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{MixmodResults}}] class
##'
##' This is a class to contain results from MIXMOD library.
##'  
##' \describe{
##'   \item{nbCluster}{integer. It indicates the number of components.}
##'   \item{model}{character. Name of the model.}
##'   \item{criterion}{list of character. This option permits to select the criterion giving the best configuration of an execution.}
##'   \item{criterionValue}{numeric. Values of the criterion.}
##'   \item{parameters}{a S4 [\code{\linkS4class{Parameter}}] object. The best model parameters.}
##'   \item{likelihood}{numeric. The model likelihood.}
##'   \item{partition}{vector of integers defining the partition.}
##'   \item{proba}{a matrix of probabilities.}
##'   \item{error}{a character. The mixmod error.}
##' }
##'
##' @examples
##'   getSlots("MixmodResults")
##'
##' @name MixmodResults-class
##' @rdname MixmodResults-class
##' @exportClass MixmodResults
##'
setClass(
    Class="MixmodResults",
    representation=representation(
        nbCluster = "numeric",
        model = "character",
        criterion = "character",
        criterionValue = "numeric",
        parameters = "Parameter",
        likelihood = "numeric",
        partition = "integer",
        proba = "matrix",
        error = "character"
    ),
    prototype=prototype(
        nbCluster = numeric(0),
        model = character(0),
        criterion = character(0),
        criterionValue = numeric(0),
        likelihood = numeric(0),
        partition = integer(0),
        proba = matrix(nrow=0,ncol=0),
        error = character(0)
    )
)
###################################################################################

###################################################################################
##' @rdname print-methods
##' @aliases print print,MixmodResults-method
##'
setMethod(
  f="print",
  signature=c("MixmodResults"),
  function(x,...){ 
    cat("* nbCluster   = ", x@nbCluster,"\n")      
    cat("* model name  = ", x@model, "\n")
    if ( x@error == "No error" ){
      cat("* criterion   = ", paste(x@criterion, "(", formatC(x@criterionValue,digits=4,format="f"), ")", sep="") ); cat("\n")
      cat("* likelihood  = ", formatC(x@likelihood,digits=4,format="f"), "\n")
      print(x@parameters)
    }else{
      cat("* No results. MIXMOD library stopped with following error: ", x@error, "!\n")
    }
  }
)
###################################################################################

###################################################################################
##' @rdname show-methods
##' @aliases show show,MixmodResults-method
##'
setMethod(
  f="show",
  signature=c("MixmodResults"),
  function(object){
    cat("* nbCluster   = ", object@nbCluster,"\n")      
    cat("* model name  = ", object@model, "\n")
    if ( object@error == "No error" ){
      cat("* criterion   = ", paste(object@criterion, "(", formatC(object@criterionValue,digits=4,format="f"), ")", sep="") ); cat("\n")
      cat("* likelihood  = ", formatC(object@likelihood,digits=4,format="f"), "\n")
      show(object@parameters)
    }else{
      cat("* No results. MIXMOD library stopped with the following error: ", object@error, "!\n")
    }
  }
)
###################################################################################

###################################################################################
##' @rdname summary-methods
##' @aliases summary summary,MixmodResults-method
##'
setMethod(
  f="summary",
  signature=c("MixmodResults"),
  function(object, ...){
    cat("**************************************************************\n")
    cat("*       Number of cluster = ", object@nbCluster,"\n")
    cat("*              Model Type = ", object@model,"\n")      
    if ( object@error == "No error" ){
      cat("*               Criterion = ", paste(object@criterion, "(", formatC(object@criterionValue,digits=4,format="f"), ")", sep="") ); cat("\n")
      cat("*              Parameters =  list by cluster\n")
      summary(object@parameters)
      cat("*          Log-likelihood = ",formatC(object@likelihood,digits=4,format="f"),"\n")     
    }else{
      cat("* No results. MIXMOD library stopped with the following error: ", object@error, "!\n")
    }
    cat("**************************************************************\n")
  }
)
###################################################################################

###################################################################################
##' Define function to draw an ellipse
##'
##' @param x an object of class [\code{\linkS4class{MixmodResults}}]
##' @param i an index of one variable from data
##' @param j an index of one variable from data
##' 
##' @keywords internal
##'
ellipse<-function(x, i, j){
  # loop over the cluster
  for ( k in 1:x@nbCluster ){
    angles <- seq(0, 2*pi, length.out=200)
    ctr<-c(x@parameters@mean[k,i],x@parameters@mean[k,j])
    A<-matrix(c(x@parameters@variance[[k]][i,i],x@parameters@variance[[k]][j,i],x@parameters@variance[[k]][i,j],x@parameters@variance[[k]][j,j]), nrow=2)
    eigVal  <- eigen(A)$values
    eigVec  <- eigen(A)$vectors
    eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
    xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
    yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
    ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
    ellRot  <- eigVec %*% t(ellBase)                                          # rotated ellipse
    
    
    lines((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, type="l", lwd=2)
    matlines(xMat, yMat, col=1, lty=2, lwd=1)
    points(ctr[1], ctr[2], pch=4, lwd=3)
  }
}
###################################################################################


###################################################################################
##' Plotting of a class [\code{\linkS4class{MixmodResults}}]  
##' 
##' Biplot of two variables from a quantitative data set. Use parameters and partition from a [\code{\linkS4class{MixmodResults}}] object to distinguish the different clusters.
##'
##' Ellipsoids (i.e. linear transformations of hyperspheres) 
##' centered at the mean can be drawn using the parameters computed by MIXMOD.
##' The directions of the principal axes of the ellipsoids are given by the eigenvectors of the covariance matrix \eqn{\Sigma}. 
##' The squared relative lengths of the principal axes are given by the corresponding eigenvalues.
##'
##' @param x an object of class [\code{\linkS4class{MixmodResults}}]
##' @param data a data frame containing a quantitative data set.
##' @param variable1 index or character containing the name of the first variable. First column of data by default.
##' @param variable2 index or character containing the name of the second variable. Second column of data by default.
##' @param col a specification for the default plotting color. By default partition is used to separate clusters with different colors.
##' @param pch either an integer specifying a symbol or a single character to be used as the default in plotting points. By default partition is used to seperate clusters with different symbols.
##' @param xlab a title for the x axis. Variable1 by default.
##' @param ylab a title for the y axis. Variable2 by default.
##' @param add.ellipse a boolean. Add ellipses to graph. TRUE by default.
##' @param ... further arguments passed to or from other methods
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##' @examples
##'   data(geyser)
##'   xem1 <- mixmodCluster(geyser,3)
##'   plotCluster(xem1["bestResult"], geyser)
##'
##'   data(iris)
##'   xem2 <- mixmodCluster(iris[1:4],2:6)
##'   plotCluster(xem2["bestResult"], iris, variable1="Sepal.Length", variable2="Sepal.Width")
##'   plotCluster(xem2["bestResult"], iris, variable1=1, variable2=4)
##'
##' @seealso \code{\link{plot}}
##' @export
##'
plotCluster <- function(x, data, variable1=colnames(data)[1], variable2=colnames(data)[2], col=x@partition+1, pch=x@partition, xlab=variable1, ylab=variable2, add.ellipse=TRUE, ...){
  if ( !is(x,"MixmodResults") )
    stop("x must be a MixmodResults object!")
  if ( !is(x@parameters, "GaussianParameter") )
    stop("x must contains Gaussian parameters!")
  if ( !is.matrix(data) & !is.data.frame(data) )
    stop("data must be a data.frame or a matrix object!")

  
  # get variable1 index
  if (is.numeric(variable1)){
    if (variable1>ncol(data))
      stop("variable1 index mismatch the data frame dimension")
    else{
      index1<-variable1
      xlab<-colnames(data)[variable1]
    }
  }
  else{
    if ( !(variable1 %in% colnames(data)) )
      stop("variable1 is unknown!")
    else
      index1<-which(colnames(data) == variable1)
  }
  # get variable2 index
  if (is.numeric(variable2)){
    if (variable2>ncol(data))
      stop("variable2 index mismatch the data frame dimension")
    else{
      index2<-variable2
      ylab<-colnames(data)[variable2]
    }
  }
  else{
    if ( !(variable2 %in% colnames(data)) )
      stop("variable2 is unknown!")
    else
      index2<-which(colnames(data) == variable2)
  }

  plot(data[,index1],data[,index2],col=col,pch=pch, xlab=xlab, ylab=ylab, ...)
  if ( add.ellipse ) ellipse(x,index1,index2)
}
###################################################################################


###################################################################################
##' Histogram of a class [\code{\linkS4class{MixmodResults}}]  
##' 
##' Histograms of data object using parameters from a [\code{\linkS4class{MixmodResults}}]
##' to plot densities.
##'
##' Data with the density of each cluster and the mixture density are drawn for each variable.
##'
##' @param x an object of class [\code{\linkS4class{MixmodResults}}]
##' @param data a vector or data frame containing a quantitative data set.
##' @param variables list of variables names (or indices) to compute a histogram. All variables from data by default.
##' @param xlab a list of title for the x axis. xlab must have the same length than variables.
##' @param main a list of title for the histogram. main must have the same length than variables.
##' @param ... further arguments passed to or from other methods
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##' @examples
##'   data(geyser)
##'   xem1 <- mixmodCluster(geyser,3)
##'   \dontrun{ histCluster(xem1["bestResult"], geyser) }
##'   histCluster(xem1["bestResult"], geyser, variables=1)
##'
##' @seealso \code{\link{hist}}
##' @export
##'
histCluster <- function(x, data, variables=colnames(data), xlab=rep("",length(variables)), main=paste("Histogram of",variables), ...){
  # check the options
  if ( !is(x,"MixmodResults") )
  stop("'x' must be a MixmodResults object!")
  if ( !is.matrix(data) & !is.data.frame(data) & !is.vector(data) )
  stop("'data' must be a vector or a data.frame or a matrix object!")
  if ( (length(variables) == 0) & (ncol(data)>1))
  stop("'variables' is empty!")
  if ( length(variables)>ncol(data) )
  stop("List of variables too long!")

  # get the indices of variables
  if (is.numeric(variables)){
    if (max(variables)>ncol(data))
      stop("At least one variable index mismatch the data frame dimension")
    else{
      indices<-variables
      main=paste("Histogram of",colnames(data)[variables])
    }
  }
  else{
    if ( sum(!(variables %in% colnames(data))) )
      stop("At least one variable is unknown!") 
    else{
      if ( ncol(data)==1 ){ indices<-1 }
      else { indices<-which(colnames(data) %in% variables) }
    }
  }
  nvar<-length(indices)
  
  # get old par 
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  
  if ( is(x@parameters, "GaussianParameter") ){
    if ( isQualitative(data) )
      stop("data must contain only quantitative variables!")
    # split the layout
    if( nvar < 4 & nvar > 1) par( mfrow = c( 1, nvar ) )
    else if ( nvar >= 4 ){
      nrow<-round(sqrt(nvar))
      if (is.wholenumber(sqrt(nvar))) ncol<-sqrt(nvar)
      else ncol<-sqrt(nvar)+1
      par( mfrow = c( nrow, ncol ) ) 
    }
    i<-1
    # loop over variables
    for (j in indices ){
      
      xaxis<-seq(min(data[,j]),max(data[,j]),by=0.0001)
      density<-matrix(nrow=x@nbCluster,ncol=length(xaxis))
      
      # loop over the clusters to generate densities
      for( k in 1:x@nbCluster ){
        density[k,]<-x@parameters["proportions",k]*dnorm(xaxis,x@parameters["mean",k][j],sqrt(x@parameters["variance",k][j,j]))
      }
      # generate mixture density
      mixture<-apply(density,2,sum)
      h<-hist(data[,j], xlab=xlab[i], main=main[i], ...)
      
      ratio<-max(h$counts)/max(mixture)
      density<-density*ratio
      mixture<-mixture*ratio
      
      lines(xaxis,mixture,col="azure4",lty=1,lwd=4)
      for( k in 1:x@nbCluster ){
        lines(xaxis,density[k,],col=k+1,lty=2,lwd=2)
      }
      i<-i+1
    }
  }
  else if ( is(x@parameters, "MultinomialParameter") ){
    stop("x must contain Gaussian parameters. See barplot() to plot multinomial parameters.")
  }
  else if ( is(x@parameters, "CompositeParameter") ){
    stop("x must contain Gaussian parameters. No plot device for composite parameters.")
  }
  else{
    stop("Uknown type of parameters!")
  }
  par(op)
}
###################################################################################



###################################################################################
##' Barplot of a class [\code{\linkS4class{MixmodResults}}]  
##' 
##' Barplot of qualitative data object using parameters from a [\code{\linkS4class{MixmodResults}}]
##' to plot probablities of modalities.
##'
##' Each line corresponds to one variable. A barplot is drawn for each cluster with the probabilities for 
##' each modality to be in that cluster.
##'
##' @param x an object of class [\code{\linkS4class{MixmodResults}}]
##' @param data a vector or data frame containing a qualitative data set.
##' @param variables list of variables names (or indices) to compute a barplot. All variables from data by default.
##' @param main a list of title for the barplot. main must have the same length than variables.
##' @param ... further arguments passed to or from other methods
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##' @examples
##'   data(birds)
##'   xem <- mixmodCluster(birds,2)
##'   barplotCluster(xem["bestResult"], birds)
##'   barplotCluster(xem["bestResult"], birds, variables=c(2,3,4))
##'   barplotCluster(xem["bestResult"], birds, variables=c("eyebrow","collar"))
##'
##' @seealso \code{\link{barplot}}
##' @export
##'
barplotCluster <- function(x, data, variables=colnames(data), main=paste("Barplot of",variables), ...){
  # check the options
  if ( !is(x,"MixmodResults") )
    stop("'x' must be a MixmodResults object!")
  if ( !is.matrix(data) & !is.data.frame(data) & !is.vector(data) )
    stop("'data' must be a vector, a data.frame or a matrix object!")
  if ( (length(variables)==0) & (ncol(data)>1) )
    stop("'variables' is empty!")
  if ( length(variables)>ncol(data) )
    stop("List of variables too long!")

  # get the indices of variables
  if (is.numeric(variables)){
    if (max(variables)>ncol(data))
      stop("At least one variable index mismatch the data frame dimension")
    else{
      indices<-variables
      main=paste("Barplot of",colnames(data)[variables])
    }
  }
  else{
    if ( sum(!(variables %in% colnames(data))) )
      stop("At least one variable is unknown!") 
    else{
      if ( ncol(data)==1 ){ indices<-1 }
      else { indices<-which(colnames(data) %in% variables) }
    }
  }
  nvar<-length(indices)

  # get old par 
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  
  
  if ( is(x@parameters, "GaussianParameter") ){
    stop("x must contain multinomial parameters. See hist() to plot Gaussian parameters.")
  }
  else if ( is(x@parameters, "MultinomialParameter") ){
    if ( !isQualitative(data) )
      stop("data must contain only qualitative variables!")      
    if ( is.factor(data) ) data<-as.integer(data)
    else if ( !is.vector(data) ){
      # loop over columns to check whether type is factor
      for ( j in 1:ncol(data) ){
        if ( is.factor(data[,j]) ) data[,j] <- as.integer(data[,j])
      }
    }
    # split the layout
    #par( mfrow = c( nvar, x@nbCluster+1 ) )
    #par(mar = par("mar")*.75)
    # split the layout
    if( nvar < 4 & nvar > 1) par( mfrow = c( 1, nvar ) )
    else if ( nvar >= 4 ){
      nrow<-round(sqrt(nvar))
      if (is.wholenumber(sqrt(nvar))) ncol<-sqrt(nvar)
      else ncol<-sqrt(nvar)+1
      par( mfrow = c( nrow, ncol ) ) 
    }
    # get number of observations
    nobs<-nrow(data)
    i<-1
    # loop over variables
    for (j in indices ){
      f<-x@parameters@factor[j]
      t<-table(data[,j])
      if ( f != length(t) ){
        g<-numeric(f)
        g[which(1:f %in% names(t))]<-t
        t<-g
      }
      names(t)<-1:f
      # freq<-barplot(as.vector(t), names.arg=names(t), xlab=xlab[i], main="", ...)
      proba<-matrix(0,nrow=x@nbCluster,ncol=f)
      for ( k in 1:x@nbCluster){
        center_k<-x@parameters@center[k,j]
        proba_k<-x@parameters@scatter[[k]][j,]
        proba_k[center_k]<-1-proba_k[center_k]
        proba[k,]<-proba_k[1:f]
        #prob<-barplot(proba_k[1:f], names.arg=names(t), main="", col=k+1, ylim=c(0,1), ...)
        #title(paste("Cluster",k),cex.main=1)
        #text(prob,proba_k[1:f]+(max(proba_k)/10),ifelse(proba_k[1:f]<0.01,format(proba_k[1:f],scientific=TRUE,digits=2),round(proba_k[1:f],2)),xpd=TRUE,font=2)
      }
      # Do the barplot and save the bar midpoints
      mp<-barplot(proba, beside = TRUE, axisnames = FALSE, names.arg=names(t), main="", ylab="Conditional frequency", ylim=c(0,1))
      # add unconditional frequencies
      for ( k in 1:f){
        lines(c(min(mp[,k])-1,max(mp[,k])+1),rep(t[k]/nobs,2),lty=2)
      }
      # add unconditional frequency legend only for the first variable
      if ( i == 1 ) text(min(mp[,1])-.5,t[1]/nobs,labels="Unconditional frequency",cex=1, pos=4, adj=c(0,0), font=2) 
      # add title
      title(main[i],cex.main=1)
      # Add the individual bar labels
      mtext(1, at = mp, text = paste("C",1:x@nbCluster), line = 0, cex = 0.5)
      # Add the group labels for each pair
      #mtext(1, at = rbind(mp[1,],colMeans(mp[-1,])), text = rep(c("freq", "proba"), f), line = 1, cex = 0.75)
      # Add the labels for each group
      mtext(1, at = colMeans(mp), text = names(t), line = 2)
      i<-i+1
    }
  }
  else if ( is(x@parameters, "CompositeParameter") ){
    stop("x must contain multinomial parameters. No plot device for composite parameters.")
  }
  else{
    stop("Uknown type of parameters!")
  }
  par(op)
}
###################################################################################



###################################################################################
##' @rdname extract-methods
##' @aliases [,MixmodResults-method
##'
setMethod(
  f="[", 
  signature(x = "MixmodResults"),
  definition=function (x, i, j, drop) {
    if ( missing(j) ){
      switch(EXPR=i,
        "nbCluster"={return(x@nbCluster)},
        "criterion"={return(x@criterion)},
        "criterionValue"={return(x@criterionValue)},
        "model"={return(x@model)},
        "likelihood"={return(x@likelihood)},
        "parameters"={return(x@parameters)},
        "partition"={return(x@partition)},
        "proba"={return(x@proba)},
        "error"={return(x@error)},
        stop("This attribute doesn't exist !")
      )
    }else{
      switch(EXPR=i,
        "criterion"={return(x@criterion[j])},
        "criterionValue"={return(x@criterionValue[j])},
        stop("This attribute doesn't exist !")
      )
    }
  }
)
##################################################################################
