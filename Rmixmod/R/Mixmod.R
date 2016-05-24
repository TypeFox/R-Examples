###################################################################################
##                               Mixmod.R                                        ##
###################################################################################

###################################################################################
##' @include global.R
##' @include MixmodResults.R
##' @include GaussianModel.R
##' @include MultinomialModel.R
NULL
###################################################################################


###################################################################################
##' Constructor of [\code{\linkS4class{Mixmod}}] class
##' 
##' This is a class to run mixmod library.
##'
##' \describe{
##'   \item{data}{numeric vector or a data frame of observations. Can be qualitative,quantitative or both(heterogeneous)}
##'   \item{dataType}{character. Type of data. It defines whether data is quantitative, qualitative or composite}
##'   \item{nbCluster}{integer. It indicates the number of classes.}
##'   \item{knownLabels}{numeric. It contains the known labels.}
##'   \item{weight}{numeric vector with n (number of individuals) rows. Weight is optionnal. This option is to be used when weight is associated to the data.}
##'   \item{nbVariable}{integer. The number of variables.}
##'   \item{nbSample}{integer. The number of observations.}
##'   \item{criterion}{list of character. This option permits to select the criterion giving the best configuration of an execution.}
##'   \item{models}{a S4 [\code{\linkS4class{Model}}] object. Defining the list of models to be tested.}
##'   \item{error}{logical. Say if at least one model finished with no error in MIXMOD.}
##'   \item{results}{a list of S4 [\code{\linkS4class{MixmodResults}}] object containing all results. Results are sorted into a ascending order according to the first criterion (descending order for the CV criterion). This order can be changed by using the sortByCriterion() method.}
##' }
##'
##' @examples
##'   getSlots("Mixmod")
##'
##' @name Mixmod-class
##' @rdname Mixmod-class
##' @exportClass Mixmod
##'
setClass(
  Class="Mixmod",
  representation=representation(
    data = "matrix",
    dataType = "character",
    factor = "numeric",
    nbCluster = "numeric",
    knownLabels = "integer",
    weight = "numeric",
    nbVariable = "integer",
    nbSample = "integer",
    criterion = "character",
    models = "Model",
    error = "logical",
    results = "list",
    "VIRTUAL"
  ),
  prototype=prototype(
    data = matrix(nrow=0,ncol=0),
    dataType = character(0),
    factor = numeric(0),
    nbCluster = numeric(0),
    knownLabels = integer(0),
    weight = numeric(0),
    nbVariable = integer(0),
    nbSample = integer(0),
    criterion = character(0),
    error = TRUE,
    results = list()
  ),
  # define validity function
  validity=function(object){
    # check for missing values
    if ( sum(is.na(object@data)) ){
      stop("data contains NA. Mixmod cannot deal with missing values!")
    }
    # check qualitative data set
    if ( object@dataType == "qualitative" ){
      # check factors
      if ( length(object@factor) != ncol(object@data) )
        stop("number of factors too small!")
      if ( sum(is.na(object@factor)) )
        stop("factor contains NA!")
      if ( sum(!is.wholenumber(object@factor)) )
        stop("factor must be integer!")
      if ( sum(is.na(object@data)) )
        stop("data contains NA!")
      if ( sum(!is.wholenumber(object@data)) )
        stop("data contains real values!")
      # get max of modalities
      #max_mod<-apply(object@data,2,max)
      #if ( sum(max_mod>object@factor) )
      #  stop("At least one modality within 'data' is greater than the number of modalities in 'factor'!")
      # get min of modalities
      #min_mod<-apply(object@data,2,min)
      #if ( sum(min_mod<1) )
      #  stop("At least one modality within 'data' is lower than 1!")
    }
    # check data type
    if ( (object@dataType != "quantitative") & (object@dataType != "qualitative") & (object@dataType != "composite") ){
      stop("unknown dataType --> dataType must be 'quantitative', 'qualitative' or 'composite'!")
    }
    # check models validity
    if ( (object@dataType == "quantitative") & !is(object@models, "GaussianModel") ){
      stop("Incorrect model specified for quantitative data !\n")
    }
    if ( (object@dataType == "qualitative") & !is(object@models, "MultinomialModel")){
      stop("Incorrect model specified for qualitative data !\n")
    }
    if ( (object@dataType == "composite") & !is(object@models, "CompositeModel")){
      stop("Incorrect model specified for heterogeneous data !\n")
    }
    # check dimensions of knownLabels
    if(length(object@knownLabels)>0){
      if ( (length(object@knownLabels)!= object@nbSample) ){
        stop("the length of knownLabels is not similar to the number of observations!")
      }   
    }
    # check if weight isn't null and equal to the number of subjects
    if(length(object@weight)>0){
      if(length(object@weight) != object@nbSample){
        stop("weight must be equal to the number of sample !")            
      }
    }
    return(TRUE)
  }
)
###################################################################################


###################################################################################
##' Create an instance of the [\code{\linkS4class{Mixmod}}] class using new/initialize.
##' 
##' Initialization method. Used internally in the `Rmixmod' package.
##' 
##' @seealso \code{\link{initialize}}
##'
##' @keywords internal
##'
##' @rdname initialize-methods
##'
setMethod(
  f="initialize",
  signature=c("Mixmod"),
  definition=function(.Object,data,dataType,models,weight,knownLabels){
    # set missing parameters as NULL
    if (missing(dataType)) dataType<-NULL
    if (missing(models)) models<-NULL
    if (missing(weight)) weight<-NULL
    if (missing(knownLabels)) knownLabels<-NULL
    
    if(!missing(data)){
      if ( is.null(dataType) ){
        .Object@dataType <- is.dataType(data)
      }
      else{
        .Object@dataType <- dataType
      }
      if ( .Object@dataType != "quantitative" ){ 
        .Object@factor <- nbFactorFromData(data)
      
        if ( is.factor(data) ){
          data<-as.integer(data) 
        }
        else if ( is.data.frame(data) ){
          # loop over columns to check whether type is factor
          for ( j in 1:ncol(data) ){
            if ( is.factor(data[,j]) ) data[,j] <- as.integer(data[,j])
          }
        }
      }
      # for quantitative data
      if( .Object@dataType == "quantitative" ){
        if( is.null(models) ){
          .Object@models = new("GaussianModel", listModels="Gaussian_pk_Lk_C", family="general", free.proportions=TRUE, equal.proportions=FALSE)
        }else{
          .Object@models <- models
        }
      }
      # for qualitative data
      else if ( .Object@dataType == "qualitative" ){
        if( is.null(models) ){
          .Object@models = new("MultinomialModel", listModels="Binary_pk_Ekjh", free.proportions=TRUE, equal.proportions=FALSE)
        }else{
          .Object@models <- models
        }
      }
      #for composite data
      else if ( .Object@dataType == "composite" ){
        if( is.null(models) ){
          .Object@models = new("CompositeModel", listModels="Heterogeneous_pk_Ekjh_Lk_Bk", free.proportions=TRUE, equal.proportions=FALSE)
        }else{
          .Object@models <- models
        }
      }
      .Object@data <- as.matrix(data)
      .Object@nbSample <- nrow(.Object@data)
      .Object@nbVariable <- ncol(.Object@data)      
      if(!is.null(weight)){
        .Object@weight <- weight
      }
      if(!is.null(knownLabels)){
        knownLabels[knownLabels==0]<-NA
        .Object@knownLabels <- as.integer(as.factor(knownLabels))
        .Object@knownLabels[which(is.na(.Object@knownLabels))]<-as.integer(0)
      }
      # call validity function
      validObject(.Object)
    }
    return(.Object)
  }
)
###################################################################################



###################################################################################
##' @rdname print-methods
##' @aliases print print,Mixmod-method
##'
setMethod(
  f="print",
  signature=c("Mixmod"),
  function(x,...){
    cat("* nbCluster = ", x@nbCluster, "\n")
    cat("* criterion = ", x@criterion, "\n")
    print(x@models)
    if(length(x@weight)!=0){ 
      cat("* weight  = ", x@weight, "\n")
    }
    if(length(x@data)!=0){
      cat("* data      =\n")
      print(x@data)
    }
    else{}
    if ( length(x@knownLabels)>0 ){
      cat("* knownLabels = ", x@knownLabels, "\n")
    }
  }
)
###################################################################################

###################################################################################
##' @rdname show-methods
##' @aliases show show,Mixmod-method
##'
setMethod(
  f="show",
  signature=c("Mixmod"),
  function(object){
    cat("* nbCluster = ", object@nbCluster, "\n")
    cat("* criterion = ", object@criterion, "\n")
    show(object@models)
    if(length(object@weight)!=0){ 
      cat("* weight  = ", object@weight, "\n")
    }
    if(length(object@data)!=0){
      nrowShow <- min(10,nrow(object@data))
      ncolShow <- min(10,ncol(object@data))
      cat("* data (limited to a 10x10 matrix) =\n")
      print(formatC(object@data[1:nrowShow,1:ncolShow]),quote=FALSE)
    }else{}
    cat("* ... ...\n")  
    if ( length(object@knownLabels)>0 ){
      if ( length(object@knownLabels)>10 ){
        cat("* knownLabels = ", object@knownLabels[1:10], " ...\n")
      }else{
        cat("* knownLabels = ", object@knownLabels, "\n")
      }
    }
  }
)
###################################################################################


###################################################################################
##' @rdname summary-methods
##' @aliases summary summary,Mixmod-method
##'
setMethod(
  f="summary",
  signature=c("Mixmod"),
  function(object, ...){
    if ( !object@error ){
      cat("**************************************************************\n")
      cat("* Number of samples    = ",object@nbSample,"\n")
      cat("* Problem dimension    = ",object@nbVariable,"\n")
      summary(object@bestResult)
    }
    return(invisible())
  }
)
###################################################################################

###################################################################################
##' Plotting of a class [\code{\linkS4class{Mixmod}}]  
##' 
##' Plotting data from a [\code{\linkS4class{Mixmod}}] object using parameters and partition
##' to distinguish the different clusters.
##'
##' For quantitative case, ellipsoids (i.e. linear transformations of hyperspheres) 
##' centered at the mean are drawn using the parameters computed by MIXMOD.
##' The directions of the principal axes of the ellipsoids are given by the eigenvectors of the covariance matrix \eqn{\Sigma}. 
##' The squared relative lengths of the principal axes are given by the corresponding eigenvalues.
##' A 1-dimensional representation of variables with the densities is drawn on the diagonal.
##'
##' For qualitative case, a Multiple Correspondance Analysis is performed to get a
##' 2-dimensional representation of the data set. Bigger symbol means that observations are similar.
##'
##' @param x an object of class [\code{\linkS4class{Mixmod}}]
##' @param y a list of variables to plot (subset). Variables names or indices. Only in a quantitive case.
##' @param ... further arguments passed to or from other methods
##'
##' @importFrom graphics plot
##' @name plot
##' @aliases plot plot,Mixmod-method
##' @docType methods
##' @rdname plot-methods
##' @exportMethod plot
##'
##' @seealso \code{\link{plot}}
##' @examples
##'   ## for quantitative case
##'   data(iris)
##'   xem <- mixmodCluster(iris[1:4],3)
##'   plot(xem)
##'   plot(xem,c(1,3))
##'   plot(xem,c("Sepal.Length","Sepal.Width"))
##'
##'   ## for qualitative case
##'   data(birds)
##'   xem2 <- mixmodCluster(birds,2)
##'   plot(xem2)
##'   legend("bottomleft",c("Cluster1","Cluster2"),col=c(2,3),pch=c(1,2))
##'
setMethod(
  f="plot",
  signature=c("Mixmod"),
  function(x, y, ...){
    # for quantitative data
    if ( x@dataType == "quantitative" ){
      # create layout
      if ( x@nbVariable == 1 ){
        stop("data has only one variable. Try hist() function to get a 1D representation of x.")
      }
      else if (!missing(y)){

        if (is.numeric(y)){
          if (max(y)>ncol(x@data))
            stop("y indices mismatch the data frame dimension")
        }
        else{
          if ( sum(y %in% colnames(x@data))!= length(y) ){
            stop(cat("unknown variable: ", paste(y[which(!(y %in% colnames(x@data)))]),"\n"))
          }
        }
        # get old par 
        op <- par(no.readonly = TRUE) # the whole list of settable par's.
        # changing marging
        par(mar = rep(2.5,4))
        # decreasing font size
        par(cex = .75)
        # create layout matrix
        #par( mfrow = c(x@nbVariable, x@nbVariable) )
        split.screen(c(length(y), length(y)))
        # create histogram on diagonal
        for ( i in 1:length(y) ){
          #par(mfg=c(i,i))
          screen(i+((i-1)*length(y)))
          histCluster(x@bestResult, x@data, variables=y[i])
        }
        if (length(y)>1){
          # create biplots
          for ( i in 2:length(y) ){
            for( j in 1:(i-1) ){
              #par(mfg=c(i,j))
              screen(j+((i-1)*length(y)))
              plotCluster(x@bestResult, x@data, variable1=y[j], variable2=y[i], ...)
            }
          }
        }
        close.screen(all.screens = TRUE)
        # restore plotting parameters
        par(op)

      }
      else{ 
        # get old par 
        op <- par(no.readonly = TRUE) # the whole list of settable par's.
        # changing marging
        par(mar = rep(2.5,4))
        # decreasing font size
        par(cex = .75)
        # create layout matrix
        #par( mfrow = c(x@nbVariable, x@nbVariable) )
        split.screen(c(x@nbVariable, x@nbVariable))
        # create histogram on diagonal
        for ( i in 1:x@nbVariable ){
          #par(mfg=c(i,i))
          screen(i+((i-1)*x@nbVariable))
          histCluster(x@bestResult, x@data, variables=colnames(x@data)[i])
        }
        # create biplots
        for ( i in 2:x@nbVariable ){
          for( j in 1:(i-1) ){
            #par(mfg=c(i,j))
            screen(j+((i-1)*x@nbVariable))
            plotCluster(x@bestResult, x@data, variable1=colnames(x@data)[j], variable2=colnames(x@data)[i], ...)
          }
        }
        close.screen(all.screens = TRUE)
        # restore plotting parameters
        par(op)
      }
    }
    # for qualitative data
    else if ( x@dataType == "qualitative" ){
      # create layout
      if ( x@nbVariable == 1 ){
        stop("data has only one variable. Try barplot() function to get a 1D representation of x.")
      }else{
        # create layout matrix
        par( mfrow = c(1, 1))
        # get binary matrix from x
        matX <- matrix2binary(as.data.frame(x@data))
        # get number of observations
        n <- dim(matX)[1]
        # get number of variables
        p <- ncol(x@data)
        
        Dc <- drop((rep(1, n)) %*% matX)
        Y <- t(t(matX)/(sqrt(p * Dc)))
        Y.svd <- svd(Y)
        individuals <- Y %*% Y.svd$v[, 2:3]/p

        # get unique points
        unique.ind<-unique(individuals)
        # get number of duplication for each individuals
        point.size<-numeric(nrow(unique.ind))
        for(i in 1:nrow(unique.ind) ){ 
          for(j in 1:nrow(individuals)){
            point.size[i]<-point.size[i]+sum(unique.ind[i,]==individuals[j,])
          }
          point.size[i]<-point.size[i]/2
        }
        # plotting the first 2 axes  
        plot( unique.ind[,2] ~ unique.ind[,1], cex=point.size,
        col=x@bestResult@partition[-which(duplicated(individuals))]+1, pch=x@bestResult@partition[-which(duplicated(individuals))], xlab='Axis 1', ylab='Axis 2', main='Multiple Correspondance Analysis', ... )
      }
    }
    else if ( x@dataType == "composite" ){
      stop("visualization for heterogeneous is not available yet")
    }
    invisible()
  }
)
###################################################################################



###################################################################################
##' Histograms of a class [\code{\linkS4class{Mixmod}}]  
##'
##' Histograms of quantitative data from a [\code{\linkS4class{Mixmod}}] object using parameters
##' to plot densities.
##'
##' Data with the density of each cluster and the mixture density are drawn for each variable.
##'  
##' @param x an object of class [\code{\linkS4class{Mixmod}}]
##' @param variables list of variables names (or indices) to compute a histogram. All variables from data by default.
##' @param ... further arguments passed to or from other methods
##'
##' @importFrom graphics hist
##' @name hist
##' @aliases hist hist,Mixmod-method
##' @docType methods
##' @rdname hist-methods
##' @exportMethod hist
##'
##' @seealso \code{\link{hist}}
##' @examples
##'   data(iris)
##'   xem <- mixmodCluster(iris[1:4],3)
##'   hist(xem)
##'   hist(xem,variables=c(1,3))
##'   hist(xem,variables=c("Sepal.Length","Sepal.Width"))
##'
setMethod(
  f="hist",
  signature=c("Mixmod"),
  function(x, ...){
    histCluster(x@bestResult, x@data, ...)
    invisible()
  }
)
###################################################################################


###################################################################################
##' Barplot of a class [\code{\linkS4class{Mixmod}}]  
##'
##' Barplot of qualitative data from a [\code{\linkS4class{Mixmod}}] object using parameters
##' to plot probablities of modalities.
##'
##' Each line corresponds to one variable. Barplot is drawn for each cluster with the probabilities for 
##' each modality to be in that cluster.
##'  
##' @param x an object of class [\code{\linkS4class{Mixmod}}]
##' @param variables list of variables names (or indices) to compute a histogram. All variables from data by default.
##' @param ... further arguments passed to or from other methods
##'
##' @importFrom graphics barplot
##' @name barplot
##' @aliases barplot barplot,Mixmod-method
##' @docType methods
##' @rdname barplot-methods
##' @exportMethod barplot
##'
##' @seealso \code{\link{barplot}}
##' @examples
##'   data(birds)
##'   xem2 <- mixmodCluster(birds,2)
##'   barplot(xem2)
##'   barplot(xem2,variables=c(2,3,4))
##'   barplot(xem2,variables=c("eyebrow","collar"))
##'
setMethod(
  f="barplot",
  signature=c("Mixmod"),
  function(height, ...){
    barplotCluster(height@bestResult, height@data, ...)
    invisible()
}
)
###################################################################################


###################################################################################
##' @rdname sortByCriterion-methods
##' @aliases sortByCriterion,Mixmod,character-method
##'
setMethod(
  f="sortByCriterion",
  signature=c("Mixmod","character"),
  definition=function(object,criterion){
    # check whether the number of criterion is greater than one
    if ( length(object@criterion) == 1 ){
      stop("Mixmod object contains only one criterion")
    }
    if ( !(criterion %in% object@criterion) ){
      stop(paste("No results for the criterion",criterion,"!"))
    }
    # sort by criterion only if there are more than one result
    if ( length(object@results) == 1 ){
      stop("only one result in Mixmod object")
    }
    # get criterion index
    index<-which(object@criterion == criterion) 
    # set the first result as the best one
    best <- object@results[[1]]@criterionValue[index]
    # loop over results
    for ( i in 2:length(object@results) ){
      x <- object@results[[i]]
      if ( x@error != "No error" ) break
      crit <- x@criterionValue[index]
      j <- i
      while ( j > 1 ){
        if ( criterion == "CV" ){
          if ( object@results[[j-1]]@criterionValue[index] < crit ){
            object@results[[j]] <- object@results[[j-1]]
          }else{
            break
          }
        }else{
          if ( object@results[[j-1]]@criterionValue[index] > crit ){
            object@results[[j]] <- object@results[[j-1]]
          }else{
            break
          }
        }
        j <- j-1
      }
      object@results[[j]] <- x 
    }
    object@bestResult<-object@results[[1]]
    object@criterion <- c(criterion,object@criterion[which(!(object@criterion %in% criterion))])
    return(object)
  }
)
###################################################################################

