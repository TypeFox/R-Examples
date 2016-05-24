
#' Summary of a boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object and produces a summary.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param ... currently not used
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return a list with summary information 
#' @method summary FDboost
#' @export
### similar to summary.mboost()
summary.FDboost <- function(object, ...) {
  
  ret <- list(object = object, selprob = NULL)
  xs <- selected(object)
  nm <- variable.names(object)
  selprob <- tabulate(xs, nbins = length(nm)) / length(xs)
  names(selprob) <- names(nm)
  selprob <- sort(selprob, decreasing = TRUE)
  ret$selprob <- selprob[selprob > 0]
  class(ret) <- "summary.FDboost"
  
  ### only show one unique offset value
  #ret$object$offset <- unique(ret$object$offset)
  
  return(ret)
}

#' Print a boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object and produces a print on the console.
#' 
#' @param x a fitted \code{FDboost}-object
#' @param ... currently not used
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return a list with information on the model 
#' @method print FDboost
#' @export
### similar to print.mboost()
print.FDboost <- function(x, ...) {
  
  cat("\n")
  if(!any(class(x)=="FDboostLong")){
    cat(cat("\t Model-based Boosting with Functional Response\n"))
  }else{
    cat("\t Model-based Boosting with Irregular Functional Response\n")
  }
  cat("\n")
  if (!is.null(x$call))
    cat("Call:\n", deparse(x$call, nlines=10), "\n\n", sep = "", nlines = 10)
  show(x$family)
  cat("\n")
  cat("Number of boosting iterations: mstop =", mstop(x), "\n")
  cat("Step size: ", x$control$nu, "\n")
  
  if(length(unique(x$offset))<10){
    cat("Offset: ", round(unique(x$offset), 3), "\n")
  }else{
    cat("Offset: ", round(unique(x$offset), 3)[1:3], "..." ,
        round(unique(x$offset), 3)[length(unique(x$offset))-3+1:3], "\n")
  }

  cat("Number of baselearners: ", length(variable.names(x)), "\n")
  cat("\n")
  invisible(x)
  
}



#' Prediction for boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object produced by \code{\link{FDboost}()} and produces 
#'  predictions given a new set of values for the model covariates or the original 
#'  values used for the model fit. This is a wrapper
#'  function for \code{\link[mboost]{predict.mboost}()}
#' 
#' @param object a fitted \code{FDboost}-object
#' @param newdata a named list or a data frame containing the values of the model 
#' covariates at which predictions are required.
#' If this is not provided then predictions corresponding to the original data are returned. 
#' If \code{newdata} is provided then it should contain all the variables needed for 
#' prediction, in the format supplied to \code{FDboost}, i.e., 
#' functional predictors must be supplied as matrices with each row corresponding to 
#' one observed function.
#' @param which a subset of base-learners to take into account for computing predictions 
#' or coefficients. If which is given (as an integer vector corresponding to base-learners) 
#' a list is returned. 
#' @param toFDboost logical, defaults to \code{TRUE}. In case of regular response in wide format 
#' (i.e. response is supplied as matrix): should the predictions be returned as matrix, or list 
#' of matrices instead of vectors
#' @param ...  additional arguments passed on to \code{\link[mboost]{predict.mboost}()}.
#' 
#' @seealso \code{\link{FDboost}} for the model fit 
#' and \code{\link{plotPredicted}} for a plot of the observed values and their predictions.
#' @return a matrix or list of predictions depending on values of unlist and which 
#' @method predict FDboost
#' @export
# predict function: wrapper for predict.mboost()
## <TODO> check which
predict.FDboost <- function(object, newdata = NULL, which = NULL, toFDboost = TRUE, ...){
  
  stopifnot(any(class(object)=="FDboost")) 
  # print("Prediction FDboost") 
  dots <- list(...)
  
  # toFDboost is only meaningful for array-data
  if(any(class(object) == "FDboostScalar") |  any(class(object) == "FDboostLong")) toFDboost <- FALSE

  if(!is.null(dots$aggregate) && dots$aggregate != "sum"){
    if(length(which) > 1 ) stop("For aggregate != 'sum', only one effect, or which=NULL are possible.")
    if(toFDboost & class(object)[1]=="FDboost"){ 
      toFDboost <- FALSE
      warning("Set toFDboost to FALSE, as aggregate!='sum'. Prediction is in long vector.")
    }
  }
  
  classObject <- class(object)
  class(object) <- "mboost"
  sel <- sort(unique(selected(object))) # which effects were selected
  
  ### <FIXME> does not work, as for bl bsignal(), bhist() and bconcurrent(), the 
  ### index is not selected
  #   # which variables do you need for the prediction
  #   allVariables <- c()
  #   for(i in which){
  #     allVariables <- c(allVariables, names(object$baselearner[[i]]$model.frame()))
  #     ## FIXME: save the index variables for bsignal(), bhist() and bconcurrent()
  #     if( grepl("bsignal", object$baselearner[[i]]$get_call() ) ){
  #       allVariables <- c(allVariables, object$baselearner[[i]]$get_call() )
  #     }
  #   }
  #   
  #   ## only keep the necessary variables
  #   allVariables <- unique(allVariables, all.vars(formula(object$timeformula)))
  #   newdata <- newdata[allVariables]
  
  ### Prepare data so that the function predict.mboost() can be used
  if(!is.null(newdata)){
    
    ## save time-variable (index over response)
    nameyind <- attr(object$yind, "nameyind")
    
    ## response observed on common grid / scalar response 
    if(!is.null(object$ydim)){
      # get the number of trajectories you want to predict
      n <- NROW(newdata[[1]])  ## FIXME check that this is not the time-variable
      
      lengthYind <- length(newdata[[nameyind]])
      #assign(nameyind, newdata[[nameyind]])
      
      # try to get more reliable information on n (number of trajectories)
      # and on lengthYind (length of time)
      # using the hmatrix-objects in newdata if available
      if(is.list(newdata) | is.data.frame(newdata)){
        classes <- lapply(newdata, class)
        alln <- c()
        alllengthYind <- c()
        for(i in 1:length(classes)){
          if( any(classes[[i]] == "hmatrix" ) ){
            # number of trajectories
            n <- length(unique(newdata[[i]][,2]) ) 
            alln <- c(alln, n)
            # number of different observation points
            lengthYind <- length(unique(newdata[[i]][,1]) ) 
            alllengthYind <- c(alllengthYind, length(unique(newdata[[nameyind]])))
          }
        }
        ## in case of scalar response, set lenth of yindex to 1
        if(object$ydim[2]==1){ 
          lengthYind <- 1 
          alllengthYind <- 1
        } 
        
        if( length(unique(alln))>1 ) stop("The hmatrix-objects in newdata imply differing numbers of trajectories.")
        if( length(unique(alllengthYind))>1 ) stop("The hmatrix-objects in newdata imply differing times or do not match the time variable.")
      }
      
      if(object$ydim[2]>1 && lengthYind==0) stop("Index of response, ", nameyind, ", must be specified and have length >0.")
      
      # dummy variable to fit the intercept
      newdata[["ONEx"]] <- rep(1.0, n)
      # dummy variable for time
      newdata$ONEtime <- rep(1.0, lengthYind)
      
      #message("Predict ", n, " x ", lengthYind," observations.")

      
    }else{ #### for response observed on irregular grid
      
      n <- 1
      lengthYind <- length(newdata[[nameyind]])
      if(is.list(newdata) | is.data.frame(newdata)){
        classes <- lapply(newdata, class)
        alllengthYind <- c(lengthYind)
        for(i in 1:length(classes)){
          if( any(classes[[i]] == "hmatrix" ) ){
            # total number of observation points
            lengthYind <- length(newdata[[i]][,1])
            alllengthYind <- c(alllengthYind, lengthYind)
          }
        }
        if( length(unique(alllengthYind))>1 ) stop("The hmatrix-objects in newdata imply differing times or do not match the time variable.")
      }
      
      # dummy variable to fit the intercept
      newdata[["ONEx"]] <- rep(1.0, lengthYind)
      # dummy variable for time
      newdata$ONEtime <- rep(1.0, lengthYind)
      
      ## message("Predict ", lengthYind, " observations in total.")
      
    } ## for response observed on irregular grid

    #     # Predict effect of offset: predOffset
    #     predOffset <- object$offsetVec # offset is just an integer 
    #     if(length(object$offsetVec)>1){ # offset is a smooth function
    #       if(!any(class(object)=="FDboostLong")){ # irregular response
    #         predOffset <- rep(object$predictOffset(newdata[[nameyind]]), each=n)
    #       }else{ # regular response 
    #         predOffset <- object$predictOffset(newdata[[nameyind]])
    #       }  
    #       names(predOffset) <- NULL
    #     }
    
    ### Predict effect of offset (scalar, regular, irregular):  
    # use the function predictOffset() on the new time-variable   
    predOffset <- object$predictOffset(newdata[[nameyind]]) 
    # for regular response: repeat offset accordingly 
    if(classObject[1]=="FDboost") predOffset <- rep(predOffset, each=n)
    
    ### In the case of bsignal(), bfpc(), bconcurrent() and bhist() it is necessary 
    # to add the index of the signal-matrix as attribute
    posBsignal <- c(grep("bsignal(", names(object$baselearner), fixed = TRUE), 
                    grep("bfpc(", names(object$baselearner), fixed = TRUE))
    posBconc <- grep("bconcurrent(", names(object$baselearner), fixed = TRUE)
    posBhist <- grep("bhist(", names(object$baselearner), fixed = TRUE)
    whichHelp <- which
    if(is.null(which)) whichHelp <- 1:length(object$baselearner)
    posBsignal <- whichHelp[whichHelp %in% posBsignal]
    posBconc <- whichHelp[whichHelp %in% posBconc]
    posBhist <- whichHelp[whichHelp %in% posBhist]
    
    if(length(c(posBsignal, posBconc, posBhist)) > 0){
      #if(!is.list(newdata)) newdata <- list(newdata)
      for(i in c(posBsignal, posBconc, posBhist)){ 
        xname <- object$baselearner[[i]]$get_names()[1] 
        ## if two ore more base-learners are connected by %X%, find the functional variable 
        if(grepl("%X", names(object$baselearner)[i])){
          form <- strsplit(object$baselearner[[i]]$get_call(), "%X")[[1]]
          findFun <- grepl("bhist", form) | grepl("bconcurrent", form) | grepl("bsignal", form) | grepl("bfpc", form)
          form <- form[findFun]
          if(sum(findFun)!=1){stop("Can only predict effect of one functional effect in %X% or %Xc%.")}
          xname <- object$baselearner[[i]]$get_names()[findFun][1]
        }
        # print(xname)
        indname <- attr(object$baselearner[[i]]$get_data()[[xname]], "indname") 
        #indname_all <- c(indname_all, indname)
        if(i %in% c(posBhist, posBconc)){
          indnameY <- attr(object$baselearner[[i]]$get_data()[[xname]], "indnameY")
        }else{
          indnameY <- NULL
        }

        attr(newdata[[xname]], "indname") <- indname
        attr(newdata[[xname]], "xname") <- xname
        attr(newdata[[xname]], "signalIndex") <-  if(indname!="xindDefault") newdata[[indname]] else seq(0,1,l=ncol(newdata[[xname]]))
        ## convert matrix to model matrix, so that as.data.frame(newdata[[xname]]) 
        ## retains matrix-structure if newdata is a list
        if( class(newdata[[xname]])[1] == "matrix" ) class(newdata[[xname]])[1] <- "model.matrix"
        
        if(i %in% c(posBhist, posBconc)){
          attr(newdata[[indnameY]], "indnameY") <-  indnameY
          attr(newdata[[xname]], "indexY") <-  if(indnameY!="xindDefault") newdata[[indnameY]] else seq(0,1,l=ncol(newdata[[xname]]))
          if(any(classObject=="FDboostLong")){
            id <- newdata[[ attr(object$id, "nameid") ]]
            attr(newdata[[xname]], "id") <-  id
          } 
        }       
        
        ### <FIXE> is this code still necessary?? changes necessary!
        #         ## <FIXME> quite ugly how to deal with %X%, is there a way to do this more generally?
        #         # save data of concurrent effects
        #         if(i %in% posBconc){
        #           newdataConc[[xname]] <- newdata[[xname]]
        #           if(grepl("%X%", names(object$baselearner)[i])){
        #             xname <- object$baselearner[[i]]$get_names() 
        #             xname <- xname[!xname %in% names(newdataConc)] # already there
        #             # print(xname)
        #             newdataConc[xname] <- newdata[xname]
        #           }
        #         } 
        #         
        #         # save data of historic effects
        #         if(i %in% posBhist){
        #           newdataHist[[xname]] <- newdata[[xname]]
        #           if(grepl("%X%", names(object$baselearner)[i])){
        #             xname <- object$baselearner[[i]]$get_names() 
        #             xname <- xname[!xname %in% names(newdataHist)] # already there
        #             # print(xname)
        #             newdataHist[xname] <- newdata[xname]
        #           }
        #         }
        
      } ## loop over posBsignal, ...
    } # end data setup for functional effects (adding index as attribute to the data)
    
    ###### Prediction using the function predict.mboost() 
    
    # Function to suppress the warning that user-specified 
    # offset of length>1 is not used in prediction, 
    # important when offset=NULL in FDboost() but not in mboost()
    muffleWarning1 <- function(w){
      if( any( grepl( "Offset not used for prediction when", w) ) ) 
        invokeRestart( "muffleWarning" )  
    }

    ## predict all effects together, model-inherent offset is included automatically
    if(is.null(which)){

      if(is.null(object$offsetFDboost) & !is.null(object$offsetMboost) ){ 
        # offset=NULL in FDboost, but not in mboost
        ## suppress the warning that the offset cannot be used if offset=NULL in FDboost
        ## as offset is predicted and included in prediction 
        predMboost <- withCallingHandlers(predict(object=object, newdata=newdata, which=NULL, ...), 
                                           warning = muffleWarning1)
        predMboost <- predMboost + predOffset
      }else{ # offset != NULL in FDboost -> treat offset like in mboost 
        predMboost <- predict(object=object, newdata=newdata, which=NULL, ...)
      }
      
      if(!is.null(dim(predMboost)) && dim(predMboost)[2] == 1) predMboost <- predMboost[,1]
      
      ## predict one or several effects separately   
    }else{   
      predMboost <- vector("list", length(which))
      for(i in seq_along(which)){
        if(which[i]!=0){
          # [,1] save vector instead of a matrix with one column
          predMboost[[i]] <- predict(object=object, newdata=newdata, which=which[i], ...)#[,1]
          if(!is.null(dim(predMboost[[i]])) && dim(predMboost[[i]])[2] == 1) predMboost[[i]] <- predMboost[[i]][,1]
          
          if(!which[i] %in% sel){
            predMboost[[i]] <- rep(0L, lengthYind*n)
          }
        }else{
          predMboost[[i]] <- predOffset # predict offset in case of which=0
        } 
      } 
    } # end else for !is.null(which) 
    
    attr(predMboost, "offset") <- predOffset
    
    
    ############################################
  }else{ # is.null(newdata)
    n <- object$ydim[1]
    lengthYind <- object$ydim[2]
    
    ## for response in long format
    if(is.null(object$ydim)){
      n <- 1
      lengthYind <- length(object$yind)
    }

    ## predict all effects together, include the offset into the prediction 
    if(is.null(which)){
      predMboost <- predict(object=object, newdata=NULL, which=NULL, ...)#[,1]
      if(!is.null(dim(predMboost)) && dim(predMboost)[2] == 1) predMboost <- predMboost[,1]
      # for which=NULL return prediction of all effects 
      # offset of original data-fit is included automatically 
      
      ## predict one or several effects separately   
    }else{
      predMboost <- vector("list", length(which))
      for(i in seq_along(which)){
        if(which[i]!=0){
          predMboost[[i]] <- predict(object=object, newdata=NULL, which=which[i], ...)#[,1]
          if(!is.null(dim(predMboost[[i]])) && dim(predMboost[[i]])[2] == 1) predMboost[[i]] <- predMboost[[i]][,1]
          if(!which[i] %in% sel){
            predMboost[[i]] <- rep(0L, lengthYind*n)
          }
        }else{
          predMboost[[i]] <- object$offset # return offset in case of which=0
        } 
      }
    } # end else for !is.null(which)   
    attr(predMboost, "offset") <- object$offset    
  }  ## end # is.null(newdata)
  
  ## remember the offset
  offsetTemp <- attr(predMboost, "offset") 
  
  
  ### unlist the prediction in case that just one effect is predicted
  if(length(which)==1){
    predMboost <- predMboost[[1]]
  } 
  
  ### return what mboost.predict would return, 
  ## i.e. vector for which=NULL, length(which)=1
  ## or a matrix for length(which)>1
  if(!toFDboost){
    if(length(which)<=1){
      return(predMboost)
    } else{
      ## check that all predictions have the same length
      stopifnot( all( length(predMboost[[1]])==lapply(predMboost, length) ) )
      ret <- matrix(unlist(predMboost), ncol=length(predMboost))
      colnames(ret) <- names(predMboost)
      attr(ret, "offset") <- offsetTemp
      return(ret)
    } 
  }
  
  # model-estimation in long format, just return vector/ list of vectors
  if(is.null(object$ydim)) return(predMboost)
  
  ### If response was in wide format in FDboost. i.e. response was given as matrix
  ## return matrix/ list of matrices 
  ## Reshape prediciton in vector to the prediction as matrix
  reshapeToMatrix <- function(pred){
    if(is.null(pred)) return(NULL)
    if(length(pred)==1) return(pred)
    
    predMat <- matrix(pred, nrow=n, ncol=lengthYind)
    return(predMat)
  }

  if(is.list(predMboost)){
    ret <- lapply(predMboost, reshapeToMatrix)
  }else{
    ret <- reshapeToMatrix(predMboost)
  }
  
  attr(ret, "offset") <- offsetTemp
  
  return(ret)
}



#' Fitted values of a boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object and computes the fitted values.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param toFDboost logical, defaults to \code{TRUE}. In case of regular response in wide format 
#' (i.e. response is supplied as matrix): should the predictions be returned as matrix, or list 
#' of matrices instead of vectors
#' @param ... additional arguments passed on to \code{\link{predict.FDboost}}
#' 
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return matrix of fitted values
#' @method fitted FDboost
#' @export
### similar to fitted.mboost() but returns the fitted values as matrix
fitted.FDboost <- function(object, toFDboost = TRUE, ...) {

  args <- list(...)
  
  if (length(args) == 0) {
    ## give back matrix for regular response and toFDboost==TRUE
    if(!any(class(object)=="FDboostLong") & toFDboost){
      ret <- matrix(object$fitted(), nrow=object$ydim[1])
    }else{ # give back a long vector
      ret <- object$fitted()
      names(ret) <- object$rownames
    }
  } else {
    ret <- predict(object, newdata=NULL, toFDboost=toFDboost, ...)
  }
  ret
}
 
#' Residual values of a boosted functional regression model 
#' 
#' Takes a fitted \code{FDboost}-object and computes the residuals.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param ... not used
#' 
#' @details The residual is missing if the corresponding value of the response was missing.
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return matrix of residual values
#' @method residuals FDboost
#' @export
### residuals (the current negative gradient)
residuals.FDboost <- function(object, ...){
  
  if(!any(class(object)=="FDboostLong")){
    resid <- matrix(object$resid(), nrow=object$ydim[1])
    resid[is.na(object$response)] <- NA 
  }else{
    resid <- object$resid()
    resid[is.na(object$response)] <- NA 
  }
  resid
}


#' Coefficients of boosted functional regression model 
#' 
#' Takes a fitted \code{FDboost}-object produced by \code{\link{FDboost}()} and 
#' returns estimated coefficient functions/surfaces \eqn{\beta(t), \beta(s,t)} and 
#' estimated smooth effects \eqn{f(z), f(x,z)} or \eqn{f(x, z, t)}. 
#' Not implemented for smooths in more than 3 dimensions.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param raw logical defaults to \code{FALSE}.
#' If \code{raw = FALSE} for each effect the estimated function/surface is calculated. 
#' If \code{raw = TRUE} the coefficients of the model are returned. 
#' @param which a subset of base-learners for which the coefficients
#' should be computed (numeric vector), 
#' defaults to NULL which is the same as \code{which=1:length(object$baselearner)}.
#' In the special case of \code{which=0}, only the coefficients of the offset are returned.
#' @param computeCoef defaults to \code{TRUE}, if \code{FALSE} only the names of the terms are returned
#' @param returnData return the dataset which is used to get the coefficient estimates as 
#' predictions, see Details. 
#' @param n1 see below
#' @param n2 see below
#' @param n3 n1, n2, n3 give the number of grid-points for 1-/2-/3-dimensional 
#' smooth terms used in the marginal equidistant grids over the range of the 
#' covariates at which the estimated effects are evaluated.
#' @param n4 gives the number of points for the third dimension in a 3-dimensional smooth term
#' @param ... other arguments, not used.
#' 
#' @return If \code{raw = FALSE}, a list containing 
#' \itemize{
#'  \item \code{pterms} a matrix containing the parametric / non-functional coefficients. 
#'  \item \code{smterms} a named list with one entry for each smooth term in the model. 
#'  Each entry contains
#'     \itemize{
#'          \item \code{x, y, z} the unique grid-points used to evaluate the smooth/coefficient function/coefficient surface
#'          \item \code{xlim, ylim, zlim} the extent of the x/y/z-axes
#'          \item \code{xlab, ylab, zlab} the names of the covariates for the x/y/z-axes
#'          \item \code{value} a vector/matrix/list of matrices containing the coefficient values 
#'          \item \code{dim} the dimensionality of the effect
#'          \item \code{main} the label of the smooth term (a short label)
#' }} 
#' If \code{raw = TRUE}, a list containing the estimated spline coefficients. 
#'  
#' @method coef FDboost
#' 
#' @details If \code{raw = FALSE} the function \code{coef.FDboost} generates adequate dummy data 
#' and uses the function \code{predict.FDboost} to 
#' compute the estimated coefficient functions. 
#' 
#' @export
### similar to coef.pffr() by Fabian Scheipl in package refund 
coef.FDboost <- function(object, raw = FALSE, which = NULL, 
                         computeCoef = TRUE, returnData = FALSE, 
                         n1 = 40, n2 = 40, n3 = 20, n4 = 10, ...){
  
  if(raw){
    return(object$coef(which = which))  
  } else {
    
    # delete an extra 0 in which as the offset is always returned
    if( length(which) > 1 && 0 %in% which) which <- which[which!=0]
    
    # List to be returned
    ret <- list()
    
    ## <FIXME> all linear terms should be part of pterms?
    # ret$pterms <- NULL 
    
    ## offset as first element
    ret$offset$x <- seq( min(object$yind), max(object$yind), l=n1)
    ret$offset$xlab <- attr(object$yind, "nameyind")
    ret$offset$xlim <- range(object$yind)
    ret$offset$value <- object$predictOffset(ret$offset$x)
    if(length(ret$offset$value)==1) ret$offset$value <- rep(ret$offset$value, n1)
    ret$offset$dim <- 1
    ret$offset$main <- "offset"
    
    # For the special case of which=0, only return the coefficients of the offset
    if(!is.null(which) & length(which)==1 && which==0){
      if(computeCoef){
        return(ret)
      }else{
        return("offset")
      }  
    }
    
    if(is.null(which)) which <- 1:length(object$baselearner)
    
    ## special case of ~1 intercept specification with scalar response
    if( object$ydim[[2]] == 1 && 
        any(which == 1) && 
        length(object$coef(which = 1)[[1]]) == 1 ){
      ret$intercept <- object$coef(which = 1)[[1]]
      which <- which[which != 1]
      if(length(which) == 0) return(ret)
    }
    
    getCoefs <- function(i){
      ## this constructs a grid over the range of the covariates
      ## and returns estimated values on this grid, with 
      ## by-variables set to 1
      ## cf. mgcv:::plots.R (plot.mgcv.smooth etc..) for original code
      
      safeRange <- function(x){
        if(is.factor(x)) return(c(NA, NA))
        return(range(x, na.rm=TRUE))
      }
      
      makeDataGrid <- function(trm){
        
        # variable for number of levels in bl2 for an effect bl1 %X% bl2
        numberLevels <- 1  

        ### generate data in the case of an bhistx()-bl
        if(grepl("bhistx", trm$get_call())){
          ng <- n2
          # get hmatrix-object
          temp <- trm$model.frame()[[1]]
          svals <- getArgvals(temp)
          svals <- seq(min(svals), max(svals),length=ng)
          tvals <- getTime(temp)
          tvals <- seq(min(tvals), max(tvals),length=ng)
          tvals <- rep(tvals, each=ng)

          if( grepl("%X", trm$get_call()) ){
            if(length(trm$get_names()) == 2){ # one %X%
              myargsHist <- environment(environment(trm$dpp)$Xfun)$args1
            }else{  # two ore more %X% (currently only works for two)
              myargsHist <- environment(environment(environment(
                environment(trm$dpp)$Xfun)$bl1$dpp)$Xfun)$args1
              if(is.null(myargsHist)) myargsHist <- environment(environment(trm$dpp)$Xfun)$args1
            }
          }else{
            myargsHist <- environment(trm$dpp)$args
          }
          
          ## use the same function for the numerical integration weights as in the model call 
          intFun <- myargsHist$intFun
          
          ## should only occur for more than two %X%
          if(is.null(intFun)){ 
            intFun <- integrationWeightsLeft
            warning("As integration function 'integrationWeightsLeft()' is used,", 
                    " which is the default in bhistx().")
          }
          
          ## generate a dummy functional variable I / integration weights
          dummyX <- I(diag(ng)) / intFun(diag(ng), svals)
          
          temph <- hmatrix(time=tvals, id=rep(1:ng, ng), x=dummyX, argvals=svals)
          
          d <- data.frame(z=I(temph))
          names(d) <- trm$get_names()[1]
          d[[ getTimeLab(temp) ]] <- tvals
          
          attr(d, "varnms") <- c(getArgvalsLab(temp), getTimeLab(temp))
          attr(d, "xm") <- seq(min(svals), max(svals),length=ng)
          attr(d, "ym") <- seq(min(tvals), max(tvals),length=ng)
          
          ## for a tensor product term: add the scalar factors to d
          if( grepl("%X", trm$get_call()) ){
            z <- trm$model.frame()[[trm$get_names()[2]]]
            if(is.factor(z)) {
              numberLevels <- length(unique(unique(z)))  
              zg <- sort(unique(z))[1] #sort(unique(z)) # use first possibility
            }else{
              zg <- 1 # use z=1 as neutral possibility
            } 
            
            d[[ trm$get_names()[2] ]] <- zg
            attr(d, "varnms") <- c(getArgvalsLab(temp), getTimeLab(temp), trm$get_names()[2])
            attr(d, "zm") <- zg
            
            ## add second factor variable to the dataset if necessary, because of two %X%
            z1 <- NULL
            if(length(trm$get_names()) > 2){
              z1 <- trm$model.frame()[[trm$get_names()[3]]]
              if(is.factor(z1)) {
                ## multiply the number of levels of both factors 
                numberLevels <- numberLevels * length(unique(unique(z1)))  
                z1g <- sort(unique(z1))[1] 
              }else{
                z1g <- 1 # use z1=1 as neutral possibility
              } 
              
              d[[ trm$get_names()[3] ]] <- z1g
              attr(d, "varnms") <- c(getArgvalsLab(temp), getTimeLab(temp), 
                                     trm$get_names()[2], trm$get_names()[3])
              attr(d, "z1m") <- z1g
            }
            
            ## make a list of data-frames with all combinations
            if(numberLevels > 1){
              
              dlist <- vector("list", numberLevels)
              zlevels <- sort(unique(z))  ##  sort(unique(z1))
              
              ## loop over all factor combinations 
              if(is.null(z1)){  # one %X%'
                dlist[[1]] <- d
                for(j in 1:numberLevels){
                  d[[ trm$get_names()[2] ]] <- zlevels[j] # use j-th factor level
                  attr(d, "add_main") <- paste0(trm$get_names()[2], "=", zlevels[j])
                  dlist[[j]] <- d
                }
              }else{ # two %X%
                z1levels <- sort(unique(z1))
                temp_d <- 1
                for(j in 1:length(zlevels)){ # loop over z 
                  d[[ trm$get_names()[2] ]] <- zlevels[j] # use j-th factor level of z
                  for(k in 1:length(z1levels)){ # loop over z1
                    d[[ trm$get_names()[3] ]] <- z1levels[k] # use k-th factor level of z1
                    attr(d, "add_main") <- paste0(trm$get_names()[2], "=", zlevels[j], ", ", 
                                                 trm$get_names()[3], "=", z1levels[k])
                    dlist[[temp_d]] <- d 
                    temp_d <- temp_d + 1
                  }
                  
                }
              } # end else !is.null(z1)
              d <- dlist  ## give back the list of data.frames
              attr(d, "numberLevels") <- numberLevels
            } ## end if(numberLevels > 1)
            
            #attr(d, "limits") <- limits
            #attr(d, "stand") <- stand
            
          }

          ## test <- predict(object, newdata=d, which=2)
          return(d)
        }
        
        ### <TODO> delete attr(d, "xm") <- xg and change code accordingly
        varnms <- trm$get_names()
        yListPlace <- NULL
        zListPlace <- NULL

        # generate grid of values in range of original data
        if(trm$dim==1){
          ng <- n1
          varnms <- varnms[!varnms %in% c("ONEx", "ONEtime")] 
          # Extra setup of dataframe in the case of a functional covariate
          if(grepl("bsignal", trm$get_call()) | grepl("bfpc", trm$get_call()) | grepl("bconcurrent", trm$get_call())){
            x <- attr(trm$model.frame()[[1]], "signalIndex")
            xg <- seq(min(x), max(x),length=ng) 
            varnms[1] <- attr(trm$model.frame()[[1]], "indname")            
          }else{
            # if(length(varnms) == 0)
            x <- trm$model.frame()[[varnms]]
            xg <- if(is.factor(x)) {
              sort(unique(x))
            } else seq(min(x), max(x), length=ng)
          }
          d <- list(xg)  # data.fame
          names(d) <- varnms
          attr(d, "xm") <- xg
          attr(d, "xname") <- varnms
          # For effect constant over index of response: add dummy-index so that length in clear
          if(attr(object$yind, "nameyind") != varnms){
            d[[attr(object$yind, "nameyind")]] <- seq(min(object$yind), max(object$yind),length=ng) 
          }           
        }        
        
        if(trm$dim > 1){          
          ng <- ifelse(trm$dim==2, n2, n3)          
          
          ### get variables for x, y and eventually z direction

          # Extra setup of dataframe in the case of a functional covariate
          if(grepl("bsignal", trm$get_call()) | grepl("bfpc", trm$get_call()) | 
             grepl("bhist", trm$get_call())){
            x <- attr(trm$model.frame()[[1]], "signalIndex")
            xg <- seq(min(x), max(x),length=ng) 
            varnms[1] <- attr(trm$model.frame()[[1]], "indname")   
          }else{ # scalar covariate
            x <- trm$model.frame()[[1]]
            xg <- if(is.factor(x)) {
              sort(unique(x))
            } else seq(min(x), max(x),length=ng)            
          }
          yListPlace <- 2
          if(grepl("by", trm$get_call()) | 
             ( !any(class(object)=="FDboostLong") &&  grepl("%X", trm$get_call())) ){
            yListPlace <- 3
          }
          if(!grepl("bhist", trm$get_call())){
            y <- trm$model.frame()[[yListPlace]] # index of response or second scalar covariate
          }else{
            y <- object$yind  ## attr(trm$model.frame()[[1]], "indexY")
            varnms[2] <- attr(object$yind, "nameyind") ## attr(trm$model.frame()[[1]], "indnameY")
            if(varnms[1]==varnms[2]) varnms[1] <- paste(varnms[2], "_cov", sep="")
            if(attr(object$yind, "nameyind") != attr(trm$model.frame()[[1]], "indnameY")){
              stop("coef.FDboost works for bhist only if time variable is the same in timeformula and bhist.")
            }
          }
          yg <- if(is.factor(y)) {
            sort(unique(y))
          } else seq(min(y), max(y),length=ng)
          if(length(varnms)==2){
            d <- list(xg, yg)  # data.frame
            attr(d, "xm") <- xg
            attr(d, "ym") <- yg    
          } else {
            zListPlace <- ifelse(yListPlace == 2, 3, 2)
            z <- trm$model.frame()[[zListPlace]]
            zg <- if(is.factor(z)) {
              sort(unique(z))
            }else{
              if(grepl("by", trm$get_call())){ 1 }else{ seq(min(z), max(z), length=n4) }
            } 
            d <- list(xg, yg, zg)  # data.frame
            ## special case of factor by-variable 
            #if(grepl("by", trm$get_call()) && grepl("bols", trm$get_call()) && is.factor(z)){
            #  d <- list(rep(xg, length(yg)), rep(yg, each=length(xg)), zg)
            #}
            # special case of factor by-variable 
            if(grepl("by", trm$get_call()) && grepl("bols", trm$get_call()) && is.factor(z)){
              d <- list(rep(xg, length(zg)), yg, rep(zg, each=length(zg)))
            }
            
            attr(d, "xm") <- d[[1]]
            attr(d, "ym") <- d[[2]]
            attr(d, "zm") <- d[[3]]
          }
          names(d) <- varnms[c(1, yListPlace, zListPlace)] # colnames(d) <- varnms
          attr(d, "varnms") <- varnms[c(1, yListPlace, zListPlace)] 
        }
        
        
        ## add dummy signal to data for bsignal()
        if(grepl("bsignal", trm$get_call()) | grepl("bfpc", trm$get_call()) ){
          d[[ trm$get_names()[1] ]] <- I(diag(ng)/integrationWeights(diag(ng), d[[varnms[1]]] ))
        }
        
        ## <FIXME> is this above dummy-matrix correct for bfpc?

        ## add dummy signal to data for bhist()
        ## standardisation weights depending on t must be multiplied to the final \beta(s,t)
        ## as they cannot be included into the variable x(s)
        ## use intFun() to compute the integration weights
        # ls(environment(trm$dpp))
        if(grepl("bhist", trm$get_call()) ){
          ## temp <- I(diag(ng)/integrationWeightsLeft(diag(ng), d[[varnms[1]]]))
          ## use intFun() of the bl to compute the integration weights
          temp <- environment(trm$dpp)$args$intFun(diag(ng), d[[varnms[1]]])
          ##if(environment(trm$dpp)$args$stand=="transform"){
          ##  xindStand <- (d[[varnms[1]]] - min(d[[varnms[1]]])) / (max(d[[varnms[1]]]) - min(d[[varnms[1]]]))
          ##  temp <- environment(trm$dpp)$args$intFun(diag(ng), xindStand)
          ##}
          d[[attr(trm$model.frame()[[1]], "xname")]] <- I(diag(ng)/temp)
          limits <- environment(trm$dpp)$args$limits
          stand <- environment(trm$dpp)$args$stand
          attr(d, "limits") <- limits
          attr(d, "stand") <- stand
        }
        
        ## add dummy signal to data for bconcurrent()
        if(grepl("bconcurrent", trm$get_call())){
          d[[ trm$get_names()[1] ]] <- I(matrix(rep(1.0, ng^2), ncol=ng))
        }
        
        if(trm$get_vary() != ""){
          d$by <- 1
          colnames(d) <- c(head(colnames(d),-1), trm$get_vary())
        } 
        
        # set time-variable to 1, if response is a scalar
        if(length(object$yind)==1){
          if( all(attr(d, "zm") == d[[attr(object$yind, "nameyind")]]) ){
            attr(d, "zm") <- 1
          }
          d[[attr(object$yind, "nameyind")]] <- 1 
        }
        
        
        # if %X% was used in combination with factor variables make a list of data-frames
        if(!any(class(object) == "FDboostLong") && grepl("%X", trm$get_call())){
          dlist <- NULL

          ## if %X% was used in combination with factor variables make a list of data-frames
          if(is.factor(x) & is.factor(z)){ ## both variables are factors 
            numberLevels <-  nlevels(x) * nlevels(z)
            xlevels <- sort(unique(x))
            zlevels <- sort(unique(z))
            
            dlist <- vector("list", numberLevels)
            temp_d <- 1
            for(j in 1:length(xlevels)){ # loop over x 
              d[[1]] <- xlevels[j] # use j-th factor level of x
              for(k in 1:length(zlevels)){ # loop over z
                d[[3]] <- zlevels[k] # use k-th factor level of z
                attr(d, "xm") <- d[[1]]
                attr(d, "zm") <- d[[3]]
                attr(d, "add_main") <- paste0(names(d)[1], "=", xlevels[j], ", ",  
                                              names(d)[3], "=", zlevels[j])
                dlist[[temp_d]] <- d 
                temp_d <- temp_d + 1
              }
            }
          }else{
            if(is.factor(z)){ ## only first variable is a factor
              numberLevels <- nlevels(z)
              zlevels <- sort(unique(z))
              dlist <- vector("list", numberLevels)
              for(j in 1:length(zlevels)){ # loop over z 
                d[[3]] <- rep(zlevels[j], length(d[[1]])) # use j-th factor level of z
                attr(d, "xm") <- d[[3]]
                attr(d, "add_main") <- paste0(names(d)[3], "=", zlevels[j])
                dlist[[j]] <- d 
              }
            }else{
              if(is.factor(x)){ ## only first variable is a factor
                numberLevels <- nlevels(x)
                xlevels <- sort(unique(x))
                dlist <- vector("list", numberLevels)
                for(j in 1:length(xlevels)){ # loop over x 
                  d[[1]] <- rep(xlevels[j], length(d[[3]])) # use j-th factor level of x
                  attr(d, "xm") <- d[[1]]
                  attr(d, "add_main") <- paste0(names(d)[1], "=", xlevels[j])
                  dlist[[j]] <- d 
                }
              }else{ ## both variables are metric
                # <FIXME> allow something more flexible, depending on bolsc / bbsc?
                ## trm$get_call()
                numberLevels <- n4
                zlevels <- seq(min(z), max(z), l = n4)
                dlist <- vector("list", numberLevels)
                for(j in 1:length(zlevels)){ # loop over x 
                  d[[3]] <- rep(zlevels[j], length(d[[1]])) # use j-th quantile of x
                  attr(d, "zm") <- d[[1]]
                  attr(d, "add_main") <- paste0(names(d)[3], "=", round(zlevels[j], 2))
                  dlist[[j]] <- d 
                }
              } 
            }
            
          }
          if(!is.null(dlist)) d <- dlist
          attr(d, "numberLevels") <- numberLevels
        }  ## end if(grepl("%X", trm$get_call()))
        
        return(d)
      } ## end of function makeDataGrid()
      
      getP <- function(trm, d){
        #return an object similar to what plot.mgcv.smooth etc. returns 
        if(trm$dim==1){
          predHelp <- predict(object, which=i, newdata=d)
          if(!is.matrix(predHelp)){ 
            X <- predHelp
          }else{
            X <- if(any(trm$get_names() %in% c("ONEtime")) | 
                    any(class(object)=="FDboostScalar")){ # effect constant in t 
              predHelp[,1]
            }else{ 
              predHelp[1,] # smooth intercept/ concurrent effect                
            }            
          } 
          P <- list(x=attr(d, "xm"), xlab=attr(d, "xname"), xlim=safeRange(attr(d, "xm")))
          ## trm$dim > 1
        }else{
          varnms <- attr(d, "varnms")
          if(trm$dim==2){
            X <- predict(object, newdata=d, which=i)
            attr(X, "offset") <- NULL
            vecStand <- NULL
           
            ## for bhist(), multiply with standardisation weights if necessary
            ## you need the args$vecStand from the prediction of X, constructed here
            if(grepl("bhist", trm$get_call())){

              ## get the args of bhist/bhistx, environment depends on use of %X%
              if( grepl("%X", trm$get_call()) ){
                if(length(trm$get_names()) == 2){ # one %X%
                  myargsHist <- environment(environment(trm$dpp)$Xfun)$args1
                }else{  # two ore more %X% (currently only works for two)
                  myargsHist <- environment(environment(environment(
                    environment(trm$dpp)$Xfun)$bl1$dpp)$Xfun)$args1
                  if(is.null(myargsHist)) myargsHist <- environment(environment(trm$dpp)$Xfun)$args1
                }
              }else{ # just bhist() or bhistx() without %X%
                myargsHist <- environment(trm$dpp)$args
              }
              
              ## this should only occur for more than two %X%
              if(is.null(myargsHist$stand)){
                warning("No standardization is used, i.e. stand = 'no',", 
                        " which is the default in bhistx().", 
                        "As intFun() integrationWeightsLeft() is used. ", 
                        "No limits are used.")
                myargsHist$stand <- "no"
                myargsHist$intFun <- integrationWeightsLeft
                myargsHist$limits <- function(s, t) TRUE
              }
              
              if(myargsHist$stand %in% c("length","time")){
                Lnew <- myargsHist$intFun(diag(length(attr(d, "ym"))), attr(d, "ym") )
                ## Standardize with exact length of integration interval
                ##  (1/t-t0) \int_{t0}^t f(s) ds
                if(myargsHist$stand == "length"){
                  ind0 <- !t(outer( attr(d, "ym"), attr(d, "xm"), myargsHist$limits) )
                  Lnew[ind0] <- 0
                  ## integration weights in s-direction always sum exactly to 1, 
                  vecStand <- rowSums(Lnew)
                }           
                ## use time of current observation for standardization
                ##  (1/t) \int_{t0}^t f(s) ds
                if(myargsHist$stand == "time"){
                  yindHelp <- attr(d, "ym")
                  yindHelp[yindHelp == 0] <- Lnew[1,1] 
                  vecStand <- yindHelp
                }
                X <- t(t(X)*vecStand)
              }
            } ## end stand for bhist / bhistx
            
            P <- list(x=attr(d, "xm"), y=attr(d, "ym"), xlab=varnms[1], ylab=varnms[2],
                      ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")),
                      z=attr(d, "zm"), zlab=varnms[3], vecStand=vecStand)
            
            ## include the second scalar covariate called z1 into the output
            if( grepl("bhistx", trm$get_call()) & length(trm$get_names()) > 2){
              extra_output <- list(z1=attr(d, "z1m"), z1lab=varnms[4])
              P <- c(P, extra_output)
            }
            
            ## save the arguments of stand and limits as part of returned object
            if(grepl("bhist", trm$get_call())){
              P$stand <- myargsHist$stand
              P$limits <- myargsHist$limits
            }
            
          }else{
            if(trm$dim==3){
              values3 <- seq(min(d[[varnms[3]]]), max(d[[varnms[3]]]), l=n4)
              xygrid <- expand.grid(d[[varnms[1]]], d[[varnms[2]]] )
              X <- lapply(values3, function(x){
                d1 <- list(xygrid[,1], xygrid[,2], x)
                names(d1) <- varnms
                matrix(predict(object, newdata=d1, which=i), ncol=length(d[[varnms[1]]]))
              }) 
              P <- list(x=attr(d, "xm"), y=attr(d, "ym"), z=values3, 
                   xlab=varnms[1], ylab=varnms[2], zlab=varnms[3],
                   ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")), zlim=safeRange(attr(d, "zm")))
            }
          }
        }
        if(!is.null(object$ydim)){
          P$value <- X
        }else{
          #### <FIXME> is dimension, byrow correct??
          P$value <- matrix(X, nrow=length(attr(d, "xm")) )
        }
        #P$coef <- cbind(d, "value"=P$value)    
        P$dim <- trm$dim
        P$main <- shrtlbls[i]
        P$add_main <- attr(d, "add_main")
        return(P)
      } ## end of function getP()
      
      trm <- object$baselearner[[i]]
      trm$dim <- length(trm$get_names())
      if(any(grepl("ONEx", trm$get_names()), 
             grepl("ONEtime", trm$get_names()))) trm$dim <- trm$dim - 1
      
      ### give error for bl1 %X% bl2 %X% bl3
      #if( grepl("bhistx", trm$get_call()) & trm$dim > 2){
      #  stop("coef.FDboost() does not work for tensor products %X% with more than two base-learners.")
      #}
          
      ## add 1 to dimension of bhist and bhistx, otherwise dim is only 1  
      if( grepl("bhist", trm$get_call()) ){
        trm$dim <- trm$dim + 1
      }
      
      # If a by-variable was specified, reduce number of dimensions
      # as smooth linear effect in several groups can be plotted in one plot 
      if( grepl("by =", trm$get_call()) && grepl("bols", trm$get_call()) || 
            grepl("by =", trm$get_call()) && grepl("bbs", trm$get_call()) ) trm$dim <- trm$dim - 1
      
      # <FIXME> what to do with bbs(..., by=factor)?

      if(trm$dim > 3 & !grepl("bhistx", trm$get_call()) ){
        warning("Can't deal with smooths with more than 3 dimensions, returning NULL for ", 
                shrtlbls[i], ".")
        return(NULL)
      }
      
      d <- makeDataGrid(trm)
      
      ### <FIXME> better solution for %X% in base-learner!!!
      if(!is.null(object$ydim) && any(grepl("%X", trm$get_call())) 
         && !any(grepl("bhistx", trm$get_call())) ) trm$dim <- trm$dim - 1
      
      ## it is necessary to expand the dataframe!
      if(!grepl("bhistx(", trm$get_call(), fixed=TRUE) && 
         is.null(object$ydim) && !grepl("bconcurrent", trm$get_call())){
        #print(attr(d, "varnms"))
        vari <- names(d)[1]
        if(is.factor(d[[vari]])){
          d[[vari]] <- d[[vari]][ rep(1:NROW(d[[vari]]), times=length(d[[attr(object$yind ,"nameyind")]]) ) ]
          if(trm$dim>1) d[[attr(object$yind ,"nameyind")]] <- rep(d[[attr(object$yind ,"nameyind")]], 
                                                                  each=length(unique(d[[vari]]))  )
          }else{
          # expand signal variable
          if( grepl("bhist(", trm$get_call(), fixed = TRUE) | 
              grepl("bsignal", trm$get_call()) | grepl("bfpc", trm$get_call()) ){
            vari <- names(d)[!names(d) %in% attr(d, "varnms")]
            d[[vari]] <- d[[vari]][ rep(1:NROW(d[[vari]]), times=NROW(d[[vari]])), ]
            
          }else{ # expand scalar variable
            vari <- names(d)[1]
            if(vari!=attr(object$yind ,"nameyind")) d[[vari]] <- d[[vari]][ rep(1:NROW(d[[vari]]), times=NROW(d[[vari]])) ]
          } 
          # expand yind 
          if(trm$dim>1) d[[attr(object$yind ,"nameyind")]] <- rep(d[[attr(object$yind ,"nameyind")]], 
                                                                  each=length(d[[attr(object$yind ,"nameyind")]])) 
        }
      }
      
      ###### just return the data, that is used for the prediction
      if(returnData){
        if(grepl("bhist", trm$get_call())){
          message("If argument stand is specified !=\"no\", the standardization will be part of the predicted coefficient.")
        } 
        return(d)
      }

      if( !is.null(attr(d, "numberLevels")) && attr(d, "numberLevels") > 1){
        if( grepl("bhistx", trm$get_call()) ) trm$dim <- 2
        ## get smooth coefficient estimates for several factor levels
        P <- lapply(d, getP, trm=trm) 
        P$numberLevels <- attr(d, "numberLevels") 
      }else{
        ## get smooth coefficient estimates
        P <- getP(trm, d)
      }
        
      # get proper labeling
      P$main <- shrtlbls[i]
      
      return(P)
    } # end of function getCoefs()
      

    ## Function to obtain nice short names for labeling of plots
    shortnames <- function(x){
      if(substr(x,1,1)=="\"") x <- substr(x, 2, nchar(x)-1)

      ## split at expressions %.% with 1 or 2 characters between % and %
      xpart <- unlist(strsplit(x, split = "%.{1,2}%"))
      ## find the expressions at which the split is done 
      operator <- gregexpr(pattern = "%.{1,2}%", text = x)[[1]]
      operator <- sapply(1:length(operator), 
                      function(i) substr(x, operator[i], operator[i] + attr(operator, "match.length")[i] -1 ) )
      
      for(i in 1:length(xpart)){
        xpart[i] <- gsub(pattern = "\\\"", replacement = "", x = xpart[i], fixed=TRUE)
        xpart[i] <- gsub(pattern = "\\", replacement = "", x = xpart[i], fixed=TRUE)
        nvar <- length(all.vars(formula(paste("Y~", xpart[i])))[-1])
        commaSep <- unlist(strsplit(xpart[i], ","))  
        
        # shorten the name to first variable and delete x= if present
        if(grepl("=", commaSep[1])){
          temp <- unlist(strsplit(commaSep[1], "="))
          temp[1] <- unlist(strsplit(temp[1], "(", fixed=TRUE))[1]
          if(substr(temp[2], 1, 1)==" ") temp[2] <- substr(temp[2], 2, nchar(temp[2]))
          if(length(commaSep) == 1){ 
            xpart[i] <- paste(temp[1], "(", temp[2], sep="")
          }else{
            xpart[i] <- paste(temp[1], "(", temp[2], ")", sep="")
          }
        }else{
          if(length(commaSep) > 1){ xpart[i] <- paste(commaSep[1], ")", sep="")}
        }
        #xpart[i] <- if(length(commaSep)==1){
        #  paste(paste(commaSep[1:nvar], collapse=","), sep="")
        #}else paste(paste(commaSep[1:nvar], collapse=","), ")", sep="") 
        #if(substr(xpart[i], nchar(xpart[i])-2, nchar(xpart[i])-1)==") "){
        #  xpart[i] <- substr(xpart[i], 1, nchar(xpart[i])-2)
        #}
      }
      ret <- xpart
      if(length(xpart)>1) ret <- paste(xpart, collapse=paste("", operator))
      if(length(xpart)>1){
        ## name is connatecated with the scheme bl1 %.% bl2 %.% bl3
        temp <- rep(NA, length(xpart) + length(operator))
        temp[seq(1, by=2, length.out = length(xpart))] <- xpart
        temp[seq(2, by=2, length.out = length(operator))] <- paste("", operator)
        ret <- paste(temp, collapse="")
      } 
      ret
    }

    ## short names for the terms, if shortnames() does not work, use the original names
    shrtlbls <- try(unlist(lapply(names(object$baselearner), shortnames)))
    if(class(shrtlbls)=="try-error") shrtlbls <- names(object$baselearner) 
    
    ###### just return the data that is used for the prediction
    if(returnData){
      ret <- list()
      ret <- lapply(which, getCoefs)
      if(length(which)==1) ret <- ret[[1]] # unlist for only one effect
      return(ret)
    }
    
    if(computeCoef){
      ### smooth terms
      ret$smterms <- lapply(which, getCoefs)         
      names(ret$smterms) <- sapply(seq_along(ret$smterms), function(i){
        ret$smterms[[i]]$main
      })
      return(ret)
    }else{
      return(shrtlbls[which])
    }  
  }
}

# help function to color perspective plots - col1 positive values, col2 negative values
getColPersp <- function(z, col1 = "tomato", col2 = "lightblue"){
  nrz <- nrow(z)
  ncz <- ncol(z)
  
  # Compute the z-value at the facet centres
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  
  # use the colors col1 and col2 for negative and positive values
  colfacet <- matrix(nrow=nrow(zfacet), ncol=ncol(zfacet))
  colfacet[zfacet < 0] <- col2
  colfacet[zfacet > 0] <- col1
  colfacet[zfacet == 0] <- "white"
  
  return(colfacet) 
}


#' Plot the fit or the coefficients of a boosted functional regression model 
#' 
#' Takes a fitted \code{FDboost}-object produced by \code{\link{FDboost}()} and 
#' plots the fitted effects or the coefficient-functions/surfaces.
#' 
#' @param x a fitted \code{FDboost}-object
#' @param raw  logical defaults to \code{FALSE}.
#' If \code{raw = FALSE} for each effect the estimated function/surface is calculated. 
#' If \code{raw = TRUE} the coefficients of the model are returned. 
#' @param rug when \code{TRUE} (default) then the covariate to which the plot applies is 
#' displayed as a rug plot at the foot of each plot of a 1-d smooth, 
#' and the locations of the covariates are plotted as points on the contour plot 
#' representing a 2-d smooth.
#' @param which a subset of base-learners to take into account for plotting. 
#' @param includeOffset logical, defaults to \code{TRUE}. Should the offset be included in 
#' the plot of the intercept (default) or should it be plotted separately.
#' @param ask logical, defaults to \code{TRUE}, if several effects are plotted the user
#' has to hit Return to see next plot.
#' @param n1 see below
#' @param n2 see below
#' @param n3 n1, n2, n3 give the number of grid-points for 1-/2-/3-dimensional 
#' smooth terms used in the marginal equidistant grids over the range of the 
#' covariates at which the estimated effects are evaluated.
#' @param n4 gives the number of points for the third dimension in a 3-dimensional smooth term
#' @param onlySelected, logical, defaults to \code{TRUE}. Only plot effects that where 
#' selected in at least one boosting iteration.
#' @param pers logical, defaults to \code{FALSE}, 
#' If \code{TRUE}, perspective plots (\code{\link[graphics]{persp}}) for 
#' 2- and 3-dimensional effects are drawn.
#' If \code{FALSE}, image/contour-plots (\code{\link[graphics]{image}}, 
#' \code{\link[graphics]{contour}}) are drawn for 2- and 3-dimensional effects. 
#' @param commonRange logical, defaults to \code{FALSE}, 
#' if \code{TRUE} the range over all effects is the same (so far not implemented).
#' 
#' @param subset subset of the observed response curves and their predictions that is plotted. 
#' Per default all observations are plotted.
#' @param posLegend location of the legend, if a legend is drawn automatically 
#' (only used in plotPredicted). The default is "topleft".
#' @param lwdObs lwd of observed curves (only used in plotPredicted)
#' @param lwdPred lwd of predicted curves (only used in plotPredicted)
#' @param ... other arguments, passed to \code{funplot} (only used in plotPredicted)
#' 
#' @aliases plotPredicted plotResiduals
#' 
#' @seealso \code{\link{FDboost}} for the model fit and 
#' \code{\link{coef.FDboost}} for the calculation of the coefficient functions. 
#' 
#' @method plot FDboost
#' @export
### function to plot raw values or coefficient-functions/surfaces of a model 
plot.FDboost <- function(x, raw = FALSE, rug = TRUE, which = NULL, 
                         includeOffset = TRUE, ask = TRUE,
                         n1 = 40, n2 = 40, n3 = 20, n4 = 11,
                         onlySelected = TRUE, pers = FALSE, commonRange = FALSE, ...){
  
  ### Get further arguments passed to the different plot-functions
  dots <- list(...)
  
  getArguments <- function(x, dots=dots){
    if(any(names(dots) %in% names(x))){
      dots[names(dots) %in% names(x)]
    }else list()
  }
  
  #argsPlot <- getArguments(x=formals(graphics::plot.default), dots=dots)
  argsPlot <- getArguments(x=c(formals(graphics::plot.default), par()), dots=dots)
  argsMatplot  <- getArguments(x=c(formals(graphics::matplot), par(), 
                                   formals(graphics::plot.default)), dots=dots)
  argsFunplot  <- getArguments(x=c(formals(funplot), par(), 
                                   formals(graphics::plot.default)), dots=dots)

  argsImage <- getArguments(x=c(formals(graphics::plot.default), 
                                formals(graphics::image.default)), dots=dots)
  argsContour <- getArguments(x=formals(graphics::contour.default), dots=dots)
  argsPersp <- getArguments(x=formals(getS3method("persp", "default")), dots=dots)
  
  plotWithArgs <- function(plotFun, args, myargs){        
    args <- c(myargs[!names(myargs) %in% names(args)], args)        
    do.call(plotFun, args)            
  }
  
  #     #### function by Fabian Scheipl
  #     # colors in rgb
  #     alpha <- function(x, alpha=25){
  #       tmp <- sapply(x, col2rgb)
  #       tmp <- rbind(tmp, rep(alpha, length(x)))/255
  #       return(apply(tmp, 2, function(x) do.call(rgb, as.list(x))))
  #     }      
  #     clrs <- alpha( rainbow(x$ydim[1]), 125) 
  
  ### get the effects to be plotted
  whichSpecified <- which
  if(is.null(which)) which <- 1:length(x$baselearner) 
  
  if(onlySelected){
    #     sel <- selected(x)
    #     if( !1 %in% sel  && length(x$offsetVec) > 1 && grepl("ONEx", names(x$baselearner)[[1]])){
    #       sel <- c(1, sel) # plot the offset as first effect
    #     } 
    which <- intersect(which, c(0, selected(x)))
  }
  
  # In the case that intercept and offset should be plotted and the intercept was never selected
  # plot the offset
  if( (1 %in% whichSpecified | is.null(whichSpecified))  & !1 %in% which & length(x$yind)>1) which <- c(0, which)
  
  if(length(which)==0){
    warning("Nothing selected for plotting.")
    return(NULL)
  } 

  ### plot coefficients of model (smooth curves and surfaces)
  if(!raw){ 
        
    # compute the coefficients of the smooth terms that should be plotted
    coefMod <- coef(x, which=which, n1=n1, n2=n2, n3=n3, n4=n4)
    terms <- coefMod$smterms
    offsetTerms <- coefMod$offset
    bl_data <- lapply(x$baselearner[which], function(x) x[["get_data"]]()) 
    
    # plot nothing but the offset
    if(length(which)==1 && which==0){
      terms <- list(offsetTerms)
      bl_data <- c(offset = list( list(x$yind) ), bl_data)
      names(bl_data[[1]]) <- attr(x$yind, "nameyind")
    }
    
    # include the offset in the plot of the intercept
    if( includeOffset && 1 %in% which && grepl("ONEx", names(terms)[1]) ){
      terms[[1]]$value <- terms[[1]]$value + matrix(offsetTerms$value, ncol=1, nrow=n1)
      terms[[1]]$main <- paste("offset", "+", terms[[1]]$main)
    }
    
    # plot the offset as extra effect
    # case 1: the offset should be included as extra plot
    # case 2: the whole model is plotted, but the intercept-base-learner was never selected
    if( (!includeOffset | (includeOffset & !1 %in% which)) & 
         is.null(whichSpecified)){
      terms <- c(offset = list(offsetTerms), terms)
      bl_data <- c(offset = list( list(x$yind) ), bl_data)
      names(bl_data[[1]]) <- attr(x$yind, "nameyind")     
    } 
   
    if((length(terms)>1 || is.null(terms[[1]]$dim) || terms[[1]]$dim==3) & ask) par(ask=TRUE)
    
    #     ### <TODO> implement common range
    #     if(commonRange){
    #       #range <- range(fit)
    #       #range[1] <- range[1]-0.02*diff(range)
    #     }else range <- NULL
    
    for(i in 1:length(terms)){
      
      trm <- terms[[i]] 
      
      myplot <- function(trm){
        
        if(grepl("bhist", trm$main)){
          ## set 0 to NA so that beta only has values in its domain
          # get the limits- function
          #limits <- get("args", (environment(x$baselearner[[which[i]]]$dpp)))$limits
          limits <- trm$limits
          if(is.null(limits)){
            warning("limits is NULL, the default limits 's<=t' are used for plotting.")
            limits <- function(s, t) {
              (s < t) | (s == t)
            }
          }
          trm$value[!outer(trm$x, trm$y, limits)] <- NA
        }
        
        # plot for 1-dim effects
        if(trm$dim==1){
          if(length(trm$value)==1) trm$value <- rep(trm$value, l=length(trm$x)) 
          
          if(!"add" %in% names(dots)){
            plotWithArgs(plot, args=argsPlot, 
                         myargs=list(x=trm$x, y=trm$value, xlab=trm$xlab, main=trm$main, 
                                     ylab="coef", type="l"))
          }else{
            plotWithArgs(lines, args=argsPlot, 
                         myargs=list(x=trm$x, y=trm$value, xlab=trm$xlab, main=trm$main, 
                                     ylab="coef", type="l"))          
          }
          
          if(rug & !is.factor(x=trm$x)){
            if(grepl("bconcurrent", trm$main) | grepl("bsignal", trm$main) | grepl("bfpc", trm$main) ){
              rug(attr(bl_data[[i]][[1]], "signalIndex"), ticksize = 0.02)
            }else ifelse(length(unique(bl_data[[i]][[1]]))!=1,
                         rug(bl_data[[i]][[1]], ticksize = 0.02),
                         rug(bl_data[[i]][[2]], ticksize = 0.02))
          } 
        } 
        
        # plot with factor variable
        if( (!grepl("bhistx", trm$main)) && trm$dim==2 &&
           ((is.factor(trm$x) | is.factor(trm$y)) | is.factor(trm$z)) ){
          if(is.factor(trm$y)){ # effect with by-variable (by-variable is factor)
            plotWithArgs(matplot, args=argsMatplot, 
                         myargs=list( x=trm$z, y=t(trm$value), xlab=trm$ylab, main=trm$main, 
                                      ylab="coef", type="l", col=as.numeric(trm$y) ) )
            if(rug){
              #rug(bl_data[[i]][[3]], ticksize = 0.02) 
              rug(x$yind, ticksize = 0.02)
            }
          }else{ # effect of factor variable
            plotWithArgs(matplot, args=argsMatplot, 
                         myargs=list(x=trm$y, y=t(trm$value), xlab=trm$ylab, main=trm$main, 
                                     ylab="coef", type="l"))
            if(rug){
              #rug(bl_data[[i]][[2]], ticksize = 0.02) 
              rug(x$yind, ticksize = 0.02)
            }
          }
          
        }else{
          # persp-plot for 2-dim effects
          if(trm$dim==2 & pers){
            if(length(unique(as.vector(trm$value)))==1){
              # persp() gives error if only a flat plane should be drawn
              plot(y=trm$value[1,], x=trm$x, main=trm$main, type="l", xlab=trm$ylab, 
                   ylab="coef")
            }else{  
              range <- range(trm$value, na.rm = TRUE)
              if(range[1]==range[2]) range <- range(0, range)
              zlim <- c(range[1] - 0.05*(range[2] - range[1]), 
                        range[2] + 0.05*(range[2] - range[1]))
              plotWithArgs(persp, args=argsPersp,
                           myargs=list(x=trm$x, y=trm$y, z=trm$value, xlab=paste("\n", trm$xlab), 
                                       ylab=paste("\n", trm$ylab), zlab=paste("\n", "coef"), 
                                       main=trm$main, theta=30, zlim=zlim,
                                       phi=30, ticktype="detailed",
                                       col=getColPersp(trm$value)))
            } 
          }
          # image for 2-dim effects
          if(trm$dim==2 & !pers){        
            plotWithArgs(image, args=argsImage,
                         myargs=list(x=trm$y, y=trm$x, z=t(trm$value), xlab=trm$ylab, ylab=trm$xlab, 
                                     main=trm$main, col = heat.colors(length(trm$x)^2)))          
            plotWithArgs(contour, args=argsContour,
                         myargs=list(trm$y, trm$x, z=t(trm$value), add = TRUE))
            
            if(rug){
              ##points(expand.grid(bl_data[[i]][[1]], bl_data[[i]][[2]]))
              if(grepl("bhist", trm$main)){
                rug(x$yind, ticksize = 0.02)
              }else{
                ifelse(grepl("by", trm$main) | ( !any(class(x)=="FDboostLong") && grepl("%X", trm$main) ) ,
                       rug(bl_data[[i]][[3]], ticksize = 0.02),
                       rug(bl_data[[i]][[2]], ticksize = 0.02))
              }
              ifelse(grepl("bsignal", trm$main) | grepl("bfpc", trm$main) | grepl("bhist", trm$main),
                     rug(attr(bl_data[[i]][[1]], "signalIndex"), ticksize = 0.02, side=2),
                     rug(bl_data[[i]][[1]], ticksize = 0.02, side=2))
            }
          }        
        }
        ### 3 dim plots
        # persp-plot for 3-dim effects
        if(trm$dim==3 & pers){
          for(j in 1:length(trm$z)){
            plotWithArgs(persp, args=argsPersp,
                         myargs=list(x=trm$x, y=trm$y, z=trm$value[[j]], xlab=paste("\n", trm$xlab), 
                                     ylab=paste("\n", trm$ylab), zlab=paste("\n", "coef"), 
                                     theta=30, phi=30, ticktype="detailed", 
                                     zlim=range(trm$value), col=getColPersp(trm$value[[j]]), 
                                     main= paste(trm$zlab ,"=", round(trm$z[j],2), ": ", trm$main, sep=""))
            )         
          }
        }
        # image for 3-dim effects
        if(trm$dim==3 & !pers){
          for(j in 1:length(trm$z)){
            plotWithArgs(image, args=argsImage,
                         myargs=list(x=trm$x, y=trm$y, z=trm$value[[j]], xlab=trm$xlab, ylab=trm$ylab,
                                     col = heat.colors(length(trm$x)^2), zlim=range(trm$value),
                                     main= paste(trm$zlab ,"=", round(trm$z[j],2), ": ", trm$main, sep="")))
            plotWithArgs(contour, args=argsContour,
                         myargs=list(trm$x, trm$y, trm$value[[j]], xlab=trm$xlab, add = TRUE))
            if(rug){
              points(bl_data[[i]][[1]], bl_data[[i]][[2]]) 
            }
            #           plotWithArgs(filled.contour, args=list(),
            #                        myargs=list(x=trm$x, y=trm$y, z=trm$value[[j]], xlab=trm$xlab, ylab=trm$ylab,
            #                                    zlim=range(trm$value), color.palette=heat.colors,
            #                                    main= paste(trm$zlab ,"=", round(trm$z[j],2), ": ", trm$main, sep="")))
          }        
        }
        
      }  ## end function myplot()
      
      ### call to function myplot()
      if(is.null(trm$numberLevels)){
        myplot(trm=trm)
      }else{  # several levels of bl1 %X% bl2
        lapply(trm[1:trm$numberLevels], myplot)
      }
      
    } # end for-loop
    
    if(length(terms)>1 & ask) par(ask=FALSE) 
    
  ### plot smooth effects as they are estimated for the original data
  }else{
        
    ################################          
    # predict the effects using the original data
    terms <- predict(x, which=which)
    offset <- attr(terms, "offset")
    
    # convert matrix into a list, each list entry for one effect
    if(is.null(x$ydim) & !is.null(dim(terms))){
      temp <- list()
      for(i in 1:ncol(terms)){
        temp[[i]] <- terms[,i]
      }
      names(temp) <- colnames(terms)
      terms <- temp
      rm(temp)
    }
    
    if(class(terms)!="list") terms <- list(terms)
    if(length(which)==1 && which==0) terms[[1]] <- offset
    if(length(which)==1 && length(terms[[1]]) == 1 && terms[[1]]==0){ terms[[1]] <- rep(0, l=length(x$yind)) }
    
    #if(length(which)==1 && !any(class(x)=="FDboostLong")) terms <- list(terms) 
    
    shrtlbls <- try(coef(x, which=which, computeCoef=FALSE))# get short names
    if(class(shrtlbls)=="try-error"){
      shrtlbls <- names(x$baselearner)[which[which!=0]]
      if(0 %in% which) shrtlbls <- c("offset", which)
    }
    if(is.null(shrtlbls)) shrtlbls <- "offset" 
    time <- x$yind
    
    # include the offset in the plot of the intercept
    if( includeOffset && 1 %in% which && grepl("ONEx", shrtlbls[1]) ){
      terms[[1]] <- terms[[1]] + x$offset
      shrtlbls[1] <- paste("offset", "+", shrtlbls[1])
    }   
    if(length(which)>1 & ask) par(ask=TRUE)
    if(commonRange){
      range <- range(terms)
      range[1] <- range[1]-0.02*diff(range)
    }else range <- NULL
    
    for(i in 1:length(terms)){
      # set values of predicted effect to missing if response is missing
      if(sum(is.na(x$response))>0) terms[[i]][is.na(x$response)] <- NA
      if(is.null(dots$ylim)){
        range <- range(terms[[i]], na.rm=TRUE)
      }else{
        range <- dots$ylim
      }
      
      if(length(time)>1){
        
        plotWithArgs(funplot, args=argsFunplot, 
                     myargs=list(x=time, y=terms[[i]], id=x$id, type="l", ylab="effect", lty=1, rug=FALSE,
                                 xlab=attr(time, "nameyind"), ylim=range, main=shrtlbls[i]))
        if(rug) rug(time)
      }else{
        plotWithArgs(plot, args=argsPlot, 
                     myargs=list(x=x$response-x$offset, y=terms[[i]], type="p", ylab="effect", 
                                 xlab="response-offset", ylim=range, main=shrtlbls[i])) 
      }
    }
    if(length(which)>1 & ask) par(ask=FALSE) 
    }       
}



########################################################################################

#' Extract information of a base-learner
#' 
#' Takes a base-learner and extracts information.
#' 
#' @param object a base-learner
#' @param what a character specifying the quantities to extract.
#' This can be a subset of "design" (default; design matrix), 
#' "penalty" (penalty matrix) and "index" (index of ties used to expand 
#' the design matrix)
#' @param asmatrix a logical indicating whether the the returned matrix should be 
#' coerced to a matrix (default) or if the returned object stays as it is 
#' (i.e., potentially a sparse matrix). This option is only applicable if \code{extract} 
#' returns matrices, i.e., \code{what = "design"} or \code{what = "penalty"}.
#' @param expand a logical indicating whether the design matrix should be expanded 
#' (default: \code{FALSE}). This is useful if ties were taken into account either manually 
#' (via argument \code{index} in a base-learner) or automatically for data sets with many 
#' observations. \code{expand = TRUE} is equivalent to \code{extract(B)[extract(B, what = "index"),]} 
#' for a base-learner \code{B}.
#' @param ... currently not used
#' @method extract blg
#' @seealso \code{\link[mboost]{extract}} for the \code{extract} function of the package mboost
extract.blg <- function(object, what = c("design", "penalty", "index"),
                        asmatrix = FALSE, expand = FALSE, ...){
  what <- match.arg(what)
  
  if(grepl("%O%", object$get_call()) | grepl("%Oz%", object$get_call())){
    object <- object$dpp( rep(1, NROW(object$model.frame()[[1]])) )    
  }else{
    object <- object$dpp(rep(1, nrow(object$model.frame())))
  }  
  return(extract(object, what = what,
                 asmatrix = asmatrix, expand = expand))
}
