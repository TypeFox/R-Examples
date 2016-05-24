##  ##############################################################################

lgnpp <- 
  function(
    object,
    newdata, locations,
    is.block = FALSE, all.pred = NULL, 
    extended.output = FALSE
  )
{
  
  ## Function post-processes predictions obtained from predict.georob
  ## for lognormal point and block kriging
  
  ## Arguments:
  
  ## object            object of class "predict.georob" (output of predict.georob)
  ## what              character scalar controlling what kind of 
  ##                   post-processing is done. Possible values are
  ##                   "point"     backtransforms point kriging predictions 
  ##                               and kriging standard errors to the 
  ##                               original scale of measurements
  ##                   "perman"    backtransforms block kriging predictions 
  ##                               and the respective standard errors to the 
  ##                               original scale of measurements using 
  ##                               the assumptions of permanence of 
  ##                               lognormality
  ##                   "optimal"   given the prediction results 
  ##                               (predictions, full covariance matrices of 
  ##                               predictions, prediction errors, etc) for a set 
  ##                               of points in a block , the function computes the optimal 
  ##                               lognormal block prediction and the associated 
  ##                               mean squared error
  ## newdata    applies for what == "perman" only: a data frame with coordinates and 
  ##                   explanatory variables for a grid of points that discretize all the 
  ##                   blocks for which block kriging predictions are of the log-transformed
  ##                   response variable are contained in 'object'
  ## locations         applies for what == "perman" only: a one-sided formula that defines the 
  ##                   coordinates in the data frame 'newdata'
  ## all.pred   applies for what == "optimal" only: a data.frame with the 
  ##                   point prediction results of the lognormally if 'object' contains 
  ###                  back-transformed response variable for all the points in the block 
  ##                   the prediction results only for a subset of these points
  ## extended.output   applies for what == "optimal" only: if TRUE the covariance matrix
  ##                   of the back-transformed point prediction errors are returned as an attribute
  ##                   of the result

  ## Value:
  
  ## For what == "point" an amended object of class "predict.georob" with
  ## the additional variables "lgn.pred", "lgn.se", "lgn.lower",
  ## "lgn.upper", containing the predictions, the prediction standard
  ## error, and the bounds of the prediction interval on the original
  ## scale of the observations.  For what == "perman" and what ==
  ## "optimal" a named vector with approximated prediction and the
  ## prediction mean squared error of the block mean on the original scale
  ## of the measurements.
  
  ## References:
  
  ## @article{Cressie-2006,
  ##     Author = {Cressie, N.},
  ##     Journal = {Mathematical Geology},
  ##     Keywords = {Lognormal Block Kriging},
  ##     Number = {4},
  ##     Pages = {413--443},
  ##     Separatanr = {274},
  ##     Title = {Block Kriging for Lognormal Spatial Processes},
  ##     Volume = {38},
  ##     Year = {2006}
  ## }
  
  
  ## 2011-12-22 A. Papritz 
  ## 2012-05-07 AP backtransformation of block predictions under permanence 
  ## of lognormality assumption
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-06-23 AP modifications for missing compontents 'lower' and 'upper'
  
  ## auxiliary function to backtransform the point predictions and 
  ## optionally the mean squared errors and the prediction intervals
  
 f.backtrf <- function( object, btMb = 0., pred.only ){
    
    ## object: dataframe with the prediction results of the 
    ##         log-transformed predictions computed with extended.output = TRUE 
    ## btMb  : vector with the terms t( coefficients ) M(B) coefficients for back-transformation
    ##         of block kriging results under assumption of permanence of lognormality
    ## pred.only: logical scalar controlling wether only predictions are computed
    
    t.result <- exp( 
      object[, "pred"] + 0.5 * ( object[, "var.target"] + btMb - object[, "var.pred"] )
    )
    
    if( !pred.only ){
      t.mu <- exp( object[, "trend"] + 0.5 * object[, "var.target"] )
      t.result <- data.frame(
        lgn.pred = t.result,
        lgn.se = t.mu * sqrt( 
          exp( object[, "var.target"] ) + exp( object[, "var.pred"] ) - 
          2 * exp( object[, "cov.pred.target"] )                       
        )
      )
      if( all( c( "lower", "upper" ) %in% names( object ) ) ){
        t.result <- cbind( t.result,
          data.frame(
            lgn.lower = exp( object[, "lower"] ),
            lgn.upper = exp( object[, "upper"] )
          )
        )
        
      }
    }
    
    return( t.result )
    
  }
  
  ## check whether back-transformation for object has been already done
  
  if( all( c( "lgn.pred", "lgn.se" ) %in% names( object ) ) ){
    return( object )  
  }
  
  ## main body of function
  
  ## find out what object contains
  
  what <- switch(
    class( object ),
    "data.frame" =,
    "SpatialPointsDataFrame" =,
    "SpatialPixelsDataFrame" =,
    "SpatialGridDataFrame" = "point",
    "SpatialPolygonsDataFrame" = "block",
    "list" = switch(
      class( object[["pred"]] ),
      "data.frame" =,
      "SpatialPointsDataFrame" =,
      "SpatialPixelsDataFrame" =,
      "SpatialGridDataFrame" = "point",
      "SpatialPolygonsDataFrame" = "block",
      stop( "'object' not generated by predict method for class 'georob'" )
    ),
    stop( "'object' not generated by predict method for class 'georob'" )
  )
  
  if( what == "point" && ( is.block || !is.null( all.pred ) ) ) what <- "optimal"
  
  ## extract data frame with prediction results
  
  t.object <- switch(
    class( object ),
    "data.frame" = object,
    "SpatialPointsDataFrame" =,
    "SpatialPixelsDataFrame" =,
    "SpatialGridDataFrame" =,
    "SpatialPolygonsDataFrame" = object@data,
    "list" = switch(
      class( object[["pred"]] ),
      "data.frame" = object[["pred"]],
      "SpatialPointsDataFrame" =,
      "SpatialPixelsDataFrame" =,
      "SpatialGridDataFrame" =,
      "SpatialPolygonsDataFrame" = object[["pred"]]@data
    )
  )
    
  ## extract variogram and type of predictions
  
  variogram.model   <- attr( t.object, "variogram.model" )
  param             <- attr( t.object, "param" )
  type              <- attr( t.object, "type" )
  
  if( is.null( variogram.model ) || is.null( param ) || is.null( type ) ) stop(
    "attributes 'variogram.model', 'param' or 'type' missing in 'object'"  
  )
  
  if( variogram.model %in% control.georob()[["irf.models"]] ) stop(
    "lognormal kriging requires a weakly stationary variogram model"
  )
  
  ## check whether object is complete
  
  if( is.null( t.object[["pred"]]) || is.null( t.object[["trend"]] ) || 
    is.null( t.object[["var.pred"]] ) || is.null( t.object[["cov.pred.target"]] ) ||
    is.null( t.object[["var.target"]] )
  ) stop( "some required items are missing, ", 
    "re-run 'predict' with argument 'extended.output = TRUE'" 
  )
  
  ## and now do the backtransformation for ...
  
  result <- switch(
    
    what,
    
    ## back-transform point predictions
    
    point = {
      
      t.result <- cbind(
        t.object,
        f.backtrf( t.object, pred.only = FALSE )
      )
      
      attr( t.result, "variogram.model" ) <- variogram.model
      attr( t.result, "param" )           <- param
      attr( t.result, "type" )            <- type
      
      result <- switch(
        class( object ),
        "data.frame" = t.result,
        "SpatialPointsDataFrame" =,
        "SpatialPixelsDataFrame" =,
        "SpatialGridDataFrame" =,
        "SpatialPolygonsDataFrame" = {
          result <- object
          result@data <- t.result
          result        
        },
        "list" = switch(
          class( object[["pred"]] ),
          "data.frame" = {
            result <- object
            result[["pred"]] <- t.result
            result          
          },
          "SpatialPointsDataFrame" =,
          "SpatialPixelsDataFrame" =,
          "SpatialGridDataFrame" =,
          "SpatialPolygonsDataFrame" = {
            result <- object
            result[["pred"]]@data <- t.result
            result          
          }
        )
      )
      
      result
      
    },
    
    block = {
      
      ## back-transform block predictions under the assumption of permanence
      ## of lognormality
      
      if( missing( newdata ) ) stop( 
        "'newdata' must be provided for back-transformation under",
        " assumption of permance of lognormality "
      )
      
      ## extract terms, regression coefficients and locations formula
      
      tt           <- attr( t.object, "terms" )
      coefficients <- attr( t.object, "coefficients" )
      locations    <- attr( t.object, "locations" )
      
      if( is.null( tt ) || is.null( coefficients ) || is.null( locations ) ) stop(
        "attributes 'terms', 'coefficients' or 'locations' missing in object"  
      )
      
      ## compute spatial covariance matrices of explanatory covariates
      
      if( !missing( locations ) ){
        locations <- as.formula( 
          paste( deparse( locations ), "-1" ), env = parent.frame() 
        )  
      }
      
      ## extract modelframe and coordinates of prediction points
      
      Terms <- delete.response( tt )
      Terms.loc <- locations
      
      ## get the model frame for newdata
      
      mf.newdata <- switch( 
        class( newdata ),
        "data.frame" = model.frame( 
          Terms, newdata, na.action = na.pass, xlev = object[["xlevels"]] 
        ),
        "SpatialPointsDataFrame" = model.frame(
          Terms, slot( newdata, "data" ), na.action = na.pass, 
          xlev = object[["xlevels"]] 
        ),
        "SpatialPixelsDataFrame" = model.frame(
          Terms, slot( newdata, "data" ), na.action = na.pass, 
          xlev = object[["xlevels"]] 
        ),
        "SpatialGridDataFrame" = model.frame(
          Terms, slot( newdata, "data" ), na.action = na.pass, 
          xlev = object[["xlevels"]] 
        ),
        "SpatialPolygonsDataFrame" = model.frame(
          Terms, slot( newdata, "data" ), na.action = na.pass, 
          xlev = object[["xlevels"]] 
        ),
        stop(
          "cannot construct model frame for class(newdata) ='", 
          class( newdata ) 
        )
      )
      
      ## check whether variables that will be used to compute the
      ## predictions agree with those in object
      
      if( !is.null( cl <- attr(Terms, "dataClasses" ) ) )
        .checkMFClasses( cl, mf.newdata )
      
      ## get fixed effects design matrix for newdata
      
      pred.X <- model.matrix( Terms, mf.newdata, contrasts.arg = object[["contrasts"]] )
      
      
      #       ## deal with non-NULL offset
      #       
      #       offset <- rep( 0, NROW(pred.X) )
      #       if( !is.null( off.num <- attr( tt, "offset" ) ) ){
      #         warning( "prediction with non-zero offset not yet debugged" )
      #         for( i in off.num ) {
      #           offset <- offset + eval( attr( tt, "variables" )[[i + 1]], newdata )
      #         }
      #       }
      #       if( !is.null( object[["call"]][["offset"]] ) ){
      #         offset <- offset + eval( object[["call"]][["offset"]], newdata )
      #       }
      
      ## get matrix of coordinates of newdata for point kriging
      
      pred.coords <- switch( 
        class( newdata ),
        "data.frame" = model.matrix(
          Terms.loc,
          model.frame( 
            Terms.loc, newdata, na.action = na.pass 
          )
        ),
        "SpatialPointsDataFrame" = model.matrix(
          Terms.loc,
          model.frame(  
            Terms.loc, as.data.frame( coordinates( newdata ) ),
            na.action = na.pass
          )
        ),
        "SpatialPixelsDataFrame" = model.matrix(
          Terms.loc,
          model.frame(  
            Terms.loc, as.data.frame( coordinates( newdata ) ),
            na.action = na.pass
          )
        ),
        "SpatialGridDataFrame" = model.matrix(
          Terms.loc,
          model.frame(  
            Terms.loc, as.data.frame( coordinates( newdata ) ),
            na.action = na.pass
          )
        ),
        stop( 
          "newdata must be of class 'data.frame', 'SpatialPointsDataFrame', 'SpatialPixelsDataFrame'\n",
          "  or 'SpatialGridDataFrame'"
        )
      )
      
      ## indices of polygons to which prediction points belong to
      
      t.ip <- over( 
        t.points <- SpatialPoints( pred.coords ), 
        t.polygons <- if( is.list( object ) ){
          geometry( object[["pred"]] )
        } else {
          geometry( object )
        }
      )
      
      ## compute t(coefficients) %*% Cov( covariates ) %*% coefficients
      ## cf Cressie, 2006, Eq. 18 & Appendix C
      
      btMb <- as.vector( 
        tapply(
          1:length( t.ip ),
          factor( t.ip ),
          function( i, x, b ){
            drop( b %*% cov( x[i, , drop = FALSE ] ) %*% b )
          },
          x = pred.X,
          b = coefficients
        )
      )
      
      ## compute back-transformed block predictions
      
      t.result <- cbind(
        t.object,
        f.backtrf( t.object, btMb, pred.only = FALSE )
      )
      
      attr( t.result, "variogram.model" ) <- variogram.model
      attr( t.result, "param" )           <- param
      attr( t.result, "type" )            <- type
      attr( t.result, "terms" )           <- tt
      attr( t.result, "coefficients" )    <- coefficients
      attr( t.result, "locations" )       <- locations
      
      result <- object
      if( class( object ) == "SpatialPolygonsDataFrame" ){
        result@data <- t.result
      } else {
        result[["pred"]]@data <- t.result
      }
      
      result
        
    },
    
    optimal = {
      
      ## number of prediction items in object
      
      n.sample <- sum( complete.cases( t.object ) )
      
      if( is.null( all.pred ) ){
        
        ## object contains the predictions of all the points in the block
        
        n <- n.sample
        
      } else if( is.numeric( all.pred ) && length( all.pred == 1 ) ){
        
        ## complete predictions contains total number of points in block
        
        n <- all.pred
        
      } else if( 
        is.data.frame( all.pred ) ||
        class( all.pred ) == "SpatialPointsDataFrame" ||
        class( all.pred ) == "SpatialPixelsDataFrame" ||
        class( all.pred ) == "SpatialGridDataFrame" 
      ){
        
        ## all.pred contains back-transformed point predictions of
        ## all the points in the block, extract the respective data frame
        
        all.pred <- switch(
          class( all.pred ),
          data.frame = all.pred,
          SpatialPointsDataFrame =,
          SpatialPixelsDataFrame =,
          SpatialGridDataFrame = all.pred@data
        )
                
        ## check whether back-transformation were done
        
        if( !all( c( "lgn.pred", "lgn.se", "lgn.lower", "lgn.upper" ) %in% names( all.pred ) ) ){ 
          all.pred <- f.backtrf( all.pred, pred.only = T )
        }
                
        ## check wether variogram models in all.pred and object match
        
        if( !all(
            c(
              identical( variogram.model, attr( all.pred, "variogram.model" ) ),
              identical( param, attr( all.pred, "param" ) ),
              identical( type, attr( all.pred, "type" ) )
            )
          )
        ) stop( "variogram or prediction type of 'object' and 'all.pred' differs" )
        
        n <- sum( complete.cases(all.pred ) )
        
      } else {
        
        stop( 
          "'all.pred' does neither contain predictions results", 
          "nor total number of points in block" 
        )
        
      }
      
      if( n < n.sample ) stop( "nrow(all.pred) < nrow(object$pred)" )
      
      ## optimal lognormal backtransformation for block prediction
      
      if( is.null( object[["pred"]] ) || is.null( object[["mse.pred"]] ) || 
        is.null( object[["var.pred"]] ) || is.null( object[["cov.pred.target"]] )
      ) stop( "some required items are missing in 'object'\n", 
        "re-run 'predict' with argument 'full.covmat = TRUE'" 
      )
      
      ## back-transform trend and kriging predictions for the points in the block
      
      t.mu <- exp( t.object[, "trend"] + 0.5 * t.object[, "var.target"] )
      t.pred <- f.backtrf( t.object, pred.only = TRUE )
      
      ## covariance matrix of log(observations)
      
      t.cov.log.obs <- object[["mse.pred"]] - object[["var.pred"]] + 
        object[["cov.pred.target"]] + t( object[["cov.pred.target"]] )
      
      
      ## covariance matrix of backtransformed prediction errors
      
      t.aux <- exp( t.cov.log.obs ) + exp( object[["var.pred"]] ) - 
        exp( object[["cov.pred.target"]] ) - exp( t( object[["cov.pred.target"]] ) )
      
      t.cov.errors <- t( t.mu * t( t.mu * t.aux ) )
      
      ## compute weighted mean of diagonal and off-diagonal elements of
      ## covariance matrix of back-transformed point kriging predictions
      
      t.sum.diag <- sum( diag( t.cov.errors ), na.rm = TRUE )
      t.sum.off.diag <- sum( t.cov.errors, na.rm = TRUE ) - t.sum.diag
      
      if( !is.data.frame( all.pred ) ){
        
        t.mean <- mean( t.pred, na.rm = TRUE )
        
        t.mse <- ( 
          t.sum.diag / n.sample * n + 
          t.sum.off.diag / ( n.sample * (n.sample-1) ) * n * (n-1)
        ) / n^2
        
        
      } else {
        
        t.mean <- mean( all.pred[["lgn.pred"]], na.rm = TRUE )
        
        t.sum.diag <- sum( all.pred[["lgn.se"]]^2, na.rm = TRUE )
        
        t.mse <- ( 
          t.sum.diag + t.sum.off.diag / ( n.sample * (n.sample-1) ) * n * (n-1)
        ) / n^2
        
      }
      
      t.result <- c( pred = t.mean, mse = t.mse )
      
      if( extended.output ){
        attr( t.result, "mse.lgn.pred" ) <- t.cov.errors
      }
      
      t.result
    }
  )
    
  invisible( result )
  
}
