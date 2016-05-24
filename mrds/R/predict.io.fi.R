#' @export
# see predict.ds for documentation
predict.io.fi <- function(object,newdata=NULL,compute=FALSE, int.range=NULL,
                          integrate=FALSE,...){
  # Functions Used: pdot.dsr.integrate.logistic, is.linear.logistic,
  #                 predict.glm (could also use predict.gam eventually)
  model <- object
  width <- model$meta.data$width
  point <- model$meta.data$point

  # if no new data supplied, use the data from the model
  if(is.null(newdata)){
    newdata <- model$mr$data
  }else{
    if(!("observer" %in% names(newdata))){
      stop("newdata does not contain a column named \"observer\"")
    }
  }

  newdata$offsetvalue <- 0

  GAM <- FALSE
  if("gam" %in% class(model$mr)){
    GAM <- TRUE
  }

  # integrate=FALSE -- predict p(y)
  if(!integrate){
    # predict here is predict.glm
    # since model here is model$mr, so model$mr is model$mr$mr
    fitted <- predict(model$mr,newdata,type="response")

    p1 <- fitted[newdata$observer==1]
    p2 <- fitted[newdata$observer==2]
    fitted <- p1+p2-p1*p2

    names(fitted) <- newdata$object[newdata$observer==1]

    return(list(fitted = fitted,
                p1     = p1,
                p2     = p2))

  }else{
    # integrate=TRUE -- get the average probability of detection
    left <- model$meta.data$left
    formula <- paste("~",as.character(model$mr$formula)[3],collapse="")

    if("gam" %in% class(model$mr)){
      integral.numeric <- TRUE
    }else{
      integral.numeric <- is.linear.logistic(newdata,formula,
                                             length(coef(model$mr)),width)
    }

    models <- list(g0model        = formula,
                   scalemodel     = NULL,
                   fullscalemodel = NULL)

    # now int.range is a vector with lower and upper bounds
    if(is.null(int.range)){
      pdot.list <- pdot.dsr.integrate.logistic(width, width, model$mr$coef,
                     newdata, integral.numeric, FALSE, models,GAM, point=point)
    }else{
      pdot.list <- pdot.dsr.integrate.logistic(int.range,width, model$mr$coef,
                      newdata, integral.numeric, FALSE, models,GAM, point=point)
    }

    # if there is left truncation, take that off the integral
    if(left !=0){
      pdot.list$pdot <- pdot.list$pdot -
                    pdot.dsr.integrate.logistic(left, width, model$mr$coef,
                                   newdata, integral.numeric, FALSE, models,
                                   GAM, point=point)$pdot
    }

    fitted <- pdot.list$pdot
    names(fitted) <- newdata$object[newdata$observer==1]

    return(list(fitted=fitted))
  }
}
