#' @export
#  see predict.ds for documentation
predict.trial.fi <- function(object, newdata=NULL,compute=FALSE, int.range=NULL,
                             integrate=FALSE,...){
  # Functions Used: pdot.dsr.integrate.logistic, is.linear.logistic,
  #                 predict.glm (could also use predict.gam eventually)
  model <- object
  point <- model$meta.data$point
  width <- model$meta.data$width

  # set newdata to be the model data if we don't have any
  if(is.null(newdata)){
    newdata <- model$data
    newdata <- newdata[newdata$object %in% as.numeric(names(model$fitted)),]
  }

  if(integrate){
    left <- model$meta.data$left
    formula <- paste("~",as.character(model$mr$formula)[3],collapse="")

    if (class(model$mr)[1]=="gam"){
      integral.numeric <- TRUE
    }else{
      integral.numeric <- is.linear.logistic(newdata,formula,
                                             length(coef(model$mr)),width)
    }

    models<- list(g0model        = formula,
                  scalemodel     = NULL,
                  fullscalemodel = NULL)

    if(is.null(int.range)){
      pdot.list <- pdot.dsr.integrate.logistic(width, width, coef(model$mr),
                                               newdata,integral.numeric, TRUE,
                                               models, point=point)
    }else{
      pdot.list <- pdot.dsr.integrate.logistic(int.range, width, coef(model$mr),
                                               newdata,integral.numeric, TRUE,
                                               models, point=point)
    }

    if(left !=0){
      pdot.list$pdot <- pdot.list$pdot -
                        pdot.dsr.integrate.logistic(left, width, coef(model$mr),
                                                    newdata,integral.numeric,
                                                    TRUE, models,
                                                    point=point)$pdot
    }

    fitted <- pdot.list$pdot
  }else{
    p1 <- predict(model$mr,newdata,type="response")
    p2 <- 0
    fitted <- p1

    if(is.null(newdata)){
      names(fitted) <- model$mr$data$object[model$mr$data$observer==1]
    }else{
      names(fitted) <- newdata$object[newdata$observer==1]
    }

    return(list(p1     = p1,
                p2     = 0,
                fitted = fitted))
  }

  if(is.null(newdata)){
    names(fitted) <- model$mr$data$object[model$mr$data$observer==1]
  }else{
    names(fitted) <- newdata$object[newdata$observer==1]
  }

  return(list(fitted=fitted))
}
