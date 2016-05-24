# Methods for s$ generic classifyUnsupv


setGeneric(
  name= "classifySupv",
  def=function(model, newdata, ...){standardGeneric("classifySupv")}
)

#' Method for old (S3-style) objects, replaces ordinary function
setMethod(
  f = "classifySupv",
  signature = c("RecLinkClassif", "RecLinkData"),
  definition = function (model, newdata, convert.na=TRUE, ...)
  {

    # type checks from previous version omitted, now enforced by 
    # method dispatching  
        
  	ret=newdata
  
    x=newdata$pairs[,-c(1,2,ncol(newdata$pairs))]
    if(any(colnames(x)!=model$attrNames))
    {
      warning("Attribute names in newdata differ from training set!")
      colnames(x)=model$attrNames
    }
    if(convert.na) x[is.na(x)]=0
  
      predict=switch(model$method,
    		svm=predict(model$model, newdata=x,...),       
  	 	 rpart=predict(model$model, newdata=x,type="class",...),       
  		  ada=predict(model$model, newdata=x,type="vector",...),       
  		  bagging=predict(model$model, newdata=x,type="class",...),
  		  nnet=predict(model$model, newdata=x,type="class",...),
        stop("Illegal classification method!"))
      # refactor to ensure uniform order of levels
      ret$prediction=factor(predict,levels=c("N","P","L"))
      class(ret)="RecLinkResult"
      return(ret)
  }
)

# Methods for big data sets
setMethod(
  f = "classifySupv",
  signature = c("RecLinkClassif", "RLBigData"),
  definition = function(model, newdata, convert.na = TRUE, withProgressBar = (sink.number()==0), ...)
  {
    prediction=switch(model$method,
  		svm=.ffpredict(model$model, newdata=newdata@pairs, withProgressBar, convert.na, ...),
	 	 rpart=.ffpredict(model$model, newdata=newdata@pairs, withProgressBar, convert.na, type="class",...),
		  ada=.ffpredict(model$model, newdata=newdata@pairs, withProgressBar, convert.na, type="vector",...),
		  bagging=.ffpredict(model$model, newdata=newdata@pairs, withProgressBar, convert.na, type="class",...),
		  nnet=.ffpredict(model$model, newdata=newdata@pairs, withProgressBar, convert.na, type="class",...),
      stop("Illegal classification method!"))

    result <- new("RLResult", data = newdata, prediction = prediction)
  }
)

if(getRversion() >= "2.15.1")  utils::globalVariables(c("i1", "i2"))
.ffpredict <- function(model, newdata, withProgressBar, convert.na, ...)
{
    if (withProgressBar)
    {
      pgb <- txtProgressBar(0, nrow(newdata))
    }

    prediction <- ff("N", length = nrow(newdata), levels = c("N", "P", "L"))
    ffrowapply(
      {
        slice <- newdata[i1:i2,,drop = FALSE]
        if (convert.na) slice[is.na(slice)] <- 0
        prediction[i1:i2] <- predict(model, slice, ...)
        if (withProgressBar) setTxtProgressBar(pgb, i2)
      },
      X = newdata)
    if (withProgressBar) close(pgb)
    prediction
}
