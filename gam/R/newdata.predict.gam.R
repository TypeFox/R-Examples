"newdata.predict.gam" <-
  function(object, newdata, type = c("link", "response", "terms"), dispersion=NULL, se.fit = FALSE, na.action=na.pass,terms=labels(object), ...)
{
  out.attrs <- attr(newdata, "out.attrs")
  is.gam<-inherits(object, "gam") && !is.null(object$smooth)
 if(is.gam) {
   if(se.fit){
     se.fit<-FALSE
     warning("No standard errors (currently) for gam predictions with newdata")
   }
   ##First get the linear predictions
   type <- match.arg(type)
   local.type<-type
   if(type=="response")local.type<-"link"
   pred<-predict.glm(object,newdata,type=local.type,dispersion=dispersion,se.fit=FALSE,terms=terms)
   ##Build up the smooth.frame for the new data
   tt <- terms(object)
    Terms <- delete.response(tt)
    smooth.frame <- model.frame(Terms, newdata, na.action = na.action, 
                     xlev = object$xlevels)
   nrows<-nrow(smooth.frame)
   old.smooth<-object$smooth
   data<-object$smooth.frame # this was the old smooth frame
   smooth.labels<-names(data)
   n.smooths<-length(smooth.labels)
   if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, smooth.frame)
    out.attrs <- attr(newdata, "out.attrs")
  

   w <- object$weights
   pred.s <- array(0, c(nrows, n.smooths), list(row.names(smooth.frame), 
                                                 smooth.labels))
   smooth.wanted <- smooth.labels[match(smooth.labels, terms,
                                         0) > 0]
   pred.s<-pred.s[,smooth.wanted,drop=FALSE]
    residuals <- object$residuals
    for(TT in smooth.wanted) {
      Call <- attr(data[[TT]], "call")
      Call$xeval <- substitute(smooth.frame[[TT]], list(TT = TT))
      z <- residuals + object$smooth[, TT]
       pred.s[, TT] <- eval(Call)
    }
    if(type == "terms")
      pred[, smooth.wanted] <- pred[, smooth.wanted] + pred.s[
                                                              , smooth.wanted]
    else pred <- pred + rowSums(pred.s)
   if(type == "response") {
     famob <- family(object)
     pred <- famob$linkinv(pred)
   }
  }
  else {
    pred<-predict.glm(object,newdata,type=type,dispersion=dispersion,se.fit=se.fit,terms=terms)
  }
  if(type != "terms" && !is.null(out.attrs)) {
    if(!is.null(out.attrs)) {
      if(se.fit) {
        attributes(pred$fit) <- out.attrs
        attributes(pred$se.fit) <- out.attrs
      }
      else attributes(pred) <- out.attrs
    }
  }
pred
}
