checkcalibration<-function (Xs, d, total, g, EPS=1e-6) 
{
 if(is.null(g)) stop("the g-weight vector is null")
 if(!is.matrix(Xs)) Xs<-as.matrix(Xs)
 tr <- crossprod(Xs, g*d)
 if(max(abs(tr-total)/total)<EPS) {result<-TRUE
                           message<-"the calibration is done"
                           value<-EPS  
                           }
else {message<-cat("the calibration cannot be done. The max EPS value is given by 'value'.\n")
      value<-max(abs(tr-total)/total)
      result<-FALSE}
list(message=message,result=result,value=value)
}
