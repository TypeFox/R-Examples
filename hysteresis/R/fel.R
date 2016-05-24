fel <-
function(x,y=NULL,method="harmonic2",period=NULL,subjects=NULL,times="unknown",subset=NULL,na.action=getOption("na.action"),control=nls.control(),boot=FALSE,...) {
  if (boot==TRUE) return(summary(fel(x,y,method,period,subjects,times,subset,na.action,control),...))
  felcall <- match.call()
  if (ncol(matrix(x)) > 2)
    times <- x[,3]
  dat <- xy.coords(x,y)
  #g
  if (!is.null(subset)) {
    dat$x<-dat$x[subset]; dat$y<-dat$y[subset];
    if (!is.null(subjects)) {
      if (!is.list(subjects)) {
        subjects<-factor(subjects[subset])
      }
        else subjects <- lapply(subjects,function (x) factor(x[subset]))}
    if (is.numeric(times))
      times<-times[subset]}
  if (!is.null(subjects)) {
    dat <- cbind("x"=dat$x,"y"=dat$y)
    if (is.numeric(times))
    ans <- by(cbind(dat,times),subjects,fel,method=method,period=period,na.action=na.action,control=control)
    else
      ans <- by(dat,subjects,fel,method=method,period=period,times=times,na.action=na.action,control=control)
    if (!is.list(subjects)) names(ans) <- levels(subjects)
    
    values <- t(sapply(ans,function (x) x["values"]$values))
    Std.Errors <- t(sapply(ans,function (x) x["Std.Errors"]$Std.Error))
    
    if (is.list(subjects)){
      subjectmat <- matrix(NA,nrow=prod(dim(ans)),ncol=length(subjects))
    if (length(subjects) > 2)  {
    for (i in 2:(length(subjects)-1)){
      subjectmat[,i] <- rep(dimnames(ans)[[i]],times=prod(dim(ans)[(i+1):length(subjects)]),each=prod(dim(ans)[1:(i-1)]))
    }}
      subjectmat[,1] <- rep(dimnames(ans)[[1]],times=prod(dim(ans)[2:length(subjects)]))
      subjectmat[,length(subjects)] <- rep(dimnames(ans)[[length(subjects)]],each=prod(dim(ans)[1:(length(subjects)-1)]))
      
      if (!is.null(names(subjects))) colnames(subjectmat) <- names(subjects)
      values <- data.frame(subjectmat,values)
      if (is.null(names(subjects))) colnames(subjectmat) <- colnames(values)[1:length(subjects)]
        
      if (method!="harmonic2" & method!="geometric") values <- values[,c(colnames(subjectmat),"b.x","b.y",
                                "cx","cy","retention","coercion","area",
                                "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
                                "semi.major","semi.minor","focus.x","focus.y","eccentricity")]
      else values <- values[,c(colnames(subjectmat),"b.x","b.y","phase.angle",
                               "cx","cy","retention","coercion","area",
                               "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
                               "semi.major","semi.minor","focus.x","focus.y","eccentricity")]
      Std.Errors <- data.frame(subjectmat,Std.Errors)
}
    ans <- list("models"=ans,"Estimates"=values,"Std.Errors"=Std.Errors)
    class(ans) <- "ellipsefitlist" 
    attr(ans,"call") <- felcall
    return(ans)
  }
  if (is.null(period))
    period <- length(dat$x)
 suppressWarnings(if (times=="unknown") {
    pred.method <- "find.times"
    dat <- data.frame(do.call(na.action,list(cbind("x"=dat$x,"y"=dat$y))))
    dat$times <-2*(0:(length(dat$x)-1))/period*pi }
  else if (is.numeric(times)){
    pred.method <- "times"
  dat <- data.frame(do.call(na.action,list(cbind("x"=dat$x,"y"=dat$y,times))))
    dat$times <-2*dat$times/period*pi}
  else {
    pred.method <- "times"
    dat <- data.frame(do.call(na.action,list(cbind("x"=dat$x,"y"=dat$y))))
    dat$times <-2*(0:(length(dat$x)-1))/period*pi })
 
  if (method=="harmonic2")
    ans <- fel.harmonic2(dat$x,dat$y,dat$times,period,pred.method)
    else if (method=="nls")
      ans <- fel.nls(dat$x,dat$y,dat$times,control,period,pred.method)
        else if (method=="direct")
        ans <- fel.direct(dat$x,dat$y,dat$times,period,pred.method)
          else if (method=="lm")
            ans <- fel.lm(dat$x,dat$y,dat$times,period,pred.method)
            else if (method=="geometric") {
              if (!is.numeric(control)) control <- 1.001
               ans <- fel.geometric(dat$x,dat$y,control=control,period=period)}
        else stop("method must be 'harmonic2', 'direct', 'geometric', 'lm' or 'nls'") 
  ans$call <- felcall
  if (method!="direct") ans$Std.Errors <- delta.error(ans)
  ans$Estimates <- ans$values
  if (method!="harmonic2" & method!="geometric")     ans$Estimates<- ans$Estimates[c("b.x","b.y",
                                              "cx","cy","retention","coercion","area",
                                              "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg",
                                              "semi.major","semi.minor","focus.x","focus.y","eccentricity")]
  else     ans$Estimates <- ans$Estimates[c("b.x","b.y","phase.angle",
                          "cx","cy","retention","coercion","area",
                          "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg",
                          "semi.major","semi.minor","focus.x","focus.y","eccentricity")]
        class(ans) <- "ellipsefit"
        ans}
