summary.frailtyHL <-
function (object,...) {
    ans<-list(0)
    names(ans)[1]<-"Model"
    if (object$Model == "gamma frailty model") ans$Model<-"Results from the gamma frailty model"
    if (object$Model == "log-normal frailty model") ans$Model<-"Results from the log-normal frailty model"
    ans$formula<-object$formula
    if (object$Method == "HL(0,1)") ans$method<-"Method : HL(0,1)"
    if (object$Method == "HL(0,2)") ans$method<-"Method : HL(0,2)"
    if (object$Method == "HL(1,1)") ans$method<-"Method : HL(1,1)"
    if (object$Method == "HL(1,2)") ans$method<-"Method : HL(1,2)"
    ans$FixCoef<-object$FixCoef
    ans$RandCoef<-object$RandCoef
    ans$likelihood<-object$likelihood
    ans$iter<-object$iter
    ans$convergence<-object$convergence
    ans$aic<-object$aic
    class(ans) <- "summary.frialtyHL"
    ans
}

