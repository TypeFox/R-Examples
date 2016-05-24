`wlr_trafo` <- function(x,...){
    UseMethod("wlr_trafo")
}

`wlr_trafo.data.frame`<-function(x,...){
    ## the main purpose of this function is for ease of use of wlr_trafo in the coin package
    ## if x[[1]] is numeric this treats the time to event for the ith observation as happening at exactly x[[1]][i] 
    if (dim(x)[2]>1 | (class(x[[1]])!="Surv" & class(x[[1]])!="numeric")) stop("data.frame method must have one variable with class either 'Surv' or 'numeric' ")
    wlr_trafo(x[[1]],...)
}



`wlr_trafo.Surv`<-function(x,...){
    ## SurvLR changes Surv object to data.frame with L and R columns
    ## type=right, left or interval are allowed, type=counting is not
    LR<-SurvLR(x)
    wlr_trafo.default(LR$L,R=LR$R,...)
}

`wlr_trafo.default` <-
function(x, R=NULL, 
    scores =c("logrank1","logrank2","wmw","normal","general"), 
    icFIT=NULL,
    initfit=NULL, 
    control=icfitControl(),
    Lin=NULL,
    Rin=NULL,
    dqfunc=NULL,...){
    L<-x
    scores<-match.arg(scores)
    if (scores!="logrank1" & scores!="logrank2" & scores!="wmw" & scores!="normal" & scores!="general") stop("scores must equal 'logrank1' or 'logrank2' or 'wmw' or 'normal' or 'general' ")
    if (scores=="general" & is.null(dqfunc)) stop("when scores='general' then dqfunc must be supplied")
    ## we allow exact event times for all subjects to be input as L numeric and R=null, then R<-L get to usual
    if (is.null(R)) R<-L
    if (is.null(icFIT)){ 
        icFIT<-icfit(L,R,initfit,control,Lin,Rin)
        if (icFIT$message!="normal convergence") warning("icFIT does not have normal convergence")   
    }  

    A<-icFIT$A
    k<-dim(A)[[2]]
    n<-dim(A)[[1]]
    if (length(icFIT$pf)!=k) stop("icFIT$pf not proper length")

    cc<-scoresFromFit(icFIT,scores, dqfunc)   
    cc
}

