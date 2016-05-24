correct.trees<-function(x,TP,FN,TN,FP,study, data){

if (!is.element("imperfect.trees", class(x))) 
        stop("Argument 'x' must be an object of class \"imperfect.trees\".")
  sy<-x$Seny
  sx<-x$Senx
  ex<-x$Espx
  ey<-x$Espy
  pp<-x$Prevl

    if (missing(data)) 
        data <- NULL
    no.data <- is.null(data)
    if (is.null(data)) {
             
    if (length(TP)!=length(x$Prevl)) stop("Classification data don't match the multinomial tree model's output")
    if (is.element(FALSE,(study==x$study))) stop("Data aren't in the same order as those used in fitting the multinomial tree model, please sort the data first.")

TPnew <-TP *pp *sy*sx/(pp *sy*sx+(1-pp )*(1-ey)*(1-ex))+FP *pp *(1-sy)*sx/(pp *(1-sy)*sx+(1-pp )*ey*(1-ex))
FNnew <-FN *pp *sy*(1-sx)/(pp *sy*(1-sx)+(1-pp )*(1-ey)*ex)+TN *pp *(1-sy)*(1-sx)/(pp *(1-sy)*(1-sx)+(1-pp )*ey*ex)
FPnew <-FP *(1-pp )*ey*(1-ex)/(pp *(1-sy)*sx+(1-pp )*ey*(1-ex))+ TP *(1-pp )*(1-ey)*(1-ex)/(pp *sy*sx+(1-pp )*(1-ey)*(1-ex))
TNnew <-TN *(1-pp )*ey*ex/(pp *(1-sy)*(1-sx)+(1-pp )*ey*ex)+FN *(1-pp )*(1-ey)*ex/(pp *sy*(1-sx)+(1-pp )*(1-ey)*ex)
         
output<-data.frame(TPnew,FNnew,FPnew,TNnew,study)
    }

    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        } 
    mf <- match.call()
   
    mf.TP <- mf[[match("TP", names(mf))]]
    mf.FN <- mf[[match("FN", names(mf))]]
    mf.TN <- mf[[match("TN", names(mf))]]
    mf.FP <- mf[[match("FP", names(mf))]]
    mf.study <- mf[[match("study", names(mf))]]
    
    TP <- eval(mf.TP, data, enclos = sys.frame(sys.parent()))
    FN <- eval(mf.FN, data, enclos = sys.frame(sys.parent()))
    TN <- eval(mf.TN, data, enclos = sys.frame(sys.parent()))
    FP <- eval(mf.FP, data, enclos = sys.frame(sys.parent()))
    study <- eval(mf.study, data, enclos = sys.frame(sys.parent()))

if (length(TP)!=length(x$Prevl)) stop("Classification data don't match the multinomial tree model's output")
    if (is.element(FALSE,(study==x$study))) stop("Data aren't in the same order as those used in fitting the multinomial tree model, please sort the data first.")


     data$TPnew <-TP *pp *sy*sx/(pp *sy*sx+(1-pp )*(1-ey)*(1-ex))+FP *pp *(1-sy)*sx/(pp *(1-sy)*sx+(1-pp )*ey*(1-ex))
data$FNnew <-FN *pp *sy*(1-sx)/(pp *sy*(1-sx)+(1-pp )*(1-ey)*ex)+TN *pp *(1-sy)*(1-sx)/(pp *(1-sy)*(1-sx)+(1-pp )*ey*ex)
data$FPnew <-FP *(1-pp )*ey*(1-ex)/(pp *(1-sy)*sx+(1-pp )*ey*(1-ex))+ TP *(1-pp )*(1-ey)*(1-ex)/(pp *sy*sx+(1-pp )*(1-ey)*(1-ex))
data$TNnew <-TN *(1-pp )*ey*ex/(pp *(1-sy)*(1-sx)+(1-pp )*ey*ex)+FN *(1-pp )*(1-ey)*ex/(pp *sy*(1-sx)+(1-pp )*(1-ey)*ex)
 
   output<-as.data.frame(data) 
  }
   class(output)<-c("corrected.freq","data.frame")
   return(output)
}