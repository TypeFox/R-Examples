getcutoff <-
function(stat,alpha,reverse) {
    nr<-length(stat)
    if (alpha<0 || alpha>1) stop("alpha must be a fraction between 0 and 1")
    if (nr<100) warning("number of replications is less than 100")
    if (nr<=0) stop("number of replications is not a positive integer")
    nr<-ceiling(nr)
    cutpos<-round(nr*alpha)
    if (cutpos<=0) cutpos=1
    else if (cutpos>nr) cutpos=nr
    if (reverse) out<-rev(sort(stat))[cutpos]
    else out<-sort(stat)[cutpos]
    return(out)
  }
