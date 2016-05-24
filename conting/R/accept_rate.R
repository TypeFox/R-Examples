accept_rate <-
function(object){

rj<-object$rj_acc
mh<-object$mh_acc

est<-list(rj_ar=100*mean(rj),mh_ar=100*mean(mh))

class(est)<-"acceptrate"

est}
