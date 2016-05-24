fun.check.gld.multi <-
function(l1,l2,l3,l4,param){
    
    n<-length(l1)
    out <- as.integer(rep(0,n));
    param <- switch(param, FKML= , fkml = , freimer = , frm = , FMKL = , fmkl = {
        ret <- .C("mult_check_gld",as.double(l1),as.double(l2),as.double(l3),
        as.double(l4),as.integer(1),as.integer(n),out)}
    , ramberg = , ram = , RS = , rs = {
        ret <- .C("mult_check_gld",as.double(l1),as.double(l2),as.double(l3),
        as.double(l4),as.integer(2),as.integer(n),out)})
    return(as.logical(ret[[7]]))
}
