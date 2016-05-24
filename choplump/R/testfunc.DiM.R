`testfunc.DiM` <-
function(choplist,SM){
    ZM<-choplist$ZM
    a<-choplist$a
    b<-choplist$b
    Z<-c(rep(0,a),rep(1,b),ZM)
    S<-c(rep(0,a+b),SM)
    ## in order to have the direction be the same 
    ## as wilcox.test use -TDiM
    return( -TDiM(S,Z) ) 
}

