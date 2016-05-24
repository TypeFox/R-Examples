ps <-
function(x,monotone=0,lambda=0,pdiff=3,ndx=NULL,deg=3,var.pen=NULL){
#se ndx NULL the empirical rule of Ruppert (2001) is used: min(n/4,40)
#se lambda<0 allora viene stimato, altrimenti deve essere un valore numerico
#var.pen: una stringa del tipo "1:k" per varying penalty
#monotone: 0: unconstrained; +1: non-descreasing; -1= non-increasing (NB sign(T)=1 and sign(F)=0)
#opzione penalty?
    r<-x
    attr(r,"ndx")<-ndx
    attr(r,"deg")<-deg
    attr(r,"pdiff")<-pdiff
    attr(r,"monot")<-monotone #isTRUE(monotone)
    attr(r,"lambda")<-lambda
    attr(r,"nomeX")<-deparse(substitute(x))
    attr(r,"var.pen")<-var.pen
    r
    }
