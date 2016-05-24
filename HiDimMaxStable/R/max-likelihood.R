
density.memory<-new.env()

clear.memory<-function(nparam=2) {
    density.memory$mem<-matrix(nrow=0,ncol=nparam+1)
}

remember<-function(params) {
    i<-which(apply(density.memory$mem[,1:length(params),drop=FALSE],1,function(p) all(p==params)))
    if(length(i)>0) {
        density.memory$mem[min(i),length(params)+1]
    } else {
        NULL
    }
}

set.memory<-function(params,value)  {
    density.memory$mem<-rbind(density.memory$mem,c(params,value))
}

maxlik.f<-function(f,data,
                  params,
                  start,
                  p.min=rep(-Inf,length(start)),
                  p.max=rep(Inf,length(start)),
                  method="NM",
                  iterlim=80,tol=0.000001,
                  trace=TRUE,
                  ln=TRUE,
                  memory=TRUE,
                  ...) {
    # Builds constraints matrices A and B
    A<-rbind(diag(length(start)),-diag(length(start)))
    B<-c(-p.min,p.max)
    A<-A[is.finite(B),]
    B<-B[is.finite(B)]
    constr<-list(ineqA=A,ineqB=B)

    if(memory) clear.memory(nparam=length(start))

    maxLik(
        function(p,...) {
            pp<-params
            pp[is.na(pp)]<-p
            v<-NULL
            if(memory) v<-remember(p)
            if(is.null(v)) {
                if(trace) cat("computes density at",pp,"\n")
                v<-tryCatch(f(data,ln=ln,params=pp,...),
                            error=function(e) -Inf)
                if(memory) set.memory(p,v)
                if(trace) cat("--> ",v,"\n")
            }
            v
        },
        start=start,
        method=method,constraints=constr,iterlim=iterlim,tol=tol,
        ...)
}
