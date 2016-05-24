#library(compiler)

F2norm<-cmpfun(function(M)sqrt(sum(M^2)))

thr.op<-cmpfun(function(x,thr){sign(x)*pmax(abs(x)-thr,0)})

#La.svd.cmp<-cmpfun(La.svd)

thresh.nuclear<-cmpfun(function(M,thr){
    s<-La.svd(M)
    dd<-thr.op(s$d,thr)
    id<-which(dd!=0)
    s$d<-dd[id]
    s$u<-s$u[,id,drop=FALSE]
    s$vt<-s$vt[id,,drop=FALSE]
    s$L<-s$u%*%(s$d*s$vt)
    s
})

thresh.l1<-cmpfun(function(x,thr){
    thr.op(x,thr)
})


rpca<-function(M,
             lambda = 1/sqrt(max(dim(M))),
             mu = prod(dim(M))/(4*sum(abs(M))),
             term.delta=10^(-7),
             max.iter=5000,
             trace=FALSE,
             thresh.nuclear.fun=thresh.nuclear,
             thresh.l1.fun=thresh.l1,
             F2norm.fun=F2norm){
    dm<-dim(M)
    term.norm<-term.delta*F2norm.fun(M)
    S<-Yimu<-matrix(0,nrow=dm[1],ncol=dm[2])

    imu<-1/mu
    limu<-lambda/mu

    i<-0
    stats<-c()
    converged<-FALSE
    while(TRUE){
        i<-i+1
        L.svd<-thresh.nuclear.fun(M-S+Yimu,imu)
        L<-L.svd$L
        S<-thresh.l1.fun(M-L+Yimu,limu)

        MLS = M-L-S

        resid.norm<-F2norm.fun(MLS)
        stats<-c(stats,resid.norm)
        if (trace) 
            print(c(iter=i,resid.norm=resid.norm))
        converged<-resid.norm<term.norm
        if ((i>max.iter)||converged)
            break;

        Yimu = Yimu + MLS
    }
    final.delta<-resid.norm*term.delta/term.norm
    if (!converged)
        warning(paste("rpca did not converge after",i,"iterations, final.delta=",final.delta))
    list(L=L,S=S,L.svd=L.svd,
         convergence=list(converged=converged,
                          iterations=i,
                          final.delta=final.delta,
                          all.delta=stats*(term.delta/term.norm)))
}

# #require(gputools)
# 
# F2norm.gpu<-cmpfun(function(M)sqrt(gputools::gpuCrossprod(as.vector(M))))
# 
# #Below function could be in principle be gpu improved, 
# #but gputools do not expose required simple multiplication functions
# thr.op.gpu<-cmpfun(function(x,thr){sign(x)*pmax(abs(x)-thr,0)})
# 
# thresh.nuclear.gpu<-cmpfun(function(M,thr){
#     s<-gputools::gpuSvd(M)
#     dd<-thr.op.gpu(s$d,thr)
#     id<-which(dd!=0)
#     s$d<-dd[id]
#     s$u<-s$u[,id,drop=FALSE]
#     s$v<-s$v[,id,drop=FALSE]
#     s$vt<-t(s$v)
#     s$L<-gputools::gpuMatMult(s$u,s$d*s$vt)
#     s
# })
# 
# 
# thresh.l1.gpu<-cmpfun(function(M,thr){
#     thr.op.gpu(M,thr)
# })
# 
# rpca.gpu<-function(M,
#              lambda = 1/sqrt(max(dim(M))),
#              mu = prod(dim(M))/(4*sum(abs(M))),
#              term.delta=10^(-5),
#              max.iter=5000,
#              trace=FALSE,
#              gpu.to.choose=NULL){
#     if (!requireNamespace("gputools", quietly = TRUE)) {
#         stop("package 'gputools' version 0.26 is needed for this function to work. Please install it with CULA library: install.packages(\"gputools_0.26.tar.gz\", configure.args=\"--with-cuda-home=/opt/cuda --with-cula-home=/opt/cula\")", call. = FALSE)
#     }
#     try(attachNamespace("gputools"),silent=TRUE)
#     if (!exists("gpuSvd")) 
#         stop("library(gputools) needs to be compiled with CULA library, and its version (0.26) needs to implement gpuSvd: install.packages(\"gputools_0.26.tar.gz\",configure.args=\"--with-cula-home=/opt/cula\")")
#     if (!is.null(gpu.to.choose))
#         gputools::chooseGpu(gpu.to.choose)    
#     rpca(M,lambda=lambda,mu=mu,term.delta=term.delta,max.iter=max.iter,trace=trace,
#          thresh.nuclear.fun=thresh.nuclear.gpu,thresh.l1.fun=thresh.l1.gpu,F2norm.fun=F2norm.gpu)
# }


