fdata2fd=function(fdataobj,type.basis=NULL,nbasis=NULL,nderiv=0,lambda=NULL,...) {
if (is.fdata(fdataobj)) DATA=fdataobj[["data"]]
else stop("No fdata object")
np=ncol(DATA)
tt =fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
if (is.null(lambda)) lambda=3e-08/diff(rtt)
if (is.null(nbasis)) {
       nbasis=ifelse(floor(np/3) > floor((np - nderiv - 4)/2),
       floor((np - nderiv - 4)/2), floor(np/3))
       }
   as <- list()
   as[[1]] <-rtt
   names(as)[[1]]<-"rangeval"
   as[[2]] <- nbasis
   names(as)[[2]]<-"nbasis"
   C <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("DATA","type.basis","nbasis","nderiv"),names(mf),0L)
   imetri <- m[2]
   if (imetri == 0) {
        type.basis1="bspline"
        a1 <- create.bspline.basis
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric))    }
   else {  a1 <- paste('create.',type.basis,'.basis',sep="")
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric)) }
  ii <- imetri + 1
  if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metric) {
            aa <- any(names(C) == names(formals(a1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
                ii <- ii + 1
                as[[ind.m]] <- C[[vv[ind.m]]]
                names(as)[[ind.m]]<-names(formals(a1)[ind.m])            }
#            else {                 as[[ind.m]] <- formals(a1)[[ind.m]]   }
            ind.m <- ind.m + 1            }
  }
 b1.1<- do.call(a1, as)
 class(DATA)="matrix"

 fd1.1 <- Data2fd(argvals=tt,y=t(DATA),basisobj=b1.1,lambda=lambda,...) ######
 if (nderiv>0) fd1.1=deriv.fd(fd1.1,int2Lfd(nderiv)) #######
 fd1.1
}

