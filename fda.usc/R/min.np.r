min.np<-function (fdataobj, h =NULL, W = NULL, Ker = Ker.norm,type.CV = GCV.S,
type.S=S.NW,par.CV=list(trim=0,draw=FALSE), verbose = FALSE,...)
{           
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")

x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names<-fdataobj[["names"]]
    nc <- nrow(fdataobj)
    np <- ncol(fdataobj)
   if (is.null(h)) {
   h<-seq(rtt[2]/80,rtt[2]/6,len=100)}
    lenh <- length(h)
    gcv <- array(NA, dim = c(lenh))
    df <- array(NA, dim = c(lenh))
    for (i in 1:lenh) {
        S2 <- type.S(tt, h[i], Ker)
        if (is.null(par.CV$trim)) par.CV$trim<-0
        gcv[i] <- type.CV(fdataobj,S=S2, W=W,trim=par.CV$trim,draw=par.CV$draw, ...)#####
        df[i]<-traza(S2)
    }
    l = which.min(gcv)
    h.opt <- h[l]
    gcv.opt <- gcv[l]
    S.opt <- type.S(tt, h.opt, Ker)
    fdata.est <- t(S.opt %*% t(x))

    if (lenh == 1) {dimnames(gcv)[1] <- list(h)   }
    else { dimnames(gcv)[[1]] <- as.numeric(round(h, 4))    }
        if (verbose){    cat("\n The minimum GCV (GCV.OPT=", round(gcv.opt, 4), sep = "",
        ") is achieved with\n the h value (h.opt=", round(h.opt,
            4), ")\n\n")
    if (h.opt< rtt[2]/80) print("Warning: h value too small")
    if (h.opt>rtt[2]/6)  print("Warning: h value too large")

if (lenh>1) {
  if (h.opt==min(h))  cat(" Warning: h.opt is the minimum bandwidth vaulue provided, range(h)=",range(h),"\n")
  else  if (h.opt==max(h))  cat(" Warning: h.opt is the maximum bandwidth vaulue provided, range(h)=",range(h),"\n")
}
}
fdata.est=fdata(fdata.est,tt,rtt,names)
output <- list(gcv = gcv,h=h,df=df,fdataobj = fdataobj, fdata.est = fdata.est,
gcv.opt = gcv.opt,h.opt = h.opt, S.opt = S.opt)
}

