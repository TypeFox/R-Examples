help.func <- function(i,j,l,Dv,x.seq,l.seq,j.seq,K,N,Index.basis.D,base,q) {
 Yt <- cbind(rep(Dv[[l.seq[l]]][[j.seq[l-1]]]$u[i],N),x.seq) #2,1
 func <- cond.cop(data=Yt,coef=Dv[[l.seq[l-1]]][[j.seq[l-1]]]$v,K=K,diff="u1",Index.basis.D,base=base,q=q) #3,1
 return(spline(x=func,y=Yt[,2],xout=Dv[[j.seq[1]]][[l.seq[1]]]$u[i])$y)
}
