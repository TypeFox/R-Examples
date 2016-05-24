mmpiter <- function(x,y,n,jlo,jhi,klo,khi,sig0,tau0,rr,ss,vv,ww){
xi <- matrix(0,10,4)
xi[1,3] <- sig0
xi[1,4] <- tau0

for (iter in 2:10){
  a <- mmp(x,y,n,jlo,jhi,klo,khi,xi[iter-1,3],xi[iter-1,4],rr,ss,vv,ww)
  xi[iter,1] <- a$jnew
  xi[iter,2] <- a$knew
  xi[iter,3] <- a$signew
  xi[iter,4] <- a$taunew
  diff1 <- xi[iter,1]-xi[iter-1,1]
  diff2 <- xi[iter,2]-xi[iter-1,2]
  if (abs(diff1) < .1 && abs(diff2) < .1) break
}
list(jhat=a$jnew,khat=a$knew,sigma2hat=a$signew,tau2hat=a$taunew)
}
