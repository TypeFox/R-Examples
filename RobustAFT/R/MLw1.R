"MLw1" <-
function(yr,th0,v0,iv,n,tp,gamm,maxit,tol){
th <- th0; v <- v0; nit <- 1; dv <-0; delta <- 0; tn <- length(yr)
beta <- tp$beta; tn <- tp$tn
repeat {
  rr  <- yr-th; vo <- v
  # scale step
  if (iv==1) v <- Scalew(vo,rr,tn-1,beta,tol/10,maxit)
  dv <- v-vo
  # location step
  rs  <- wi <- rr/v; cnd <- rs!=0
  wi[cnd] <- tPsiw(rs[cnd],-Inf, Inf)/rs[cnd]
  delta   <- sum(wi[cnd]*rs[cnd])/sum(wi[cnd])
  th <- th+gamm*delta
  if (nit==maxit) cat("MLw1: nit=maxit","\n")
  if (nit==maxit | (abs(delta)<tol & abs(dv)<tol)) break
  nit <- nit+1}
list(th1=th,v1=v,nit=nit)}

