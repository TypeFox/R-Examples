"MLwp" <-
function(Xr,yr,th0,v0,iv,n,tp,gamm,maxit,tol,nitmon){
# Maximum likelihood
th   <- th0; v <- v0; nit <- 1; dv <-0; np <- ncol(Xr); delta <- rep(0,np); tn <- length(yr)
beta <- tp$beta; tn <- tp$tn
repeat {
  rr  <- as.vector(yr-Xr%*%th)
# scale step
  vo <- v 
  if (iv==1) v <- Scalew(vo,rr,tn-np,beta, tol/10,maxit)
  dv <- v-vo
# coefficient step
  rs  <- wi <- rr/v; cnd <- rs!= 0
  wi[cnd] <- tPsiw(rs[cnd],-Inf,Inf)/rs[cnd]
  sqw    <- sqrt(wi)
  rs     <- sqw*rs
  WX     <- sqw*Xr
  XW2X   <- t(WX)%*%WX
  XWr    <- t(WX)%*%rs
  delta  <- solve(XW2X,XWr)
  th     <- th+gamm*as.vector(delta)
  if (nitmon) cat("nit,v,th : ",nit,round(v,4),"\n",round(th,4),"\n")
  if (nit==maxit) cat("MLwp: nit=maxit","\n")
  if (nit==maxit | (all(abs(delta)<tol) & abs(dv)<tol)) break
  nit    <- nit+1}
list(th1=th,v1=v,nit=nit)}

