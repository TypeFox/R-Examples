"E.Smat" <-
function(XtX,xbar) {
#xk    <- 1.717817
#i0.11 <- integrate(Psi2phi.w, lower=-xk,upper=xk)$value
#i0.12 <- integrate(PsiChiphi.w, lower=-xk,upper=xk)$value
#i0.22 <- integrate(Chi2phi.w, lower=-25,upper=4)$value
i0.11 <- 0.4096672; i0.12 <- -0.01102064; i0.22 <- 0.1458827
np    <- length(xbar)
E.S <- matrix(0,ncol=np+1,nrow=np+1)
E.S[1:np,1:np] <- i0.11*XtX;  E.S[1:np,np+1] <- i0.12*xbar
E.S[np+1,1:np] <- i0.12*xbar; E.S[np+1,np+1] <- i0.22
E.S}

