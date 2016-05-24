"E.Smat.n" <-
function(XtX,xbar) {
#xk    <- 1.5477
#i0.11 <- integrate(Psi2phi.n, lower=-xk,upper=xk)$value
#i0.12 <- integrate(PsiChiphi.n, lower=-xk,upper=xk)$value
#i0.22 <- integrate(Chi2phi.n, lower=-25,upper=4)$value
i0.11 <- 0.5284873; i0.12 <- 0; i0.22 <- 0.1406566
np    <- length(xbar)
E.S <- matrix(0,ncol=np+1,nrow=np+1)
E.S[1:np,1:np] <- i0.11*XtX;  E.S[1:np,np+1] <- i0.12*xbar
E.S[np+1,1:np] <- i0.12*xbar; E.S[np+1,np+1] <- i0.22
E.S}

