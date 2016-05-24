eq.2sp <- function(A,r){
# is determinant positive 
if(A[1,1]*A[2,2]>A[1,2]*A[2,1]) detpos <- TRUE else detpos <- FALSE
# is there a coexistence equilibrium
if(((-r[1]*A[2,2]+r[2]*A[1,2]) >0) 
&& ((-r[2]*A[1,1]+r[1]*A[2,1])>0)) coex <- TRUE else coex <- FALSE

if(coex==T){
 # find equilibrium
 xeq <- -solve(A,r)
 # community matrix
 C <- diag(xeq)%*%A
 # det, Tr, and Disc
 delta <- det(C); Tr<- sum(diag(C));D <- Tr^2-4*det(C)
 Det.Tra.Disc <- c(delta, Tr, D)
 # eigenvalues
 eval <- eigen(C)$values
 out <- list(detpos=detpos,coex=coex,xeq=round(xeq,3),
             Det.Tra.Disc=round(Det.Tra.Disc,3),eval=round(eval,3))
} else out <- list(detpos=detpos,coex=coex)

return(out)
}
