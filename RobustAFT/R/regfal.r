regfal <-
function(f,cc,lower,upper,nint=50,tol=0.001,maxit=20,...){
# to solve f(x,...)=cc in the interval [lower,upper]
as  <- seq(from=lower,to=upper,length=nint)
ck  <- NULL; nit <- 0; i <- 0
while(i < (nint-1)) {
 i  <- i+1
 fa <- f(as[i]  ,...)-cc
 fb <- f(as[i+1],...)-cc
# cat(i,as[i],as[i+1],fa,fb,"\n")
 if (fa*fb <= 0) break }
if (i==(nint-1) & fa*fb>0) {cat("no solution ", "i=", i, "\n"); return(list(solution=ck,nit=-1))}
ak <- as[i]
bk <- as[i+1]
fa <- f(ak,...)-cc
fb <- f(bk,...)-cc
nit <- 0; conv <- FALSE
while(!conv) {
 fak <- f(ak,...)-cc
 fbk <- f(bk,...)-cc
 ck  <- ak-(ak-bk)/(fak-fbk)*fak
 if (is.nan(ck)) {cat("solution non assigned \n"); return(list(solution=NULL,nit=nit))}
 fck <- f(ck,...)-cc
 if (fak*fck > 0) ak <- ck else bk <- ck
 conv <- max(abs(fck)) < tol
# cat(nit,ak,bk,ck,fak,fbk,fck,"\n")
 nit <- nit+1
 if (nit==maxit)  break}
return(list(solution=ck,nit=nit))}

