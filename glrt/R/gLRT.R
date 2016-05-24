gLRT <-
function(A, k=2, method=c("glrt1", "glrt2", "glrt3", "glrt4", "score"), M=50, rho=0, gamma=0, EMstep=TRUE, ICMstep=TRUE, tol=1e-06, maxiter=1000, inf=Inf)
{
if(method == "glrt1")
gLRT1(A, k, M, EMstep, ICMstep, tol, maxiter, inf)
else if(method == "glrt2")
gLRT2(A, k, rho, gamma, EMstep, ICMstep, tol, maxiter, inf)
else if(method == "glrt3")
gLRT3(A, k, rho, gamma, EMstep, ICMstep, tol, maxiter, inf)
else if(method == "glrt4")
gLRT4(A, k, rho, gamma, EMstep, ICMstep, tol, maxiter, inf)
else if(method == "score")
ScoreTest(A, k, EMstep, ICMstep, tol, maxiter, inf)
else
stop("Invalid method was specified!")
}

