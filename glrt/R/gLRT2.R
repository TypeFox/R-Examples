gLRT2 <-
function(A, k=2, rho=0, gamma=0, EMstep=TRUE, ICMstep=TRUE, tol=1e-06, maxiter=1000, inf=Inf)
{
A[A[,2]==inf,2] = Inf
if(ncol(A) == 3 && all(A[,2] > A[,1]) && length(unique(A[,3])) == k && all(A[,3]>=0) && all(A[,3]< k) )
{
AA = A[,-3] 
trt = A[,3]
est = ModifiedEMICM(AA, EMstep=EMstep, ICMstep=ICMstep, tol=tol, maxiter=maxiter)
tiny = .Machine$double.eps*100
est$sigma = ifelse(abs(est$sigma) < tiny, 0, est$sigma)
est$sigma = ifelse(abs(1.0 - est$sigma) < tiny, 1.0, est$sigma)
cens = CensorType(AA, inf=inf)
u = Teststat2(trt, k, cens, est, rho=rho, gamma=gamma, c0=1)
v = Var2(trt, k, cens, est, rho=rho, gamma=gamma, c0=1)
chisq = u[1:k-1] %*% solve(v[1:k-1, 1:k-1]) %*% u[1:k-1]
p = 1-pchisq(chisq, k-1)
}
else
{
if(any(A[,1] == A[,2]))
stop("Exact observations are not allowed. Please use method 'glrt3' instead!")
else
stop("Please Verify data format, # of samples, and treatment indicator!")
}
out = data.frame()
class(out) = "glrt2"
out$method = "Generalized log-rank test (Sun, Zhao, and Zhao, 2005)"
out$u = u
out$var = v
out$chisq = chisq
out$p = p

out
}

