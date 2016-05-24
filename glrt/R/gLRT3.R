gLRT3 <-
function(A, k=2, rho=0, gamma=0, EMstep=TRUE, ICMstep=TRUE, tol=1e-06, maxiter=1000, inf=Inf)
{
A[A[,2]==inf,2] = Inf
if(ncol(A) == 3 && all(A[,2] >= A[,1]) && length(unique(A[,3])) == k && all(A[,3]>=0) && all(A[,3]< k) )
{
AA = A[,-3] 
trt = A[,3]
est = ModifiedEMICM(AA, EMstep=EMstep, ICMstep=ICMstep, tol=tol, maxiter=maxiter)
tiny = .Machine$double.eps*100
est$sigma = ifelse(abs(est$sigma) < tiny, 0, est$sigma)
est$sigma = ifelse(abs(1.0 - est$sigma) < tiny, 1.0, est$sigma)
cens = CensorType(AA, inf=inf)
temp = cens
temp[temp != 4] = 5  
counts = table(c(0, 0, trt), c(5, 4, temp))  
counts[1,] = counts[1,] - 1  
u = Teststat3(trt, k, cens, counts, est, rho=rho, gamma=gamma, c0=1)
v = Var3(trt, k, cens, counts, est, rho=rho, gamma=gamma, c0=1)
chisq = Chisqstat3(u, v, counts)
chisqstat = chisq[1]
df = chisq[2]
p = 1-pchisq(chisqstat, df)
}
else
{
stop("Please Verify data format, # of samples, and treatment indicator!")
}
out = data.frame()
class(out) = "glrt3"
out$method = "Generalized log-rank test (Zhao, Zhao, Sun, and Kim, 2008)"
out$u = u
out$var = v
out$chisq = chisqstat
out$df = df
out$p = p
out
}

