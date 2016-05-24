gLRT4 <- 
function(A, k = 2, rho = 0, gamma = 0, EMstep = TRUE, ICMstep = TRUE, tol = 1e-06, maxiter = 1000, inf = Inf)
{
    A[A[,2]==inf,2] = Inf
    if (ncol(A) == 3 && all(A[, 2] > A[, 1]) && (k == 2) && all(A[, 3] >= 0) && all(A[, 3] < k)) 
    {
        AA = A[, -3]
        trt = A[, 3]
	
	AA0 = AA[trt == 0,]
	AA1 = AA[trt == 1,]
	
	n1 = sum(trt)
	n0 = length(trt) - n1

	if(n0 < 1 || n1 < 1)
	   stop("Not a two-sample comparision! So quit.")
	else 
	{
 	   c0 = 1
       est0 = ModifiedEMICM(AA0, EMstep = EMstep, ICMstep = ICMstep, tol = tol, maxiter = maxiter)
       est1 = ModifiedEMICM(AA1, EMstep = EMstep, ICMstep = ICMstep, tol = tol, maxiter = maxiter)
	   tiny = .Machine$double.eps*100
	   est0$sigma = ifelse(abs(est0$sigma) < tiny, 0, est0$sigma)
	   est0$sigma = ifelse(abs(1.0 - est0$sigma) < tiny, 1.0, est0$sigma)
	   est1$sigma = ifelse(abs(est1$sigma) < tiny, 0, est1$sigma)
	   est1$sigma = ifelse(abs(1.0 - est1$sigma) < tiny, 1.0, est1$sigma)
       cens0 = CensorType(AA0)
       cens1 = CensorType(AA1)
       u = Teststat4(AA0, n0, est0, cens0, AA1, n1, est1, cens1, rho, gamma, c0)
	   if(abs(u[2]) > tiny)  # u[2] != 0
		S0 = (u[1]/u[2])^2
	   if(S0 > 1)
               p = 2*(1 - pf(S0, 1, 1))
	   else
		p = 2*pf(S0, 1, 1)
	}
    }
    else 
    {
        if (any(A[, 1] == A[, 2])) 
            stop("Exact observations are not allowed. Please use method 'glrt1', 'glrt3', or 'score' instead!")
        else 
	    stop("Please verify data format, # of samples, and treatment indicator!")
    }
    out = data.frame()
    class(out) = "glrt4"
    out$method = "Generalized log-rank test (Zhao, Duan, Zhao, and Sun, 2013)"
    out$u = u
    out$v = NA
    out$fstat = S0
    out$df = c(1, 1)
    out$p = p
    out
}

