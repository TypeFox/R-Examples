ModifiedEMICM <-
function (A, EMstep = TRUE, ICMstep = TRUE, keepiter = FALSE, tol = 1e-06, maxiter = 1000) 
{
    if (ncol(A) == 2 && all(A[, 2] >= A[, 1])) 
    {
	temp = Aintmap(A[,1],A[,2]) 
	A = t(temp$A)
	intmap = temp$intmap
	startend = StartEnd(temp$A)
    }
    else
    { 
	stop("More than 2 columns or Li > Ri, so quit!")
    }
    temp <- ModifiedEMICMmac(A, EMstep = EMstep, ICMstep = ICMstep, keepiter = keepiter, tol = tol, maxiter = maxiter)
    if (is.null(temp)) 
    {
        return(NULL)
    }
    class(temp) <- "icsurv"
    temp$intmap <- intmap
    temp$ppairs <- startend
    temp$alpha <- A
    temp
}

