# Calculate Goodman-Kruskal Gamma
# Original from: Laura Thompson,
# https://home.comcast.net/~lthompson221/Splusdiscrete2.pdf

GKgamma<-function(x, level=0.95)
{
    # x is a matrix of counts.  You can use output of crosstabs or xtabs in R.    
    # Confidence interval calculation and output from Greg Rodd
    
    # Check for using S-PLUS and output is from crosstabs (needs >= S-PLUS 6.0)
    if(is.null(version$language) && inherits(x, "crosstabs")) { oldClass(x)<-NULL; attr(x, "marginals")<-NULL}
        
	## TODO: add tests for matrix or table
	
    n <- nrow(x)
    m <- ncol(x)
    pi.c<-pi.d<-matrix(0, nrow=n, ncol=m)
        
    row.x<-row(x)
    col.x<-col(x)
    
    for(i in 1:(n)){
        for(j in 1:(m)){
            pi.c[i, j]<-sum(x[row.x<i & col.x<j]) + sum(x[row.x>i & col.x>j])
            pi.d[i, j]<-sum(x[row.x<i & col.x>j]) + sum(x[row.x>i & col.x<j])
        }
    }

    C <- sum(pi.c*x)/2
    D <- sum(pi.d*x)/2
    
    psi<-2*(D*pi.c-C*pi.d)/(C+D)^2
    sigma2<-sum(x*psi^2)-sum(x*psi)^2
    
    gamma <- (C - D)/(C + D)
    pr2 <- 1 - (1 - level)/2
    CI <- qnorm(pr2) * sqrt(sigma2) * c(-1, 1) + gamma

    result <- list(gamma = gamma, C = C, D = D, sigma = sqrt(sigma2), 
        CIlevel = level, 
        CI = c(max(CI[1], -1), min(CI[2], 1))
        )  
    class(result) <- "GKgamma"
    result   
}

print.GKgamma <-
function (x, digits = 3, ...) 
{
	cat("gamma        :", round(x$gamma, digits = digits), 
			"\n")
	cat("std. error   :", round(x$sigma, digits = digits), 
			"\n")
	cat("CI           :", round(x$CI, digits = digits), 
			"\n")
	invisible(x)
}