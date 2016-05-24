T_stationary <- function(sample, L = 49, J = 500, MC_rep=1000, cumulative_var = .90, Ker1 = FALSE, Ker2 = TRUE, h = ncol(sample)^.5, pivotal=FALSE)
{
    xrefine = N = ncol(sample)
    trefine = nrow(sample)
    if(Ker1)
    {
        K=function(x)
        {
            output = min(1, max(1.1-abs(x),0))
            return(output)
        }
    }
    if(Ker2)
    {
        K=function(x)
        {
            output = min(1, max(2-2*abs(x),0))
            return(output)
        }
    }
    basis = create.fourier.basis(c(0,1),L)
    ld = length(cumulative_var)
    h1 = h
    X1_bar = rowMeans(sample)
    mean_subtracted_X1 = sample - X1_bar
    gamma_hat = list()
    for(i in 0:(N-1))
    {
        temp_matrix = matrix(rep(0, trefine^2), trefine,trefine)
        for(j in (i+1):N)
        {
            temp_matrix = temp_matrix + (mean_subtracted_X1[,j] %*% t(mean_subtracted_X1[,j-i]))
        }
        gamma_hat[i+1]=list(temp_matrix/N)
    }
    cov_sample1 = gamma_hat[[1]]
    for(index in 1:(N-1))
    {
        cov_sample1 = cov_sample1 + K(index/h1) *(gamma_hat[[index+1]] + t(gamma_hat[[index+1]]))
    }
    Z_matrix = cov_sample1

    e1 = list()
    for(index in 1:L)
    {
        e1[index] = list(as.matrix(eval.basis(evalarg =(1:trefine)/trefine, basisobj = basis,Lfdobj=0)[,index]))
    }
    eigenvalues1 = (eigen(Z_matrix)$values)/trefine
    D = matrix(0,L,L)
    for(k in 1:L)
    {
        for(ell in 1:L)
        {
            Integrand_matrix = Z_matrix * (e1[[k]] %*% t(e1[[ell]]))
            D[k,ell] = 1/(trefine^2)*sum(Integrand_matrix)
        }
    }
    eigenpairs = eigen(D)
    eigenvectors = eigenpairs$vec
	eigenvectors2 = eigen(Z_matrix)$vectors
    eigenvalues = eigenpairs$val
    evals = eigenvalues
  	if(pivotal)
  	{
		d = c(1:ld)
        switch = 0
        stoper = 1
        spot = 1
        while(switch==0)
        {
            while((sum(eigenvalues[c(1:spot)])/sum(eigenvalues)) < cumulative_var[stoper])
            {
              spot = spot+1
            }
            d[stoper] = spot
            stoper = stoper+1
            if(stoper == (length(d)+1))
            {
              switch = 1
            }
        }
    	T_N0=1:ld
		for(r in 1:ld)
		{
        	ds=d[r]
        	inp.matrix=matrix(0,ds,N)
        	eig.v.norm=((trefine)^.5)*eigenvectors2
			for(j in (1:ds))
			{
				for(k in (1:N))
				{
					inp.matrix[j,k]=t(sample[,k])%*%(eig.v.norm[,j])/trefine
				}
			}
      	  	T_Nsum=rep(0,ds)
        	for(j in (1:ds))
        	{
			    s.0=sum(inp.matrix[j,(1:xrefine)])
			    for(x in (1:xrefine))
			    {
				    T_Nsum[j]=T_Nsum[j]+(1/xrefine)*((1/N^.5)*(sum(inp.matrix[j,(1:x)])-(x/xrefine)*s.0))^2
			    }
			}
          	T_N0[r]=sum(T_Nsum/eigenvalues[1:ds])
    	}
        T = vector(, MC_rep)
        T_array = matrix(0, length(d), MC_rep)
        lambda = eigenvalues
        for(dd in 1:length(d))
        {
            for(k in 1:MC_rep)
            {
                z=rnorm(d[dd]*J)
                tot=0
                for(n in c(1:d[dd]))
                {
                  sum1 = 0
                  sum1 =sum((z[c(((n-1)*d[dd]+1):((n-1)*d[dd]+J))]/(pi*c(1:J)))^2)
                  tot = tot+sum1
                }
                T_array[dd,k] = T[k] = tot
            }
        }
        p_values=(1:ld)
        for(dd in 1:length(d))
        {
            p_values[dd] = round(1 - ecdf(T_array[dd,])(T_N0[dd]),4)
        }
        cat("\n")
        cat("Pivotal test of stationarity for a functional time series\n")
        cat("\n")
        cat("null hypothesis: the series is stationary\n")
        cat("\n")

        cat(paste("p-value = ", p_values, sep=""),"\n")
        cat(paste("N (number of functions) = ", N, sep=""),"\n")
        cat(paste("number of MC replications = ", MC_rep, sep=""))
        cat("\n")
    }
  	if(pivotal==FALSE)
  	{
        T = vector(, MC_rep)
        T_array = (1: MC_rep)
        lambda = eigenvalues
	    d = min(c(length(which(lambda>0)),15))       
        for(k in 1:MC_rep)
        {
            z=rnorm(d*J)
            tot=0
            for(n in c(1:d))
            {
                sum1 = 0
                sum1 =sum((z[c(((n-1)*d+1):((n-1)*d+J))]/(pi*c(1:J)))^2)
                tot = tot+lambda[n]*sum1
            }
            T_array[k] = T[k] = tot
        }
    	int = sum(((1/sqrt(N)) *((sample[,1])-(1/xrefine)*rowSums(sample)))^2/(xrefine * trefine))
	    for(x in 2:xrefine)
	    {
    	    int = int + sum(((1/sqrt(N)) * (rowSums(sample[,1:x]) -(x/xrefine) * rowSums(sample)))^2/(xrefine * trefine))
    	}
	    statT_N = int
	    p_values=(1:ld)
    	for(dd in 1:length(d))
    	{
        	p_values = round(1 - ecdf(T_array)(statT_N),4)
	    }
	    cat("\n")
    	cat("Monte Carlo test of stationarity of a functional time series\n")
	    cat("\n")
	    cat("null hypothesis: the series is stationary\n")
    	cat("\n")

	    cat(paste("p-value = ", p_values, sep=""),"\n")
	    cat(paste("N (number of functions) = ", N, sep=""),"\n")
	    cat(paste("number of MC replications = ", MC_rep, sep=""))
    	cat("\n")    
    }
}
