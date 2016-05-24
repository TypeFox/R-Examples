cvpvs.logreg <- function(X,Y,tau.o=10, find.tau=FALSE, delta=2, tau.max=80, tau.min=1,
	pen.method=c("vectors","simple","none"),progress=TRUE)
{
	pen.method <- match.arg(pen.method)
	X <- as.matrix(X)
        n <- dim(X)[1]
	d <- dim(X)[2]
	Y <- factor(Y)
	Y <- unclass(Y)
	L <- max(Y)
        Ylevels <- levels(Y)
        # Stop if lengths of X[,1] and Y do not match
        if(length(Y) != length(X[,1])) {
          stop('length(Y) != length(X[,1])')
        }
	
	if (progress)
	{
		print('Computation of cross-validated p-values',
			quote=FALSE)
		print(paste('for ',as.character(n),
			' training observations.'),
			quote=FALSE)
		print('Preliminary log. regression:',
			quote=FALSE)
	}
	tmp <- penlogreg(X,Y,tau.o,
		pen.method=pen.method,progress=progress)
	a0 <- tmp$a
	b0 <- tmp$b
	
	PV <- matrix(1,nrow=n,ncol=L)
        tau.opt <- PV
	for (i in 1:n)
	{
		if (progress)
		{
			print(paste('Observation no. ',as.character(i),' ...'),
				quote=FALSE)
		}
		NewX <- X[i,]
		Xr <- X[(1:n)!=i,]
		Yr <- Y[(1:n)!=i]
                tmp <- pvs.logreg(NewX,Xr,Yr,tau.o=tau.o,
                                            find.tau=find.tau, delta=delta,
                                            tau.max=tau.max, tau.min=tau.min,
                                            pen.method=pen.method,a0=a0,b0=b0)
                PV[i,] <- tmp
                if(find.tau==TRUE)
                  {
                    tau.opt[i,] <- attributes(tmp)$tau.opt
                  }
              }
        if(find.tau==TRUE)
          {
            attributes(PV)$tau.opt <- tau.opt
          }
        dimnames(PV)[[2]] <- Ylevels
	return(PV)
}
