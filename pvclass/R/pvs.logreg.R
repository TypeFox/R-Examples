pvs.logreg <- function(NewX, X, Y, tau.o=10, find.tau=FALSE, delta=2, tau.max=80, tau.min=1,
	a0=NULL, b0=NULL,
	pen.method=c("vectors","simple","none"),
	progress=FALSE)
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
	
	if (is.null(a0) || is.null(b0))
	{
		if (progress)
		{
			print('Preliminary log. regression:',quote=FALSE)
		}
		tmp <- penlogreg(X,Y,tau.o,
			pen.method=pen.method,progress=progress)
		a0 <- tmp$a
		b0 <- tmp$b
	}

        NewX <- as.matrix(NewX)
  
        if(d > 1 & NCOL(NewX) == 1) {
          NewX <- t(NewX)
        }
  
        if(d == 1 & NCOL(NewX) > 1) {
          NewX <- t(NewX)
        }
  
        nr <- NROW(NewX)
        s <- NCOL(NewX)

        # Stop if dimensions of NewX[i,] and X[j,] do not match
        if(s != d) {
          stop('dimensions of NewX[i,] and X[j,] do not match!')
        } 
        
	PV <- matrix(1,nrow=dim(NewX)[1],ncol=L)
	Xa <- rbind(rep(0,d),X)
	Ya <- c(0,Y)

        if(find.tau==TRUE) {
          tau.opt <- PV
          for (i in 1:dim(NewX)[1])
            {
              Xa[1,] <- NewX[i,]
              for (theta in 1:L)
		{
                  if (progress)
                    {
                      txt <- paste('Observation no. ',as.character(i),
                                   ', class ',as.character(theta),' ...')
                      print(txt,quote=FALSE)
                    }
                  Ya[1] <- theta
                  
                  tau.opt[i,theta] <- find.tau.opt(Xa, Ya, theta,
                                        tau.o, delta, tau.max, tau.min,
                                        pen.method, a0, b0)   
                  tmp <- penlogreg(Xa,Ya,tau.o=tau.opt[i,theta],
                                             a0=a0,b0=b0)$PM[Ya==theta,theta]
                  PV[i,theta] <- mean(tmp <= tmp[1])
		}
              attributes(PV)$tau.opt <- tau.opt
            }
        } else {
          for (i in 1:dim(NewX)[1])
            {
              Xa[1,] <- NewX[i,]
              for (theta in 1:L)
		{
                  if (progress)
                    {
                      txt <- paste('Observation no. ',as.character(i),
                                   ', class ',as.character(theta),' ...')
                      print(txt,quote=FALSE)
                    }
                  Ya[1] <- theta
                  tmp <- penlogreg(Xa,Ya,tau.o,
                                             a0=a0,b0=b0)$PM[Ya==theta,theta]
                  PV[i,theta] <- mean(tmp <= tmp[1])
		}
            }
        }
        dimnames(PV)[[2]] <- Ylevels
	return(PV)
}
