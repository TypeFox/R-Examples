chow.test <- function(y1,x1,y2,x2,x=NULL)
{
	mlr <- function(xy)
	{
		N <- nrow(xy)
		P <- ncol(xy)
		R <- cor(xy)
		b <- solve(R[-P,-P], R[,P][-P])
		variance <- var(xy)
		vr <- diag(variance)
		vr <- (vr/vr[P])[-P]
		b <- b/sqrt(vr)
		sse <- (variance[P, P]-(var(xy)[,P][-P])%*%b)*(N-1)
		sse
	}
        xy1<-cbind(x1,y1)
        xy2<-cbind(x2,y2)
	sse12 <- mlr(xy1)+mlr(xy2)
# in case of pooled x is known
        if(!is.null(x))
        {
          xy <- cbind(x,c(y1,y2))
          sse <- mlr(xy)
        }
        else 
	sse <- mlr(rbind(xy1, xy2))
        df1 <- ncol(xy1)
	df2 <- nrow(xy1)+nrow(xy2)-2*(df1)
	f <- (sse-sse12)*df2/(df1*sse12)
	p <- pf(f, df1, df2, lower.tail=FALSE)
	z <- c(f, df1, df2, p)
	names(z) <- c("F value", "d.f.1", "d.f.2", "P value")
	z
}
