"sim.yo" <- 
function(x, coord=NULL, method="(2*a)/((2*a) + b + c)", dn=NULL, normalize = FALSE, listin = FALSE, listout = FALSE, ...)
{	
	if (listin) {
		x <- mama(x)
		x <- as.matrix(x)
		}
	x <- ifelse(x>0, 1, 0)
	df <- as.matrix(x)
	zeina <- row.names(df)
    anz <- nrow(df)
	a <- df %*% t(df) 
	b <- df %*% (1 - t(df)) 
	c <- (1 - df) %*% t(df) 
	d <- ncol(df) - a - b - c
	if (normalize) {
		an <- a/(a+b+c)
		bn <- b/(a+b+c)
		cn <- c/(a+b+c)
		a <- an
		b <- bn
		c <- cn
		}
	dis <- eval(parse(text = method))	
	dis <- as.dist(dis)
	attr(dis, "Size") <- anz
    attr(dis, "Labels") <- zeina
    attr(dis, "method") <- method
    attr(dis, "call") <- match.call()
    class(dis) <- "dist"
    if (listout) {
        dis <- liste(dis, entry=method)
        dis$a <- a[row(a) > col(a)]
	    dis$b <- b[row(b) > col(b)]
	    dis$c <- c[row(c) > col(c)]
	    dis$d <- c[row(d) > col(d)]
    	}
    if (!is.null(coord)){
	   xydist <- liste(dist(coord), entry="distance")
	   dis <- cbind(xydist, as.vector(dis))
	   names(dis)[4] <- method
	   X <- (outer(coord[,1], coord[,1], FUN="+"))*0.5
	   Y <- (outer(coord[,2], coord[,2], FUN="+"))*0.5	   
	   dis$X <- X[row(X) > col(X)]
	   dis$Y <- Y[row(Y) > col(Y)]
	   dis$xdist <- dist(coord[,1])
	   dis$ydist <- dist(coord[,2])
	   dis$a <- a[row(a) > col(a)]
	   dis$b <- b[row(b) > col(b)]
	   dis$c <- c[row(c) > col(c)]
	   dis$d <- c[row(d) > col(d)]
	   if (!is.null(dn)) {
	       if(length(dn)==1){
	           dis <- dis[(dis$distance <= dn), ]
	       }
	       else{
	           dis <- dis[((dis$distance >= min(dn)) & (dis$distance <= max(dn))), ]
	       }
	   }
    }
    return(dis)
}