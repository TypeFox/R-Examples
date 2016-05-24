
bvbox <- function(a, d = 7, mtitle = "Bivariate Boxplot",
 method = "robust",xlab="X",ylab="Y", add = FALSE, ...)
{
#
#a is data matrix
#d is constant(usually 7)
#
	p <- length(a[1,  ])
	if(method == "robust") {
		param <- biweight(a[, 1:2])
		m1 <- param[1]
		m2 <- param[2]
		s1 <- param[3]
		s2 <- param[4]
		r <- param[5]
	}
	else {
		m1 <- mean(a[, 1])
		m2 <- mean(a[, 2])
		s1 <- sqrt(var(a[, 1]))
		s2 <- sqrt(var(a[, 2]))
		r <- cor(a[, 1:2])[1, 2]
	}
	x <- (a[, 1] - m1)/s1
	y <- (a[, 2] - m2)/s2
	e <- sqrt((x * x + y * y - 2 * r * x * y)/(1 - r * r))
	e2 <- e * e
	em <- median(e)
	emax <- max(e[e2 < d * em * em])
	r1 <- em * sqrt((1 + r)/2)
	r2 <- em * sqrt((1 - r)/2)
	theta <- ((2 * pi)/360) * seq(0, 360, 3)
	xp <- m1 + (r1 * cos(theta) + r2 * sin(theta)) * s1
	yp <- m2 + (r1 * cos(theta) - r2 * sin(theta)) * s2
	r1 <- emax * sqrt((1 + r)/2)
	r2 <- emax * sqrt((1 - r)/2)
	theta <- ((2 * pi)/360) * seq(0, 360, 3)
	xpp <- m1 + (r1 * cos(theta) + r2 * sin(theta)) * s1
	ypp <- m2 + (r1 * cos(theta) - r2 * sin(theta)) * s2
	maxxl <- max(xpp)
	minxl <- min(xpp)
	maxyl <- max(ypp)
	minyl <- min(ypp)
	b1 <- (r * s2)/s1
	a1 <- m2 - b1 * m1
	y1 <- a1 + b1 * minxl
	y2 <- a1 + b1 * maxxl
	b2 <- (r * s1)/s2
	a2 <- m1 - b2 * m2
	x1 <- a2 + b2 * minyl
	x2 <- a2 + b2 * maxyl
	maxx <- max(c(a[, 1], xp, xpp, x1, x2))
	minx <- min(c(a[, 1], xp, xpp, x1, x2))
	maxy <- max(c(a[, 2], yp, ypp, y1, y2))
	miny <- min(c(a[, 2], yp, ypp, y1, y2))
        if (add) {
            points(a[,1], a[,2])
        } else {
	    plot(a[, 1], a[, 2], xlim = c(minx, maxx), ylim = c(miny, maxy), xlab =xlab, ylab =ylab, ...)
        }
	lines(xp, yp)
	lines(xpp, ypp, lty = 2)
	segments(minxl, y1, maxxl, y2, lty = 3)
	segments(x1, minyl, x2, maxyl, lty = 4)
}
#
biweight<-function(a,const1=9,const2=36,err=0.0001) {
#
#a is data matrix with two cols.
#const1=common tuning constant
#const2=bivariate tuning constant
#err=convergence criterion.
#
              x<-a[,1]
              y<-a[,2]
              n<-length(x)
              mx<-median(x)
              my<-median(y)
              madx<-median(abs(x-mx))
              mady<-median(abs(y-my))
              if(madx != 0) { ux<-(x-mx)/(const1*madx)
                              ux1<-ux[abs(ux)<1]
                              tx<-mx+(sum((x[abs(ux)<1]-mx)*(1-ux1*ux1)^2)/
                                     sum((1-ux1^2)^2))
                              sx<- sqrt(n)*sqrt(sum((x[abs(ux)<1]-mx)^2*
                                   (1-ux1*ux1)^4))/abs(sum((1-ux1*ux1)*
                                    (1-5*ux1*ux1)))
                             }
                  else { tx<-mx
                         sx<-sum(abs(x-mx))/n
                       }
              if(mady != 0) { uy<-(y-my)/(const1*mady)
                              uy1<-uy[abs(uy)<1]
                              ty<-my+(sum((y[abs(uy)<1]-my)*(1-uy1*uy1)^2)/
                                     sum((1-uy1^2)^2))
                              sy<- sqrt(n)*sqrt(sum((y[abs(uy)<1]-my)^2*
                                   (1-uy1*uy1)^4))/abs(sum((1-uy1*uy1)*
                                    (1-5*uy1*uy1)))
                             }
                  else { ty<-my
                         sy<-sum(abs(y-my))/n
                       }
               z1<-(y-ty)/sy+(x-tx)/sx
               z2<-(y-ty)/sy-(x-tx)/sx
              mz1<-median(z1)
              mz2<-median(z2)
              madz1<-median(abs(z1-mz1))
              madz2<-median(abs(z2-mz2))
              if(madz1 != 0) { uz1<-(z1-mz1)/(const1*madz1)
                              uz11<-uz1[abs(uz1)<1]
                           tz1<-mz1+(sum((z1[abs(uz1)<1]-mz1)*(1-uz11*uz11)^2)/
                                     sum((1-uz11^2)^2))
                              sz1<- sqrt(n)*sqrt(sum((z1[abs(uz1)<1]-mz1)^2*
                                   (1-uz11*uz11)^4))/abs(sum((1-uz11*uz11)*
                                    (1-5*uz11*uz11)))
                             }
                  else { tz1<-mz1
                         sz1<-sum(abs(z1-mz1))/n
                       }
              if(mady != 0) { uz2<-(z2-mz2)/(const1*madz2)
                              uz21<-uz2[abs(uz2)<1]
                          tz2<-mz2+(sum((z2[abs(uz2)<1]-mz2)*(1-uz21*uz21)^2)/
                                     sum((1-uz21^2)^2))
                              sz2<- sqrt(n)*sqrt(sum((z2[abs(uz2)<1]-mz2)^2*
                                   (1-uz21*uz21)^4))/abs(sum((1-uz21*uz21)*
                                    (1-5*uz21*uz21)))
                             }
                  else { tz2<-mz2
                         sz2<-sum(abs(z2-mz2))/n
                       }
              esq<-((z1-tz1)/sz1)^2+((z2-tz2)/sz2)^2
              w<-numeric(length=n)
              c2<-const2
              for(i in 1:10) {
              w[esq<const2]<-(1-esq[esq<const2]/const2)^2
              w[esq>=const2]<-0
              l<-length(w[w==0])
              if(l<0.5*n) break
                  else const2<-2*const2
                            }
               tx<-sum(w*x)/sum(w)
               sx<-sqrt(sum(w*(x-tx)^2)/sum(w))
               ty<-sum(w*y)/sum(w)
               sy<-sqrt(sum(w*(y-ty)^2)/sum(w))
               r<-sum(w*(x-tx)*(y-ty))/(sx*sy*sum(w))
               const2<-c2
               wold<-w
            for(i in 1:100) {
                     z1<-((y-ty)/sy+(x-tx)/sx)/sqrt(2*(1+r))
                     z2<-((y-ty)/sy-(x-tx)/sx)/sqrt(2*(1+r))
                     esq<-z1*z1+z2*z2
                     for(j in 1:10) {
                                    w[esq<const2]<-(1-esq[esq<const2]/const2)^2
                                    w[esq>=const2]<-0
                                    l<-length(w[w==0])
                                    if(l<0.5*n) break
                                         else const2<-2*const2
                                     }
               tx<-sum(w*x)/sum(w)
               sx<-sqrt(sum(w*(x-tx)^2)/sum(w))
               ty<-sum(w*y)/sum(w)
               sy<-sqrt(sum(w*(y-ty)^2)/sum(w))
               r<-sum(w*(x-tx)*(y-ty))/(sx*sy*sum(w))
               term<-sum((w-wold)^2)/(sum(w)/n)^2
               if(term<-err) break
                    else {wold<-w
                          const2<-c2
                         }
                             }
            param<-c(tx,ty,sx,sy,r)
            param
}


mahal <- function(x, index) {
    if (!is.matrix(x)) stop("x is not a matrix")
    xbar <- apply(x[index,], 2, mean)
    S <- var(x[index,])
    S <- solve(S)
    xcent <- t(t(x) - xbar)
    apply(xcent, 1, function(x) x %*% S %*% x)          
}
 
stalac <- function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    rn <- rownames(x)
    if (is.null(rn)) rn <- 1:nrow(x)
    n <- length(x[,1])
    p <- length(x[1,])
    s <- 1:n
    ind <- matrix(0, n-p, n)
    ind1 <- 0
    thresh <- qchisq((n-0.5)/n,p)

    index<-1:(p+1)
    for(i in (p+1):n) {
        ind1<-ind1+1
        if(i==(p+1)) D<-mahal(x,index)
        index<-order(D)
        index1<-sort(index[1:i])
        D<-mahal(x,index1) 
        index2<-s[D>thresh]
        ind[ind1,index2]<-ind[ind1,index2]+1
    }
    y<-rep(1:(n-p),rep(n,(n-p)))
    x<-rep(1:n,(n-p))
    par(mai = par("mai") * c(1, 1, 2, 1))
    plot(x, y, type = "n", axes = FALSE, xlab = "",
         ylab = "Number of observations used for estimation")
    axis(3, at = 1:n, labels = rn, las = 2)
    axis(2, at = 1:(n-p), labels = n:(p+1))
    txt <- c("", "*")[as.vector(t(ind[(n-p):1,])) + 1]
    text(x, y, labels = txt)
}

chiplot <- function(x, y, ...) {

    n <- length(x)
    ind <- numeric(n)

    for (i in 1:n) {
	for (j in (1:n)[-i])	
            if(x[i] > x[j] & y[i] > y[j]) ind[i] <- ind[i] + 1
    }
    ind <- ind/(n - 1)

    ind1 <- rank(x) - 1
    ind1 <- ind1/(n - 1)

    ind2 <- rank(y) - 1
    ind2 <- ind2/(n - 1)

    s <- sign((ind1 - 0.5) * (ind2 - 0.5))

    chi <- (ind - ind1  *ind2)/sqrt(ind1 * (1 - ind1) * ind2 *(1 - ind2))

    lambda <- 4 * s * pmax((ind1 - 0.5)^2, (ind2 - 0.5)^2)
    thresh <- 4 * (1/(n - 1) - 0.5)^2

    plot(lambda[abs(lambda) < thresh],
         chi[abs(lambda) < thresh], ylim = c(-1, 1), xlim = c(-1, 1),
         xlab = expression(lambda),
         ylab = expression(chi), ...)
    abline(h = c(-1, 1) * 1.78 / sqrt(n), lty = 2)
}
