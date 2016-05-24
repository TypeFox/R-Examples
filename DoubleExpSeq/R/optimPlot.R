optimPlot <-
function (y,m,groups,contrast=c(1,2),use.all.groups=TRUE,...)
                {
                  if( !use.all.groups )
                   {
                      cols1 <- groups == unique(groups)[contrast[1]]
                      cols2 <- groups == unique(groups)[contrast[2]]
                      wh.cols <- c(which(cols1),which(cols2))
                      y <- y[,wh.cols]
                      m <- m[,wh.cols]
                      groups <- as.character(groups[ wh.cols ])
                      groups <- as.factor(groups)
                      cols1 <- groups == unique(groups)[1]
                      cols2 <- groups == unique(groups)[2]
                   }

                   if( is.null(groups) ) groups = rep(1,ncol(y))
                   targets <- rownames(y)
                   groups <- as.factor(groups)
                   K = length(unique(groups))

                   z <- y/m
                   neff <- apply(z,1,FUN=function(vec) sum(vec!=1 & vec!=0,na.rm=TRUE))
                   neff[neff<=K] <- K+1

                   S <- .DBS( y=y , m=m , groups=groups )
                    names(S) <- as.character(neff)
                   a = (neff-K)/2
                   par = optim( par = 0.5 , fn = .nll.gbp.delta2 , a = a , S = S , method = "Brent" , lower = 0 , upper = 1 )$par

                   gammas <- seq( 0.01 , .99 , by = 0.01 )
                   vals <- sapply( gammas , .nll.gbp.delta2 , a=a , S=S )
                   plot( gammas , vals , xlim=c(0,1) , type="l" , lwd=1.25 , xlab="gamma" , ylab="Negative Log Likelihood" , ...)
                    abline( h=.nll.gbp.delta2( par , a=a , S=S ) , col = "red" , lty="dashed" , lwd = .75 )
                    abline( v=par , col = "red" , lty="dashed" , lwd = .75 )
                    points( par , .nll.gbp.delta2( par , a=a , S=S ) , pch=8 , col ="red" , cex=1.25 )
                }
