madsim <-
function(mdata=NULL,n=10000,ratio=0,
         fparams=data.frame(m1=7,m2=7,shape2=4,lb=4,ub=14,pde=0.02,sym=0.5),
         dparams=data.frame(lambda1=0.13,lambda2=2,muminde=1.0,sdde=0.5),
         sdn=0.4, rseed=50) {
    set.seed(rseed);
    m1 <- fparams$m1;
    m2 <- fparams$m2;
    m <- m1+m2+1;  #first column will be reference if ratio=1
    n1 <- length(mdata);
    if (n1>100) {
        if (n == 0) {
           n <- n1;
           x2 <- mdata;
        } else {
                  if (n < n1) {
                     x2 <- sample(mdata,n,replace=FALSE);
                  } else {
                     x2 <- sample(mdata,n,replace=TRUE);
                  }
           }
    } else {
           x <- rbeta(n,2,fparams$shape2);
           x2 <- fparams$lb+(fparams$ub*x);
    }

    xdat <- matrix(c(rep(0,n*m)),ncol=m);
    xid <- matrix(c(rep(0,n)),ncol=1);
    for (i in 1:n) {
        alpha <- dparams$lambda1*exp(-dparams$lambda1*x2[i]);
        xi_val <- runif(m,min=(1-alpha)*x2[i],max=(1+alpha)*x2[i]);
        if (sample(1:100,1)>(100*fparams$pde)) { # case of non DE gene
           xdat[i,] <- xi_val;
        } else { # case of DE gene
               xi1 <- xi_val[1:(m1+1)];
               mude <- dparams$muminde+rexp(1,dparams$lambda2);
               if (sample(1:100,1)>(100*fparams$sym)) { # up regulated gene
                  xi2 <- xi_val[(m1+2):m]+rnorm(m2,mean=mude,sd=dparams$sdde);
                  xid[i] <- 1;
               } else { # down regulated gene
                  xi2 <- xi_val[(m1+2):m]-rnorm(m2,mean=mude,sd=dparams$sdde);
                  xid[i] <- -1;
               }
               xdat[i,] <- c(xi1,xi2);
        }
    }
    xsd <- sd(xdat[,1]);
    if (sdn>0) {
       ndat <- matrix(c(rnorm(n*m,mean=0,sd=sdn)),ncol=m);
       xdat <- xdat+ndat;
    }
    xdata <- matrix(c(rep(0,n*(m-1))),ncol=(m-1));
    if (ratio) {  # ratio expression values
       xdata <- xdat[,2:m]-xdat[,1];
    } else {    # absolute expression values
           xdata <- xdat[,2:m];
    }
    list(xdata=xdata, xid=xid, xsd=xsd)
}
