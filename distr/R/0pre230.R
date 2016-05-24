### for working under R < 2.3.0
if(getRversion()<'2.3.0') 
{ ## ignore ncp 

 ###beta Distribution
 
   qbeta <- function(p, shape1, shape2, ncp = 0, lower.tail = TRUE, 
                     log.p = FALSE)
           {if(isTRUE(all.equal(ncp,0)))
               stats::qbeta(p, shape1, shape2, lower.tail, log.p)
            else
              {x <- c(0.0,1.0)
               pfun <- function(x)
                       { pbeta(x, shape1=shape1, shape2=shape2, ncp=ncp)}
               qfun <- P2Q(pfun,x)
               p <- ifelse(log.p,exp(p),p)
               p <- ifelse(lower.tail,p,1-p)
               qfun(p)
              }   
            }
  
  rbeta <- function(n, shape1, shape2, ncp = 0)
           {if(isTRUE(all.equal(ncp,0)))
               stats::rbeta(n, shape1, shape2)
            else
               {X <- rchisq(n,df=2*shape1,ncp=ncp)
                Y <- rchisq(n,df=2*shape2,ncp=0)
                X/(X+Y)}
           }
 
 ###F Distribution
 
  qf    <- function(p, df1, df2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
           {if(isTRUE(all.equal(ncp,0)))
               stats::qf(p, df1, df2, lower.tail, log.p)
            else
              {TQ <- getdistrOption("TruncQuantile")
               xz <- qchisq(1-TQ,df=df1,ncp=ncp); xn<-qchisq(TQ,df=df2,ncp=0)
               x <- c(0,xz/xn*df2/df1)
               pfun <- function(x){pf(x, df1=df1, df2=df2, ncp=ncp)}
               qfun <- P2Q(pfun,x)
               p <- ifelse(log.p,exp(p),p)
               p <- ifelse(lower.tail,p,1-p)
               qfun(p)
              }   
           }
  rf    <- function(n, df1, df2, ncp = 0)
           {if(isTRUE(all.equal(ncp,0)))
                 stats::rf(n, df1, df2)
            else df2*rchisq(n, df=df1, ncp=ncp)/rchisq(n, df=df2, ncp=0)/df1
           }
 
  ###T Distribution
  
  qt    <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
           {if(isTRUE(all.equal(ncp,0)))
               stats::qt(p, df, lower.tail, log.p)
            else
              {TQ <- getdistrOption("TruncQuantile")*2
               xz <- qnorm(1-TQ,mean=df); xn<-sqrt(qchisq(TQ,df=df,ncp=0)/df)
               x <- c(-xz/xn,xz/xn)
               pfun <- function(x){pt(x, df=df, ncp=ncp)}
               qfun <- P2Q(pfun,x)
               p <- ifelse(log.p,exp(p),p)
               p <- ifelse(lower.tail,p,1-p)
               qfun(p)
              }   
           }
  rt    <- function(n, df, ncp = 0)
           {if(isTRUE(all.equal(ncp,0)))
                 stats::rt(n, df)
            else rnorm(n,mean=ncp)/sqrt(rchisq(n,df=df)/df)
           }
}
