## to be exported: berechnet Konfidenzbänder, simultan und punktweise
qqbounds <- function(x,D,alpha,n,withConf.pw, withConf.sim,
                     exact.sCI=(n<100),exact.pCI=(n<100),
                     nosym.pCI = FALSE, debug = FALSE){
   x <- sort(unique(x))
   if("gaps" %in% names(getSlots(class(D))))
       {if(!is.null(gaps(D)))
            x <- sort(unique(c(x, gaps(D))))
       }
   c.c <- matrix(NA,nrow=length(x),ncol=4)
   colnames(c.c) <- c("sim.left","sim.right","pw.left","pw.right")

   SI <- .SingleDiscrete(x,D)
   SI.in <- SI<4
   SIi <- SI[SI.in]
   x.in <- x[SI.in]
   p.r <- p(D)(x.in)
   p.l <- p.l(D)(x.in)
   l.x <- length(x.in)
   if(debug){
     cat("the partition into discrete mass points and cont. points gives\n")
     print(SI)
     cat("x.in is\n")
     print(x.in)
     cat("number of mass points:\n")
     print(sum(SI.in))
     cat("p.r and p.l\n")
     print(cbind(p.r,p.l))
     cat("length of l.x\n")
     print(l.x)
     print(c(alpha=alpha,n=n,exact.sCI=exact.sCI))
   }

   silent0 <- !debug

   c.crit <- if(withConf.sim) try(.q2kolmogorov(alpha,n,exact.sCI, silent0),
                                   silent=silent0) else NULL
   c.crit.i <- if(withConf.pw) try(.q2pw(x.in,p.r,D,n,alpha,exact.pCI,nosym.pCI,
                                         silent0),silent=silent0) else NULL
   #print(cbind(c.crit,c.crit.i))
   if(debug){
      cat("returned c.crit is\n")
      print(str(c.crit))
      cat("returned c.crit.i is\n")
      print(str(c.crit.i))
   }
   te.i <- withConf.pw  & !is(c.crit.i,"try-error")
   te.s <- withConf.sim & !is(c.crit,  "try-error")

   if(te.s){
      c.crit.r <- q.r(D)(pmax(1-p.r-c.crit/sqrt(n),
                         # alternative: pmax(1-(1:l.x)/l.x-c.crit/sqrt(n),
                         getdistrOption("DistrResolution")),lower.tail=FALSE)
      c.crit.l <- q(D)(pmax(p.l-c.crit/sqrt(n),
                       # alternative: pmax(((1:l.x)-1)/l.x-c.crit/sqrt(n),
                       getdistrOption("DistrResolution")))
      c.crit.l[abs(c.crit.l)==Inf] <- NA
      c.crit.r[abs(c.crit.r)==Inf] <- NA
      c.crit.l[SIi == 2 | SIi == 3] <- NA
      c.crit.r[SIi == 2 | SIi == 3] <- NA
      c.c[SI.in,1:2] <- cbind(c.crit.l,c.crit.r)
   }
   if(te.i){
      c.crit.i <- x.in + c.crit.i/sqrt(n)
      c.crit.i[SIi == 2 | SIi == 3] <- NA
      c.c[SI.in,3:4] <- c.crit.i
      c.c[SI.in & abs(c.crit.i[,1])==Inf,3] <- NA
      c.c[SI.in & abs(c.crit.i[,2])==Inf,4] <- NA
   }
   return(list(crit = c.c, err=c(sim=te.s,pw=te.i)))
}
# returnlevelplot(xex,datax=FALSE,GEVFamilyMuUnknown(loc=es[1],shape=es[3],scale=es[2]))
