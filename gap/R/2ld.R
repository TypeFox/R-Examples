# worked 28/6/03
# note that tables are symmetric do not fix, see kbyl below
LD22<-function(h,n) 
{
   D<-VarD<-Dmax<-VarDmax<-Dprime<-VarDprime<-x2<-lor<-vlor<-00

   z<-.C("tbyt",h=as.vector(h), haplotypes=as.double(n),
          D=as.double(D), VarD=as.double(VarD),
          Dmax=as.double(Dmax), VarDmax=as.double(VarDmax),
          Dprime=as.double(Dprime), VarDprime=as.double(VarDprime),
          x2=as.double(x2),lor=as.double(lor),vlor=as.double(vlor),PACKAGE="gap")

    invisible(list(h=h,n=n,D=z$D,VarD=z$VarD,
         Dmax=z$Dmax,VarDmax=z$VarDmax,Dprime=z$Dprime,
         VarDprime=z$VarDprime,x2=z$x2,lor=z$lor,vlor=z$vlor))
}

# refine on 17/4/2005
# verbose, default values, etc.
LDkl<-function(n1=2,n2=2,h,n,optrho=2,verbose=FALSE)
{
   Dp<-x2<-seX2<-rho<-seR<-klinfo<-0
   VarDp<-0
   Dijtable<-Dmaxtable<-Dijptable<-VarDijtable<-X2table<-VarDijptable<-matrix(rep(0,n1*n2),nrow=n1)

   z<-.C("kbyl",nalleles1=as.integer(n1), nalleles2=as.integer(n2),
          h=as.double(h), haplotypes=as.double(n),
          Dp=as.double(Dp),VarDp=as.double(VarDp),
          Dijtable=matrix(Dijtable,nrow=n1),
          VarDijtable=matrix(VarDijtable,nrow=n1),
          X2table=matrix(X2table,nrow=n1),
          Dmaxtable=matrix(Dmaxtable,nrow=n1),
          Dijptable=matrix(Dijptable,nrow=n1),
          VarDijptable=matrix(VarDijptable,nrow=n1),
          x2=as.double(x2), seX2=as.double(seX2),
          rho=as.double(rho), seR=as.double(seR), optrho=as.integer(optrho),
          klinfo=as.double(klinfo),verbose=as.integer(verbose),PACKAGE="gap")

   n1 <- z$nalleles1
   n2 <- z$nalleles2
   h <- t(matrix(z$h,nrow=n2))
   n <- z$haplotypes
   Dp <- z$Dp
   VarDp <- z$VarDp
   Dijtable <- t(matrix(z$Dijtable,nrow=n2))
   VarDijtable <- t(matrix(z$VarDijtable,nrow=n2))
   X2table <- t(matrix(z$X2table,nrow=n2))
   Dmaxtable <- t(matrix(z$Dmaxtable,nrow=n2))
   Dijptable <- t(matrix(z$Dijptable,nrow=n2))
   VarDijptable <- t(matrix(z$VarDijptable,nrow=n2))
   ptable <- 1-pchisq(X2table,1)
   x2 <- z$x2
   seX2 <- z$seX2
   rho <- z$rho
   seR <- z$seR
   optrho <- z$optrho
   klinfo <- z$klinfo
   df <- (n1-1)*(n2-1)
   if(verbose)
   {
      cat("\nEstimated haplotype frequencies\n\n")
      print(h)
      cat("\nTable of D\n\n")
      print(Dijtable)
      cat("\nTable of Dmax\n\n")
      print(Dmaxtable)
      cat("\nTable of D'\n\n")
      print(Dijptable)
      cat("\nTable of SE(D')\n\n")
      print(sqrt(VarDijptable))
      cat("\nTable of Chi-squares (based on D)\n\n")
      print(X2table)
      cat("\nTable of p values\n\n")
      print(ptable)
      cat("\nChi-squared statistic=",x2,"df=",df,"p=",1-pchisq(x2,df),"\n\n")
      cat("\nGlobal disequilibrium statistics and their standard errors\n\n")
      cat("D' coefficient=",Dp,'SE=',sqrt(VarDp),"\n\n")
      cat("Kullback-Leibler information",klinfo,"\n")
   }
   invisible(list(n1=n1, n2=n2, h=h, n=n,
   Dp=Dp,VarDp=VarDp,Dijtable=Dijtable, VarDijtable=VarDijtable, 
   Dmaxtable=Dmaxtable,
   Dijptable=Dijptable, VarDijptable=VarDijptable,
   X2table=X2table, ptable=ptable,
   x2=x2, seX2=seX2, rho=rho, seR=seR, optrho=optrho, klinfo=klinfo))
}

klem <- function(obs,k=2,l=2)
{
  if(length(obs)!=k*l*(k+1)*(l+1)/4) stop("incorrect length of genotype table")
  h <- rep(0,k*l)
  l0 <- l1 <- 0
  z <- .C("kbylem",obs=as.double(obs),nalleles1=as.integer(k),nalleles2=as.integer(l),
        Rh=as.double(h),l0=as.double(l0),l1=as.double(l1))
  invisible(list(h=z$Rh,l0=z$l0,l1=z$l1))
}

