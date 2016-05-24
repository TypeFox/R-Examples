ppfcoeffab <-
function(ppf,nsamp,printflg,datastr){

     T=nrow(ppf)
     nc=ncol(ppf)
     maxn=nc-1
     nout= floor((T+1)/2)
     ppf[is.nan(ppf)]=0

      pkab=phth2ab(ppf)
      pkab=pkab$a
      pkab=as.matrix(pkab)

      r=nrow(pkab)
      c=ncol(pkab)

      pkpv=matrix(1,r,c)
      colsamp=matlab::sum(nsamp,FALSE)
      sigma0=1/sqrt(colsamp)
      pkpv[1,]=2*(1-pnorm(abs(pkab[1,]),0,sigma0))

     sigma=sqrt(2)*sigma0
       for (j in 2:T)
           {if (j==T & matlab::rem(T,2)==0)
              {sigma=sigma0}
               pkpv[j,]=2*(1-pnorm(abs(pkab[j,]),0,sigma))
           }
      if (printflg)
         {  for (k in 1:nc)
            {cat(paste('\n'))
            cat(paste('pi coeffs in ab form for',datastr, ',n= ',k-1,'\n'))
            detail <- matrix(c(pkab[,k],Re(pkpv[,k])),ncol=2)
                 colnames(detail) <- c("pihat_k", "pv")
                  row.names(detail)<-paste("k=",seq(0,T-1), sep="")
                 print(detail)
        }
      }

      result = list(pkab=pkab,pkpv=pkpv)
      class(result) = "ppfcoeffab"
      result
}
