
freq_binom_two_bryantday_twostage=function(r0=0.2,r1=0.35,t0=0.3,t1=0.1,alpha.r,power,nrange,alpha.t=alpha.r){

  # toxicity is written as an event is no toxicity!
  t1=1-t1
  t0=1-t0

  output=matrix(0,1,14)
  imax=1000

  # generate a matrix of all designs which fit criteria
  # note this is slow since using data.frame - data.table would be much faster
  cat("\n Searching for solutions with sample size (final : interim)")
  for(n in nrange){
    cat("\n",n,": ")
    irange=floor(n/4):ceiling(n*3/4)

    for(n1 in irange){
      cat(n1,"")
      c1r_range=floor(r0*n1):floor(r1*n1)
      for(c1r in c1r_range){
        c2r_range=max(c1r,floor(r0*n)):ceiling(r1*n)
        for(c2r in c1r:n){

          ar0=sum(dbinom(c1r:n1,n1,r0)*(1-pbinom(c2r-c1r:n1-1,n-n1,r0)))
          ar1=sum(dbinom(c1r:n1,n1,r1)*(1-pbinom(c2r-c1r:n1-1,n-n1,r1)))

          c1t_range=floor(t0*n1):ceiling(t1*n1)
          for(c1t in c1t_range){
            c2t_range=max(c1t,floor(t0*n)):floor(t1*n)
            for(c2t in c2t_range){

              at1=sum(dbinom(c1t:n1,n1,t1)*(1-pbinom(c2t-c1t:n1-1,n-n1,t1)))
              at0=sum(dbinom(c1t:n1,n1,t0)*(1-pbinom(c2t-c1t:n1-1,n-n1,t0)))

              if(ar0*at1<alpha.r & ar1*at0<alpha.t & ar1*at1>power){

                samp00=n1+(n-n1)*sum(dbinom(c1r:n1,n1,r0))*sum(dbinom(c1t:n1,n1,t0))
                samp10=n1+(n-n1)*sum(dbinom(c1r:n1,n1,r1))*sum(dbinom(c1t:n1,n1,t0))
                samp01=n1+(n-n1)*sum(dbinom(c1r:n1,n1,r0))*sum(dbinom(c1t:n1,n1,t1))
                samp11=n1+(n-n1)*sum(dbinom(c1r:n1,n1,r1))*sum(dbinom(c1t:n1,n1,t1))

                a=dim(output)[1]
                if(n==output[a,4]){
                  if(max(samp00,samp10,samp01)<output[a,13]){
                    output[a,]=c(n1,c1r,c1t,n,c2r,c2t,ar0*at1,ar1*at0,ar1*at1,samp00,samp10,samp01,max(samp00,samp10,samp01),samp11)
                  }
                } else {
                  output=rbind(output,c(n1,c1r,c1t,n,c2r,c2t,ar0*at1,ar1*at0,ar1*at1,samp00,samp10,samp01,max(samp00,samp10,samp01),samp11))
                }
                imax=min(imax,max(samp00,samp10,samp01))
              }
            }
          }
        }
      }
    }
  }

  output=output[-1,]
  output=data.frame(output)
  names(output)=c("Interim","Iresponse","Itoxicity","final","response","toxicity","alpha.r","alpha.t","power","ssr0t1","ssr1t1","ssr0t0","ssanynull","ssalt")

  min=which.min(output[,13])

  # return optimal minmax and all possible designs.
  return(binom_two_bryantday(optimal=output[min,],minmax=output[1,],all.fit=output))
}


# ended
