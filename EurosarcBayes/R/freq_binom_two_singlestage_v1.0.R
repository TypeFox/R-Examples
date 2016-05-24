freq_binom_two_singlestage=function(r0,r1,t0,t1,power,alpha.r,nmax=100,alpha.t=alpha.r,nmin=1,adjust=TRUE){

  output=matrix(0,0,6)
  for(i in nmin:nmax){
    for(r in 1:i){
      alphar=1-pbinom(r-1,i,r0)
      powerr=1-pbinom(r-1,i,r1)
      for(t in 1:i){
        alphat=pbinom(t,i,t0)
        powert=pbinom(t,i,t1)
        if(adjust==TRUE){
          if(powerr*powert>power & alphar*powert<alpha.r & alphat*powerr<alpha.t){
            output=rbind(output,c(i,r,alphar,t, alphat,powerr*powert))
          }
        } else {
          if(powerr*powert>power & alphar<alpha.r & alphat<alpha.t){
            output=rbind(output,c(i,r,alphar,t, alphat,powerr*powert))
          }
        }
      }
    }
  }
  output=data.frame(output)
  names(output)=c("n","response","alphar","toxicity","alphat","power")

  return(binom_two_singlestage(optimal=output[1,],output=output))
}

# ended
