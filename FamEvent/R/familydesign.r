familyDesign <-
function(n=1000, affectnum=0,  m.carrier= 0, dominant.m = TRUE, dominant.s = TRUE, depend=1, 
		  base.dist="Weibull", frailty.dist="gamma", 
          vbeta= c(-1.126, 2.55, 1.6), parms=c(0.016, 3), variation="none", 
          allelefreq=c(0.02, 0.2), mrate=0.1, age1=c(65,2.5), age2=c(45,2.5), 
          agemin=20) 
{
data<-numeric()
cumind<-0
i<- 1
j<- 0
while (i <= n) {
  j <- j + 1
  dat<-try(familyStructure(i,cumind=cumind, m.carrier=m.carrier, 
  		base.dist=base.dist, frailty.dist=frailty.dist,
	    depend=depend, parms=parms, vbeta=vbeta, dominant.m=dominant.m, dominant.s=dominant.s,
	    variation=variation, allelefreq=allelefreq, mrate=mrate, age1=age1, age2=age2,
	    agemin=agemin))
	    
  if(is.null(attr(dat, "class"))){
	   # At least one parent in first gen and two sibs in the second gen should be affected
	if(affectnum==3) until<- ifelse(sum(dat[dat[,7]==1,13]) >= 1 & sum(dat[dat[,7]==2,13]) > 1, TRUE, FALSE) 
	#[,7]=generation, [,13]=status
	else until <- TRUE
	
		if(!is.null(dim(dat))){
	    if(nrow(dat)>0 ){
	    	if(until){
		   	    data<-rbind(data, dat)
        		cumind<-cumind+nrow(dat)
        		i<-i+1
	   		}    
		  }
	  }
} # close "is.null(attr(dat, "class"))"
} # close while

data
}
