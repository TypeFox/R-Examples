#' @keywords internal
llspps<-function(counts, sizes, s){
  #note counts and sizes refer to the counts of each size in the full population.  
  # s is the sample.
  if(min(s)<1)print("Error - s values must be positive")
  totalsize<-sum(counts*sizes)
  obs_size<-cumsum(s) #cumulative sizes observed
  obs_size2<-c(0,obs_size)[1:length(obs_size)]
  obscounts<-tabulate(s)[sizes]
  result<-sum(lfactorial(counts)-lfactorial(counts-obscounts))+sum(log(s/(totalsize-obs_size2)))
  result
}


#' @keywords internal
lllspps<-function(lcountdif,sizes,s){
  #couputes the log likelihood.  This is used to feed to optim.  The first argument is the log of the 
  # difference between the population counts and the sample counts.  This formulation forces this difference 
  # to be non-negative.
  countdif<-exp(lcountdif)
  counts<-countdif+tabulate(s)[sizes]
  llspps(counts,sizes,s)  
}

#' @keywords internal
sppssample<-function(counts,sizes,n){
   rstuff<-counts
   ss<-rep(0,n)
   for(i in 1:n){
     ss[i] <- sample(sizes, size=1,prob=rstuff*(sizes))
     rstuff[which(sizes==ss[i])] <- rstuff[which(sizes==ss[i])] - 1 
   }
   ss
}

#' @keywords internal
sppssamplei<-function(pop,n){
  sample(x=(1:length(pop)),size=n,prob=pop,replace=FALSE)
}

#' @keywords internal
roundstoc<-function(vec){# takes a vector and makes it integers, keeping the total the same
	#vec<-c(1.2,1.3,2.5)
	target<-sum(vec)
	temp<-floor(vec)
	needed<-target-sum(temp)
	away<-vec-temp
	while(needed>.5){
		toget<-sample(c(1:length(vec)),size=1,prob=away)
		temp[toget]<-temp[toget]+1
		away<-vec-temp
		away[away<0]<-0
		needed<-needed-1
		}
	temp
	}

#' @keywords internal
Zee<-function(lam,n,nk,sizes,D){
(1/n)*
sum(1/(sum((nk*sizes/n)/(1-exp(-lam*sizes)))-(D/n)))-lam
}


#' @keywords internal
bnwnest<-function(lst,toplam=10,bottomlam=1e-10){
  s<-lst
  sizes<-sort(unique(s))
  nk<-tabulate(s)[tabulate(s)>0]
  D2<-cumsum(s)
  D<-c(0,D2[-length(D2)])
  n<-length(s)
  if(sign(Zee(bottomlam,n=n,nk=nk,sizes=sizes,D=D))!=
    sign(Zee(toplam,n=n,nk=nk,sizes=sizes,D=D))){
  lest<-uniroot(Zee,n=n,nk=nk,sizes=sizes,D=D,
    interval=c(bottomlam,toplam))$root
  bbnest<-nk/(1-exp(-lest*sizes))
  out<-bbnest
  }else{out<-rep(0,length(sizes))}
  out
}

#' @keywords internal
mleNest<-function(lst,start=NULL){
  if(is.null(start)){start<-bnwnest(lst)}
  s<-lst
  sizes<-sort(unique(s))
  nk<-tabulate(s)[tabulate(s)>0]
  D2<-cumsum(s)
  D<-c(0,D2[-length(D2)])
  n<-length(s)
  if(sum(start)>0){
  xstart<-log(start-tabulate(s)[sizes])
  }else{xstart<-log(tabulate(s)[sizes])}
  out<-optim(par=xstart,
    fn=lllspps,
    sizes=sizes, s=s,
    control=list(maxit=100000,fnscale=-1))
  Nest<-exp(out$par)+tabulate(s)[sizes]
  list(Nk=Nest,conv=out$conv)
}


#' @keywords internal
bins<-function(vec,breaks,maxval, below=TRUE){
#returns the vec re-coded to assign the values of vec to bins
#defined by the values of breaks, with ranges of [ ) structure
# unless below=FALSE.  New values are bin means.
  newvec<-vec
  mini<-min(c(0,vec))
  maxi<-breaks[1]
  midi<-(mini+maxi)/2
  newvec[vec<maxi]<-midi
  for(i in 2:length(breaks)){
      mini<-breaks[i-1]
      maxi<-breaks[i]
      midi<-(mini+maxi)/2    
      if(below){newvec[vec>=mini & vec<maxi]<-midi
       }else{newvec[vec>mini & vec<=maxi]<-midi}
  }
  if(below){newvec[vec>=max(breaks)]<-maxval
    }else{newvec[vec>max(breaks)]<-maxval}
  newvec
}




