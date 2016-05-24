`protein.con` <-
function (D0,D,c,a,d.D,d.c, d.a, data.dilutes,r=1.2,minimal.err=5) {
x.weighted.mean= rep(NA,nrow(data.dilutes))
x.err = x.weighted.mean
xflag = x.weighted.mean		# takes values of 0,1,2, which means under detection, OK, saturated
K = ncol(data.dilutes) 		# number of total dilution steps for a sample
igamma = log(D0)/log(D)	#where gamma is 1/gamma, a parameter in Sips model
M =min(1e9,1/c+a)			#when M is too large, take 1e9.
x.saturation.level=   D0^(K-1)/((1/( M/r - a)- 1/(M-a)))^igamma 
x.nodetection.level = D0^(1-1)/((1/( r*a - a)- 1/(M-a)))^igamma
for (Np in 1:nrow(data.dilutes)){ # for each sample
   x=rep(NA,K); w=x; xL=x; xH=x;  #initialization
   y = data.dilutes[Np,]
   if((y[K]> M/r) && length(y[y<M/r])<2) {#condition to call saturation
		xflag[Np] = 2 
		x.weighted.mean[Np] = x.saturation.level # Use M/r value
		x.err[Np] = NA
   } else {
   if((y[1]<r*a) & length(y[y>r*a])<2) {#condition to call undetected
		xflag[Np] = 0
		x.weighted.mean[Np] = x.nodetection.level # Use r*a value
		x.err[Np] = NA
   } else {
   y[y>M/r] = NA # for removing signals near saturation 
   y[y<a*r] = NA # for removing signals near bg noise
   
   for (k in 1:K){# for each signal in a dilution series
      y[k] =max(min(M/1.01,y[k]), a+minimal.err) # limit y[k] to be within a+minimal.err and M/1.01
     	x[k] =   D0^(k-1) /(1/(y[k]-a)- c)^igamma #estimated protein concentration prior dilution
	de.x.over.de.a = igamma * D0^(k-1)*(1/(y[k]-a)- c)^(-igamma-1)/(y[k]-a)^2
      de.x.over.de.c = igamma * D0^(k-1)*(1/(y[k]-a)- c)^(-igamma-1)
      de.x.over.de.D = x[k] *log(1/(y[k]-a)- c) * igamma/D/log(D)/D0^(k-1)
	w[k] = (de.x.over.de.a * d.a)^2 + ( de.x.over.de.c * d.c)^2 + (de.x.over.de.D * d.D)^2
	}
   w = w[!is.na(x)] # removing signals near saturation or bg noise
   x = x[!is.na(x)] # removing signals near saturation or bg noise
   if(length(x) > 0 ) {
	x.range = 3* max(1, median(x)*0.01,mad(x)) 
	x.f = (abs(x-median(x)) < x.range) # removing outliers
	x=x[x.f]
	w=w[x.f] # removing outliers
	w= 1/w
	x.weighted.mean[Np] = sum (x*w) /sum(w)
	x.err[Np]=1/sqrt(sum(w))
	}
  }#end of else saturation
  }# end of else below detection 
}#end of for each sample Np
ret  <-  cbind(x.weighted.mean, x.err,xflag)
rownames(ret) <- rownames(data.dilutes)
return(ret)
}#end of function

