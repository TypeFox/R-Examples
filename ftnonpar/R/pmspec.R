"pmspec"<-function(x,pks=0,alpha=0.9,sqzf=0.9,mult=0,lcl=0,ln=0,fig=0,pow=10^-2)  
{
        m<-length(x)
	kk<-1	
	while (2**kk <m ) {
			 kk<-kk+1
			 }
	n<-2**kk
	edf<-double(n)
	edf[1:n]<-0		
	edf[1:m]<-x-mean(x)
	edf<-fft(edf)
	edf<-Re(edf)**2+Im(edf)**2
	edf<-edf/(4*pi*m)	
	n<-n/2
	edf<-edf[1:n]
	if(n<=256){
		mult<-1
		}
	edf<-edf+pow*mean(edf)
	kk<-kk-1
#
#
	sqzf<-min(0.95,sqzf)	
#	calculate cutoff value
#	beta2<-1-(1-alpha)/(2*n)
#	beta1<-1-beta2
	if (mult > 0){
		kk<-n-1
		}
	qxtrm<-double(2*kk+2)
#
                ### fortran subroutine###
#	
    	tmp <- .Fortran(
		"npspcdn",
		as.double(edf),
		double(n+1),
		double(n+1),
		double(n+1),
		double(n),
		double(n),
		double(n+1),
		as.double(qxtrm),
		integer(n+1),
		integer(n+1),
		integer(n+1),
		as.integer(n),
		integer(1),
		as.double(sqzf),
		as.integer(pks),
		as.integer(kk),
		as.integer(mult),
		as.integer(lcl),
		as.double(alpha),
                PACKAGE="ftnonpar"
		)
#	
              ###spectral density###
#
	df<- tmp[[5]]
	edf<-edf
	sy<-tmp[[2]]
	ll<-tmp[[3]]
	uu<-tmp[[4]]
	nkn<-tmp[[12]]
	str<-double(2*nkn)
	dim(str)<-c(nkn,2)
	dm<-tmp[[9]]
	str[,1]<-dm[1:nkn]
	dm<-tmp[[7]]
	str[,2]<-dm[1:nkn]
	pks<-tmp[15]
#
              ###number of peaks##
#		
#              ###plot scale###
#
	if(fig==0){	
                xx<-0:(n-1)
		xx<-pi*xx/n
	if(ln==0) {
		  low<- min(log(edf[2:(n-1)]),log(df[2:(n-1)]))		
	          upp<- max(log(edf[2:(n-1)]),log(df[2:(n-1)]))
		  rng<-upp-low
		  low<-low-0.02*rng
		  upp<-upp+0.02*rng
	          plot(xx[2:(n-1)],log(edf[2:(n-1)]),ylim=range(low,upp),
	              xlab="Frequency ",ylab="Log(spectral density)",
			col=2) 
	          lines(xx[2:(n-1)],log(df[2:(n-1)]))
	          }
	else      {
	          upp<- 1.02*max(edf,df)
	          plot(xx,df,type="l",ylim=range(0,upp),
	              xlab="Frequency ",ylab="Spectral density") 	
	          points(xx,edf,col=2)
	          }
	}
#
	list(edf=edf,df=df,pks=pks,ll=ll,uu=uu,str=str,qxt=qxtrm)
	
#  
#
#
                       ###INPUTS###
#
# x	     = data
#
# alpha      = level for scaled observations
#
# sqzf	     = squeeze factor for tube
#
	               ###OUTPUTS###
#
# scl	     = scale
#
# pks	     = number of peaks
#							  	
}


