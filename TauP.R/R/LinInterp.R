LinInterp <-
function(xin,yin,xout,mode='data')
{
yout=NULL
for(i in xout){

	if(i %in% xin){
		yout=c(yout,yin[xin==i])
	}else{
		f=switch(mode,jump='ordered',data=mean,all='ordered')
		yout=c(yout,approx(xin,yin,i,ties=f)$y)
	}
}
return(yout)

}

