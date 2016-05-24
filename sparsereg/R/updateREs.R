#Update Random Effects
	
	
	updateREs<-function(y, fits.main,ran.effect1,ran.effect2,ran.effect3,id,id2,id3,
	sigma.sq.b,sigma.sq.b2,sigma.sq.b3,sigma.sq,fix.eff=FALSE,EM0=FALSE){

	var.zero<-1
	if(EM0==TRUE) var.zero<-0
	
	if(length(id)>0){
	v.curr<-as.vector(y-as.vector(fits.main)-ran.effect2-ran.effect3)
	mean.vec<-tapply(v.curr,id,mean)
	length.vec<-tapply(y,id,length)
	ran.prec<-1/(length.vec/sigma.sq+1/sigma.sq.b)
	ran.effect.mean<-mean.vec*ran.prec*length.vec/sigma.sq
	ran.effect.samp<-rnorm(length(ran.effect.mean),mean=ran.effect.mean,sd=1/ran.prec^.5*var.zero)
	
	ran.effect<-mean.vec[id]
	ran.effect<-ran.effect-mean(ran.effect)
	sigma.sq.b.shape<-(mean(length.vec)-1)/2
	sigma.sq.b.scale<-sum(ran.effect.samp^2)/2
	sigma.sq.b<-rinvgamma(1,shape=sigma.sq.b.shape,scale=sigma.sq.b.scale)	
	if(EM0) sigma.sq.b<-sigma.sq.b.scale/(max(c(sigma.sq.b.shape-1,1)))
	if(fix.eff) sigma.sq.b<-1e5
	} 

	if(length(id2)>0){
	v.curr2<-as.vector(y-as.vector(fits.main)-ran.effect-ran.effect3)
	mean.vec2<-tapply(v.curr2,id2,mean)
	length.vec2<-tapply(y,id2,length)
	ran.prec2<-1/(length.vec2/sigma.sq+1/sigma.sq.b2)
	ran.effect.mean2<-mean.vec2*ran.prec2*length.vec2/sigma.sq
	ran.effect.samp2<-rnorm(length(ran.effect.mean2),mean=ran.effect.mean2,sd=1/ran.prec2^.5*var.zero)
	
	ran.effect2<-mean.vec2[id2]
	ran.effect2<-ran.effect2-mean(ran.effect2)

	sigma.sq.b.shape<-(mean(length.vec2)-1)/2
	sigma.sq.b.scale<-sum(ran.effect.samp2^2)/2
	sigma.sq.b2<-rinvgamma(1,shape=sigma.sq.b.shape,scale=sigma.sq.b.scale)
	if(EM0) sigma.sq.b2<-sigma.sq.b.scale/(max(c(sigma.sq.b.shape-1,1)))

	if(fix.eff) sigma.sq.b2<-1e5
	} 

	if(length(id3)>0){
	v.curr3<-as.vector(y-as.vector(fits.main)-ran.effect-ran.effect2)
	mean.vec3<-tapply(v.curr3,id3,mean)
	length.vec3<-tapply(y,id3,length)
	ran.prec3<-1/(length.vec3/sigma.sq+1/sigma.sq.b3)
	ran.effect.mean3<-mean.vec3*ran.prec3*length.vec3/sigma.sq
	ran.effect.samp3<-rnorm(length(ran.effect.mean3),mean=ran.effect.mean3,sd=1/ran.prec3^.5*var.zero)
	
	ran.effect3<-mean.vec3[id3]
	ran.effect3<-ran.effect3-mean(ran.effect3)
	sigma.sq.b.shape<-(mean(length.vec3)-1)/2
	sigma.sq.b.scale<-sum(ran.effect.samp3^2)/2
	sigma.sq.b3<-rinvgamma(1,shape=sigma.sq.b.shape,scale=sigma.sq.b.scale)	
	if(EM0) sigma.sq.b3<-sigma.sq.b.scale/(max(c(sigma.sq.b.shape-1,1)))

	if(fix.eff) sigma.sq.b3<-1e5

	} 
	
	output<-list("fits"=fits.main+ran.effect+ran.effect2+ran.effect2,"REs"=cbind(ran.effect,ran.effect2,ran.effect3),"sigmas"=c(sigma.sq.b,sigma.sq.b2,sigma.sq.b3))
	
	
		
	}##Close out update REs
