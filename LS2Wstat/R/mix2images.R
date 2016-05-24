mix2images <-function(imageA,imageB,prop=.25,pos="e"){

dima<-ncol(imageA)
dimb<-ncol(imageB)

dimprop<-prop*dimb 	# dimprop is how big the subimage proportion is

if(dimprop>dima){	
	side<-(dimb-dima)/2
	imageB[(side+1):(dimb-side),(side+1):(dimb-side)]<-imageA
}
else{
	side<-(dimb-dimprop)

	# get central part of imageA
	Aside<-(dima-dimprop)/2
	Asection<-imageA[(Aside+1):(dima-Aside),(Aside+1):(dima-Aside)]

	if(!is.numeric(pos)){
		switch(pos,
		c={
			a<-aa<-0
			b<-bb<-dimprop
		},
		d={
			a<-side
			b<-dimb
			aa<-0	
			bb<-dimprop
		},
		a={
			a<-0
			b<-dimprop
			aa<-side
			bb<-dimb
		},
		b={
			a<-aa<-side
			b<-bb<-dimb
		},
		e={
			a<-aa<-side/2
			b<-bb<-a+dimprop
		}
		)
	}
	else{	
		# specify the top-left corner pixel of the subimage
		a<-pos[2]-1
		aa<-pos[1]-1
		b<-a+dimprop
		bb<-aa+dimprop
	}

	imageB[(a+1):b,(aa+1):bb]<-Asection

}

return(imageB)

}

