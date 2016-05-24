cropimage <-function(image,newsize=NULL,pos="e"){

d<-dim(image)

# work out how much to remove from edges of original image:

ld1<-log2(d[1])
ld2<-log2(d[2])

if(is.null(newsize)){
	newsizea<-2^floor(ld1)
	newsizeb<-2^floor(ld2)
}
else{
	newsizea<-newsizeb<-newsize
}

h<-d[1]-newsizea	
v<-d[2]-newsizeb	

if(!is.numeric(pos)){
	switch(pos,
	c={
		a<-aa<-0
		b<-h
		bb<-v
	},
	d={
		a<-bb<-0
		aa<-v	
		b<-h
	},
	a={
		a<-h
		b<-aa<-0
		bb<-v
	},
	b={
		a<-h
		b<-0
		aa<-v
		bb<-0
	},
	e={
		if((h%%2)==1){
			a<-h/2
			b<-h-a
		}
		else{
			a<-b<-h/2
		}
		if((v%%2)==1){
			aa<-v/2
			bb<-v-aa
		}
		else{
			aa<-bb<-v/2
		}
	}
	)

}
else{		
	# specify the position of the subimage
	# by giving the position of the top-left 
	# corner pixel

	a<-pos[1]-1
	b<-a+newsizea
	aa<-pos[2]-1
	bb<-aa+newsizeb
}

subim<-image[(a+1):(d[1]-b),(aa+1):(d[2]-bb)]

return(subim)
}

