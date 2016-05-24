print_uncertainty_nd <-
function(model,T,type="pn",lower=NULL,upper=NULL,
			resolution=20, nintegpoints=400,main="",
			cex.main=1,cex.lab=1,cex.contourlab=1,cex.axis=1,
			nlevels=10,levels=NULL,
			xdecal=3,ydecal=3, option="mean"){
	
	d <- model@d
  mynames <- colnames(model@X)
  if (resolution>40) resolution <- 40 #otherwise calculation time too long
  
	if(is.null(lower)) lower <- rep(0,times=d)
	if(is.null(upper)) upper <- rep(1,times=d)
	
	if (d==1){
		print("Error in print_uncertainty_nd, number of dimension is equal to 1. Please use function print_uncertainty_1d instead")
		return(0)
	}
	if (d==2){
		print("Error in print_uncertainty_nd, number of dimension is equal to 2. Please use function print_uncertainty_2d instead")
		return(0)
	}
	
	#d variables. d (d-1) / 2 pairs
	#number of points for the predict:
	# d(d-1)/2 . resolution^2 . nintegpoints
	#ex: d = 4, resolution = 50, nintegpoints = 1000 ----> 15,000,000
	
	sub.d <- d-2
	integration.tmp <- matrix(sobol(n=nintegpoints, dim = sub.d),ncol=sub.d)
	sbis <- c(1:resolution^2)
	numrow <- resolution^2 * nintegpoints
	integration.base <- matrix(c(0),nrow=numrow,ncol=sub.d)
	for(i in 1:sub.d){
		col.i <- as.numeric(integration.tmp[,i])
		my.mat <- expand.grid(col.i,sbis)
		integration.base[,i] <- my.mat[,1]	
	}
	
	
	s <- seq(from=0,to=1,length=resolution)
	sbis <- c(1:nintegpoints)
	s.base <- expand.grid(sbis,s,s)
	s.base <- s.base[,c(2,3)]
	
	prediction.points <- matrix(c(0),nrow=numrow,ncol=d)
  if(type=="vorob"){
    print("Vorob'ev plot not available in n dimensions.")
    print("We switch to a pn plot.")
    type <- "pn"
  }

	par(mfrow=c(d-1,d-1))
	for (d1 in 1:(d-1)){
		for (d2 in 1:d){
			
      if(d2!=d1){
			  if(d2 < d1){
          plot.new()
			  }else{
          #build prediction points
    			myindex <- 1
    			for (ind in 1:d){
    				if(ind==d1){
    					prediction.points[,ind] <- lower[ind] + s.base[,1] * ( upper[ind] - lower[ind] )
    					scale.x <- lower[ind] + s * ( upper[ind] - lower[ind] )
    					name.x <- mynames[d1]
    				}else if(ind==d2){
    					prediction.points[,ind] <- lower[ind] + s.base[,2] * ( upper[ind] - lower[ind] )
    					scale.y <- lower[ind] + s * ( upper[ind] - lower[ind] )
    					name.y <- mynames[d2] 
    				}else{
    					prediction.points[,ind] <- lower[ind] + integration.base[,myindex] * ( upper[ind] - lower[ind] )
    					myindex <- myindex + 1
    				}				
    			}
    			#prediction.points built !
          #nrow(prediction.points)
    			pred <- predict_nobias_km(object=model,newdata=prediction.points,type="UK",low.memory=TRUE)
          pn <- pnorm((pred$mean - T)/pred$sd)
    			
    			if(type=="pn") { 
            myvect <- pn
            zlim <- c(0,1)
    			}else if(type=="sur"){ 
            myvect <- pn * (1-pn)
            zlim <- c(0,0.25)
    			}else if(type=="timse"){
    				Wn <- 1/(sqrt(2*pi)*pred$sd)*exp(-1/2*((pred$mean-T)/pred$sd)^2)
    				myvect <- (pred$sd)^2 * Wn
            zlim <- NULL
    			}else if(type=="imse"){ 
            myvect <- (pred$sd)^2
            zlim <- NULL
    			}else{ 
            myvect <- pn
            zlim <- c(0,1)
          }
          
          #now we have to calculate resolution^2 integrals...
          myvect <- matrix(myvect,nrow=nintegpoints)
    			
          if (option == "mean") {myvect <- colMeans(myvect)
          }else if(option == "max"){myvect <- apply(X=myvect,MARGIN=2,FUN=max)
          }else if(option == "min"){myvect <- apply(X=myvect,MARGIN=2,FUN=min)
          }else{myvect <- colMeans(myvect)}
          
    			mymatrix <- matrix(myvect, nrow=resolution,ncol=resolution)
          if(is.null(zlim)) zlim <- c(min(mymatrix),max(mymatrix))
          
    			image(x=scale.x,y=scale.y,z=mymatrix,zlim=zlim,col=grey.colors(10),
    					xlab="",ylab="",cex.axis=cex.axis,main=main,
    					cex.main=cex.main,cex.lab=cex.lab,axes=TRUE)
    			
    			mtext(name.x, side=1, line=xdecal,cex=cex.lab ) 
    			mtext(name.y, side=2, line=ydecal,cex=cex.lab )
    			
    			if(!is.null(levels)){
    				contour(x=scale.x,y=scale.y,z=mymatrix,add=TRUE,labcex=cex.contourlab,levels=levels)	
    			} else {
    				contour(x=scale.x,y=scale.y,z=mymatrix,add=TRUE,labcex=cex.contourlab,nlevels=nlevels)
    			}
			  }        
      }
	    
      
		}
	}
	
	
		
	return(mean(myvect))
}
