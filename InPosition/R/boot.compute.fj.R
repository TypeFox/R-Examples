

####The DATA that comes in is dependent on method (e.g., PCA vs. CA)
	#In the case of PCA, send in the already normalized matrix?
	#In the case of CA, just let it do its thing.

#for now, this has to take in the full Fixed.Data results.
boot.compute.fj <- function(DATA,res,DESIGN=NULL,constrained=FALSE){
	boot.index <- boot.samples(DATA=DATA,DESIGN=DESIGN,constrained=constrained)	

	output.types <- c("expoOutput")#,"texpoOutput","mexpoOutput")
	data.types <- c("ExPosition.Data")#,"TExPosition.Data","MExPosition.Data")	
	mds.types <- c('epMDS')#can add DiSTATIS to this.
	pca.types <- c('epPCA','epGPCA')#,'tepBADA')
	ca.types <- c('epCA','epMCA')#,'tepDICA')
	
	#I can probably trim this back...
	if(class(res)[1] %in% output.types){
		indicator <- which(class(res)[1] %in% output.types)
		if(sum(names(res) %in% data.types)==1 && length(names(res))==2){
			if(output.types[indicator]=="expoOutput"){
				if(!(class(res$ExPosition.Data)[1] %in% "epCA")){
					res$ExPosition.Data$fi <- res$ExPosition.Data$fi[boot.index,]
				}
				res <- res$ExPosition.Data
			}else{
				stop(paste("class of expoOutput required.", class(res)[1], "entered as res",sep=" "))	
			}				
		}else{
			stop(paste("res class type is unknown:",names(res),sep=" "))
		}
	}
	
	###some recognition needs to happen here...
	if((class(res)[1] %in% "epCA")){		
		boot.sup.data <- contingency.data.break(DATA,boot=TRUE)
		boot.index <- NULL
	}
	else{	
		boot.sup.data <- DATA[boot.index,]
	}

	#this should catch that things are NULL and skip them.
	if((class(res)[1] %in% c(pca.types))){
		fjj <- supplementaryCols(boot.sup.data,res,center=res$center,scale=res$scale)$fjj
	}else if((class(res)[1] %in% c(ca.types))){
		fjj <- supplementaryCols(boot.sup.data,res)$fjj
	}else if((class(res)[1] %in% c(mds.types))){ #this is the same as rows. 
		#fjj <- supplementaryCols(boot.sup.data,res)$fjj
		stop("epMDS data cannot be bootstrapped at this time. Please use epPCA.inference.battery()")
	}
	fjj<-replace(fjj,is.nan(fjj),0)
	return(fjj)
}