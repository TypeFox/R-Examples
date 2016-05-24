discrepancyCriteria <- function(design,type='all'){
	#---------------------------------------
	# source code by Jessica FRANCO (2006.10.05)
  # modified by Bertrand Iooss (2013.26.12)
	#---------------------------------------
	# inputs
	# - design of experiments
	# - type of dicrepancies to be computed
	
	X <- as.matrix(design)
	dimension 	<- dim(X)[2]	# dimension 
	n 			<- dim(X)[1]	# number of points
	
	if ( n < dimension ){
	   stop('Warning : the number of points is lower than the dimension.')
	}
	# To check the experimental region
	if ( min(X)<0 || max(X)>1 ){
		warning("The design is rescaling into the unit cube [0,1]^d.")
		M <- apply(X,2,max)
		m <- apply(X,2,min)
		for (j in 1:dim(X)[2]){
			X[,j] <- (X[,j]-m[j])/(M[j]-m[j])
		}	
	}

	R <- list()
	DisC2 <- FALSE
	DisL2 <- FALSE
	DisL2star <- FALSE
	DisM2 <- FALSE
	DisS2 <- FALSE
	DisW2 <- FALSE
	
	if (length(type)==1 && type=='all'){
		type <- c('C2','L2','L2star','M2','S2','W2')	
	}
	for(i in 1:length(type)){
		type_ <- type[i]
	  	switch(type_,
        		C2 = {DisC2 <- TRUE},
	  	      L2 = {DisL2 <- TRUE},
	  	      L2star = {DisL2star <- TRUE},
        		M2 = {DisM2 <- TRUE},
        		S2 = {DisS2 <- TRUE},
        		W2 = {DisW2 <- TRUE})
	}

	# centered L2-discrepancy
	#------------------------
	if(DisC2==TRUE){
		s1 <- 0; s2 <- 0
		for (i in 1:n){
			p  <- prod((1+0.5*abs(X[i,]-0.5)-0.5*((abs(X[i,]-0.5))^2)))
			s1 <- s1+p
  			for (k in 1:n){
				q  <- prod((1+0.5*abs(X[i,]-0.5)+0.5*abs(X[k,]-0.5)-0.5*abs(X[i,]-X[k,])))
      			s2 <- s2+q
  			}
		}
		R <- c(R,DisC2 = sqrt(((13/12)^dimension)-((2/n)*s1) + ((1/n^2)*s2)))
	}
	# L2-discrepancy
	#------------------------
	if(DisL2==TRUE){
	  s1 <- 0; s2 <- 0
	  for (i in 1:n){
	    p  <- prod(X[i,]*(1-X[i,]))
	    s1 <- s1+p
	    for (k in 1:n){
	      q <- 1
	      for (j in 1:dimension){
	        q <- q*(1-max(X[i,j],X[k,j]))*min(X[i,j],X[k,j])
	      }
	      s2 <- s2+q
	    }
	  }
	  R <- c(R,DisL2 = sqrt(12^(-dimension) - (((2^(1-dimension))/n)*s1) + ((1/n^2)*s2)))
	}
	
	# L2star-discrepancy
	#------------------------
	if(DisL2star==TRUE){
	  dL2<-0
	  for (j in 1:n){
      for (i in 1:n){
	    if(i!=j){
        t<-c()
	      for (l in 1:dimension) t<-c(t,1-max(X[i,l],X[j,l]))
	      t<-(prod(t))/(n^2)
	    }
	    else{
        t1<-1-X[i,]
        t1<-prod(t1)
        t2<-1-X[i,]^2
        t2<-prod(t2)
        t<-t1/(n^2)-((2^(1-dimension))/n)*t2
	    }
      dL2<-dL2+t}
	  }  
	  R <- c(R,DisL2star = sqrt(3^(-dimension)+dL2))
	}
	
	# modified L2-discrepancy
	#------------------------
	if(DisM2 == TRUE){
		s1 <- 0; s2 <- 0
		for (i in 1:n){
   			p  <- 1
			p  <- prod((3-(X[i,]*X[i,])))
			s1 <- s1+p
			for (k in 1:n){
				q <- 1
     				for (j in 1:dimension){
					q <- q*(2-max(X[i,j],X[k,j]))
     				}
	 			s2 <- s2+q
  			}
		}
		R <- c(R,DisM2 = sqrt(((4/3)^dimension) - (((2^(1-dimension))/n)*s1) + ((1/n^2)*s2)))
	}

	# symmetric L2-discrepancy
	#------------------------
	if(DisS2 == TRUE){
		s1 <- 0; s2 <- 0
		for (i in 1:n){
			p <- prod((1+(2*X[i,])-(2*X[i,]*X[i,])))
			s1 <- s1+p
			for (k in 1:n){
      			q <- prod((1-abs(X[i,]-X[k,])))
     				s2 <- s2+q
  			}
		}
		R <- c(R,DisS2 = sqrt(((4/3)^dimension) - ((2/n)*s1) + ((2^dimension/n^2)*s2)))
	}
	
	# wrap-around L2-discrepancy
	#------------------------
	if(DisW2 == TRUE){
		s1 <- 0
		for (i in 1:n){
			for (k in 1:n){
      			p  <- prod((1.5-((abs(X[i,]-X[k,]))*(1-abs(X[i,]-X[k,])))))
     				s1 <- s1+p
  			}
		}
		R <- c(R , DisW2 = sqrt((((4/3)^dimension) + ((1/n^2)*s1))))
	}
	
	return(R)
}
