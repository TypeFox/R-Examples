dmaxDesign <- function(n,dimension,range,niter_max=1000,seed=NULL){
	# if no seed is provided in argument, choice of the seed for 'runif'
	if (is.null(seed)){
		seed <- as.numeric(Sys.time())
	}
	set.seed(seed)
	
	# The spherical variogram model
	variosphe <- function(h,a){
		if (h <= a){
			gamma <- 1.5*(h/a)- 0.5*(h/a)^3
		}
		else{
			gamma <- 1
		}
		return(gamma)
	}

	p <- matrix(runif(n*dimension,0,1),n,dimension)
	p_init <- p
	
	# Hypothesis validation
	if (n < dimension){
		stop('Warning : the number of points is lower than the dimension.')
	}

 	if (range == 0){
		# To choose a large range in order to fill the region experimental
		range <- sqrt(dimension)/2 
 	}

	Distance <- as.matrix(dist(p))

	# Compute Correlation Matrix
	C <- diag(1,n,n)
	for (i in 1:(n-1)){
		for (k in (i+1):n){
			if (Distance[i,k] <= range){	 
		 		C[i,k] <- 1-variosphe(Distance[i,k],range)
		 	}
			C[k,i] <- C[i,k]
		}
	}

	# Compute determinant of the matrix C
 	dinit <- det(C)
 
	compt <- 1
	Deter <- rep(0,niter_max)
	Deter[compt] <- dinit

	while(compt < niter_max){
		compt <- compt+1
		# Choose an index 
		iunif <- floor(runif(1,1,n+1))
 		# Choose a potential new point 
		u <- runif(dimension)
		pu <- p
		pu[iunif,] <- u
		# Compute the potential new determinant of the matrix C (compute only row number iunif)
		Distu <- matrix(0,n,n)
       	for (ii in 1:n){
            	for (j in 1:dimension){          
            		Distu[iunif,ii]  <- Distu[iunif,ii] + (pu[iunif,j] - pu[ii,j])*(pu[iunif,j] - pu[ii,j])
               	}		   
       	}	
		Distu[iunif,] <- sqrt(Distu[iunif,])
		Cu <- C
		for (k in 1:n){
			if (k != iunif){
				if (Distu[iunif,k] <= range){	
					Cu[iunif,k] <- 1-variosphe(Distu[iunif,k],range)
				}
				Cu[k,iunif] <- Cu[iunif,k]
			}
		}
		du <- det(Cu)
		Deter[compt] <- dinit
		if((du - dinit) >= 0){
			Deter[compt] <- du
			p <- pu
			C <- Cu
			Distance <- Distu
			dinit <- du
		}
	
	 }# End while

	# Outputs:
 	return(list(n=n,dimension=dimension,range=range,niter_max=niter_max,design_init=p_init,design=p,det_init=Deter[1],det_end=dinit,seed=seed))
 }

