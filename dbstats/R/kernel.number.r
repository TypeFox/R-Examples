
 #######################
 #### kernel.number ####
 #######################

 ## description: Internal function. Kernel weights 
 ##
 ##        methods: 1 = Epanechnikov
 ##                 2 = Biweight
 ##                 3 = Triweight
 ##                 4 = Gaussiano
 ##                 5 = triangular
 ##                 6 = uniforme
 ##

kernel.number <- function(u,j=1){
	ku <- 0*u
	Iu <- which(abs(u)<1)
	aux <- switch(j,
		#
		# 1 = Epanechnikov
		.75*(1-u[Iu]^2),
		#
		# 2 = Biweight
		(15/16)*(1-u[Iu]^2)^2,
		#
		# 3 = Triweight
		(35/32)*(1-u[Iu]^2)^3,
		#
		# 4 = Gaussiano
		(exp(-.5* u^2))/sqrt(2*pi),
		#
		# 5 = triangular
		1-abs(u[Iu]),
		#
		# 6 = uniforme
		.5,
		#
		# Otherwise: Epanechnikov
		#.75*(1-u[Iu]^2)
	)
  if (j==4){
		ku <- aux
	}else{
		ku[Iu] <- aux
	}
	return(ku)
}