setClass("EBS",
         representation(model = "character", data= "numeric", length="numeric", Kmax="numeric", HyperParameters="numeric", Variance="numeric", overdispersion="numeric", Li="matrix", Col="matrix", matProba="matrix", unif="logical"),
         prototype(model = "Poisson", Kmax=15),
)

setMethod("show", "EBS",
	function (object){
		cat ( "Object of class EBS \n " )
		cat("\n Model used for the segmentation: \n")
		print(object@model)
		cat("\n Length of data: \n")
		print(object@length)
		cat("\n data: \n")
		str(object@data)
		cat("\n Maximum number of segments considered for the segmentation \n")
		print(object@Kmax)
		cat("\n Hyper-parameters used for prior on data distribution: \n")
		if(object@model=="Poisson")
		{
		  p=c(object@HyperParameters[1],object@HyperParameters[2])
		  cat("Gamma(alpha,beta), with (alpha,beta)= ") 
		  print(p)
		} else if (object@model=="Negative Binomial")
		{
		  p=c(object@HyperParameters[1],object@HyperParameters[2])
		  cat("Beta(alpha, beta) with (alpha,beta) = ") 
		  print(p)
		} else if (object@model=="Normal Homoscedastic")
		{
		  p=c(object@HyperParameters[1],object@HyperParameters[2])
		  cat("Normal(mu, sigma^2) with (mu,sigma^2)= ") 
		  print(p)
		} else 
		{
			p=c(object@HyperParameters[1],object@HyperParameters[2])
			q=c(object@HyperParameters[3],object@HyperParameters[4])
		  cat(" for mean: Normal(mu,sigma) with (mu,variance/sigma)= ") 
		  print(p)
		  cat("for inverse of variance: Gamma(alpha, beta) with 2(alpha,beta) = ") 
		  print(q)
		}
		if(object@model=="Negative Binomial")
		{
		  cat("\n used value for inverse of overdispersion \n")
		  print(object@overdispersion)
		}  
		if(object@model=="Normal Homoscedastic")
		{
		  cat("\n used value for variance\n")
		  print(object@Variance)
		}  
		cat("\n Log-proba [1,i[ in j segments: (getLi)")
		str(object@Li)
		cat("\n Log-proba [i,n+1[ in j segments: (getCol)")
		str(object@Col)
		cat("\n Log-proba [i,j[: (getP)")
		str(object@matProba)
})

setGeneric ("getLength",
	function(object){ standardGeneric ("getLength" )}
)
setMethod("getLength", "EBS",
	function (object){
	return ( object@length )
	}
)

setGeneric ("getModel",
	function(object){ standardGeneric ("getModel" )}
)
setMethod("getModel", "EBS",
	function (object){
	return ( object@model )
	}
)

setGeneric ("getData",
	function(object){ standardGeneric ("getData" )}
)
setMethod("getData", "EBS",
	function (object){
	return ( object@data )
	}
)

setGeneric ("getKmax",
	function(object){ standardGeneric ("getKmax" )}
)
setMethod("getKmax", "EBS",
	function (object){
	return ( object@Kmax )
	}
)

setGeneric ("getHyperParameters",
	function(object){ standardGeneric ("getHyperParameters" )}
)
setMethod("getHyperParameters", "EBS",
	function (object){
	return ( object@HyperParameters )
	}
)

setGeneric ("getVariance",
	function(object){ standardGeneric ("getVariance" )}
)
setMethod("getVariance", "EBS",
	function (object){
	return ( object@Variance )
	}
)

setGeneric ("getOverdispersion",
	function(object){ standardGeneric ("getOverdispersion" )}
)
setMethod("getOverdispersion", "EBS",
	function (object){
	return ( object@overdispersion )
	}
)

setGeneric ("getLi",
	function(object){ standardGeneric ("getLi" )}
)
setMethod("getLi", "EBS",
	function (object){
	return ( object@Li )
	}
)

setGeneric ("getCol",
	function(object){ standardGeneric ("getCol" )}
)
setMethod("getCol", "EBS",
	function (object){
	return ( object@Col )
	}
)

setGeneric ("getP",
	function(object){ standardGeneric ("getP" )}
)
setMethod("getP", "EBS",
	function (object){
	return ( object@matProba )
	}
)

setGeneric ("getPriorm",
	function(object){ standardGeneric ("getPriorm" )}
)
setMethod("getPriorm", "EBS",
	function (object){
	return ( object@unif )
	}
)  

############################################################

setClass("EBSProfiles",
         representation(model = "character", data= "matrix", length="numeric", NbConditions="numeric",  K="numeric", HyperParameters="numeric", Variance="numeric", overdispersion="numeric", Li="list", Col="list", P="list", unif="logical"),
         prototype(model = "Poisson"),
)

      
setMethod("show", "EBSProfiles",
	function (object){
		cat ( "Object of class EBSProfiles \n " )
		cat("\n Model used for the segmentation: \n")
		print(object@model)
		cat("\n Number of profiles considered: \n")
		print(object@NbConditions)
		cat("\n Length of each profile: \n")
		print(object@length)
		cat("\n data: (Data) \n")
		str(object@data)
		cat("\n Maximum number of segments considered for each profile \n")
		print(object@K)
		cat("\n Hyper-parameters used for prior on data distribution: \n")
		if(object@model=="Poisson")
		{
			cat('for each profile, Gamma(alpha,beta), with \n')
			for (i in 1:object@NbConditions)
			{
				p = c(object@HyperParameters[2*(i-1)+1],object@HyperParameters[2*(i-1)+2])
				cat("Profile ")
				print(i)
				cat("  (alpha,beta)= ") 
				print(p)
		  }
		} else if (object@model=="Negative Binomial")
		{
			cat('for each profile, Beta(alpha,beta), with \n')
			for (i in 1:object@NbConditions)
			{
				p = c(object@HyperParameters[2*(i-1)+1], object@HyperParameters[2*(i-1)+2])
				cat("Profile ")
				print(i)
				cat(" (alpha, beta)= ") 
				print(p)
		  }
		} else if (object@model=="Normal Homoscedastic")
		{
			cat('for each profile, Normal(mu, sigma^2), with \n')
			for (i in 1:object@NbConditions)
			{
				p = c(object@HyperParameters[2*(i-1)+1], object@HyperParameters[2*(i-1)+2])
				cat("Profile ")
				print(i)
				cat("  (mu, sigma^2)= ") 
				print(p)
		  }
		} else 
		{
			cat('for each profile, for mean: Normal(mu, sigma), for inverse of variance: Gamma(alpha, beta), with \n')
		  for (i in 1:object@NbConditions)
			{
				p = c(object@HyperParameters[4*(i-1)+1],object@HyperParameters[4*(i-1)+2])
				q = c(object@HyperParameters[4*(i-1)+3],object@HyperParameters[4*(i-1)+4])
				cat("Profile ")
				print(i)
				cat("  (mu,variance/sigma)=") 
				print(p)
				cat("  Gamma(2*alpha, 2* beta): \n 2 *alpha = ") 
				print(q)
		  }
  	}		
  	
  	if(object@model=="Negative Binomial")
		{
		  cat("\n used values for inverse of overdispersion \n")
		  print(object@overdispersion)
		}  
		if(object@model=="Normal Homoscedastic")
		{
		  cat("\n used values for variance\n")
		  print(object@Variance)
		}  
		cat("\n For each profile l: ")
		cat("\n   Log-proba [1,j[ in i segments: (Li(x)[[l]])")
		str(object@Li)
		cat("\n   Log-proba [i,n+1[ in j segments: (Col(x)[[l]])")
		str(object@Col)
		cat("\n   Log-proba [i,j[: (matProba(x)[[l]])")
		str(object@P)
})


setGeneric ("Length",
	function(object){ standardGeneric ("Length" )}
)
setMethod("Length", "EBSProfiles",
	function (object){
	return ( object@length )
	}
)

setGeneric ("Model",
	function(object){ standardGeneric ("Model" )}
)
setMethod("Model", "EBSProfiles",
	function (object){
	return ( object@model )
	}
)

setGeneric ("Data",
	function(object){ standardGeneric ("Data" )}
)
setMethod("Data", "EBSProfiles",
	function (object){
	return ( object@data )
	}
)

setGeneric ("Kmax",
	function(object){ standardGeneric ("Kmax" )}
)
setMethod("Kmax", "EBSProfiles",
	function (object){
	return ( object@K )
	}
)

setGeneric ("HyperParameters",
	function(object){ standardGeneric ("HyperParameters" )}
)
setMethod("HyperParameters", "EBSProfiles",
	function (object){
	return ( object@HyperParameters )
	}
)

setGeneric ("Variance",
	function(object){ standardGeneric ("Variance" )}
)
setMethod("Variance", "EBSProfiles",
	function (object){
	return ( object@Variance )
	}
)

setGeneric ("Overdispersion",
	function(object){ standardGeneric ("Overdispersion" )}
)
setMethod("Overdispersion", "EBSProfiles",
	function (object){
	return ( object@overdispersion )
	}
)

setGeneric ("Li",
	function(object){ standardGeneric ("Li" )}
)
setMethod("Li", "EBSProfiles",
	function (object){
	return ( object@Li )
	}
)

setGeneric ("Col",
	function(object){ standardGeneric ("Col" )}
)
setMethod("Col", "EBSProfiles",
	function (object){
	return ( object@Col )
	}
)

setGeneric ("matProba",
	function(object){ standardGeneric ("matProba" )}
)
setMethod("matProba", "EBSProfiles",
	function (object){
	return ( object@P )
	}
)

setGeneric ("NbConditions",
	function(object){ standardGeneric ("NbConditions" )}
)
setMethod("NbConditions", "EBSProfiles",
	function (object){
	return ( object@NbConditions )
	}
)

setGeneric ("Priorm",
	function(object){ standardGeneric ("Priorm" )}
)
setMethod("Priorm", "EBSProfiles",
	function (object){
	return ( object@unif )
	}
) 


