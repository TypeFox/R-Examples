# Used from Alain Zuur Support package : MCMCSupportHighstat.R
# Zuur, Hilbe, Ieno (2013), Beginner's Guide to GLM and GLMM with R, Highlands
# http://www.highstat.com/
#############################################################

MyBUGSOutput <- function(xx,vars){
	 x <- xx$sims.matrix
    OUT <- matrix(nrow = length(vars), ncol=4)
    j<-1
	for(i in vars){
	  xi <- x[,i]
   	  OUT[j,3:4] <- quantile(xi, probs = c(0.025, 0.975))
   	  OUT[j,1] <- mean(xi)
   	  OUT[j,2] <- sd(xi)
   	  j <- j + 1
	}
	colnames(OUT) <- c("mean", "se", "2.5%", "97.5%")
	rownames(OUT) <- vars
	OUT
}

################################################################

uNames <- function(k,Q){
  #Function to make a string of variables names of the form:
  #c("u[1]","u[2]", etc, "u[50]")
  #Q=50 knots were used
  String<-NULL
  for (j in 1:Q){String <- c(String, paste(k,"[",j,"]",sep = ""))}
  String
}
