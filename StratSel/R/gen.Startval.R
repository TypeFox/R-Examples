gen.Startval <-
function(Startval, user.supplied.startval, corr, ys, xs11, xs14, xs24, dim.x11, dim.x14, dim.x24){
	if (user.supplied.startval==0) {
		beta1 <- coef(glm(as.numeric(ys!=1) ~ xs11[,-1] + xs14[,-1], family = binomial(link = "probit")))
		if (dim.x24[2]>=3) beta2 <- coef(glm(as.numeric(ys[which(ys!=1)]==4) ~ xs24[which(ys!=1),-1], family = binomial(link = "probit")))
		if (dim.x24[2]==2) beta2 <- coef(glm(as.numeric(ys[which(ys!=1)]==4) ~ xs24[which(ys!=1),-1], family = binomial(link = "probit")))
		Beta1 <- c(beta1[1:dim.x11[2]],0,beta1[(1+dim.x11[2]):(dim.x11[2]+dim.x14[2]-1)])
		Startval <- c(Beta1,beta2)
		#print(summary(glm(as.numeric(ys!=1) ~ xs11[,-1] + xs14[,-1], family = binomial(link = "probit"))))
	} 
	
	if(corr==TRUE & user.supplied.startval==0)  Startval <- c(Startval,0)
	
	if (NA %in% Startval) Startval[is.na(Startval)==TRUE] <- 0
	return(Startval)

}
