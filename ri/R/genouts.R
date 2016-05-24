genouts <-
function(Y,Z,ate=0) {
	Y0 <- Y1 <- Y
	Y0[Z==1] <- Y[Z==1] - ate
	Y1[Z==0] <- Y[Z==0] + ate
	return(Ys=list(Y0=Y0,Y1=Y1))
	}
