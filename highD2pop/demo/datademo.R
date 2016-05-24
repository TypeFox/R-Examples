
data(chr1qseg)

impute<-function(x) { 	x[which(is.na(x))] <- mean(x,na.rm=TRUE)
						return(x)
					}

X <- apply(chr1qseg$X,2,impute)
Y <- apply(chr1qseg$Y,2,impute)

ChenQin.test(X,Y)

GCT.test(X,Y,r=20,smoother="parzen")
