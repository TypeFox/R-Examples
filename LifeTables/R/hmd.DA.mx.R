hmd.DA.mx <-
function(data, age, sex="female"){
	age.groups <- rep(NA, length(age))
	for(i in 1:length(age)){
	  age.groups[i] <- which(age[i]==c(0,1,seq(5,110,5)))
	}
# 	if(is.matrix(data)){
# 		age.groups <- ncol(data)
# 		} else age.groups <- length(data)
	if(!is.matrix(data)){
	  stop("'data' should be in matrix form with ages in columns")
	}
	
	if(sex=="male"){
		lt.m0.xp <- t(log(get("mlt.mx",envir=.GlobalEnv))[age.groups,])
		}	
		
	if(sex=="female"){
		lt.m0.xp <- t(log(get("flt.mx",envir=.GlobalEnv))[age.groups,])
		}
		
	hmd.m0.train <- MclustDA(data=exp(lt.m0.xp), class=get("class5",envir=.GlobalEnv))
	hmd.m0.test <- predict(hmd.m0.train, newdata=data)
	out.dens <- hmd.m0.test$z
	classification <- hmd.m0.test$classification
	
	return(list(train=hmd.m0.train, out.dens=out.dens, classification=classification))
	
	}

