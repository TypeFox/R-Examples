hmd.DA <-
function(x, child.mort=4, sex="female", adult.mort=NULL){
	#data(MLTobs)
	if(is.numeric(adult.mort)){
		x <- as.matrix(cbind(x, adult.mort))
		}

	if(sex=="male"){
		if(!is.numeric(adult.mort)){
			if(child.mort==1){mod.train <- get("hmd.1m0.train.m",envir=.GlobalEnv)}
			if(child.mort==2){mod.train <- get("hmd.5m0.train.m",envir=.GlobalEnv)}
			if(child.mort==3){mod.train <- get("hmd.1q0.train.m",envir=.GlobalEnv)}
			if(child.mort==4){mod.train <- get("hmd.5q0.train.m",envir=.GlobalEnv)}
			}
			
		if(is.numeric(adult.mort)){
			if(child.mort==1){mod.train <- get("hmd.1m0a.train.m",envir=.GlobalEnv)}
			if(child.mort==2){mod.train <- get("hmd.5m0a.train.m",envir=.GlobalEnv)}
			if(child.mort==3){mod.train <- get("hmd.1q0a.train.m",envir=.GlobalEnv)}
			if(child.mort==4){mod.train <- get("hmd.5q0a.train.m",envir=.GlobalEnv)}
			}
		}	
		
	if(sex=="female"){
		if(!is.numeric(adult.mort)){
			if(child.mort==1){mod.train <- get("hmd.1m0.train.f",envir=.GlobalEnv)}
			if(child.mort==2){mod.train <- get("hmd.5m0.train.f",envir=.GlobalEnv)}
			if(child.mort==3){mod.train <- get("hmd.1q0.train.f",envir=.GlobalEnv)}
			if(child.mort==4){mod.train <- get("hmd.5q0.train.f",envir=.GlobalEnv)}
			}
			
		if(is.numeric(adult.mort)){
			if(child.mort==1){mod.train <- get("hmd.1m0a.train.f",envir=.GlobalEnv)}
			if(child.mort==2){mod.train <- get("hmd.5m0a.train.f",envir=.GlobalEnv)}
			if(child.mort==3){mod.train <- get("hmd.1q0a.train.f",envir=.GlobalEnv)}
			if(child.mort==4){mod.train <- get("hmd.5q0a.train.f",envir=.GlobalEnv)}
			}
		}
		

	test.mod <- predict(mod.train, newdata=x)
	out.dens <- test.mod$z
	classification <- test.mod$classification

	
	return(list(train=mod.train, out.dens=out.dens, classification=classification))
	}

