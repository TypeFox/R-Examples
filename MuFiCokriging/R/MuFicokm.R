MuFicokm <-
function (	formula , 
		MuFidesign, 
		response, 
		nlevel , 
		formula.rho = ~1, 
		covtype = "matern5_2", 
		coef.trend = NULL, 
		coef.rho = NULL,
		coef.cov = NULL, 
		coef.var = NULL, 
		nugget = NULL, 
   	 	nugget.estim = FALSE, 
		noise.var = NULL, 
		estim.method="MLE",
		penalty = NULL, 
		optim.method = "BFGS", 
    		lower = NULL, 
		upper = NULL, 
		parinit = NULL, 
		control = NULL, 
    		gr = TRUE, 
		iso = FALSE, 
		scaling = FALSE, 
		knots = NULL)
{

	model <- list()
for(i in 2:nlevel){
	if(is.null(coef.rho[[i-1]]) | is.null(coef.trend[[i]])){
		if(!is.null(coef.trend[[i]]) & !is.null(coef.rho[[i-1]])){
			stop("coef.trend and coef.rho must be both NULL or both non-NULL")
		}
	}
}

for(i in 2:nlevel){
	if(!is.null(coef.rho[[i-1]]) & !is.null(coef.trend[[i]])){
		coef.trend[[i]] <- c(coef.rho[[i-1]],coef.trend[[i]])
	}
}

	##initialisation
if(as.numeric(length(covtype))==1){
	covtype1 <- covtype 
}else{
	covtype1 <- covtype[[1]]
}
if(as.numeric(length(estim.method))==1){
	estim.method1 <- estim.method
}else{
	estim.method1 <- estim.method[[1]]
}
if(as.numeric(length(coef.trend))==1){
	coef.trend1 <- coef.trend
}else{
	coef.trend1 <- coef.trend[[1]]
}
if(as.numeric(length(coef.cov))==1){
	coef.cov1 <- coef.cov
}else{
	coef.cov1 <- coef.cov[[1]]
}
if(as.numeric(length(coef.var))==1){
	coef.var1 <- coef.var
}else{
	coef.var1 <- coef.var[[1]]
}
if(as.numeric(length(nugget))==1){
	nugget1 <- nugget
}else{
	nugget1 <- nugget[[1]]
}
if(as.numeric(length(nugget.estim))==1){
	nugget.estim1 <- nugget.estim
}else{
	nugget.estim1 <- nugget.estim[[1]]
}
if(as.numeric(length(noise.var))==1){
	noise.var1 <- noise.var
}else{
	noise.var1 <- noise.var[[1]]
}
if(as.numeric(length(penalty))==1){
	penalty1 <- penalty
}else{
	penalty1 <- penalty[[1]]
}
if(as.numeric(length(optim.method))==1){
	optim.method1 <- optim.method
}else{
	optim.method1 <- optim.method[[1]]
}
if(as.numeric(length(lower))==1){
	lower1 <- lower
}else{
	lower1 <- lower[[1]]
}
if(as.numeric(length(upper))==1){
	upper1 <- upper
}else{
	upper1 <- upper[[1]]
}
if(as.numeric(length(parinit))==1){
	parinit1 <- parinit
}else{
	parinit1 <- parinit[[1]]
}
if(as.numeric(length(control))==1){
	control1 <- control
}else{
	control1 <- control[[1]]
}
if(as.numeric(length(gr))==1){
	gr1 <- gr
}else{
	gr1 <- gr[[1]]
}
if(as.numeric(length(iso))==1){
	iso1 <- iso
}else{
	iso1 <- iso[[1]]
}
if(as.numeric(length(scaling))==1){
	scaling1 <- scaling
}else{
	scaling1 <- scaling[[1]]
}
if(as.numeric(length(knots))==1){
	knots1 <- knots
}else{
	knots1 <- knots[[1]]
}


if(dim(as.matrix(MuFidesign$PX))[2]==1){
	PX <-  data.frame(MuFidesign$PX)
	names(PX) <- "X1"
}else{
	PX <-  data.frame(MuFidesign$PX)
}


km.Z1 <- km(formula[[1]], 
		design =PX, 
		response = data.frame(response[[1]]),
		covtype=covtype1,
   		coef.trend = coef.trend1, 
		coef.cov = coef.cov1, 
		coef.var = coef.var1,
  		nugget = nugget1, 
		nugget.estim=nugget.estim1, 
		noise.var=noise.var1, 
		estim.method = estim.method1, 
		penalty = penalty1, 
   		optim.method =  optim.method1, 
		lower = lower1, 
		upper = upper1, 
		parinit = parinit1, 
   		control = control1, 
		gr = gr1, 
		iso=iso1, 
		scaling=scaling1, 
		knots=knots1)



zD <- list()
zD[[1]] <- list()
for(i in 2:nlevel){		
	zD[[1]][[i]] <- predict(
				object = km.Z1,
				newdata = data.frame(ExtractNestDesign(MuFidesign,i)),
				type = "UK",
				checkNames = FALSE)
}

km.Z <- list()
km.Z[[1]] <- km.Z1

for(k in 2:nlevel){

	##initialisation
	if(identical(formula.rho ,~1)){
		formulari <- formula.rho 
	}else{
		formulari <- formula.rho[[k-1]]
	}
	if(as.numeric(length(estim.method))==1){
		estim.method1 <- estim.method
	}else{
		estim.method1 <- estim.method[[k]]
	}
	if(as.numeric(length(covtype))==1){
		covtype1 <- covtype 
	}else{
		covtype1 <- covtype[[k]]
	}
	if(as.numeric(length(coef.trend))==1){
		coef.trend1 <- coef.trend
	}else{
		coef.trend1 <- coef.trend[[k]]
	}
	if(as.numeric(length(coef.cov))==1){
		coef.cov1 <- coef.cov
	}else{
		coef.cov1 <- coef.cov[[k]]
	}
	if(as.numeric(length(coef.var))==1){
		coef.var1 <- coef.var
	}else{
		coef.var1 <- coef.var[[k]]
	}
	if(as.numeric(length(nugget))==1){
		nugget1 <- nugget
	}else{
		nugget1 <- nugget[[k]]
	}
	if(as.numeric(length(nugget.estim))==1){
		nugget.estim1 <- nugget.estim
	}else{
		nugget.estim1 <- nugget.estim[[k]]
	}
	if(as.numeric(length(noise.var))==1){
		noise.var1 <- noise.var
	}else{
		noise.var1 <- noise.var[[k]]
	}
	if(as.numeric(length(penalty))==1){
		penalty1 <- penalty
	}else{
		penalty1 <- penalty[[k]]
	}
	if(as.numeric(length(optim.method))==1){
		optim.method1 <- optim.method
	}else{
		optim.method1 <- optim.method[[k]]
	}
	if(as.numeric(length(lower))==1){
		lower1 <- lower
	}else{
		lower1 <- lower[[k]]
	}
	if(as.numeric(length(upper))==1){
		upper1 <- upper
	}else{
		upper1 <- upper[[k]]
	}
	if(as.numeric(length(parinit))==1){
		parinit1 <- parinit
	}else{
		parinit1 <- parinit[[k]]
	}
	if(as.numeric(length(control))==1){
		control1 <- control
	}else{
		control1 <- control[[k]]
	}
	if(as.numeric(length(gr))==1){
		gr1 <- gr
	}else{
		gr1 <- gr[[k]]
	}
	if(as.numeric(length(iso))==1){
		iso1 <- iso
	}else{
		iso1 <- iso[[k]]
	}
	if(as.numeric(length(scaling))==1){
		scaling1 <- scaling
	}else{
		scaling1 <- scaling[[k]]
	}
	if(as.numeric(length(knots))==1){
		knots1 <- knots
	}else{
		knots1 <- knots[[k]]
	}

	if(dim(as.matrix(MuFidesign$PX))[2]==1){
		PX <-  data.frame(ExtractNestDesign(MuFidesign,k))
		names(PX) <- "X1"
	}else{
		PX <-  data.frame(ExtractNestDesign(MuFidesign,k))
	}
	
	km.Zi <- kmCok(formula[[k]], 
			design = PX , 
			response = data.frame(response[[k]]),
			formula.rho  = formulari, 
			Z = zD[[k-1]][[k]]$mean, 
			covtype=covtype1,
  			coef.trend = coef.trend1, 
			coef.cov = coef.cov1, 
			coef.var = coef.var1,
  			nugget = nugget1, 
			nugget.estim=nugget.estim1, 
			noise.var=noise.var1, 
			estim.method = estim.method1, 
			penalty = penalty1, 
   			optim.method =  optim.method1, 
			lower = lower1, 
			upper = upper1, 
			parinit = parinit1, 
   			control = control1, 	
			gr = gr1, iso=iso1, 
			scaling=scaling1, 
			knots=knots1)
	
	if((k+1)<=nlevel){
		zD[[k]] <- list()
		for(i in (k+1):nlevel){	
			if(identical(coef.trend1,NULL)){
				typepred <- "UK"
			}else{
				typepred <- "SK"
			}

			zD[[k]][[i]] <- predict.kmCok(
				object = km.Zi,
				newdata = data.frame(ExtractNestDesign(MuFidesign,i)),
				newZ = zD[[k-1]][[i]]$mean,
				type = "UK")
		}
	}

	km.Z[[k]] <- km.Zi
}

if(is.null(nugget)){
	nuggetout <- rep(0,nlevel)
}else{
	if(as.numeric(length(nugget))==1){
		nuggetout <- rep(nugget,nlevel)
	}else{
		nuggetout <- nugget[[1]]
		for(inug in 2:nlevel){
			nuggetout <- c(nuggetout, nugget[[inug ]])
		}
	}
}

model$cok <- km.Z
model$ZD <- zD
model$response <- response 
model$nlevel <- nlevel
model$Dnest <- MuFidesign
model$nuggets <- nuggetout

class(model) <- "MuFicokm"

return(model)

}
