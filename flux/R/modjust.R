modjust <- function(models, alpha = 0.1, minimum = 0.8, prmtrs = list(t1 = 20), ...){
	# save class of models
	smc <- class(models)
	# define p.value extract function
	lmp <- function(mf){
		f <- summary(mf)$fstatistic
		if(is.null(f)){f <- c(70, 1, 30)}
		p <- pf(f[1],f[2],f[3],lower.tail=F)
		attributes(p) <- NULL
		return(p)
	}
	# define adjustment function
	checkfunc <- function(x){
		if(names(x$mod)[1]=="neu"){return(x)}
		else{
			mf <- lm(resid(x$mod[[1]]) ~ x$mod[[1]]$model$R)
			x$mod$n.out.adj <- 0
			min.n <- minimum*length(x$mod[[1]]$model$R)
			while((lmp(mf) > alpha) & (min.n < length(x$mod[[1]]$model$R))){
				selout <- which.max(x$mod[[1]]$model$R)
				newdata <- data.frame(R = x$mod[[1]]$model$R[-selout], Temp = x$mod[[1]]$model$Temp[-selout])
				reconew <- reco(newdata$R, newdata$Temp, Tref=x$mod[[1]]$model$Tref, T0=x$mod[[1]]$model$T0, ...)
				reconew <- reconew[names(reconew)==names(x$mod)[1]]
				if(class(reconew[[1]])=="try-error"){break}
				x$mod[1] <- reconew
				class(x$mod[[1]]) <- c("reco", class(reconew[[1]]))
				x$mod$n.out.adj <- x$mod$n.out.adj + 1
				mf <- lm(resid(reconew[[1]]) ~ reconew[[1]]$model$R)
			}
		return(x)
		}
	}
	# use the checkfunc
	models <- lapply(models, checkfunc)
	class(models) <- smc
	# if wanted skip models with parameters that go astray
	if(!is.null(prmtrs)){
		for(i in seq(length(prmtrs))){
			models <- models[tbl8(models)[,names(prmtrs)[i]] < prmtrs[[i]]]
		}
	}
	# reclass models
	class(models) <- smc
	return(models)
}