pim.plot <- function(object,predict=NULL,fun=exp,...){

	tapply.mi <- function(formula, data, index, family,...){
	
	if(!is.factor(index)) stop("Index must be a factor.")
	
	if(family=="coxph"){
		fits <- lapply(levels(index),function(grp){
			coxph(formula,data[index==grp,])
		})
	}
	else{
		fits <- lapply(levels(index),function(grp){
			glm(formula,data[index==grp,],family=family)
		})		
	}
	
	names(fits) <- levels(index)

fits
}

	get.group <- function(x,G=10){

	brks <- quantile(x,seq(0,1,length=(G+1)))

	list(
		 group = cut(x,br=brks,lab=1:G,include.lowest=TRUE),
		 levels = (brks[-1]+brks[-(G+1)])/2
		 )
	}

	get.estimates <- function(fit,trt,fun){
		fun(c(fit$coef[length(fit$coef)],confint(fit,trt)))
	}


        if(is.null(predict)){
          if(object@formula@family=="coxph"){
		predict <- predict(coxph(object@formula@prognostic,object@data))
              }
          else{
		predict <- predict(glm(object@formula@prognostic,object@data,
									family=object@formula@family))
              }

        }
        
	risk.index <- get.group(predict)
	f <- as.character(object@formula@formula)
	f[3] <- object@formula@trt
	f <- make.formula(f)
	
	if(object@formula@family=="coxph"){
		overall <- coxph(f,object@data)
	}
	else{
		overall <- glm(f,object@data,family=object@formula@family)
	}


	ci.effect <- get.estimates(overall,trt=object@formula@trt,fun=fun)	
	fits <- tapply.mi(f,object@data,risk.index$group,object@formula@family)
	effects <- mapply(get.estimates,fit=fits,
				MoreArgs=list(trt=object@formula@trt,fun=fun))
	
	ylim <- range(c(effects,ci.effect))
	
	arg.list <- list(...)
	
	if(length(arg.list)==0) arg.list <- list(xlab="Baseline prognosis")
	if(length(grep("xlab",names(arg.list)))==0) arg.list$xlab <- "Baseline prognosis"
	if(length(grep("ylab",names(arg.list)))==0) arg.list$ylab <- "Treatment effect"
	if(length(grep("ylim",names(arg.list)))==0) arg.list$ylim <- ylim
	
	arg.list$x <- risk.index$levels
	arg.list$y <- effects[1,]
	
	do.call(plot,arg.list)
	
	x <- par("usr")
			
	rect(xleft=x[1],xright=x[2],ybottom=ci.effect[2],ytop=ci.effect[3],
					col="grey80",border="transparent")
	
	points(x=risk.index$levels,y=effects[1,])
	
	for(i in 1:ncol(effects)){
		lines(x=rep(risk.index$levels[i],2),y=effects[2:3,i])
	}
}

setMethod("plot","anoint",
	function(x,...){
		pim.plot(object=x,...)
	}
)
