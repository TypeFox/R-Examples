MDR.stage.forward <-
function(x,y,order=NULL, s1.rsquared = NULL, s1.pvalue = NULL, s2.rsquared = NULL, s2.pvalue = NULL, max.step = NULL, trace = NULL,...){

	if(is.null(s1.rsquared))s1.rsquared <- 0.02 
	if(is.null(s1.pvalue))s1.pvalue <- 0.01
	if(is.null(s2.rsquared))s2.rsquared <- 0.02 
	if(is.null(s2.pvalue))s2.pvalue <- 0.01
	if(is.null(max.step))max.step <- 100
	if(is.null(trace))trace <- FALSE
	if(is.null(order)) order <- 2
	res <- MDR.stage.forward1(x=x,y=y,order=order, s1.rsquared = s1.rsquared, s1.pvalue = s1.pvalue, s2.rsquared = s2.rsquared, s2.pvalue = s2.pvalue, max.step = max.step, trace = trace)
	return(res)


}



MDR.stage.forward1 <-
function(x,y,order, s1.rsquared, s1.pvalue, s2.rsquared, s2.pvalue, max.step, trace,...){

	n <- ncol(x)
	index <- t(combn(n,order))
	


	c = ncol(index)
	r = nrow(index)

	#### Step 1: One-Way ANOVA ####
	index <- sing.mod(x,y,order,alpha = s1.pvalue,beta = s1.rsquared,delete = TRUE)
	index.step1 <- as.matrix(index[,1:c])
		
	#### Step 2: Stagewise Forward Selection ####
	id <- which(index[,c+1] > s2.pvalue)
	if(length(id) > 0)index <- index[-id,]
	r2 <- max(index[,c+3])
	id <- which(index[,c+3] == r2)
	index <- as.matrix(index[,1:c])	
	selection.index <- index[id[1],]				
	index <- IndexModify(index,selection.index)

	## Modify y ##	
	x.ba <- x[,selection.index]

	################################################
	ok <- complete.cases(x.ba) ## missing 3/3/2014
	x <- x[ok,]
	x.ba <- x[,selection.index]
	y <- as.matrix(y,,1)
	y <- y[ok,] 
	################################################

	if(c > 1){		
		for(i in 1:(c-1)){
			x.ba[,i+1]<- paste(x.ba[,i],x.ba[,i+1],sep=":")
		}
		x.ba <- x.ba[,c]
	}
	x.ba <- factor(x.ba)
	reg <- lm(y~x.ba)
	y <- resid(reg) ## y <- y - y.hat

	selection <- TRUE	
	step.num <- 1
	if(trace)cat("Step:",step.num,"\n")
	R <- NULL
	R[step.num] <- r2

	while(selection){#start(1)
	if(trace)cat("Step:",step.num + 1,"\n")
		r <- nrow(index)
		r2max <- 0
		################################################
		y <- as.matrix(y,,1) ## missing 3/3/2014
		################################################	
		for(i in 1:r){#start(2)
			x.ba <- x[,index[i,]]
			################################################
			ok <- complete.cases(x.ba) ## missing 3/3/2014
			################################################
			if(c > 1){		
				for(j in 1:(c-1)){
					x.ba[,j+1]<- paste(x.ba[,j],x.ba[,j+1],sep=":")
				}
				x.ba <- x.ba[,c]
			}
			x.ba <- factor(x.ba)
			################################################			
			reg <- lm(y[ok,]~x.ba[ok]) ## missing 3/3/2014
			################################################
			#reg <- lm(y~x.ba)
			summ <- summary(reg)
			F <- summ[[10]][[1]]
			df1 <- summ[[10]][[2]]
			df2 <- summ[[10]][[3]] 
			Pvalue <- 1-pf(F,df1,df2)
			if(!is.nan(Pvalue)){
				if(Pvalue < s2.pvalue){
					nR = length(y)
					r2 <- summ$r.squared
					r2adjust <- 1-(1-r2)*(nR-1)/(nR-2) # Adjusted R square
					r2adjust <- r2adjust*(1-sum(R))
					if(r2max < r2adjust){
						r2max <- r2adjust
						ind <- index[i,]
					}
				}
			}		
		}#end(2)

		step.num <- step.num + 1
	
		if(r2max >= s2.rsquared){
			index <- IndexModify(index,ind)
			r <- nrow(index)
			#index <- as.matrix(index,,c)
			selection.index <- rbind(selection.index,ind)
			R[step.num] <- r2max
			if(r>0){
				## modify y
				x.ba <- x[,ind]
				################################################
				ok <- complete.cases(x.ba) ## missing 3/3/2014
				x <- x[ok,]
				x.ba <- x[,selection.index]
				y <- as.matrix(y,,1)
				y <- y[ok,] 
				################################################
				if(c > 1){		
					for(i in 1:(c-1)){
						x.ba[,i+1]<- paste(x.ba[,i],x.ba[,i+1],sep=":")
					}
					x.ba <- x.ba[,c]
				}
				x.ba <- factor(x.ba)
				reg <- lm(y~x.ba)
				y <- resid(reg) ## y <- y - y.hat
			}else if(r<=0){
				selection <- FALSE
			}

		}else{
			selection <- FALSE
		}

		if(max.step == step.num)selection <- FALSE

	}#end(1)

	rownames(selection.index) <- NULL
	#res <- list(selection.index <- selection.index, R2 <- R, index.step1 <- index.step1)
	res <- list(selection.index <- selection.index, R2 <- R)
	return(res)
}
