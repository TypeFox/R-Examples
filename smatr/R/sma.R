
sma <- function(formula, data, subset, na.action, log='',
	 method=c("SMA","MA","OLS"), type=c("elevation","shift"), alpha=0.05, 
	 slope.test=NA, elev.test=NA, multcomp=FALSE, multcompmethod=c("default","adjusted"),
	 robust=FALSE,V=matrix(0,2,2),
	 ...)
{
	method <- match.arg(method)
	type <- match.arg(type)
	multcompmethod <- match.arg(multcompmethod)
	
	# Model frame (code borrowed from 'lm')
	mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	if(ncol(mf) == 3)mf[,3] <- as.factor(mf[,3])
	
	# Throw out groups with <3 observations, if a grouping variable is used.
	# (Otherwise, there are in total 2 rows of data, surely no user is that daft).
	if(ncol(mf) == 3){
	tab <- table(mf[,3])
	if(any(tab < 3)){
		
		thislevel <- levels(mf[,3])[which(tab < 3)]
		mf <- mf[!(mf[,3] %in% thislevel),]
		mf <- droplevels(mf)
		attributes(mf)$row.names <- rownames(mf)
		cat("Warning: dropped level of grouping variable (sample size < 3) :",names(mf)[3]," = ",thislevel,"\n")
	}
	}
	
	# Log-transform 
	log <- tolower(log)
	
	# Throw out zeroes and negative values.
	xdel <- FALSE
	ydel <- FALSE
	if(log == "xy" || log == "y"){
		if(any(mf[,1] <= 0)){
			ind <- which(mf[,1] <= 0)
			ydel <- TRUE
			mf <- droplevels(mf[-ind,])
		}
	}
	if(log == "xy" || log == "x"){
		if(any(mf[,2] <= 0)){
			ind <- which(mf[,2] <= 0)
			xdel <- TRUE
			mf <- droplevels(mf[-ind,])
		}
	}
	
	if(log == "x")mf[,2] <- log10(mf[,2])
	if(log == "y")mf[,1] <- log10(mf[,1])
	if(log == "xy"){
		mf[,2] <- log10(mf[,2])
		mf[,1] <- log10(mf[,1])
	}
	if(!(log %in% c("","x","y","xy")))
		warning("Log transformation ignored! Use one of '', 'x','y' or 'xy'")
	
	# Prepare testing by group.
	mt <- attr(mf, "terms")
	vn <- names(mf)  # variable names
	tn <- attr(mt, "term.labels")   # term names
	
	# Check if there is group testing, and check that formula is compatible.
	if(length(tn) == 1)grouptest <- "none"
	if(length(tn) == 2){
		formulacheck <- all(tn == c(vn[2], vn[3]))
		if(formulacheck){
			if(type == "elevation")grouptest <- "elevcom"
			if(type == "shift")grouptest <- "shiftcom"
		} else grouptest <- "malformed"
	}
	# Slope test implied by formula (z+y+z:y = z*y)
	if(length(tn) == 3){
		formulacheck <- all(tn == c(vn[2], vn[3], paste(vn[2],":",vn[3],sep="")))
		if(formulacheck)grouptest <- "slopecom" else grouptest <- "malformed"
	}
	if(length(tn) > 3 || grouptest == "malformed"){
		warning("Formula not supported by sma() and/or ma(), and is ignored.")
		grouptest <- "none"
	}

	# Check for intercept. Also note that it is not allowed to drop the intercept for some group tests,
	# but this is not yet tested here!
	intercept <- if(attr(mt, "intercept") == 0) FALSE else TRUE
	
	# Halt execution when robust=T and no intercept.
	if(robust & !intercept)stop("Cannot perform robust estimation when no intercept included.")
	
	#Determine grouping 
	if(grouptest %in% c("elevcom","shiftcom","slopecom")){
		ngroups <- nlevels(mf[,3])
		
		# Fix the V matrix, when multiple grouping (cf. email david on Feb 4)
		if(length(dim(V)) == 2)V2 <- array( V, c(dim(V), ngroups) )
		grps<-mf[,3]
		lv <- levels(grps)
		
		
		if(method == "OLS"){
			commonslopetest <- NA
			commonslopetestval <- NA
			grouptestresult <- ""
		} else {
		
		# Whenever there are groups, do test for common slope.
		commonslopetest <- slope.com(mf[,1], mf[,2], mf[,3], alpha=alpha, intercept=intercept, method=method, group.names=lv, V=V2, robust=robust)
		
		# Test the common slope against hypthesized value, if this option is set.
		if(!is.na(slope.test)){
		commonslopetestval <- slope.com(mf[,1], mf[,2], mf[,3], alpha=alpha, slope.test=slope.test, intercept=intercept, method=method, group.names=lv, V=V2, robust=robust)
		} else {
		commonslopetestval <- NA
		}
		
		#run group tests
		if(grouptest == "elevcom"){
			if(!intercept)stop("Cannot perform elevation test without fitted intercept.")
			grouptestresult <- elev.com(mf[,1], mf[,2], mf[,3], alpha=alpha, method=method, group.names=lv, V=V2, robust=robust)
		}
		if(grouptest == "shiftcom"){
			grouptestresult <- shift.com(mf[,1], mf[,2], mf[,3], intercept=intercept, method=method, group.names=lv, V=V2, robust=robust)
		}
		if(grouptest == "slopecom")grouptestresult <- ""  #<-- already stored in commonslopetest
		
		}
	 }
	else{
	
	# single group
		ngroups<-1	
		grps <- as.factor(rep("all", length(mf[,1])))
		lv <- levels(grps)
		commonslopetest <- NA
		commonslopetestval <- NA
		grouptestresult <- ""
	}
	 
	#Calculate stuff for each group. Get the sma coefficients.
	coeff <- list(); n<- list(); r2<- list(); pval <- list(); from <- list()
	to<-list(); slopetest <-list(); elevtest <-list()
		
	for(i in 1:ngroups){
		X <- mf[grps == lv[i],2]
		Y <- mf[grps == lv[i],1]
		
		#groupsize
		n[[i]] <- length(X)
	
		#sma coefficients 
		coeff[[i]] <- line.cis(Y,X,intercept=intercept, method=method, 
			alpha=alpha, robust=robust, ...)   

		# correlation and P-value
		if(intercept){
			r2[[i]]<- cor(X, Y)^2
			pval[[i]] <- cor.test(X, Y, method = "pearson")$p.value
		} else {
			r2[[i]] <- sum(X*Y)^2/(sum(X^2) * sum(Y^2))
			pval[[i]] <- 1 - pf(r2[[i]]/(1-r2[[i]])*(n[[i]]-1),1,n[[i]]-1)  
		}
	  
      	# Test slope against some specified value
     	if(!is.na(slope.test)){
			slopetest[[i]] <- slope.test(Y,X,  test.value=slope.test, method=method, 
				alpha=alpha, intercept=intercept, robust=robust)
			slopetestdone <- TRUE
		} else {
			slopetest[[i]] <- slope.test(Y,X,  test.value=NA, method=method, 
				alpha=alpha, intercept=intercept, robust=robust)
			slopetestdone <- FALSE
		}
	
		# Test elevation against some specified value
		if(!is.na(elev.test)){
			if(!intercept)stop("Cannot perform elevation test without fitted intercept.")
				elevtest[[i]] <- elev.test( Y,X, test.value=elev.test, method=method, alpha=alpha, robust=robust)
				elevtestdone <- TRUE
		} else {
				elevtest[[i]] <- elev.test( Y,X, test.value=NA, method=method, alpha=alpha, robust=robust)
				elevtestdone <- FALSE
		}
      	
      	#determine range of fitted values (as X value)
        B <- coeff[[i]][2,1]
		a <- coeff[[i]][1,1]
		
        if(method=="SMA"){
  	    	from[[i]] <- (min(Y+B*X) - a)/(2.0*B) 
  	    	to[[i]] <-(max(Y+B*X)-a)/(2.0*B)
  	    } else if (method =="MA"){
  	    	from[[i]] <- (min(X+B*Y) - B*a)/(1+B^2)
  	    	to[[i]] <-(max(X+B*Y) - B*a)/(1+B^2)
   	    } else if (method == "OLS"){
			from[[i]] <- min(X)
			to[[i]] <- max(X)
		}
		
  	   	if(log %in% c("x","xy")){
			from[[i]] <- 10^from[[i]]
			to[[i]] <- 10^to[[i]]
		}    
	}
	
	# apply group names to new variables, if more than one group.
	if(ngroups > 1){
		names(coeff) <- lv
		names(n) <- lv
		names(r2) <- lv
		names(pval) <- lv
		names(from) <- lv
		names(to) <- lv
	}

	# coefficients under H0 (nullcoef)
	nullcoef <- NA
	if(grouptest == "none" & slopetestdone){   # sm2
		b <- slope.test
		a <- mean(Y) - b*mean(X)
		nullcoef <- c(a,b)
	}
	if(grouptest == "none" & elevtestdone){    # sm3
		a <- elev.test
		b <- NA
		nullcoef <- c(a,b)
	}
	if(grouptest == "slopecom" & !slopetestdone & !elevtestdone){ #sm4
		
		if(method == "OLS")
			nullcoef <- NA
		else {
		bcom <- commonslopetest$b
		nullcoef <- vector("list",ngroups)
		for(i in 1:ngroups){
			X <- mf[grps == lv[i],2]
			Y <- mf[grps == lv[i],1]
			a <- mean(Y) - bcom*mean(X)
			b <- bcom
			nullcoef[[i]] <- matrix(rep(NA,6),nrow=2)
			nullcoef[[i]][,1] <- c(a,b)
		}
		}
	}
	if(grouptest == "slopecom" & slopetestdone & !elevtestdone){ # sm5

		if(method == "OLS")
			nullcoef <- NA
		else {
		nullcoef <- vector("list",ngroups)
		for(i in 1:ngroups){
			X <- mf[grps == lv[i],2]
			Y <- mf[grps == lv[i],1]
			b <- slope.test
			a <- mean(Y) - b*mean(X)
			
			nullcoef[[i]] <- matrix(rep(NA,6),nrow=2)
			nullcoef[[i]][,1] <- c(a,b)
		}
		}
	}
	if(grouptest %in% c("elevcom","shiftcom") & !slopetestdone & !elevtestdone){
		
		# Overwrite the coefficients : these are assuming a common slope.
		bcom <- commonslopetest$b
		
		for(i in 1:ngroups){
			X <- mf[grps == lv[i],2]
			Y <- mf[grps == lv[i],1]
			a <- mean(Y) - bcom*mean(X)
			b <- bcom
			coeff[[i]][,1] <- c(a,b)
			coeff[[i]][2,2:3] <- commonslopetest$ci
			
			# Fix January 2012. New elev.com now reports CIs by group.
			if(grouptest=="elevcom")coeff[[i]][1,2:3] <- grouptestresult$as[i,2:3]
		}
		
		if(grouptest == "elevcom"){
			nullcoef <- vector("list",ngroups)
			for(i in 1:ngroups){
				a <- grouptestresult$a
				b <- bcom
				nullcoef[[i]] <- matrix(rep(NA,6),nrow=2)
				nullcoef[[i]][,1] <- c(a,b)
			}
		} else {
			nullcoef <- coeff
		}
	}
	
	l <- list()
	l$coef <- coeff
	l$nullcoef <- nullcoef
	l$commoncoef <- commonslopetest
	l$commonslopetestval <- commonslopetestval
	l$alpha <- alpha
	l$method <- method
	l$intercept <- intercept
	l$call <- match.call()
	l$data <- mf
	l$log <- log
	l$variables <- names(mf)
	l$origvariables <- all.vars(match.call()$formula)
  l$formula <- formula
	l$groups <- lv
	l$groupvarname <- if(ncol(mf) == 3)names(mf)[3] else NA
	l$gt <- grouptest
	l$gtr <- grouptestresult
	l$slopetest <- slopetest
	l$elevtest <- elevtest
	l$slopetestdone <- slopetestdone
	l$elevtestdone <- elevtestdone
	l$n <- n
	l$r2 <- r2
	l$pval <- pval
	l$from <- from
	l$to <- to
	sm <- list()
	f <- function(x)unname(unlist(x))
	
	for(i in 1:ngroups){
		sm[[i]] <- list(group=l$groups[i], n=f(l$n[i]), r2=f(l$r2[i]), pval=f(l$pval[i]),
		              Slope=l$coef[[i]][2,1], Slope_lowCI = l$coef[[i]][2,2], Slope_highCI = l$coef[[i]][2,3],  
					        Int = l$coef[[i]][1,1], Int_lowCI = l$coef[[i]][1,2], Int_highCI = l$coef[[i]][1,3],
					        Slope_test = if(all(is.na(slopetest)))NA else f(l$slopetest[[i]]$test.value), 
					        Slope_test_p= if(all(is.na(slopetest)))NA else f(l$slopetest[[i]]$p), 
					        Elev_test = if(all(is.na(elevtest)))NA else f(l$elevtest[[i]]$test.value),
					        Elev_test_p= if(all(is.na(elevtest)))NA else f(l$elevtest[[i]]$p))
	}
	tmp <- as.data.frame(do.call("rbind",sm))
	l$groupsummary <- as.data.frame(lapply(tmp, unlist))
	
	class(l) <- "sma"
	
	
	if(multcomp & ngroups == 2)
		warning("Multiple comparisons not performed, because there are only two groups (there is only one comparison to do)!")
	
	if(multcomp & ngroups > 2){
	
		# Pairwide group combinations:
		pairmat <- combn(l$groups, 2)
		npairs <- ncol(pairmat)
		paircall <- list()
		
		# Call sma with a data subsets:
		for(k in 1:npairs){
			
			# List of arguments as 'sma' was called.
			thiscall <- as.list(match.call(expand.dots = TRUE))[-1]
			
			# Set multcomp to FALSE
			thiscall$multcomp <- FALSE
			
			# Make data subset.
			#datasubs <- mf[mf[,3] %in% pairmat[,k],]
			datasubs <- data[data[,l$groupvarname] %in% pairmat[,k],]
			datasubs <- droplevels(datasubs)  # drop empty levels.
			thiscall$data <- datasubs
			
			# Set dummy argument, to keep track of what type of call this is:
			thiscall$dummy <- TRUE
			
			# Re-Call sma, using the data subset.
			paircall[[k]] <- do.call(what="sma", args=thiscall)
		}
		
		# Slope p values.
		multres <- as.data.frame(cbind(t(pairmat)))
		
		# Slope test (default, stored in commoncoef element).
		if(l$gt == "" | l$gt == "slopecom"){
			pvalsSlope <- sapply(paircall, function(x)x$commoncoef$p)
			multres$Pval <- pvalsSlope
			multcompdone <- "slope"
			
			# Test stat., df. 
			multres$TestStat <- sapply(paircall, function(x)x$commoncoef$LR)
			multres$df <- sapply(paircall, function(x)x$commoncoef$df)
			
			# Slope values
			multres$Slope1 <- sapply(paircall, function(x)x$commoncoef$bs[1,1])
			multres$Slope2 <- sapply(paircall, function(x)x$commoncoef$bs[1,2])
			
		}
		# elevation test, if performed.
		if(l$gt == "elevcom"){
			pvalsElev <- sapply(paircall, function(x)x$gtr$p)
			multres$Pval <- pvalsElev
			multcompdone <- "elevation"
			
			# Test stat., df. 
			multres$TestStat <- sapply(paircall, function(x)x$gtr$stat)
			multres$df <- sapply(paircall, function(x)x$gtr$df)
			
			# Elevation values.
			multres$Elev1 <- sapply(paircall, function(x)x$gtr$as[1])
			multres$Elev2 <- sapply(paircall, function(x)x$gtr$as[2])
			
		}
		if(l$gt == "shiftcom"){
			pvalsShift <- sapply(paircall, function(x)x$gtr$p)
			multres$Pval <- pvalsShift
			multcompdone <- "shift"
			
			# Test stat., df. 
			multres$TestStat <- sapply(paircall, function(x)x$gtr$stat)
			multres$df <- sapply(paircall, function(x)x$gtr$df)
			
			# Shift values.
			multres$Shift1 <- sapply(paircall, function(x)x$gtr$f.mean[1])
			multres$Shift2 <- sapply(paircall, function(x)x$gtr$f.mean[2])
		}
		
		# Label first two variable names in multcomp result.
		names(multres)[1:2] <- paste(l$groupvarname, "_", 1:2, sep="")
		
		# P-value Bonferroni adjustment
		if(multcompmethod == "adjusted"){
			adjusted.p.multi <- 1-(1-multres$Pval)^sum(multres$Pval >= 0,na.rm=TRUE)
			multres$Pval <- adjusted.p.multi
		}
		l$multcompresult <- multres
		l$multcompdone <- multcompdone
		l$multcompmethod <- multcompmethod

	}
	if(!multcomp){
		l$multcompresult <- NA
		l$multcompdone <- "none"
		l$multcompmethod <- NA
	}
	
	# Now print warning that X or Y values were deleted (complicated here, because only print once).
	thiscall <- as.list(match.call(expand.dots = TRUE))[-1]
	if(!("dummy" %in% names(thiscall)) && xdel)warning("Deleted negative and/or zero values in X variable.")
	if(!("dummy" %in% names(thiscall)) && ydel)warning("Deleted negative and/or zero values in Y variable.")
	
    return(l)
}
