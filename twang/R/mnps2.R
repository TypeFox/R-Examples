## require(twang); data(AOD); mnps.AOD4 <- mnps2(treat ~ illact + crimjust + subprob + subdep + white, data = AOD, estimand = "ATT", stop.method = "es.max", n.trees = 1000, treatATT = "community")
mnps2 <- function(formula = formula(data), data, n.trees = 10000,
interaction.depth = 3,
shrinkage = 0.01,
bag.fraction = 1.0,
perm.test.iters = 0,
print.level = 2,
iterlim = 1000,
verbose =TRUE,
estimand = "ATE",
stop.method = "es.max",
sampw = NULL,
treatATT = NULL, ...){
	stop.method <- levels(as.factor(stop.method))  ## alphbetizes for consitency in ordering of plots
	origStopMeth <- stop.method
	
	multinom <- TRUE
	
	if(is.null(sampw)) sampw <- rep(1, nrow(data))
	
	data <- data.frame(data, sampw=sampw)
	
		terms <- match.call()
   # all this is just to extract the variable names
   mf <- match.call(expand.dots = FALSE)
   
  
   m <- match(c("formula", "data"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf[[1]] <- as.name("model.frame")
   mf$na.action <- na.pass
#   mf$na.action <- NULL   
#   mf$na.action <- "na.pass"
   mf$subset <- rep(FALSE, nrow(data)) # drop all the data

   mf <- eval(mf, parent.frame())
   Terms <- attr(mf, "terms")
   resp <- attr(mf, "response")
   var.names <- attributes(Terms)$term.labels
   treat.var <- as.character(formula[[2]])		
		
	allowableStopMethods <- c("ks.mean", "es.mean","ks.max", "es.max")

		
	multFit <- gbm(formula = formula, data = data, weights = sampw, n.trees = n.trees, 
		interaction.depth = interaction.depth, shrinkage = shrinkage, bag.fraction=bag.fraction,
		verbose = verbose, distribution = "multinomial")
	
	nMethod <- length(stop.method)
	methodList <- vector(mode="list", length = nMethod)
	for(i in 1:nMethod){
		if(is.character(stop.method[i])){
			if (!(stop.method[i] %in% allowableStopMethods)){
				print(allowableStopMethods)
				stop("Each element of stop.method must either be one of \nthe above character strings, or an object of the stop.method class.")	
			}
		
		methodList[[i]] <- get(stop.method[i])
		methodName <- paste(stop.method[i], ".", estimand, sep="")
		methodList[[i]]$name <- methodName
		}
		else {
			if (class(stop.method[i]) != "stop.method"){
				print(allowableStopMethods)
				stop("Each element of stop.method must either be one of \nthe above character strings, or an object of the stop.method class.")	
			}
			methodList[[i]] <- stop.method[i]	
		}
	}

stop.method <- methodList	
	
	w <- ps <- data.frame(matrix(NA, nrow = nrow(data), ncol = nMethod))
	
	stop.method.names <- sapply(stop.method,function(x){x$name})
	
	names(w) <- names(ps) <- stop.method.names
	
	balance <- matrix(NA, ncol = nMethod, nrow = 25)
	
	respAll <- model.frame(formula, data = data, na.action = na.pass)[,1]
	
	if(!is.factor(respAll)) stop("The treatment variable must be a factor variable with at least 3 levels.")
	
	if(length(levels(respAll))<3) stop("The treatment variable must be a factor variable with at least 3 levels.")	
	
	respLev <- levels(respAll)
	M <- length(respLev)
	
	levExceptTreatATT <- NULL
	
	if(estimand == "ATT"){
		if(!(treatATT %in% respLev)) stop("'treatATT' must be one of the levels of the treatment variable.")
		else levExceptTreatATT <- respLev[respLev != treatATT]
	}
	
	dummyPS <- list(gbm.obj = NA,treat = NULL, treat.var = treat.var,
                 	desc = NULL, ps = NULL, w = NULL, sampw = NULL, estimand = estimand,
                  	datestamp  = date(), parameters = NULL, alerts = NULL,
                  	iters = NA, balance = NULL, n.trees = NA, data = NA)




	
	desc <- vector(mode = 'list', length = nMethod + 1)
	
	dummyPS$desc <- desc
		
	nFits <- ifelse(estimand == "ATT", M-1, M)
		
	mnpsObj <- list(psList = vector(mode = 'list', length = nFits), nFits = nFits, estimand = estimand, 
					treatATT = treatATT, treatLev = respLev,
					levExceptTreatATT = levExceptTreatATT, data = data, treatVar = respAll, treat.var = treat.var, stopMethods = stop.method,
					sampw = sampw)
					
	dummyPS$ps <- dummyPS$w <- matrix(NA, nrow = nrow(data), ncol = nMethod)	
	
	class(dummyPS) <- "ps"			
					
	for(i in 1:nFits) mnpsObj$psList[[i]] <- dummyPS

	
	for(i.tp in 1:nMethod){
		tp <- stop.method.names[i.tp]
	
	if(tp %in% c("es.max.ATT","es.mean.ATT","es.max.ATE","es.mean.ATE")) ES <- TRUE
	else ES <- FALSE
	
	if(tp %in% c("es.max.ATT","ks.max.ATT","es.max.ATE","ks.max.ATE")) summaryFcn <- max
	else summaryFcn <- mean
	
	compFcn <- max
   
   if(estimand == "ATT"){
   	if(is.null(treatATT)) stop("Must specify the 'treated' condition via the treatATT argument \n
   	when the estimand is set equal to ATT")
   }
   
					
		iters <- round(seq(1,multFit$n.trees,length=25))
	
	bal <- rep(0, length(iters))
	
	bal[1] <- pairwiseComp(wtVec = rep(1,nrow(data)), data = data, treat.var = treat.var, 
		vars = var.names, sampW = sampw, estimand = estimand, 
		compFcn = match.fun(compFcn), summaryFcn = match.fun(summaryFcn), treatATT = treatATT)[[1+ES]]
	
	for(j in 2:length(iters)){
		hldComp <- pairwiseComp(gbm1 = multFit, data = data, treat.var = treat.var, 
		vars = var.names, sampW = sampw, i = iters[j], estimand = estimand, 
		compFcn = match.fun(compFcn), summaryFcn = match.fun(summaryFcn), treatATT = treatATT)
		bal[j] <- hldComp[[1+ES]]
		#bal[i] <- max(apply(hldComp[[1+ES]], 1, match.fun(summaryFcn)))		
	}
	
	balance[,i.tp] <- bal
	
      interval <- which.min(bal) +c(-1,1)
      interval[1] <- max(1,interval[1])
      interval[2] <- min(length(iters),interval[2])
	
      opt<-optimize(MetricI,
                    interval=iters[interval],
                    maximum   = FALSE,
                    tol       = 1,
                    fun       = match.fun(pairwiseComp),
                    compFcn = match.fun(compFcn),
                    summaryFcn = match.fun(summaryFcn), 
                    vars      = var.names,
                    treat.var = treat.var,
                    data      = data,
                    sampW     = sampw,
#                    na.action = stop.method[[i.tp]]$na.action,
                    na.action = "level",
                    gbm1      = multFit,
                    estimand    = estimand,
                    treatATT = treatATT,
                    onlyES = (ES == 1), onlyKS = (ES == 0))
		
		
	if(verbose) cat("   Optimized at",round(opt$minimum),"\n")
	if(multFit$n.trees-opt$minimum < 100) warning("Optimal number of iterations is close to the specified n.trees. n.trees is likely set too small and better balance might be obtainable by setting n.trees to be larger.")
	
	fittedProbs <- predict(multFit, n.trees = opt$minimum, type = "response")[,,1]
	
	ps[,i.tp] <- fittedProbs[cbind(1:nrow(fittedProbs), as.numeric(respAll))]
	
	if(estimand == "ATT") w[,i.tp] <- makeWeightsMNPS(ps[,i.tp], estimand = "ATT", sampW = sampw, 
		treatATT = (respAll == treatATT), treatATTps = fittedProbs[,treatATT])
	else w[,i.tp] <- makeWeightsMNPS(ps[,i.tp], estimand = "ATE")
	
	if(estimand == "ATE")
		for(i in 1:M){
			tempDt <- data
			tempDt[,treat.var] <- data[,treat.var] == respLev[i]
			unwDesc <- desc.wts(tempDt,c(treat.var, var.names), 
									treat.var = treat.var,
									w = 1 + 0 * w[,i.tp],
									sampw = sampw,
									tp = "Make descriptions",
									na.action    = stop.method[[i.tp]]$na.action,
									verbose=verbose,
									alerts.stack = NULL,
									estimand       = estimand,
									multinom = TRUE)
			unwDesc$n.trees <- NA

			hldDesc <- desc.wts(tempDt,c(treat.var, var.names), 
									treat.var = treat.var,
									w = w[,i.tp],
									sampw = sampw,
									tp = "Make descriptions",
									na.action    = stop.method[[i.tp]]$na.action,
									verbose=verbose,
									alerts.stack = NULL,
									estimand       = estimand,
									multinom = TRUE)
			hldDesc$n.trees <- round(opt$minimum)						
			
			mnpsObj$psList[[i]]$gbm.obj <- NULL
			mnpsObj$psList[[i]]$treat <- tempDt[,treat.var]
			mnpsObj$psList[[i]]$treat.var <- treat.var
			mnpsObj$psList[[i]]$desc[[1]] <- unwDesc
			mnpsObj$psList[[i]]$desc[[i.tp + 1]] <- hldDesc		
			mnpsObj$psList[[i]]$ps[,i.tp] = fittedProbs[, i]
			mnpsObj$psList[[i]]$w[,i.tp] = w[,i.tp]
			mnpsObj$psList[[i]]$sampw = sampw
			mnpsObj$psList[[i]]$iters = iters
			mnpsObj$psList[[i]]$balance = balance
			mnpsObj$psList[[i]]$n.trees = round(opt$minimum)
			mnpsObj$psList[[i]]$data = NULL			
			names(mnpsObj$psList[[i]]$desc) <- c("unw", stop.method.names)	
			
	}
	
	
		if(estimand == "ATT")
		for(i in 1:(M-1)){
			tempDt <- data
			tempDt[,treat.var] <- data[,treat.var] == treatATT
			currentSubset <- data[,treat.var] %in% c(treatATT, levExceptTreatATT[i])
			tempDt <- subset(tempDt, currentSubset)
			unwDesc <- desc.wts(tempDt,c(var.names), 
									treat.var = treat.var,
									w = 1 + 0 * w[currentSubset,i.tp],
									sampw = sampw[currentSubset],
									tp = "Make descriptions",
									na.action    = stop.method[[i.tp]]$na.action,
									verbose=verbose,
									alerts.stack = NULL,
									estimand       = estimand,
									multinom = TRUE)
			unwDesc$n.trees <- NA

			hldDesc <- desc.wts(tempDt,c(var.names), 
									treat.var = treat.var,
									w = w[currentSubset,i.tp],
									sampw = sampw[currentSubset],
									tp = "Make descriptions",
									na.action    = stop.method[[i.tp]]$na.action,
									verbose=verbose,
									alerts.stack = NULL,
									estimand       = estimand,
									multinom = TRUE)
			hldDesc$n.trees <- round(opt$minimum)	
			
			if(i.tp == 1){
				mnpsObj$psList[[i]]$ps <- mnpsObj$psList[[i]]$w <-  data.frame(matrix(NA, nrow = sum(currentSubset), ncol = nMethod))
				names(mnpsObj$psList[[i]]$ps) <- names(mnpsObj$psList[[i]]$w) <- stop.method.names
				mnpsObj$psList[[i]]$sampw <- rep(NA, sum(currentSubset))
			}					
			
			mnpsObj$psList[[i]]$gbm.obj <- NULL
			mnpsObj$psList[[i]]$treat <- tempDt[,treat.var]
			mnpsObj$psList[[i]]$treat.var <- treat.var
			mnpsObj$psList[[i]]$desc[[1]] <- unwDesc
			mnpsObj$psList[[i]]$desc[[i.tp + 1]] <- hldDesc		
			mnpsObj$psList[[i]]$ps[,i.tp] <- fittedProbs[currentSubset, i]
			mnpsObj$psList[[i]]$w[,i.tp] <- w[currentSubset,i.tp]
			mnpsObj$psList[[i]]$sampw <- sampw[currentSubset]
			mnpsObj$psList[[i]]$iters <- iters
			mnpsObj$psList[[i]]$balance <- balance
			mnpsObj$psList[[i]]$n.trees <- round(opt$minimum)
			mnpsObj$psList[[i]]$data <- data			
			names(mnpsObj$psList[[i]]$desc) <- c("unw", stop.method.names)	

			
	}

	
   
   	
   	
   }

	origStopMethLong <- origStopMeth
	for(i in 1:length(origStopMeth)) origStopMethLong[i] <- paste(origStopMeth[i], ".", estimand, sep = "")
	for(i in 1:length(mnpsObj$psList)){
		mnpsObj$psList[[i]]$ps <- data.frame(mnpsObj$psList[[i]]$ps)
		names(mnpsObj$psList[[i]]$ps) <- origStopMethLong
	}

		
	returnObj <- mnpsObj
	returnObj$stopMethods <- origStopMeth
	returnObj$gbm.obj <- multFit
	
	class(returnObj) <- "mnps"


	return(returnObj)


}

#mnps(formula = treat ~ age + educ, data = lalonde2, estimand = "ATT", treatATT = "0")
#ft1 <- mnps(formula = treat ~ age + educ, data = lalonde2, estimand = "ATT", treatATT = "0")
