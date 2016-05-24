
mnps <- function(formula = formula(data), data, n.trees = 10000,
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
	multinom <- TRUE
	
	if(is.null(sampw)) sampw <- rep(1, nrow(data))
	
	
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
   
   if(estimand == "ATT"){
   	if(is.null(treatATT)) stop("Must specify the 'treated' condition via the treatATT argument \n
   	when the estimand is set equal to ATT")
   }
   


	respAll <- model.frame(formula, data = data, na.action = na.pass)[,1]
	
	if(!is.factor(respAll)) stop("The treatment variable must be a factor variable with at least 3 levels.")
	
	respLev <- levels(respAll)
	M <- length(respLev)
	if(M < 3) stop("The treatment variable must be a factor variable with at least 3 levels.")
	
	levExceptTreatATT <- NULL
	
	if(estimand == "ATT"){
		if(!(treatATT %in% respLev)) stop("'treatATT' must be one of the levels of the treatment variable.")
		else levExceptTreatATT <- respLev[respLev != treatATT]
	}
	
	

	if(estimand == "ATE"){
		nFits <- M
		
		if(length(n.trees) == 1) n.trees <- rep(n.trees, nFits)
		
		hldFts <- vector(mode = "list", length = nFits)
		
		for(i in 1:nFits){
			## fit GBMs, etc
			currResp <- as.numeric(respAll == respLev[i])
			currDat <- data.frame(currResp = currResp, data)
			currFormula <- update(formula, currResp ~ .)
			currPs <- ps(formula = currFormula, data = currDat, n.trees = n.trees[i], interaction.depth = interaction.depth,
			shrinkage = shrinkage, bag.fraction = bag.fraction, perm.test.iters = perm.test.iters, print.level = print.level, 
			iterlim = iterlim,
			verbose = verbose, estimand = "ATE", stop.method = stop.method, sampw = sampw, multinom = TRUE)
			
			hldFts[[i]] <- currPs
				
		}
		names(hldFts) <- respLev	
	}
	
	if(estimand == "ATT"){
		nFits <- M - 1
		if(length(n.trees) == 1) n.trees <- rep(n.trees, nFits)
		
		hldFts <- vector(mode = "list", length = nFits)		
		for(i in 1:nFits){
			## subset data, do 
			currDat <- data[respAll == treatATT | respAll == levExceptTreatATT[i], ]
			currResp <- respAll[respAll == treatATT | respAll == levExceptTreatATT[i]]
			sampwCurr <- sampw[respAll == treatATT | respAll == levExceptTreatATT[i]]
#			currResp <- currResp == levExceptTreatATT[i]
			currResp <- currResp == treatATT			
			currDat <- data.frame(currResp = currResp, currDat)
			currFormula <- update(formula, currResp ~ .)
			currPs <- ps(formula = currFormula, data = currDat, n.trees = n.trees[i], interaction.depth = interaction.depth,
			shrinkage = shrinkage, bag.fraction = bag.fraction, perm.test.iters = perm.test.iters, print.level = print.level, 
			iterlim = iterlim,
			verbose = verbose, estimand = "ATT", stop.method = stop.method, sampw = sampwCurr, multinom = TRUE)
			
			hldFts[[i]] <- currPs

			
		}
		names(hldFts) <- levExceptTreatATT
		}
	
	returnObj <- list(psList = hldFts, nFits = nFits, estimand = estimand, treatATT = treatATT, treatLev = respLev, levExceptTreatATT = levExceptTreatATT, data = data, treatVar = respAll, treat.var = treat.var, stopMethods = stop.method, sampw = sampw)
	
	class(returnObj) <- "mnps"


	return(returnObj)


}

#mnps(formula = treat ~ age + educ, data = lalonde2, estimand = "ATT", treatATT = "0")
#ft1 <- mnps(formula = treat ~ age + educ, data = lalonde2, estimand = "ATT", treatATT = "0")
