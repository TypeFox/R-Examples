wrightMap <-
function(thetas, thresholds = NULL, item.side = itemModern, person.side = personHist,  main.title = "Wright Map", min.logit.pad = 0.25, max.logit.pad = 0.25, min.l = NULL, max.l = NULL, item.prop = 0.8, return.thresholds = TRUE, new.quartz = FALSE, use.hist = NULL,...) {

	## Helper Functions


	if(is.null(thresholds))
		thresholds <- thetas
	
	thetas <- personData(thetas)
	
	thresholds <- itemData(thresholds,...)
	
	# Setting plot parameters
	thresholds <- as.matrix(thresholds)
	
	nD <- ncol(thetas)
	
	min.theta <- quantile(thetas, probs = c(0.01), na.rm = TRUE)
	max.theta <- quantile(thetas, probs = c(0.99), na.rm = TRUE)

	if (is.null(min.l)) {
		min.l <- min(c(min.theta, thresholds), na.rm = TRUE) - min.logit.pad
	}

	if (is.null(max.l)) {
		max.l <- max(c(max.theta, thresholds), na.rm = TRUE) + max.logit.pad
	}

	yRange <- c(min.l, max.l)
	
	# Generating Full Map
	
	if (new.quartz) 
		dev.new(width = 9, height = 5)
	
	op <- par("oma", "mar", "mgp")
	
	#par(oma = op$oma + c(0, 5, 0, 5))
	par(mar = c(op$mar[1], 0.2, op$mar[3], 0.1))
	#par(mgp = c(op$mar[1] - 2.4, 1, 0))

	par(mar = c(op$mar[1], 0.2, op$mar[3], 0.1))
	left.marg <- 0.05
	right.marg <- .1
	plots <- 1 - (left.marg + right.marg)
	divider <- plots * (1 - item.prop) + left.marg
	person.screen <- c(left.marg,divider,0,1)
	item.screen <- c(divider,1 - right.marg,0,1)
	screens <- matrix(c(item.screen,person.screen),ncol = 4,byrow = TRUE)
	old.screens <- split.screen()
	split.screen(screens)
	item.screen <- screen()
	dots <- list(...)
	dots[c("yRange","oma")] <- NULL
	item.params <- list(thr = thresholds,yRange = yRange,oma = c(0,0,0,0))
	do.call(item.side,c(item.params,dots))
	par(oma = c(0, 0, 3, 0))
	mtext(main.title, side = 3, line = 1, outer = TRUE, font = 2)
	close.screen(item.screen)
	if(!is.null(use.hist)) {
		message("Parameter 'use.hist' is deprecated. Please use 'person.side' parameter instead.")
		person.side <- ifelse(use.hist,personHist,personDens)
	}
	
	dots[c("axis.logits","show.axis.logits","close.on.close")] <- NULL
	person.params <- list(thetas = thetas,yRange = yRange,close.on.close = FALSE,oma = c(0,0,0,0),axis.logits = "",show.axis.logits = FALSE)
	do.call(person.side,c(person.params,dots))
	par(op)
	
	curr.screens <- split.screen()
	new.screens <- curr.screens[!(curr.screens %in% old.screens) ]
	close.screen(new.screens)
	
	if (return.thresholds) {
		return(thresholds)
	}



}
