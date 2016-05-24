calculate.accuracy.measures <-
function(data, marker, status, tag.healthy, direction = c("<", ">"), pop.prev, control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95) {
	direction <- match.arg(direction)
	
	# Validate the prevalence:
	if (is.na(pop.prev) || is.null(pop.prev)) {
		pop.prev <- calculate.sample.prev(data = data, status = status, tag.healthy = tag.healthy)
	}
	validate.prevalence(pop.prev)	
  
	cutoff <- sort(unique(data[,marker]))
	marker.healthy = data[data[,status] == tag.healthy, marker]
	marker.diseased = data[data[,status] != tag.healthy, marker]

	n = list(h = length(marker.healthy), d = length(marker.diseased))

	if(n$h == 0) {
  		stop("There are no healthy subjects in your dataset, so Specificity cannot be calculated.")
	}
	if(n$d == 0) {
  		stop("There are no diseased subjects in your dataset, so Sensitivity cannot be calculated.")
	}
	c.names <- if(ci.fit) {
		c("Value", "ll", "ul")
	} else {
		"Value"
	}
	Se <- Sp <- PPV <- NPV <- DLR.Positive <- DLR.Negative <- matrix(ncol=ifelse(ci.fit, 3, 1), nrow = length(cutoff), dimnames = list(1:length(cutoff), c.names))
	
	if(direction == "<") {
		testSe <- outer(marker.diseased,cutoff,">=")
		testSp <- outer(marker.healthy,cutoff,"<")
	} else {
		testSe <- outer(marker.diseased,cutoff,"<=")
		testSp <- outer(marker.healthy,cutoff,">")			
	}
	Se[,1] <- apply(testSe,2,sum)/(n$d)
	Sp[,1] <- apply(testSp,2,sum)/(n$h)

	PPV[,1] <- (pop.prev*Se[,1])/(pop.prev*Se[,1] + (1-pop.prev)*(1-Sp[,1]))	

	NPV[,1] <- ((1-pop.prev)*Sp[,1])/((1-pop.prev)*Sp[,1] + pop.prev*(1-Se[,1]))

	DLR.Positive[,1] <- Se[,1]/(1-Sp[,1])

	DLR.Negative[,1] <- (1-Se[,1])/Sp[,1]

	if(ci.fit == TRUE) {
  		ci <- confidence.intervals(Se[,1], Sp[,1], PPV[,1], NPV[,1], DLR.Positive[,1], DLR.Negative[,1], pop.prev, n, control, conf.level)
  		Se[,-1] <- ci$ci.Se
  		Sp[,-1] <- ci$ci.Sp
  		PPV[,-1] <- ci$ci.PPV
  		NPV[,-1] <- ci$ci.NPV
  		DLR.Positive[,-1] <- ci$ci.DLR.positive
  		DLR.Negative[,-1] <- ci$ci.DLR.negative
	}
		
	AUC <- calculate.empirical.AUC(data, marker, status, tag.healthy, direction, conf.level)

	res <- list(cutoffs = cutoff, Se = Se, Sp = Sp, PPV = PPV, NPV = NPV, DLR.Positive = DLR.Positive, DLR.Negative = DLR.Negative, AUC = AUC, pop.prev = pop.prev, n = n)
	return(res)
}
