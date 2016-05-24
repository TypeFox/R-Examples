# These are some utility functions used in XNomial


chiStat <-
function(obs, exp) {
	pr <- to.probs(exp);
	nexp <- pr * sum(obs);
	return(sum((obs - nexp)^2/nexp))
}

constrain <-
function(x,b) {
	z <- x;
	if(x < b[[1]]) {z <- b[[1]]};
	if(x > b[[2]]) {z <- b[[2]]};
	return(z);
}


logLRmultinomial <-
function(obs, exp) {
	p <- to.probs(exp);
	expected <- p * sum(obs);
	r  <- sapply(obs/expected, function(x) {if(x > 0) log(x) else 0});
	return (-sum(r * obs))
}

ntables <-
function(x) {choose(sum(x) + length(x) -1, sum(x))}


to.probs <-
function(x) {
	n  <- sum(x);
	return(x/n);
}


statHistogram <-
function(c, showCurve=T) {
	n  <-  c$nn;
	nbins  <-  c$histobins;
	bounds  <-  c$histobounds;
	switch(c$statType,
	{  # If LLR
		statObs  <-  -2*c$observedLLR;
		r <- (statObs - bounds[[1]])/(bounds[[2]] - bounds[[1]]);
		name = "-2 ln(LR)";
	},{  # If Prob, the histoData and bounds are scalled to -2 log (pr/perfectProb), so we must compensate now
		statObs  <-  -2*log(c$observedProb);
		perfectProb  <- dmultinom(as.integer(sum(c$obs) * c$expr), prob = c$expr);
		r <- (-2 * log(c$observedProb/perfectProb) - bounds[[1]])/(bounds[[2]] - bounds[[1]]);
		name = "-2 ln(probability)"
		offset = -2 * log(perfectProb);
		bounds = bounds + offset;
		showCurve  <- F; # Not the asymptotic curve for probability
	},{  # If Chi
		statObs = c$observedChi;
		r <- (statObs - bounds[[1]])/(bounds[[2]] - bounds[[1]]);
		name = "Chi Sq"
	})

	nblack = constrain(nbins*r, c(0,nbins));
	colr <- c(rep("gray40", nblack), rep("lightcoral", nbins - nblack + 5));
	xlab = paste("Statistic:",name," = ", formatC(statObs), "    Range: (", formatC(bounds[[1]]),"  -  ", formatC(bounds[[2]]), ")");
	barplot(c$histoData, col = colr, border = NA, xlab = xlab, space=0);
	if(showCurve){
		if(is.null(c$ntrials)) {ntrials <- 1} else {ntrials <- c$ntrials};  #It's NULL if this is an exact test, not Monte Carlo
		dx <- 1:c$histobins;
		dy <- dchisq(seq(from=c$histobounds[[1]], to=c$histobounds[[2]], length.out=c$histobins), c$nn-1);
		lines(dx,dy * (ntrials) * (c$histobounds[[2]] - c$histobounds[[1]])/c$histobins, col="blue", lwd=2)
	}
}
