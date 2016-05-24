assignColorBreaks <- function(rates, NCOLORS = 64, spex = "s", logcolor = FALSE, method = c("linear","quantile","jenks"), JenksSubset = NULL) {
	method = match.arg(method, c("linear", "quantile", "jenks"));
		if (mode(rates) == "numeric") {
		if (logcolor == FALSE) {
			if (method == "quantile") {
				bks <- quantile(rates, seq(0,1, length.out=(NCOLORS+1)));
			}
			if (method == 'jenks') {
				bks <- getJenksBreaks(rates, k=(NCOLORS + 1), subset = JenksSubset);
			}
			if (method == 'linear') {
				bks <- seq(min(rates), max(rates), length.out = (NCOLORS+1));
			}
		}
		else {
			if (method == "quantile") {
				bks <- quantile(log(rates), seq(0,1, length.out=(NCOLORS+1)));
			}
			if (method == 'jenks') {
				bks <- getJenksBreaks(log(rates), k=(NCOLORS + 1), subset = JenksSubset);
			}
			if (method == 'linear') {
				bks <- seq(min(log(rates)), max(log(rates)), length.out = (NCOLORS+1));
			}
		}	
	}
	else if (mode(rates) == "list") {
		if (tolower(spex) == "s") {
			if (logcolor == FALSE) {
				if (method == "quantile") {
					bks <- quantile(rates[[1]], seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(rates[[1]], k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(rates[[1]]), max(rates[[1]]), length.out = (NCOLORS+1));
				}
			}
			else {
				if (method == "quantile") {	
					bks <- quantile(log(rates[[1]]), seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(log(rates[[1]]), k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(log(rates[[1]])), max(log(rates[[1]])), length.out = (NCOLORS+1));
				}
			}
		}
		else if (tolower(spex) == "e") {
			if (logcolor == FALSE) {
				if (method == "quantile") {
					bks <- quantile(rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(rates[[2]], k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(rates[[2]]), max(rates[[2]]), length.out = (NCOLORS+1));
				}
			}
			else {
				if (method == "quantile") {	
					bks <- quantile(log(rates[[2]]), seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(log(rates[[2]]), k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(log(rates[[1]])), max(log(rates[[1]])), length.out = (NCOLORS+1));
				}
			}
		}
		else if (tolower(spex) == "netdiv") {
			if (logcolor == FALSE) {
				if (method == "quantile") {	
					bks <- quantile(rates[[1]] - rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(rates[[1]] - rates[[2]], k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(rates[[1]] - rates[[2]]), max(rates[[1]] - rates[[2]]), length.out = (NCOLORS+1));
				}
			}
			else { 
				z <- safeLog(rates[[1]] - rates[[2]]);
				if (method == "quantile") {
					bks <- quantile(z, seq(0,1, length.out=(NCOLORS+1)));
					#bks <- quantile(log(rates[[1]] - rates[[2]]), seq(0,1, length.out=(NCOLORS+1)))
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(z, k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(z), max(z), length.out = (NCOLORS+1));
					#bks <- seq(min(log(rates[[1]] - rates[[2]])), max(min(log(rates[[1]] - rates[[2]]))), length.out=(NCOLORS+1) );
				}
				attr(bks, "increment") <- attr(z, "increment");
				return (safeLog(bks, inverse = TRUE));
			}
		}
	}
	
	if (logcolor)
		return (exp(bks));
	return (bks)
}


safeLog <- function(x, inverse = FALSE) {
	
	if (inverse)
		y <- exp(x) - attr(x, "increment")
	else {
		y <- log(x + abs(min(x)) + 0.0001);
		attr(y, "increment") <- abs(min(x)) + 0.0001; 
	}	
	return (y);
}
