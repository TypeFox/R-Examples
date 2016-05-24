"emplot" <- 
function(data, alog = "x", labels = TRUE, ...)
{
    data <- sort(as.numeric(data))
    ypoints <- 1 - ppoints(data)
    plot(data, ypoints, log = alog, xlab = "", ylab = "", ...)
    if(labels) {
        xxlab <- "x"
	yylab <- "1 - F(x)"
	if(alog == "x")
	    xxlab <- paste(xxlab, "(on log scale)")
        if(alog == "y")
	    yylab <- paste(yylab, "(on log scale)")
	if(alog == "xy" || alog == "yx") {
	    yylab <- paste(yylab, "(on log scale)")
            xxlab <- paste(xxlab, "(on log scale)")
        }
	title(xlab = xxlab, ylab = yylab)
    }
    invisible(list(x = data, y = ypoints))
}

"exindex" <- 
function(data, block, start = 5, end = NA, reverse = FALSE,
    auto.scale = TRUE, labels = TRUE, ...)
{
    sorted <- rev(sort(as.numeric(data)))
    n <- length(sorted)
    if(is.character(block)) {
        times <- as.POSIXlt(attributes(data)$times)
        if(block %in% c("semester", "quarter")) {
            sem <- quart <- times$mon
            sem[sem %in% 0:5] <- quart[quart %in% 0:2] <- 0
            sem[sem %in% 6:11] <- quart[quart %in% 3:5] <- 1
            quart[quart %in% 6:8] <- 2
            quart[quart %in% 9:11] <- 3
        }
        grouping <- switch(block,
            semester = paste(times$year, sem),
            quarter = paste(times$year, quart),
            month = paste(times$year, times$mon),
            year = times$year,
            stop("unknown time period"))
        b.lengths <- as.numeric(tapply(data, grouping, length))
        b.maxima <- as.numeric(tapply(data, grouping, max))
    }
    else {
	data <- as.numeric(data)
	nblocks <- (length(data) %/% block) + 1
	grouping <- rep(1:nblocks, rep(block, nblocks))[1:length(data)]
	b.lengths <- tapply(data, grouping, length)
	b.maxima <- tapply(data, grouping, max)
    }
    b.lengths <- b.lengths[!is.na(b.lengths)]
    b.maxima <- rev(sort(b.maxima[!is.na(b.maxima)]))
    if(is.numeric(block)) r <- block
    else r <- round(mean(b.lengths[2:(length(b.lengths) - 1)]))
    k <- round(n/r)
    un <- unique(b.maxima)[-1]
    K <- match(un, b.maxima) - 1
    N <- match(un, sorted) - 1
    if(is.na(end)) end <- k
    cond <- (K < end) & (K >= start)
    un <- un[cond]
    K <- K[cond]
    N <- N[cond]
    theta2 <- K/N
    theta <- logb(1 - K/k)/(r * logb(1 - N/n))
    out <- cbind(N, K, un, theta2, theta)
    yrange <- range(theta)
    index <- K
    if(reverse)	index <-  - K
    if(auto.scale)
        plot(index, theta, ylim = yrange, type = "l", xlab = "", ylab = "",
             axes = FALSE, ...)
    else plot(index, theta, type = "l", xlab = "", ylab = "", axes =
              FALSE, ...)
    axis(1, at = index, labels = paste(K), tick = FALSE)
    axis(2)
    axis(3, at = index, labels = paste(format(signif(un, 3))), tick = FALSE)
    box()
    if(labels) {
      	ylabel <- paste("theta (", k, " blocks of size ", r, ")", sep = "")
	title(xlab = "K", ylab = ylabel)
	mtext("Threshold", side = 3, line = 3)
    }
    invisible(out)
}

"hill" <- 
function(data, option = c("alpha","xi","quantile"), start = 15, end = NA,
    reverse = FALSE, p = NA, ci = 0.95, auto.scale = TRUE, labels = TRUE, ...)
{
    data <- as.numeric(data)
    ordered <- rev(sort(data))
    ordered <- ordered[ordered > 0]
    n <- length(ordered)
    option <- match.arg(option)
    if((option == "quantile") && (is.na(p)))
        stop("Input a value for the probability p")
    if((option == "quantile") && (p < 1 - start/n)) {
	cat("Graph may look strange !! \n\n")
	cat(paste("Suggestion 1: Increase `p' above",
                  format(signif(1 - start/n, 5)), "\n"))
	cat(paste("Suggestion 2: Increase `start' above ",
                  ceiling(length(data) * (1 - p)), "\n"))
    }
    k <- 1:n
    loggs <- logb(ordered)
    avesumlog <- cumsum(loggs)/(1:n)
    xihat <- c(NA, (avesumlog - loggs)[2:n])
    alphahat <- 1/xihat
    y <- switch(option,
	    alpha = alphahat,
	    xi = xihat,
	    quantile = ordered * ((n * (1 - p))/k)^(-1/alphahat))
    ses <- y/sqrt(k)
    if(is.na(end)) end <- n
    x <- trunc(seq(from = min(end, length(data)), to = start))
    y <- y[x]
    ylabel <- option
    yrange <- range(y)
    if(ci && (option != "quantile")) {
       	qq <- qnorm(1 - (1 - ci)/2)
       	u <- y + ses[x] * qq
       	l <- y - ses[x] * qq
       	ylabel <- paste(ylabel, " (CI, p =", ci, ")", sep = "")
       	yrange <- range(u, l)
    }
    if(option == "quantile") ylabel <- paste("Quantile, p =", p)
    index <- x
    if(reverse) index <-  - x
    if(auto.scale)
        plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "",
	     axes = FALSE, ...)
    else plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
    axis(1, at = index, labels = paste(x), tick = FALSE)
    axis(2)
    threshold <- findthresh(data, x)
    axis(3, at = index, labels = paste(format(signif(threshold, 3))),
         tick = FALSE)
    box()
    if(ci && (option != "quantile")) {
       	lines(index, u, lty = 2, col = 2)
       	lines(index, l, lty = 2, col = 2)
    }
    if(labels) {
       	title(xlab = "Order Statistics", ylab = ylabel)
       	mtext("Threshold", side = 3, line = 3)
    }
    invisible(list(x = index, y = y))
}

"meplot" <- 
function(data, omit = 3, labels = TRUE, ...)
{
    data <- as.numeric(data)
    n <- length(data)
    myrank <- function(x, na.last = TRUE)
    {
        ranks <- sort.list(sort.list(x, na.last = na.last))
	if(is.na(na.last))
	     x <- x[!is.na(x)]
	for(i in unique(x[duplicated(x)])) {
	    which <- x == i & !is.na(x)
	    ranks[which] <- max(ranks[which])
	}
	ranks
    }
    data <- sort(data)
    n.excess <- unique(floor(length(data) - myrank(data)))
    points <- unique(data)
    nl <- length(points)
    n.excess <- n.excess[-nl]
    points <- points[-nl]
    excess <- cumsum(rev(data))[n.excess] - n.excess * points
    y <- excess/n.excess
    xx <- points[1:(nl-omit)] ; yy <- y[1:(nl-omit)]
    plot(xx, yy, xlab = "", ylab = "", ...)
    if(labels) title(xlab = "Threshold", ylab = "Mean Excess")
    invisible(list(x = xx, y = yy))
}

"qplot" <- 
function(data, xi = 0, trim = NA, threshold = NA, line = TRUE,
    labels = TRUE, ...)
{
    data <- as.numeric(data)
    if(!is.na(threshold)) data <- data[data >= threshold]
    if(!is.na(trim)) data <- data[data < trim]
    if(xi == 0) {
        add <- "Exponential Quantiles"
	y <- qexp(ppoints(data))
    }
    if(xi != 0) {
        add <- paste("GPD Quantiles; xi =", xi)
	y <- qgpd(ppoints(data), xi = xi)
    }
    plot(sort(data), y, xlab = "", ylab = "", ...)
    if(labels) title(xlab = "Ordered Data", ylab = add)
    if(line) abline(lsfit(sort(data), y))
    invisible(list(x = sort(data), y = y))
}

"records" <- 
function(data, do.plot = TRUE, conf.level = 0.95, ...)
{
    data <- as.numeric(data)
    record <- cummax(data)
    expected <- cumsum(1/(1:length(data)))
    se <- sqrt(expected - cumsum(1/((1:length(data))^2)))
    trial <- (1:length(data))[!duplicated(record)]
    record <- unique(record)
    number <- 1:length(record)
    expected <- expected[trial]
    se <- se[trial]
    if(do.plot) {
       	ci <- qnorm(0.5 + conf.level/2)
       	upper <- expected + ci * se
       	lower <- expected - ci * se
       	lower[lower < 1] <- 1
       	yr <- range(upper, lower, number)
       	plot(trial, number, log = "x", ylim = yr, xlab = "Trial",
             ylab = "Records", main = "Plot of Record Development", ...)
	lines(trial, expected)
	lines(trial, upper, lty = 2)
	lines(trial, lower, lty = 2)
    }
    data.frame(number, record, trial, expected, se)
}

