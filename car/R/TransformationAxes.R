# Axes for transformations (J. Fox)

# last modified 27 December 2009 by J. Fox

# functions to add untransformed axis to right or top of a plot
#  for power, Box-Cox,  or Yeo-Johnson transformations

basicPowerAxis <- function(power, base=exp(1), side=c("right", "above", "left", "below"), 
	at, start=0, lead.digits=1, n.ticks, grid=FALSE, grid.col=gray(0.50), grid.lty=2,
	axis.title="Untransformed Data", cex=1, las=par("las")) {
	side <- if(is.numeric(side)) side 
		else which(match.arg(side) == c("below", "left", "above", "right"))
	axp <- if (side %% 2 == 1) par("xaxp") else par("yaxp")
	if (missing(n.ticks)) n.ticks <- axp[3] + 1
	ticks <- nice(seq(from=axp[1], to=axp[2], length=n.ticks), lead.digits=lead.digits)
	ticks.x <- if (power != 0) nice(ticks[ticks > 0]^(1/power), lead.digits=lead.digits) 
		else nice(log(base)*exp(ticks), lead.digits=lead.digits)
	ticks.x <- if (missing(at)) ticks.x
		else at
	ticks.text <- as.character(ticks.x - start)
	ticks.trans <- if (power != 0) ticks.x^power else log(ticks.x, base)
	axis(side, labels=ticks.text, at=ticks.trans, las=las)
	if (grid && (side %% 2 == 0)) abline(h=ticks.trans, lty=grid.lty, col=grid.col)
	if (grid && (side %% 2 == 1)) abline(v=ticks.trans, lty=grid.lty, col=grid.col)
	mtext(axis.title, side=side, line=3, cex=cex)
}

bcPowerAxis <- function(power, side=c("right", "above", "left", "below"), 
	at, start=0, lead.digits=1, n.ticks, grid=FALSE, grid.col=gray(0.50), grid.lty=2,
	axis.title="Untransformed Data", cex=1, las=par("las")) {
	inverse.power <- function(x, p){
		if (p == 0) exp(x)
		else (1 + p*x)^(1/p)
	}
	side <- if (is.numeric(side)) side 
		else which(match.arg(side) == c("below", "left", "above", "right"))
	axp <- if (side %% 2 == 1) par("xaxp") else par("yaxp")
	if (missing(n.ticks)) n.ticks <- axp[3] + 1
	ticks <- nice(seq(from=axp[1], to=axp[2], length=n.ticks), lead.digits=lead.digits)
	ticks.x <- if (power != 0) nice(inverse.power(ticks[ticks > 0], power), lead.digits=lead.digits)
		else nice(inverse.power(ticks, 0), lead.digits=lead.digits)
	ticks.x <- if (missing(at)) ticks.x
		else at
	ticks.text <- as.character(ticks.x - start)
	ticks.trans <- bcPower(ticks.x, power)
	axis(side, labels=ticks.text, at=ticks.trans, las=las)
	if (grid && (side %% 2 == 0)) abline(h=ticks.trans, lty=grid.lty, col=grid.col)
	if (grid && (side %% 2 == 1)) abline(v=ticks.trans, lty=grid.lty, col=grid.col)
	mtext(axis.title, side=side, line=3, cex=cex)
}
	
yjPowerAxis <- function(power, side=c("right", "above", "left", "below"), 
	at, lead.digits=1, n.ticks, grid=FALSE, grid.col=gray(0.50), grid.lty=2,
	axis.title="Untransformed Data", cex=1, las=par("las")) {
	inverse.bc <- function(x,p){
		if (p == 0) exp(x)
		else (1 + p*x)^(1/p)
	}
	inverse.power <- function(x, p){
		ifelse(x == 0, 0, ifelse(x > 0, inverse.bc(x, p) - 1, -inverse.bc(abs(x), 2 - p) + 1))
	}
	side <- if(is.numeric(side)) side 
		else which(match.arg(side) == c("below", "left", "above", "right"))
	axp <- if (side %% 2 == 1) par("xaxp") else par("yaxp")
	if (missing(n.ticks)) n.ticks <- axp[3] + 1
	ticks <- nice(seq(from=axp[1], to=axp[2], length=n.ticks), lead.digits=lead.digits)
	ticks.x <- nice(inverse.power(ticks, power), lead.digits=lead.digits)
	ticks.x <- if (missing(at)) ticks.x
		else at
	ticks.text <- as.character(ticks.x)
	ticks.trans <- yjPower(ticks.x, power)
	axis(side, labels=ticks.text, at=ticks.trans, las=las)
	if (grid && (side %% 2 == 0)) abline(h=ticks.trans, lty=grid.lty, col=grid.col)
	if (grid && (side %% 2 == 1)) abline(v=ticks.trans, lty=grid.lty, col=grid.col)
	mtext(axis.title, side=side, line=3, cex=cex)
}


# function to add a right or top probability axis to a plot of logits or probits

probabilityAxis <- function(scale=c("logit", "probit"), side=c("right", "above", "left", "below"),
	at, lead.digits=1, grid=FALSE, grid.lty=2, grid.col=gray(0.50),
    axis.title = "Probability", interval = 0.1, cex = 1, las=par("las")){
    side <- if(is.numeric(side)) side 
        else which(match.arg(side) == c("below", "left", "above", "right"))
	scale <- match.arg(scale)
	trans <- if (scale == "logit") function(p) log(p/(1 - p)) else qnorm
	inv.trans <- if (scale == "logit") function(x) 1/(1 + exp(-x)) else pnorm
    x <- if (side %% 2 == 1) par("usr")[c(1, 2)] else par("usr")[c(3, 4)]
    fact <- 10^( - (floor(log(interval, 10))))
    p.min <- nice(inv.trans(x[1]), direction="down", lead.digits=lead.digits)
    p.max <- nice(inv.trans(x[2]), direction="up", lead.digits=lead.digits)
    tick.min <- max(interval, (floor(fact*p.min))/fact)
    tick.max <- min(1 - interval, (ceiling(fact*p.max))/fact)
    ticks.p <- seq(tick.min, tick.max, interval)
	mins <- c(.05, .01, .005, .001, .0005, .0001)
	maxs <- c(.95, .99, .995, .999, .9995, .9999)
	ticks.p <- c(mins[mins >= p.min], ticks.p)
	ticks.p <- c(ticks.p, c(maxs[maxs <= p.max]))
    ticks.p <- if (missing(at)) ticks.p else at
    ticks.text <- as.character(ticks.p)
    ticks.x <- trans(ticks.p)
    axis(side, labels=ticks.text, at=ticks.x, las=las)
    if (grid && (side %% 2 == 0)) abline(h=ticks.x, lty=grid.lty, col=grid.col)
    if (grid && (side %% 2 == 1)) abline(v=ticks.x, lty=grid.lty, col=grid.col)
    mtext(axis.title, side=side, line=3, cex=cex)
}
