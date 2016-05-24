setGeneric("t.test",
	t.test
)

setMethod("t.test", "DVH.list",
	function(x, y=NULL, ..., paired=FALSE, mu=0, conf.level=0.95, var.equal=FALSE, alternative=c("two.sided", "greater", "less"), type=c("cumulative", "differential"), dose=c("absolute", "relative"), volume=c("relative", "absolute")) {
		alternative <- match.arg(alternative)
		type <- match.arg(type)
		dose <- match.arg(dose)
		volume <- match.arg(volume)
		if (class(y) != "DVH.list") {
			y <- new("DVH.list", y)
		}
		x <- new("DVH.list", lapply(x, convert.DVH, type=type, dose=dose, volume=volume))
		x <- x[!unlist(lapply(x, is.empty))]		
		y <- new("DVH.list", lapply(y, convert.DVH, type=type, dose=dose, volume=volume))
		y <- y[!unlist(lapply(y, is.empty))]		
		N.x <- length(x)
		N.y <- length(y)
		if (N.x < 2) {
			stop("not enough 'x' observations")
		}
		if (N.y < 2) {
			stop("not enough 'y' observations")
		}
		if (paired & (N.x != N.y)) {
			stop("'x' and 'y' must have the same length")
		}
		doses <- var(c(x, y))$dose
		data.x <- matrix(NA, nrow=length(doses), ncol=N.x)
		for (i in 1:N.x) {
			data.x[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses, rule=2)$y
		}
		data.y <- matrix(NA, nrow=length(doses), ncol=N.y)
		for (i in 1:N.y) {
			data.y[,i] <- approx(y[[i]]$doses, y[[i]]$volumes, doses, rule=2)$y
		}
		ttest.x <- ttest.y <- ttest.p <- ttest.conf1 <- ttest.conf2 <- c()
		for (i in 1:length(doses)) {
			t.i <- 1
			try(t.i <- t.test(data.x[i,], data.y[i,], paired=paired, mu=mu, conf.level=conf.level, var.equal=var.equal, alternative=alternative), silent=TRUE)
			if (identical(t.i, 1)) {
				ttest.x <- c(ttest.x, mean(data.x[i,], na.rm=TRUE))
				ttest.y <- c(ttest.y, mean(data.y[i,], na.rm=TRUE))
				ttest.p <- c(ttest.p, NA)
				ttest.conf1 <- c(ttest.conf1, NA)
				ttest.conf2 <- c(ttest.conf2, NA)
			}
			else {
				if (paired) {
					ttest.x <- c(ttest.x, mean(data.x[i,], na.rm=TRUE))
					ttest.y <- c(ttest.y, mean(data.y[i,], na.rm=TRUE))
				}
				else {
					ttest.x <- c(ttest.x, t.i$estimate[1])
					ttest.y <- c(ttest.y, t.i$estimate[2])
				}
				ttest.p <- c(ttest.p, t.i$p.value)
				ttest.conf1 <- c(ttest.conf1, t.i$conf.int[1])
				ttest.conf2 <- c(ttest.conf2, t.i$conf.int[2])
			}
		}
		return(list(dose=doses, x.mean=ttest.x, y.mean=ttest.y, p=ttest.p, conf.int1=ttest.conf1, conf.int2=ttest.conf2))
	}
)
