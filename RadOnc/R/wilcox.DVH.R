setGeneric("wilcox.test",
	wilcox.test
)

setMethod("wilcox.test", "DVH.list",
	function(x, y=NULL, alternative=c("two.sided", "greater", "less"), mu=0, paired=FALSE, exact=TRUE, correct=TRUE, conf.level=0.95, ..., type=c("cumulative", "differential"), dose=c("absolute", "relative"), volume=c("relative", "absolute")) {
		type <- match.arg(type)
		dose <- match.arg(dose)
		volume <- match.arg(volume)
		alternative <- match.arg(alternative)
		if (class(y) != "DVH.list") {
			y <- new("DVH.list", y)
		}
		x <- new("DVH.list", lapply(x, convert.DVH, type=type, dose=dose, volume=volume))
		x <- x[!unlist(lapply(x, is.empty))]
		y <- new("DVH.list", lapply(y, convert.DVH, type=type, dose=dose, volume=volume))
		y <- y[!unlist(lapply(y, is.empty))]
		N.x <- length(x)
		N.y <- length(y)
		if (N.x < 1) {
			stop("not enough 'x' observations")
		}
		if (N.y < 1) {
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
		wilcox.x <- wilcox.y <- wilcox.p <- wilcox.conf1 <- wilcox.conf2 <- wilcox.est <- c()
		for (i in 1:length(doses)) {
			w.i <- 1
			suppressWarnings(try(w.i <- wilcox.test(data.x[i,], data.y[i,], paired=paired, mu=mu, conf.level=conf.level, exact=exact, alternative=alternative, conf.int=TRUE, correct=correct), silent=TRUE))
			wilcox.x <- c(wilcox.x, median(data.x[i,], na.rm=TRUE))
			wilcox.y <- c(wilcox.y, median(data.y[i,], na.rm=TRUE))
			if (identical(w.i, 1)) {
				wilcox.p <- c(wilcox.p, NA)
				wilcox.conf1 <- c(wilcox.conf1, NA)
				wilcox.conf2 <- c(wilcox.conf2, NA)
				wilcox.est <- c(wilcox.est, NA)
			}			
			else {
				wilcox.p <- c(wilcox.p, w.i$p.value)
				wilcox.conf1 <- c(wilcox.conf1, w.i$conf.int[1])
				wilcox.conf2 <- c(wilcox.conf2, w.i$conf.int[2])
				if (paired) {
					wilcox.est <- c(wilcox.est, median(data.x[i,]-data.y[i,], na.rm=TRUE))
				}
				else {
					wilcox.est <- c(wilcox.est, w.i$estimate)
				}
			}
		}
		return(list(dose=doses, x.med=wilcox.x, y.med=wilcox.y, p=wilcox.p, conf.int1=wilcox.conf1, conf.int2=wilcox.conf2, estimate=wilcox.est))
	}
)
