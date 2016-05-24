# based on .... by ....
# ported to dismo by Brian Anacker
# Date : October 2012
# Version 1.0
# Licence GPL v3


nicheEquivalency <- function (sp1, sp2, predictors, n=99, model=maxent, verbose=TRUE, ...) {

	runModel <- function(xy){
		x <- model(predictors, p=xy, a=bg, ...)
		predict(x, predictors)
	}
	

	di.enm <- function(preds){
	 	psum <- cellStats(preds, 'sum', na.rm = TRUE)
		s <- preds / psum
		d <- overlay(s, fun = function(i, j) { abs(i - j)	} )

		D <- 1 - 0.5 * cellStats(d, 'sum', na.rm = TRUE)

		h  <- overlay(s, fun = function(i, j) { (sqrt(i) - sqrt(j))^2 } )
		H2 <- cellStats(h, 'sum', na.rm = TRUE)
		
		I <- 1 - H2 * 0.5
		return(c(D = D, I = I))
	}

	sp1 <- .getMatrix(sp1)
	sp2 <- .getMatrix(sp2)

	p <- rbind(cbind(sp=1, sp1), cbind(sp=2, sp2))
    bg <- sampleRandom(predictors, size = min(ncell(predictors), 9999), na.rm = TRUE, xy = T)
    bg <- bg[, c('x', 'y')]

	pred_sp1 <- runModel(sp1)
	pred_sp2 <- runModel(sp2)
	if (verbose) { message('2 baseline models done'); flush.console() }
	
	preds <- stack(pred_sp1, pred_sp2)
	di <- di.enm(preds)

 	di.random <- sapply( 1:n, function(i) {
			s <- sample(p[,1])
			a <- runModel(p[s == 1, 2:3])
			b <- runModel(p[s == 2, 2:3])
			if (verbose) { message(paste(i, '-'), appendLF = FALSE); flush.console() }
			di.enm(stack(a, b))	
		} 
	)
	
	di.random <- t(di.random)

	if (verbose) message('done')	

    m <- colMeans(di.random)
    s <- apply(di.random, 2, sd)
    
	p.D <- pnorm(di[1], m[1], s[1], lower.tail = FALSE)
    p.I <- pnorm(di[2], m[2], s[2], lower.tail = FALSE)
	
    out <- list(method = "niche identity test",  
        null = "niche models are identical", statistic = di, 
        p.value = c(p.D, p.I), null.distribution = di.random)
    class(out) <- "ntest"
    out
}



#library(dismo)
#predictors <- stack(list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE ))  
#p = sampleRegular(predictors, 10000, asRaster=TRUE)
#occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
#occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
#spa <- occ[1:50, ]
#spb <- occ[51:116, ]
#ne <- nicheEquivalency(spa, spb, p, n=5, factors='biome')



