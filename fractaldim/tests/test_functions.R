
generate.RF1d <- function(n, alpha=1) {
	rf <- GaussRF(x = c(0, 1, 1/n), model = "stable", 
				grid = TRUE, gridtriple = TRUE,
    				param = c(mean=0, variance=1, nugget=0, scale=1, kappa=alpha))
    return(rf)				
}

generate.RF2d <- function(n1, n2=n1, alpha=1) {
	rf <- GaussRF(x = c(0, 1, 1/n1), y = c(0, 1, 1/n2), model = "stable", 
				grid = TRUE, gridtriple = TRUE,
    				param = c(mean=0, variance=1, nugget=0, scale=1, kappa=alpha))
    return(rf)				
}

test.estimate.get <- function() {
	rf <- generate.RF1d(n=50)
	fd <- fd.estimate(rf, methods=c('variation', 'hallwood'))
	stopifnot(length(fd$fd) == 2)
	summary(fd)
	fd.hw <- fd.get(fd, 'hallwood')
	stopifnot(length(fd.hw$fd) == 1)
	cat('\nTest of estimate and get FD OK.\n')
}

test.estimate.boxcount <- function() {
	rf <- generate.RF1d(n=100)
	fd <- fd.estimate(rf, methods='boxcount', window.size=50)
	stopifnot(length(fd$fd) == 2)
	summary(fd)
}	
	
test.estimate.variation <- function() {
	rf <- generate.RF1d(n=50)
	fd <- fd.estimate(rf, methods=list('madogram', 'rodogram', 'variogram',
							'hallwood', list(name='variation',p.index=1),
							list(name='variation', p.index=0.5)))
	stopifnot(length(fd$fd) == 6)
	summary(fd)
	fd.var <- fd.get(fd, 'variation')
	stopifnot(length(fd.var$fd) == 2)
	stopifnot(fd$fd[1] == fd.var$fd[1])
	stopifnot(fd$fd[2] == fd.var$fd[2])
	cat('\nTest of estimate using variation OK.\n')
}


test.estimate.sw.get <- function() {
	rf <- generate.RF1d(n=64)
	fd <- fd.estimate(rf, methods=c('genton', 'boxcount', 'dctII'),
							window.size=33, step.size=10)
	stopifnot(all(dim(fd$fd) == c(4, 3)))
	summary(fd)
	fd.hw <- fd.get(fd, 'boxcount')
	stopifnot(all(dim(fd.hw$fd) == c(4,1)))
	summary(fd.hw)
	cat('\nTest of estimate with sliding window and get FD OK.\n')
}


test.estimate2d <- function() {
	rf <- generate.RF2d(n1=40)
	fd <- fd.estimate(rf, methods = c("squareincr", 'transect.var', 'transect.incr1', 'filter1'),
						window.size=20, p.index=1)
	stopifnot(all(dim(fd$fd) == c(2,2,4)))
	summary(fd)
	cat('\nTest of estimate using 2d estimators OK.\n')
}



