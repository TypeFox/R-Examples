# author: Jean-Pierre Rossi <jean-pierre.rossi@supagro.inra.fr>
# modifications by Robert Hijmans
# messi2 function by Paulo van Breugel


.messi2 <- function(p,v){
  f<-100*findInterval(p,sort(v))/length(v)
  a <- length(p)
  maxv <- max(v)
  minv <- min(v)
  opt1 <- 100*(p-minv)/(maxv-minv)
  opt2 <- 2*f
  opt3 <- 2 * (100-f)
  opt4 <- 100*(maxv-p)/(maxv-minv)
  simi <- ifelse(f==0, opt1, ifelse(f<=50, opt2, ifelse(f<100, opt3,opt4)))
  return(simi)
}


.messi2b <- function(p,v) {
# seems slightly faster than messi2
	v <- stats::na.omit(v)
	f <- 100*findInterval(p, sort(v)) / length(v)
	minv <- min(v)
	maxv <- max(v)
	ifelse(f == 0, 100*(p-minv)/(maxv-minv), 
		ifelse(f <= 50, 2*f, 
		ifelse(f < 100, 2*(100-f),
			100*(maxv-p)/(maxv-minv)
	)))
}


.messi3 <- function(p,v) {
# seems 2-3 times faster than messi2
	v <- stats::na.omit(v)
	f <- 100*findInterval(p, sort(v)) / length(v)
	minv <- min(v)
	maxv <- max(v)
	res <- 2*f 
	f[is.na(f)] <- -99
	i <- f>50 & f<100
	res[i] <- 200-res[i]

	i <- f==0 
	res[i] <- 100*(p[i]-minv)/(maxv-minv)
	i <- f==100
	res[i] <- 100*(maxv-p[i])/(maxv-minv)
	res
}


.messi4 <- function(p,v) {
	v <- stats::na.omit(v)
	f <- findInterval(p, sort(v)) / length(v)
	r <- range(v)
	res <- 2*f 
	f[is.na(f)] <- -99
	i <- f > 0.5 & f < 1
	res[i] <- 2-res[i]
	i <- f==0 
	res[i] <- (p[i]-r[1])/(r[2]-r[1])
	i <- f==1
	res[i] <- (r[2]-p[i])/(r[2]-r[1])
	res * 100
}



.messix <- function(p,v) {
# a little bit different, no negative values.
	a <- ecdf(v)(p)
	a[a>0.5] <- 1-a[a>0.5]
	200 * a
}



#P = runif(1000)
#V=runif(25)
#system.time(x <- .messi2a(P,V))
#system.time(y <- .messi2b(P,V))
#system.time(z <- .messi3(P,V))
#system.time(a <- .messi4(P,V))

	


.messOld <- function(x, v, full=FALSE) {
	stopifnot(NCOL(v) == nlayers(x))
	out <- raster(x)
	E <- getValues(x)

	nl <- nlayers(x)
	if (nl == 1) {
		rmess <- .messi2(E, v)
		names(out) <- 'mess'
		return( setValues(out, rmess) )

	} else {
		E <- sapply(1:ncol(E), function(i) .messi2(E[,i], v[,i]))
		rmess <- apply(x, 1, min, na.rm=TRUE)
		if (full) {
			out <- brick(out, nl=nl+1)
			names(out) <- c(names(x), "mess")
			return( setValues(out, cbind(E, rmess)) )
		} else {
			names(out) <- 'mess'
			return( setValues(out, rmess) )
		}
	}	
}




mess <- function(x, v, full=FALSE, filename='', ...) {

	stopifnot(NCOL(v) == nlayers(x))
	out <- raster(x)
	nl <- nlayers(x)
	filename <- trim(filename)
	nms <- paste(names(x), '_mess', sep='')
	
	if (canProcessInMemory(x)) {
		x <- getValues(x)
		if (nl == 1) {
			rmess <- .messi3(x, v)
			names(out) <- 'mess'
			out <- setValues(out, rmess)
		} else {
			x <- sapply(1:ncol(x), function(i) .messi3(x[,i], v[,i]))
			rmess <- apply(x, 1, min, na.rm=TRUE)
			if (full) {
				out <- brick(out, nl=nl+1)
				names(out) <- c(nms, "mess")
				out <- setValues(out, cbind(x, rmess))
			} else {
				names(out) <- 'mess'
				out <- setValues(out, rmess)
			}
		}	
		if (filename != '') {
			out <- writeRaster(out, filename, ...)
		}
		return(out)
		
	} else {

		if (nl == 1) {
		
			names(out) <- "mess"
			tr <- blockSize(out)
			pb <- pbCreate(tr$n, ...)	
			out <- writeStart(out, filename, ...)
			for (i in 1:tr$n) {
				vv <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
				vv <- .messi3(vv, v)
				out <- writeValues(out, vv, tr$row[i])
				pbStep(pb) 
			}
		
		} else {
	
			if (full) {
				out <- brick(out, nl=nl+1)
				names(out) <- c(nms, "mess")
				tr <- blockSize(out)
				pb <- pbCreate(tr$n, ...)	
				out <- writeStart(out, filename, ...)
				for (i in 1:tr$n) {
					vv <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
					vv <- sapply(1:ncol(v), function(i) .messi3(vv[,i], v[,i]))
					m <- apply(vv, 1, min, na.rm=TRUE)
					out <- writeValues(out, cbind(vv, m), tr$row[i])
					pbStep(pb) 
				}
				
			} else {
			
				names(out) <- "mess"
				tr <- blockSize(out)
				pb <- pbCreate(tr$n, ...)	
				out <- writeStart(out, filename, ...)
				for (i in 1:tr$n) {
					vv <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
					vv <- sapply(1:ncol(v), function(i) .messi3(vv[,i], v[,i]))
					m <- apply(vv, 1, min, na.rm=TRUE)
					out <- writeValues(out, m, tr$row[i])
					pbStep(pb) 
				}
			}
		}
		out <- writeStop(out)
		pbClose(pb) 
	}	
	out
}
