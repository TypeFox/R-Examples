dzipfman <- function (x, s = 1, b = NULL, N = NULL, log = FALSE) 
{
	if(is.null(N)){
		stop(paste("Must specify N.", "\n"))
	}
	out<-rep(0,length(x))
	temp <- (x!=floor(x) | x<1 | ((x==floor(x))&(x>N)))
	if(sum(temp)>0){
		out[which(temp)] <- 0
	}
	if(sum(temp)!=length(x)){
    	if ((is.null(b)) & (N < Inf)) {
        	if (s <= 0) {
            	stop(paste("Invalid value for s!", "\n"))
        	}
        	out[which(!temp)] <- x[which(!temp)]^(-s)/(sum((1:N)^(-s)))
    	}
    else if (!is.null(b)) {
        if (s <= 0 | b < 0) {
            stop(paste("Invalid value for s and/or b!", "\n"))
        }
        if (N == Inf) {
            stop(paste("N must be finite!", "\n"))
        }
        out[which(!temp)] <- (x[which(!temp)] + b)^(-s)/(sum((c(1:N) + b)^(-s)))
    }
    else {
        if (s <= 1) {
            stop(paste("Invalid value for s!", "\n"))
        }
        out[which(!temp)] <- x[which(!temp)]^(-s)/zeta.fun(s)
    }
}
    if (log) out <- log(out)
	if(any(x!=floor(x))){
		ind <- which(x!=floor(x))
		for(i in 1:length(ind)){
			warning(paste("non-integer x =",x[ind[i]]), call. = TRUE)
		}
	}
	out
}

rzipfman <- function(n, s = 1, b = NULL, N = NULL)
{
	if(is.null(N)){
		stop(paste("Must specify N.", "\n"))
	}
    if ((is.null(b)) & (N < Inf)) {
        if (s <= 0) {
            stop(paste("Invalid value for s!", "\n"))
        }
        out <- qzipfman(p = runif(n), s = s, N = N)
    }
    else if (!is.null(b)) {
        if (s <= 0 | b < 0) {
            stop(paste("Invalid value for s and/or b!", "\n"))
        }
        if (N == Inf) {
            stop(paste("N must be finite!", "\n"))
        }
        out <- qzipfman(runif(n), s = s, b = b, N = N)
    }
    else {
        if (s <= 1) {
            stop(paste("Invalid value for s!", "\n"))
        }
        out <- qzipfman(runif(n), s = s, N = Inf)
    }
    out.lvl <- names(sort(table(out),decreasing = TRUE))
    y <- 1:length(unique(out))
    out <- y[match(out,out.lvl)]
    out
}


qzipfman <- function (p, s = 1, b = NULL, N = NULL, lower.tail = TRUE, log.p = FALSE)
{
	if(is.null(N)){
		stop(paste("Must specify N.", "\n"))
	}
    if (log.p) p <- exp(p)
    if (lower.tail == FALSE) p <- 1 - p
    if ((is.null(b)) & (N < Inf)) {
        if (s <= 0) {
            stop(paste("Invalid value for s!", "\n"))
        }
	all.p<-rep(1,length(p))
	temp.ind <- NULL
	if(any((p > 1) | (p < 0))){
		temp.ind <- which((p > 1) | (p < 0))
		all.p[temp.ind] <- -Inf
	}
		temp <- (pzipfman(q = 1, s = s, N = N) >= p | p>1)
		if (sum(!temp)>0 & N<1e6){
		p.temp <- cumsum(dzipfman(1:N, s=s, N=N))
		out.p <- sapply(1:length(which(!temp)), function(i) min(which(p.temp >= p[which(!temp)[i]])))
		out.p2 <- sapply(1:length(which(!temp)), function(i) suppressWarnings(max(which(p.temp == p[which(!temp)[i]]))))
		all.p[which(!temp)] <- pmax(out.p, out.p2)
            	} else if(sum(!temp)>0 & N >= 1e6){
			x <- c(1e+06, 1e+300)
			y <- c(sum(dzipfman(1:1e6, s = s, N = N)) , 1)
			which.p <- which(p >= y[1] & p < 1)
			out.lm <- lm(y~I(-1/x))
			be <- out.lm$coef
			new.q <- -1/((p[which.p]-be[1])/be[2])
			if(any(new.q >= 1e+300 | new.q == -Inf)) all.p[which(new.q >= 1e+300)] = Inf
			if(any((new.q < 1e+300) & (new.q >= 1e+06))) all.p[which.p[which((new.q < 1e+300) & (new.q >= 1e+06))]] = floor(new.q[which((new.q < 1e+300) & (new.q >= 1e+06))])
		}
    } else if (!is.null(b)) {
        if (s <= 0 | b < 0) {
            stop(paste("Invalid value for s and/or b!", "\n"))
        }
        if (N == Inf) {
            stop(paste("N must be finite!", "\n"))
        }
	all.p<-rep(1,length(p))
	temp.ind <- NULL
	if(any((p > 1) | (p < 0))){
		temp.ind <- which((p > 1) | (p < 0))
		all.p[temp.ind] <- -Inf
	}
		temp <- (pzipfman(q = 1, s = s, b =b, N = N) >= p | p > 1)
            if (sum(!temp)>0 & N < 1e6) {
		p.temp <- cumsum(dzipfman(1:N, s = s, b = b, N = N))
		out.p <- sapply(1:length(which(!temp)), function(i) min(which(p.temp >= p[which(!temp)[i]])))
		out.p2 <- sapply(1:length(which(!temp)), function(i) suppressWarnings(max(which(p.temp == p[which(!temp)[i]]))))
		all.p[which(!temp)] <- pmax(out.p, out.p2)
            } else if(sum(!temp)>0 & N >= 1e6){
			x <- c(1e+06, 1e+300)
			y <- c(sum(dzipfman(1:1e6, s = s, b = b, N = N)) , 1)
			which.p <- which(p >= y[1] & p < 1)
			out.lm <- lm(y~I(-1/x))
			be <- out.lm$coef
			new.q <- -1/((p[which.p]-be[1])/be[2])
			if(any(new.q >= 1e+300 | new.q == -Inf)) all.p[which(new.q >= 1e+300)] = Inf
			if(any((new.q < 1e+300) & (new.q >= 1e+06))) all.p[which.p[which((new.q < 1e+300) & (new.q >= 1e+06))]] = floor(new.q[which((new.q < 1e+300) & (new.q >= 1e+06))])
		}
    } else {
        if (s <= 1) {
            stop(paste("Invalid value for s!", "\n"))
        }
        all.p <- rep(NA,length(p))
	if(any((p > 1) | (p < 0))){
		all.p[which((p > 1) | (p < 0))] <- -Inf
	}
        temp <- (1/zeta.fun(s)) * cumsum((1:1e+06)^(-s))
        temp.1 <- (1/zeta.fun(s))
        temp.max <- max(temp)
	if(any(p <= temp.1)){
		all.p[which((p <= temp.1)&(p >= 0))] <- 1
	}
	if(any(p==1)) all.p[which(p==1)] <- Inf
	if(any(p > temp.max)){
		x <- c(1e+06, 1e+300)
		y <- c(temp.max, 1)
		out.lm <- lm(y~I(-1/x))
		be <- out.lm$coef
		which.p <- which((p > temp.max) & (p <= 1))
		new.q <- -1/((p[which.p]-be[1])/be[2])
		if(any(new.q >= 1e+300 | new.q == -Inf)) all.p[which(new.q >= 1e+300)] = Inf
		if(any((new.q < 1e+300) & (new.q >= 1e+06))) all.p[which.p[which((new.q < 1e+300) & (new.q >= 1e+06))]] = floor(new.q[which((new.q < 1e+300) & (new.q >= 1e+06))])
	}
	if(any(is.na(all.p))){
		all.p[which(is.na(all.p))] <- sapply(1:length(which(is.na(all.p))), function(i) min(which(temp >= p[which(is.na(all.p))[i]])))
	}
    }
	if(any(all.p==-Inf)){
		all.p[which(all.p==-Inf)] <- NaN
	}
all.p
}

pzipfman <- function (q, s = 1, b = NULL, N = NULL, lower.tail = TRUE, log.p = FALSE)
{
	if(is.null(N)){
		stop(paste("Must specify N.", "\n"))
	}
	q <- floor(q)
	temp <- vector("list",length(q))
    if ((is.null(b)) & (N < Inf)) {
        if (s <= 0) {
            stop(paste("Invalid value for s!", "\n"))
        }
            if (any(q <= 0)) {
                temp[which(q <= 0)] <- 0
            }
		if (any(q > N)) {
                temp[which(q > N)] <- 1
            }
		if (any(q>0 & q<=N)){
			ind <- which(q>0 & q<=N)
			temp[ind] <- lapply(1:length(ind), function(i) dzipfman(x = 1:q[ind[i]], s = s, N = N))
		}
		} else if (!is.null(b)) {
        if (s <= 0 | b < 0) {
            stop(paste("Invalid value for s and/or b!", "\n"))
        }
        if (N == Inf) {
            stop(paste("N must be finite!", "\n"))
        }
            if (any(q <= 0)) {
                temp[which(q <= 0)] <- 0
            }
		if (any(q > N)) {
                temp[which(q > N)] <- 1
            }
		if (any(q>0 & q<=N)){
			ind <- which(q>0 & q<=N)
			temp[ind] <- lapply(1:length(ind), function(i) dzipfman(x = 1:q[ind[i]], s = s, b = b, N = N))
		}
		} else {
        if (s <= 1) {
            stop(paste("Invalid value for s!", "\n"))
        }
            if (any(q <= 0)) {
                temp[which(q <= 0)] <- 0
            }
		if (any(q == Inf)) {
                temp[which(q == Inf)] <- 1
            }
		if (any(q>0 & q<Inf)){
			ind <- which(q>0 & q<Inf)
			temp[ind] <- lapply(1:length(ind), function(i) dzipfman(x = 1:q[ind[i]], s = s, N = Inf))
		}
		}
	if(lower.tail==FALSE){
		temp <- sapply(1:length(temp), function(i) 1 - sum(signif(temp[[i]],14)))
		if(any(temp < 0)) temp[which(temp < 0)] <- 0
		if(any(temp > 1)) temp[which(temp > 1)] <- 1
		if(log.p) temp <- log(temp)
	} else if(log.p){
		temp <- sapply(1:length(temp), function(i) log(temp[[i]][1]) + log(1 + sum(temp[[i]][-1])/temp[[i]][1]))
	} else {
		temp <- sapply(temp,sum)
		if(any(temp < 0)) temp[which(temp < 0)] <- 0
		if(any(temp > 1)) temp[which(temp > 1)] <- 1
	}	
temp
}
