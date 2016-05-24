zm.ll <- function (x, N=NULL, s = 1, b = 1, dist = c("Zipf", "Zipf-Man", "Zeta"), ...) {	
dist <- match.arg(dist)
if(class(x)!="table"){ 
	x <- table(x)
	}
	names(x)=1:length(x)
	x.labs <- as.numeric(names(x))
	N.temp <- max(x.labs)
	x <- x[order(as.numeric(names(x)))]
	if(dist == "Zeta"){
		N <- N.temp
		}
	if(is.null(N)){
		N <- N.temp
	}
	if(N<N.temp) stop(paste("N cannot be smaller than the maximum number of categories in x!",
		"\n"))
	N.seq <- 1:N
	x <- c(x, rep(0,N-length(x)))
	if(dist=="Zipf"){
		ll.zipf <- function(s) sum(x*(s*log(N.seq)+log(sum(1/(N.seq)^s))))
		fit <- suppressWarnings(stats4::mle(ll.zipf,start=list(s=s),lower=0,...))
	}
	if(dist=="Zipf-Man"){
		ll.zima <- function(s,b) sum(x*(s*log(N.seq+b)+log(sum(1/(N.seq+b)^s))))
		fit <- suppressWarnings(stats4::mle(ll.zima,start=list(s=s,b=b),lower=c(0,0),...))
	}
	if(dist=="Zeta"){
		ll.zeta <- function(s) sum(x*(s*log(N.seq)+log(zeta.fun(s))))
		fit <- suppressWarnings(stats4::mle(ll.zeta,start=list(s=s),lower=1+1e-14,...))
	}
	fit
	}