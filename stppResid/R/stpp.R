stpp <- function(x, y, t, stw)
{
	if(missing(x) || missing(y) || missing(t))
		stop("x, y, and t must be specified")
	if(missing(stw)) {
		stw <- stwin()
		cat("stw missing - using default space-time window [0, 1]x[0, 1]\n")
	}
	if(!is.stwin(stw))
		stop("stw must be an object of class stwin")
	x.out <- which((x < stw$xcoord[1]) | (x > stw$xcoord[2]))
	y.out <- which((y < stw$ycoord[1]) | (y > stw$ycoord[2]))
	t.out <- which((t < stw$tcoord[1]) | (t > stw$tcoord[2]))
	all.out <- unique(c(x.out, y.out, t.out))
	if(length(all.out) > 0) {
		cat("Warning:", length(all.out), " point(s) are outside of the space-time window and will be removed\n")
		x <- x[-all.out]
		y <- y[-all.out]
		t <- t[-all.out]
	}
	if((length(x) != length(y)) || (length(x) != length(t)) || (length(y) != length(t)))
		stop("x, y, and t must be of same length")
	if(!is.numeric(x) || !is.numeric(y) || !is.numeric(t))
		stop("x, y, and t must be numeric vectors")
  if(length(unique(t)) != length(t))
    stop("t values must be unique")
	x <- x[order(t)]
	y <- y[order(t)]
	t <- sort(t)
	X <- list(x = x, y = y, t = t)
	X <- c(X, stw)	
	class(X) <- "stpp"
	return(X)
}