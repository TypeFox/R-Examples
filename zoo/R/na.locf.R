na.locf <- function(object, na.rm = TRUE, ...)
	UseMethod("na.locf")

na.locf.default <- function(object, na.rm = TRUE, fromLast, rev, maxgap = Inf, rule = 2, ...) {

	L <- list(...)
	if ("x" %in% names(L) || "xout" %in% names(L)) {

		if (!missing(fromLast)) {
			stop("fromLast not supported if x or xout is specified")
		}
		return(na.approx(object, na.rm = na.rm, 
			maxgap = maxgap, method = "constant", rule = rule, ...))
	}

	na.locf.0 <- function(x) {
	      L <- !is.na(x)
	      idx <- if (fromLast)
	         rev(c(NA,rev(which(L)))[cumsum(rev(L))+1])
              else
	         c(NA,which(L))[cumsum(L)+1]
	      # na.index(x,i) returns x[i] except if i[j] is NA then
	      # x[i[j]] is NA too
	      na.index <- function(x, i) {
		L <- !is.na(i)
		x[!L] <- NA
		x[L] <- coredata(x)[i[L]]
	  	x
	      }
	      xf <- na.index(x, idx)
              .fill_short_gaps(x, xf, maxgap = maxgap)
	}
   	if (!missing(rev)) {
	   warning("na.locf.default: rev= deprecated. Use fromLast= instead.")
	   if (missing(fromLast)) fromLast <- rev
	} else if (missing(fromLast)) fromLast <- FALSE
	rev <- base::rev
	object[] <- if (length(dim(object)) == 0)
		na.locf.0(object)
	else
		apply(object, length(dim(object)), na.locf.0)
	if (na.rm) na.trim(object, is.na = "all") else object
}

na.contiguous.data.frame <-
na.contiguous.zoo <- function(object, ...) 
{
    if (length(dim(object)) == 2) 
        good <- apply(!is.na(object), 1, all)
    else good <- !is.na(object)
    if (!sum(good)) 
        stop("all times contain an NA")
    tt <- cumsum(!good)
    ln <- sapply(0:max(tt), function(i) sum(tt == i))
    seg <- (seq_along(ln)[ln == max(ln)])[1] - 1
    keep <- (tt == seg)
    st <- min(which(keep))
    if (!good[st]) 
        st <- st + 1
    en <- max(which(keep))
    omit <- integer(0)
    n <- NROW(object)
    if (st > 1) 
        omit <- c(omit, 1:(st - 1))
    if (en < n) 
        omit <- c(omit, (en + 1):n)
    cl <- class(object)
    if (length(omit)) {
        object <- if (length(dim(object))) 
            object[st:en, ]
        else object[st:en]
        attr(omit, "class") <- "omit"
        attr(object, "na.action") <- omit
        if (!is.null(cl)) 
            class(object) <- cl
    }
    object
}

na.contiguous.list <- function(object, ...)
	lapply(object, na.contiguous)
