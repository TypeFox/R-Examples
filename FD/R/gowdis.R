`gowdis` <-
function(x, w, asym.bin = NULL, ord = c("podani", "metric", "classic") ){

if (length(dx <- dim(x)) != 2 || !(is.data.frame(x) || is.numeric(x))) stop("x is not a dataframe or a numeric matrix\n")

# n = number of rows, p = number of variables
n <- dx[1]
p <- dx[2]

ord <- match.arg(ord)

varnames <- dimnames(x)[[2]]

# check for weight vector, add equal weights if missing
if (!missing(w)){

	# check if correct class and length
	if (length(w) != p | !is.numeric(w) ) stop("w needs to be a numeric vector of length = number of variables in x\n")
	if (all(w == 0) ) stop("Cannot have only 0's in 'w'\n")
	w <- w / sum(w)
	}
else w <- rep(1, p) / sum(rep(1, p))


if (is.data.frame(x)) {
        type <- sapply(x, data.class)
        }
    else {
        type <- rep("numeric", p)
        names(type) <- colnames(x)
	}

# replace character variables by factors
if (any(type == "character") ) for (i in 1:p) if (type[i] == "character") x[,i] <- as.factor(x[,i])

# check for binary variables
is.bin <- function(k) all(k[!is.na(k)] %in% c(0,1))

bin.var <- rep(NA,p); names(bin.var) <- varnames
for (i in 1:p) bin.var[i] <- is.bin(x[,i])

if (any(type[bin.var] != "numeric")) stop("Binary variables should be of class 'numeric'\n")

type[type %in% c("numeric", "integer")] <- 1
type[type == "ordered"] <- 2
type[type %in% c("factor", "character")] <- 3
type[bin.var] <- 4

# convert asymmetric binary variables, if present
if (!is.null(asym.bin) ){
	if (!all(bin.var[asym.bin])) stop("Asymetric binary variables must only contain 0 or 1\n")
	else type[asym.bin] <- 5
	}
	
type <- as.numeric(type)

# convert factors to their internal numeric codes
x <- data.matrix(x)

# convert ordinal variables to ranks, following Podani (1999)
# or when ord = "classic", convert to numeric
if (any(type == 2) ) {
   if (ord != "classic")   for (i in 1:p) if (type[i] == 2) x[,i] <- rank(x[,i], na.last = "keep")
   else   for (i in 1:p) if (type[i] == 2) x[,i] <- as.numeric(x[,i])
}

# compute the range of each variable (this will only be used for numeric and ordinal variables)
range.Data <- function(v){
	r.Data <- range(v, na.rm = T)
	res <- r.Data[2]-r.Data[1]
	return(res)
	}

range2<- apply(x, 2, range.Data)

# compute Timax, and Timin for each variable (these will only apply to ordinal variables, see Podani [1999], eq. 2b.)
comp.Timax <- function(v){
	Ti.max <- max(v, na.rm = T)
	no.na <- v[!is.na(v)]
	res <- length(no.na[no.na == Ti.max])
	return(res)
	}

Timax <- apply(x, 2, comp.Timax)
	
comp.Timin <- function(v){
	Ti.min <- min(v, na.rm = T)
	no.na <- v[!is.na(v)]
	res <- length(no.na[no.na == Ti.min])
	return(res)
	}

Timin <- apply(x, 2, comp.Timin)

if (ord == "podani") pod <- 1 else pod <- 2

res <- .C("gowdis", as.double(x), as.double(w), as.integer(type), as.integer(n), as.integer(p), as.double(range2), as.integer(pod), as.double(Timax), as.double(Timin), res = double(n*(n-1)/2), NAOK = T, PACKAGE = "FD")$res

type[type == 1] <- "C"
type[type == 2] <- "O"
type[type == 3] <- "N"
type[type == 4] <- "B"
type[type == 5] <- "A"

if (any(is.na(res) ) ) attr(res, "NA.message") <- "NA's in the dissimilarity matrix!"
attr(res, "Labels") <- dimnames(x)[[1]]
attr(res, "Size") <- n
attr(res, "Metric") <- "Gower"
attr(res, "Types") <- type
class(res) <- "dist"

return(res)
}

