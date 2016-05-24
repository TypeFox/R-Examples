genwwn.thpower <-
function(N=128, ar=NULL, ma=NULL, plot.it=FALSE, sigsq=1, alpha=0.05,
	away.from="standard", filter.number=10, family="DaubExPhase",
	verbose=FALSE){

if (is.na(J.N <- IsPowerOfTwo(N)))
        stop("Data has to be of length a power of two")
    if (N < 16)
        stop("Need time series to be of length 16 or greater")
#
# Away.from can be an integer (keep away from the away.from finest scales,
# or the text string "standard" where we internally set the away.from based
# on the length of the data
#
    if (!is.numeric(away.from)) {
        if (away.from != "standard")
                stop("away.from argument has to be positive integer or the string `standard'")
        if (N==16 || N==32)
                away.from <- 2
        else if (N > 32 && N <= 1024)
                away.from <- J.N - 4
        else if (N > 1024)
                away.from <- J.N - 5

        if (verbose==TRUE)
                cat("N is: ", N, " J.N is: ", J.N, " away.from is: ", away.from, "\n")
        }


if (!is.null(ar))
	ar.coef <- c(1, -ar)
else
	ar.coef <- 1

if (!is.null(ma))
	ma.coef <- c(1, ma)
else
	ma.coef <- 1

ar.fn <- as.function(polynomial(ar.coef))
ma.fn <- as.function(polynomial(ma.coef))

spec.fn <- function(w, ar.fn, ma.fn, sigsq)	{

	Cf <- exp(complex(real=0, imaginary=-w))

	top <- Mod(ma.fn(Cf))^2
	bot <- Mod(ar.fn(Cf))^2

	spec <- sigsq*top/(2*pi*bot)

	return(spec)
	}
	

N <- N/2

wp <- 2*pi*(1:N)/(2*N)

spec <- spec.fn(w=wp, ar.fn=ar.fn, ma.fn=ma.fn, sigsq=sigsq)

the.var <- 2*integrate(spec.fn, lower=0, upper=pi, ar.fn=ar.fn, ma.fn=ma.fn,
	sigsq=sigsq)$value

if (plot.it==TRUE)	{
	plot(wp, spec, type="l", ylim=c(0, max(spec)), main="ARMA spectrum",
		xlab="Frequency", ylab="Spectrum")
	scan()
	}


norspec <- 2*pi*spec/the.var

norspecwd <- wd(norspec, filter.number=filter.number, family=family)

# New
norspecvarwd <- sqwd(norspec^2, filter.number=filter.number, family=family,
	type="wavelet", m0=nlevelsWT(norspecwd))

if (plot.it==TRUE)	{
	plot(norspecwd, main="Wavelet coefficients of normalized spectrum")
	scan()
	}

#
# Collect wavelet coefficients
#

all.hwc <- NULL
# New
all.sdwc <- NULL

for(j in (nlevelsWT(norspecwd)-1-away.from):0)	{
	all.hwc <- c(all.hwc, accessD(norspecwd, level=j))
	all.sdwc <- c(all.sdwc, sqrt(accessD(norspecvarwd, level=j))) # New
	}

ncoef <- length(all.hwc)

if (ncoef != length(all.sdwc))
	stop("SDs vector is different length to hwc vector")

alpha.c <- alpha/ncoef

C.alpha.c <- qnorm( 1 - alpha.c/2)


# Old
#vec.for.prod <- pnorm(C.alpha.c-all.hwc) - pnorm(-C.alpha.c-all.hwc)
# New
vec.for.prod <- pnorm((C.alpha.c-all.hwc)/all.sdwc) - pnorm((-C.alpha.c-all.hwc)/all.sdwc)

th.power <- 1 - prod(vec.for.prod)

ll <- list(C.alpha.c=C.alpha.c, th.power=th.power, norspecwd=norspecwd,
	norspecvarwd=norspecvarwd, all.hwc=all.hwc, all.sdwc=all.sdwc)

return(ll)




}
