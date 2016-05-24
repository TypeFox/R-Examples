coeftofn <-
function (alpha, beta, n = 256, filter.number = 1,
	family = c("DaubExPhase", "DaubLeAsymm")) 
{
family <- match.arg(family)
#
# Turn alpha and beta wavelet coefficients into a function at
# resolution n using specific wavelet
#
#
# Pad coefficients to the right length 
lD <- n - 1
alpha <- c(rep(0, lD - length(alpha)), alpha)
beta <- c(rep(0, lD - length(beta)), beta)
#
# Create empty vessels for coefficients
#
bwd <- awd <- wd(rep(0, n), filter.number = filter.number, family = family)
#
# Store wavelet coefficients in vessels
#
awd$D <- alpha
bwd$D <- beta
#
# Create functions by inverting coefficients
#
alpha <- wr(awd)
beta <- wr(bwd)
#
# Build and return answer object
#
l <- list(alpha = alpha, beta = beta)
class(l) <- "csBiFunction"
return(l)
}
