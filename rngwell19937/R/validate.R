well19937.validate <- function()
{
	result <- function(x) {
		if (x) { "OK" } else { "FAILED" }
	}
	if (RNGkind(NULL)[1] != "user-supplied") {
		cat("Restoring user-supplied generator\n")
		RNGkind("user-supplied")
	}
	# SFMT 32
	set.initialization("sfmt")
	set.resolution(32)
	set.seed(123456789)
	aux <- runif(1000000)
	d1 <- 2^32 * runif(10) - c(136727666, 2709878109, 787020490, 2005653182,
		361632219, 4120149904, 383905305, 3321628803, 4029726286, 2820685373)
	r1 <- all(d1 == 0)
	cat("initialization     sfmt, resolution 32:", result(r1),"\n")
	# SFMT 53
	set.resolution(53)
	set.seed(234567891)
	aux <- runif(1000000)
	d2 <- 2^53 * runif(7) - c(5375099290856534, 4584046116196875, 2797402997351966,
		2337032758473196, 8713402114218560, 4133043146500834, 3229315589793403)
	r2 <- all(d2 == 0)
	cat("initialization     sfmt, resolution 53:", result(r2),"\n")
	# MRG32k5a 32
	set.initialization("mrg32k5a")
	set.resolution(32)
	set.seed(123456789)
	aux <- runif(1000000)
	d3 <- 2^32 * runif(10) - c(1134116397, 2841234926, 4186228405, 3779293683,
		3795008926, 1981908955, 3867072533, 2562139842, 1586001182, 25052338)
	r3 <- all(d3 == 0)
	cat("initialization mrg32k5a, resolution 32:", result(r3),"\n")
	# MRG32k5a 53
	set.resolution(53)
	set.seed(234567891)
	aux <- runif(1000000)
	d4 <- 2^53 * runif(7) - c(5235163687128757, 2097982372150467, 7765147415699790,
		3158530122713276, 6281840767500381, 5274792491989681, 5008469347891762)
	r4 <- all(d4 == 0)
	cat("initialization mrg32k5a, resolution 53:", result(r4),"\n")
	# MRG32k5a 32 vector 1
	set.resolution(32)
	set.vector.seed(123456789)
	aux <- runif(1000000)
	d51 <- 2^32 * runif(10) - c(1134116397, 2841234926, 4186228405, 3779293683,
		3795008926, 1981908955, 3867072533, 2562139842, 1586001182, 25052338)
	# MRG32k5a 32 vector 2
	set.resolution(32)
	set.vector.seed(1:5)
	aux <- runif(1000000)
	d52 <- 2^32 * runif(10) - c(3532849506, 2678408862, 2385013101, 2189132525,
		1546434356, 3862404650, 2891933366, 1569423855, 256923335, 2853845277)
	r5 <- all(d51 == 0) & all(d52 == 0)
	cat("initialization   vector, resolution 32:", result(r5),"\n")
	# MRG32k5a 53 vector 1
	set.resolution(53)
	set.vector.seed(234567891)
	aux <- runif(1000000)
	d61 <- 2^53 * runif(7) - c(5235163687128757, 2097982372150467, 7765147415699790,
		3158530122713276, 6281840767500381, 5274792491989681, 5008469347891762)
	# MRG32k5a 53 vector 2
	set.resolution(53)
	set.vector.seed(2:6)
	aux <- runif(1000000)
	d62 <- 2^53 * runif(7) - c(5093507248644797, 8528502054353349, 7732655783590473,
		329624472083424.0, 7287030050425235, 4690770001917934, 1579164290395871)
	r6 <- all(d61 == 0) & all(d62 == 0)
	cat("initialization   vector, resolution 53:", result(r6),"\n")
	if (!(r1 & r2 & r3 & r4 & r5 & r6)) stop("validation of rngwell19937 FAILED")
	invisible(NULL)
}

