expfourier_gif <-
function(data, evalpts, params, TT = NA, tplus = FALSE)
{
	# params=c(p, constant, a, b)
	# a and b must be of the same length
	P <- params[1]
	constant <- params[2]
	A <- params[3:(length(params)/2 + 1)]
	B <- params[(length(params)/2 + 2):(length(params))]	
	Z1 <- length(A)
	if(any(is.na(TT))) {
		if (is.vector(evalpts)) Time <- evalpts
			else Time <- evalpts[, "time"]
		times <- cos((pi/P) * (outer(Time, seq(2, 2 * Z1, 2))))
		timec <- times %*% A
		times <- sin((pi/P) * (outer(Time, seq(2, 2 * Z1, 2))))
		times <- times %*% B
		ci <- exp(constant + times + timec)
	}
	else {
		lambda <- function(Time, z1, A, B, constant, p)
		{
			times <- cos((pi/p) * (outer(Time, seq(2, 2 * z1, 2))))
			timec <- times %*% A
			times <- sin((pi/p) * (outer(Time, seq(2, 2 * z1, 2))))
			times <- times %*% B
			as.numeric(exp(constant + timec + times))
		}
		ci <- integrate(lambda, TT[1], TT[2], z1 = Z1, A = A, B = B, 
			constant = constant, p = P)
		if (ci$message == "OK")
		  ci <- ci$value
		else
		  stop(paste("Problems with Numerical Integration: ", ci$message))
	}
	ci <- as.vector(ci)
	return(ci)
}
attr(expfourier_gif, "rate") <- "bounded"


exppoly_gif <-
function(data, evalpts, params, TT = NA, tplus = FALSE)
{
	B <- as.matrix(params)
	if(any(is.na(TT))) {
		N <- length(B)
		if (is.vector(evalpts)) eval.times <- evalpts
			else eval.times <- evalpts[, "time"]
		ci <- exp(outer(eval.times, (0:(N - 1)), "^") %*% B)
	}
	else {
		#The trivial case (I = 0).
		if(nrow(B) == 1) {
			ci <- (TT[2] - TT[1]) * exp(B)
		}
		# Case where I = 1.
		else {
			if(nrow(B) == 2) {
				ci <- (exp(B[1,  ]) * (exp(B[2,  ] * TT[2]) - 
				  exp(B[2,  ] * TT[1])))/B[2,  ]
			}
			#If I > 1, need to numerically integrate.
			else {
				lambda <- function(times, b)
				{
				  N <- length(b) - 1
				  as.numeric(exp(outer(times, seq(0, N), "^")
					%*% as.matrix(b)))
				}
				ci <- integrate(lambda, TT[1], TT[2], b = 
				  params)
				if (ci$message == "OK")
				  ci <- ci$value
				else
				  stop(paste("Problems with Numerical Integration: ", ci$message))
			}
		}
	}
	ci <- as.vector(ci)
	return(ci)
}
attr(exppoly_gif, "rate") <- "bounded"


fourier_gif <-
function(data, evalpts, params, TT = NA, tplus = FALSE)
{
	# params=c(p, constant, a, b)
	# a and b must be of the same length
	P <- params[1]
	constant <- params[2]
	B <- matrix(params[(length(params)/2 + 2):(length(params))], ncol = 1)
	A <- matrix(params[3:(length(params)/2 + 1)], ncol = 1)
	Z <- length(A)
	if(any(is.na(TT))) {
		if (is.vector(evalpts)) Time <- evalpts
			else Time <- evalpts[, "time"]
		times <- cos((pi/P) * (outer(Time, seq(2, 2 * Z, 2))))
		timec <- times %*% A
		times <- sin((pi/P) * (outer(Time, seq(2, 2 * Z, 2))))
		times <- times %*% B
		ci <- (constant + times + timec)
	}
	else {
		times <- sin((pi/P) * (outer(TT[2], seq(2, 2 * Z, 2)))) * (P/(2 * 
			pi * seq(1, Z))) - sin((pi/P) * (outer(TT[1], seq(2, 2 * 
			Z, 2)))) * (P/(2 * pi * seq(1, Z)))
		timec <- times %*% A
		times <-  - cos((pi/P) * (outer(TT[2], seq(2, 2 * Z, 2)))) * (P/(
			2 * pi * seq(1, Z))) + cos((pi/P) * (outer(TT[1], seq(2, 
			2 * Z, 2)))) * (P/(2 * pi * seq(1, Z)))
		times <- times %*% B
		constant <- constant * (TT[2] - TT[1])
		ci <- (constant + times + timec)
		if(ci < 0)
			ci <- NA
	}
	ci <- as.vector(ci)
	return(ci)
}
attr(fourier_gif, "rate") <- "bounded"


poly_gif <-
function(data, evalpts, params, TT = NA, tplus = FALSE)
{
	B <- as.matrix(params)
	N <- length(B)
	if(any(is.na(TT))) {
		if (is.vector(evalpts)) eval.times <- evalpts
			else eval.times <- evalpts[, "time"]
		ci <- outer(eval.times, 0:(N-1), "^")%*%B
	}
	else {
		ci <- sum((TT[2]^seq(1, N) - TT[1]^seq(1, N))*B/seq(1, N))
		if (ci<0) ci <- NA
	}
	names(ci) <- NULL
	ci <- as.vector(ci)
	return(ci)
}
attr(poly_gif, "rate") <- "bounded"


simple_gif <-
function(data, evalpts, params, TT = NA, tplus = FALSE)
{
	nms <- names(params)
	if(is.null(nms)) {
		fullnames <- c("a", "b", "g")
		names(params) <- fullnames
	}
	A <- params["a"]
	B <- params["b"]
	G <- params["g"]
	if(any(is.na(TT))) {
		if (is.vector(evalpts)) eval.times <- evalpts
			else eval.times <- evalpts[, "time"]
		ci <- A + B * (eval.times^G)
	}
	else {
		ci <- A * (TT[2] - TT[1]) + 
			B/(G+1)*(TT[2]^(G + 1) - TT[1]^(G + 1))
		if (ci<0) ci <- NA
	}
	names(ci) <- NULL
	return(ci)
}
attr(simple_gif, "rate") <- "bounded"



