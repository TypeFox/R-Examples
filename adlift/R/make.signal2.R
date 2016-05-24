"make.signal2" <-
function(name, x, snr = Inf, ...)
{

n<-length(x)

z <- switch(name,
	ppoly =	{
        	y <- rep(0, n)
                xsv <- (x <= 0.5)
                y[xsv] <- -16 * x[xsv]^3 + 12. * x[xsv]^2
                xsv <- (x > 0.5) & (x <= 0.75)
                y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 40 * x[xsv] + 28))/3 - 1.5
                xsv <- x > 0.75
                y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 32 * x[xsv] + 16))/3
        	y
	},
	dirac = n * (x == floor(0.37 * n)/n),
	kronecker = (x == floor(0.37 * n)/n),
	heavisine = 4 * sin(4 * pi * x) - sign(x - 0.3) - sign(0.72 - x),
	bumps = {
		pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
		hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
		wth <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
		y <- rep(0, n)
		for(j in 1:length(pos))
			y <- y + hgt[j]/(1 + abs((x - pos[j]))/wth[j])^4
		y
	},
	blocks = {
		pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
		hgt <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
		y <- rep(0, n)
		for(j in 1:length(pos))
			y <- y + ((1 + sign(x - pos[j])) * hgt[j])/2
		y
	},
	doppler = sqrt(x * (1 - x)) * sin((2 * pi * 1.05)/(x + 0.05)),
	ramp = x - (x >= 0.37),
	cusp = sqrt(abs(x - 0.37)),
	crease = exp(-4 * abs(x - 0.5)),
	sing = 1/abs(x - (floor(n * 0.37) + 0.5)/n),
	hisine = sin(pi * n * 0.6902 * x),
	midsine = sin(pi * n * 0.3333 * x),
	losine = sin(pi * n * 0.03 * x),
	linchirp = sin(0.125 * pi * n * x^2),
	twochirp = sin(pi * n * x^2) + sin((pi/3) * n * x^2),
	quadchirp = sin((pi/3) * n * x^3),
	# QuadChirp + LinChirp + HiSine
	mishmash1 = sin((pi/3) * n * x^3) + sin(pi * n * 0.6902 * x) + sin(pi * n * 0.125 * x^2),
	# QuadChirp + LinChirp + HiSine + Bumps
	mishmash2 = {
		# wernersorrows
		y <- sin(pi * (n/2) * x^3) + sin(pi * n * 0.6902 * x) + sin(pi * n * x^2)
		pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
		hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
		wth <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
		for(j in 1:length(pos))
			y <- y + hgt[j]/(1 + abs((x - pos[j])/wth[j]))^4
		y
	},
	# QuadChirp + MidSine + LoSine + Sing/200
	mishmash3 = sin((pi/3) * n * x^3) + sin(pi * n * 0.3333 * x) + sin(pi * n * 0.03 * x) + (1/abs(x - (floor(n * 0.37) + 0.5)/n))/((200 * n)/512),
	gauss = dnorm(x, 0.3, 0.025),
	jumpsine = 10 * (sin(4 * pi * x) + as.numeric(x >= 0.625 & x < 0.875)),
	levelshift = as.numeric(x >= 0.25 & x < 0.39),
	patches = {
		if(n < 16)
			stop("n must be >= 16 to generate patches\n")
		J <- floor(log(n, 2))
		y <- rep(0, n)
		for(j in 0:(J - 4))
			y[(1:2^j) + 3 * 2^(j + 2)] <- 1
		y
	},
	linear = 2 * x - 1,
	quadratic = 4 * (1 - x) * x,
	cubic = (64 * x * (x - 1) * (x - 0.5))/3,
	stop("Unknown signal name.")
)

if(snr > 0) z <- z + (rnorm(n) * sqrt(var(z)))/snr		

z
}
