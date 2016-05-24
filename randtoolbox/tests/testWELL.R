library(randtoolbox)

params <- c(
	"512a",
	"521a",
	"521b",
	"607a",
	"607b",
	"800a",
	"800b",
	"1024a",
	"1024b",
	"19937a",
	"19937b",
	"19937c",
	"21701a",
	"23209a",
	"23209b",
	"44497a",
	"44497b")

result1 <- rep(NA, times=length(params))
result2 <- rep(NA, times=length(params))

seed <- floor(2^31*runif(1))
cat("using seed", seed, "for test of the output of WELL RNG\n")

m <- 100
cat("generating sequences of the length", m, "from each generator\n")

for (i in 1:length(params))
{
	cat(i, "")
	set.generator("WELL", version=params[i], seed=seed)
	s0 <- getWELLState()
	x <- runif(m)
	out <- rngWELLScriptR(m, s0, params[i], includeState=TRUE)
	s1 <- getWELLState()
	result1[i] <- all(x == out$x)
	result2[i] <- all(s1 == out$state)
}

cat("\n\n")
print(data.frame(params, result1, result2))

stopifnot(all(result1, result2))

