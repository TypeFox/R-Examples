library(psd)

data(magnet)
data(Tohoku)

x <- subset(Tohoku, epoch=="preseismic")$areal * 1e-6
n <- length(x)

# detrend and demean
t. <- 1L:n - (n + 1)/2
sumt2. <- n * (n**2 - 1)/12
x <- x - mean(x) - sum(x * t.) * t./sumt2.

message("\t-->\tProject Magnet example:")

pm1 <- psdcore(magnet$clean, verbose=TRUE)
pmn <- spectrum(magnet$clean, plot=FALSE)
try(with(normalize(pmn, src='spectrum'), plot(freq, dB(spec), type='l')))
try(with(pm1, lines(freq, dB(spec), col='blue')))

print(summary(dB(pmn$spec) - dB(pm1$spec)))

message("\t-->\tRandom-noise example:")

plot(p1 <- psdcore(rnorm(n), verbose=TRUE), log='dB')

message("\t-->\tTohoku example:")

pmn <- spectrum(x, plot=FALSE)
try(plot(normalize(pmn,src = 'spectrum'), log='dB'))

p1 <- psdcore(x, verbose=TRUE)
try(with(p1, lines(freq, dB(spec), col='blue')))

p2 <- psdcore(x, ntaper = 10, verbose=TRUE)
try(with(p2, lines(freq, dB(spec), col='red')))

message("\t-->\tdone.")
