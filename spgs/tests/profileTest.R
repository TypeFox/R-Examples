#ProfileTest.R

library(spgs)

#Set random number generator and initialise random seed
set.seed(331, kind="Mersenne-Twister", normal.kind="Inversion")

#Simulate a sequence of nucleotides as a Markov chain
len <- 200 #length of sequence to simulate
s <- simulateMarkovChain(len, matrix(0.25, 4, 4), states=c("a", "c", "g", "t"))

#Add some ambiguous symbols and uppercase symbols to s
idx <- sample(len, len/4)
s[idx] <- toupper(s[idx])
idx <- sample(len, len/8)
ambiguousSymbols <- c("b", "B", "d", "D", "h", "H", "v", "V", "s", "S", "w", "W", "m", "M", "k", "K", "r", "R", "y", "Y", "n", "N", "x", "X", "-")
s[idx] <- sample(ambiguousSymbols, len/8, replace=TRUE)

#Test OligoProfile class
op <- oligoProfile(s,2, name="Markov chain", plot=FALSE)
print(op)
op <- oligoProfile(s,2, case="u", name="Markov chain", plot=FALSE)
print(op)
op <- oligoProfile(s,2, case="a", name="Markov chain", plot=FALSE)
print(op)
op <- oligoProfile(s,2, disambiguate=FALSE, name="Markov chain", plot=FALSE)
print(op)
