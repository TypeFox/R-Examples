library(spgs)

#Set random number generator and initialise random seed
set.seed(331, kind="Mersenne-Twister", normal.kind="Inversion")

#Simulate a sequence of nucleotides as a Markov chain
s <- simulateMarkovChain(500, matrix(0.25, 4, 4), states=c("a", "c", "g", "t"))

#Circular tests
pc <- pair.counts(s)
tc <- triple.counts(s)
qc <- quadruple.counts(s)
cc4 <- cylinder.counts(s, 1:4)
c2l2c <- spgs:::cyl2lag2.counts(s, 10)
circ.res <- c(
	identical(apply(tc, 1:2, sum), apply(tc, 2:3, sum)),
	identical(apply(qc, 1:2, sum), apply(qc, 2:3, sum)),
	identical(apply(qc, 2:3, sum), apply(qc, 3:4, sum)),
	identical(pc, apply(cc4, 1:2, sum)),
	identical(pc, apply(cc4, 2:3, sum)),
	identical(pc, apply(cc4, 3:4, sum)),
	identical(tc, apply(cc4, 1:3, sum)),
	identical(tc, apply(cc4, 2:4, sum)),
	identical(qc, cc4),
	identical(qc, c2l2c[,,,,3])
)
print(circ.res)

#Non-circular tests
p <- pair.counts(s, circular=FALSE)
t <- triple.counts(s, circular=FALSE)
q <- quadruple.counts(s, circular=FALSE)
c4 <- cylinder.counts(s, 1:4, circular=FALSE)
c2l2 <- spgs:::cyl2lag2.counts(s, 10, circular=FALSE)

print(c(sum(p!=pc), sum(t!=tc), sum(q!=qc)))

noncirc.res <- c(
	identical(p, cylinder.counts(s, 1:2, circular=FALSE)),
	identical(t, cylinder.counts(s, 1:3, circular=FALSE)),
	identical(q, cylinder.counts(s, 1:4, circular=FALSE)),
	identical(q, c2l2[,,,,3])
)
print(noncirc.res)

