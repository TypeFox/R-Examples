suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

# Normalize phi.Obs
phi.Obs <- ex.train$phi.Obs / mean(ex.train$phi.Obs)

# Generate sequences.
da.roc <- simu.orf(length(phi.Obs), b.Init$roc,
                   phi.Obs = phi.Obs, model = "roc")
names(da.roc) <- names(phi.Obs)
write.seq(da.roc, "toy_roc.fasta")

# Read seqeuences back.
seq.roc <- read.seq("toy_roc.fasta")
seqstring.roc <- convert.seq.data.to.string(seq.roc)
phi <- data.frame(ORF = names(phi.Obs), phi.value = phi.Obs)

# Generate data structures from sequences.
aa.names <- names(b.Init$roc)
reu13.df <- gen.reu13.df(seqstring.roc, phi, aa.names = aa.names)
n <- gen.n(seqstring.roc, aa.names = aa.names)
y <- gen.y(seqstring.roc, aa.names = aa.names)

# Run codon fits.
.CF.AC$renew.iter <- 3
ret.time <- system.time({
  ret <- cubfits(reu13.df, ex.train$phi.Obs, y, n,
                 nIter = 20,
                 verbose = TRUE, report = 5,
                 model = "roc", adaptive = "simple")
})

x <- rowMeans(do.call("cbind", ret$phi.Mat)[, 11:20])
y <- ex.train$phi.Obs
plotprxy(x, y)

x <- log10(x / mean(x))
y <- log10(y / mean(y))
print(mean(x))
print(summary(lm(y ~ x))$r.squared)

