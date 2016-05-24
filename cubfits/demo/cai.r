suppressMessages(library(cubfits, quietly = TRUE))

y <- ex.train$y
y.list <- convert.y.to.list(y)
CAI <- calc_cai_values(y, y.list)$CAI
plot(CAI, log10(ex.train$phi.Obs), main = "Expression vs CAI",
     xlab = "CAI", ylab = "Expression (log10)")

### Verify with the seqinr example.
library(seqinr)
inputdatfile <- system.file("sequences/input.dat", package = "seqinr")
input <- read.fasta(file = inputdatfile, forceDNAtolower = FALSE)
names(input)[65] <- paste(names(input)[65], ".1", sep = "") # name duplicated.
input <- input[order(names(input))]

### Convert to cubfits format.
seq.string <- convert.seq.data.to.string(input)
new.y <- gen.y(seq.string)
new.y.list <- convert.y.to.list(new.y)
ret <- calc_cai_values(new.y, new.y.list)

### Rebuild w.
w <- rep(1, 64)
names(w) <- codon.low2up(rownames(caitab))
for(i in 1:64){
  id <- which(names(ret$w) == names(w)[i])
  if(length(id) == 1){
    w[i] <- ret$w[id]
  }
}
CAI.res <- sapply(input, seqinr::cai, w = w)

### Plot.
plot(CAI.res, ret$CAI,
     main = "Comparison of seqinR and cubfits results",
     xlab = "CAI from seqinR", ylab = "CAI from cubfits", las = 1)
abline(c(0, 1))
