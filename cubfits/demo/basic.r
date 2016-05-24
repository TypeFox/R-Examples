suppressMessages(library(cubfits, quietly = TRUE))

seq.data <- read.seq(get.expath("seq_200.fasta"))
phi.df <- read.phi.df(get.expath("phi_200.tsv"))
aa.names <- c("A", "C", "D")

# Read in from FASTA file.
seq.string <- convert.seq.data.to.string(seq.data)
reu13.df <- gen.reu13.df(seq.string, phi.df, aa.names)
reu13.list.new <- gen.reu13.list(seq.string, aa.names)
y <- gen.y(seq.string, aa.names)
head(y[[1]])
n <- gen.n(seq.string, aa.names)
head(n[[1]])
scuo <- gen.scuo(seq.string, aa.names)
head(scuo)

# Convert to list format.
reu13.list <- convert.reu13.df.to.list(reu13.df)
head(reu13.list[[1]])
y.list <- convert.y.to.list(y)
head(y.list[[1]])
n.list <- convert.n.to.list(n)
head(n.list[[1]])
