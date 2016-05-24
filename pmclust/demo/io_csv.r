### Setup environment.
suppressMessages(library(pmclust, quietly = TRUE))
comm.set.seed(123, diff = TRUE)
if(comm.size() != 4){
  comm.stop("4 processors are needed.")
}

### Generate an example data.
N <- 400
p <- 2
K <- 2
data.spmd <- generate.MixSim(N, p, K, MaxOmega = 0.001)

### Dump fake data to a csv file.
X.df.org <- as.data.frame(cbind(data.spmd$X.spmd, data.spmd$CLASS.spmd))
colnames(X.df.org) <- c(paste("V", 1:p, sep = ""), "ID")
comm.write.csv(X.df.org, file = "toys_org.csv", row.names = FALSE)

### Read data from the csv file.
X.df.new <- comm.read.csv("toys_org.csv")
X.spmd <- as.matrix(X.df.new[, 1:p])

### Run clustering.
PARAM.org <- set.global(K = K)          # Set global storages.
PARAM.org <- initial.em(PARAM.org)      # One initial.
PARAM.new <- apecma.step(PARAM.org)     # Run APECMa.
em.update.class()                       # Get classification.

### Dump data and clustering results to a new csv file.
X.df.new$ID.est <- .pmclustEnv$CLASS.spmd
comm.write.csv(X.df.new, file = "toys_new.csv", row.names = FALSE)

### Quit.
finalize()

