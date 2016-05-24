ReadOtuTable <- function(fn) read.table(fn, header=T, row.names=1, check.names=F)

ReadPairTable <- function(fn) {
  # read a samples table. headers are: 'name', 'before', 'after'
  # before and after are the labels of the relevant columns in the otu table
  read.table(fn, header=T, colClasses=c("character", "character", "character"))
}

ReadQuadTable <- function(fn) {
  # read a quads table. headers are: 'name', 'control', 'treatment'
  # control and treatment are the names of samples in the quad
  read.table(fn, header=T, colClasses=c("character", "character", "character"))
}

Experiment <- function(otu.table, pair.table, quad.table) {
  # read in the pair, quad, and otu tables. make the samples and fit them.

  # read in the tables if filenames are given
  if (class(pair.table) == "character") pair.table <- ReadPairTable(pair.table)
  if (class(quad.table) == "character") quad.table <- ReadQuadTable(quad.table)

  # check that the pairs in the quad table are all in the pairs table
  check.pair <- function(name) {
    if (!(name %in% pair.table$name)) {
      stop(sprintf("the pair named '%s' is in the quads list but not the pairs list", name))
    }
  }
  for (i in 1:nrow(quad.table)) {
    row <- quad.table[i,]
    check.pair(row$control)
    check.pair(row$treatment)
  }

  # read in the otu table
  if (class(otu.table) == "character") otu.table <- ReadOtuTable(otu.table)

  # check that all the pair names are in the otu table
  check <- function(name) {
    if (!(name %in% names(otu.table))) {
      stop(sprintf("the pair '%s' was in the pairs file but not the otu table"))
    }
  }
  for (i in 1:nrow(pair.table)) {
    row <- pair.table[i,]
    check(row$before)
    check(row$after)
  }

  # create a pairs list
  pairs <- list()
  for (i in 1:nrow(pair.table)) {
    row <- pair.table[i,]
    before <- Sample(otu.table[[row$before]], row$before)
    after <- Sample(otu.table[[row$after]], row$after)
    p <- SamplePair(before, after, row$name)
    p$treatment <- row$treatment
    pairs[[row$name]] <- p
  }

  # create a quads list
  quads <- list()
  for (i in 1:nrow(quad.table)) {
    row <- quad.table[i,]
    q <- SampleQuad(pairs[[row$control]], pairs[[row$treatment]], name=row$name)
    quads[[row$name]] <- q
  }

  # package samples and quads into an experiment object
  expt <- list(otu.table=otu.table, quads=quads, pairs=pairs)
  class(expt) <- "texmex.experiment"
  expt
}

ExtractQuadData <- function(expt, quad.name) {
  # extract relevant information on a quad from the otu table and
  # computed z's

  quad <- expt$quads[[quad.name]]

  ct.t0 <- Sample(quad$control$name0, expt$otu.table)
  ct.t1 <- Sample(quad$control$name1, expt$otu.table)
  tx.t0 <- Sample(quad$treatment$name0, expt$otu.table)
  tx.t1 <- Sample(quad$treatment$name1, expt$otu.table)

  # compute log fold changes
  ct.lfc <- log(ct.t1$ra) - log(ct.t0$ra)
  tx.lfc <- log(tx.t1$ra) - log(tx.t0$ra)

  # store z scores
  ct.t0.z <- quad$control$fit0$z
  ct.t1.z <- quad$control$fit1$z
  tx.t0.z <- quad$treatment$fit0$z
  tx.t1.z <- quad$treatment$fit1$z

  # compute changes in z score
  ct.dz <- ct.t1.z - ct.t0.z
  tx.dz <- tx.t1.z - tx.t0.z

  # store F scores
  ct.t0.F <- quad$control$fit0$F
  ct.t1.F <- quad$control$fit1$F
  tx.t0.F <- quad$treatment$fit0$F
  tx.t1.F <- quad$treatment$fit1$f

  # compute changes in F scores
  ct.dF <- ct.t1.F - ct.t0.F
  tx.dF <- tx.t1.F - tx.t0.F

  # construct a data frame
  data.frame(
    ct.t0.n=ct.t0$n, ct.t0.ra=ct.t0$ra, ct.t0.z=ct.t0.z, ct.t0.F=ct.t0.F,
    ct.t1.n=ct.t1$n, ct.t1.ra=ct.t1$ra, ct.t1.z=ct.t1.z, ct.t1.F=ct.t1.F,
    ct.lfc=ct.lfc, ct.dz=ct.dz, ct.dF=ct.dF,
    tx.t0.n=tx.t0$n, tx.t0.ra=tx.t0$ra, tx.t0.z=tx.t0.z, tx.t0.F=tx.t0.F,
    tx.t1.n=tx.t1$n, tx.t1.ra=tx.t1$ra, tx.t1.z=tx.t1.z, tx.t1.F=tx.t1.F,
    tx.lfc=tx.lfc, tx.dz=tx.dz, tx.dF=tx.dF
  )
}