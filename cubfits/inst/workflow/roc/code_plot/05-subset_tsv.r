### This script collect the posterior means of MCMC runs.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))
source("00-set_env.r")

### Get all cases.
for(i.case in case.names){
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    next
  }
  load(fn.in)
  # fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  # if(!file.exists(fn.in)){
  #   next
  # }
  # load(fn.in)

  ### For log.mu.
  all.names <- names(b.logmu.PM)
  id.var <- grep("log.mu", all.names)

  AA <- gsub("(.)\\.(.*)", "\\1", b.logmu.label)
  CODON <- gsub("(.)\\.(.*)", "\\2", b.logmu.label)

  tmp.PM <- b.logmu.PM[id.var]
  tmp.MED <- b.logmu.MED[id.var]
  tmp.CI <- b.logmu.ci.PM[id.var,]
  tmp.STD <- b.STD[id.var]
  ret <- data.frame(AA = AA, CODON = CODON, Mean = tmp.PM, Median = tmp.MED,
                    CI.025 = tmp.CI[, 1], CI.975 = tmp.CI[, 2], Std = tmp.STD)
  for(i.aa in unique(AA)){
    tmp <- .CF.GV$synonymous.codon.split[[i.aa]]
    tmp <- data.frame(AA = i.aa, CODON = tmp[which(!(tmp %in% CODON))],
                      Mean = 0, Median = 0, CI.025 = 0, CI.975 = 0, Std = 0)
    ret <- rbind(ret, tmp)
  }
  order.id <- order(as.character(ret$AA), as.character(ret$CODON))
  ret <- ret[order.id,] 

  fn.out <- paste(prefix$table, "logmu_", i.case, "_PM.tsv", sep = "")
  write.table(ret, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)

  ### For Delta.t.
  all.names <- names(b.negsel.PM)
  id.var <- grep("Delta.t", all.names)

  AA <- gsub("(.)\\.(.*)", "\\1", b.negsel.label)
  CODON <- gsub("(.)\\.(.*)", "\\2", b.negsel.label)

  tmp.PM <- b.negsel.PM[id.var]
  tmp.MED <- b.negsel.MED[id.var]
  tmp.CI <- b.negsel.ci.PM[id.var,]
  tmp.STD <- b.STD[id.var]
  ret <- data.frame(AA = AA, CODON = CODON, Mean = tmp.PM, Median = tmp.MED,
                    CI.025 = tmp.CI[, 1], CI.975 = tmp.CI[, 2], Std = tmp.STD)
  for(i.aa in unique(AA)){
    tmp <- .CF.GV$synonymous.codon.split[[i.aa]]
    tmp <- data.frame(AA = i.aa, CODON = tmp[which(!(tmp %in% CODON))],
                      Mean = 0, Median = 0, CI.025 = 0, CI.975 = 0, Std = 0)
    ret <- rbind(ret, tmp)
  }
  order.id <- order(as.character(ret$AA), as.character(ret$CODON))
  ret <- ret[order.id,] 

  fn.out <- paste(prefix$table, "deltat_", i.case, "_PM.tsv", sep = "")
  write.table(ret, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)

  ### For prior.
  param.name <- c("sigmaW", "m.Phi", "s.Phi", "bias.Phi")
  if(.CF.CT$type.p != "lognormal_bias"){
    param.name <- param.name[-which(param.name == "bias.Phi")]
  }
  if(length(grep("wophi", i.case)) > 0){
    param.name <- param.name[-which(param.name == "sigmaW")]
    param.name <- param.name[-which(param.name == "bias.Phi")]
  }
  ret <- data.frame(param = param.name, Mean = p.PM, Median = p.MED,
                    CI.025 = p.CI[, 1], CI.975 = p.CI[, 2], Std = p.STD)
  fn.out <- paste(prefix$table, "prior_", i.case, "_PM.tsv", sep = "")
  write.table(ret, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)

  ### For E[Phi].
  ret <- data.frame(ORF = names(phi.PM), Mean = phi.PM, Median = phi.MED,
                    CI.025 = phi.CI[, 1], CI.975 = phi.CI[, 2], Std = phi.STD)

  fn.out <- paste(prefix$table, "phi_", i.case, "_PM.tsv", sep = "")
  write.table(ret, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)
}
