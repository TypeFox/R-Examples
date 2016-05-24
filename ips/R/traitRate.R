# Interface to traitRate program (Mayrose & Otto, 2011)
# PACKAGE: ips
# CALLED BY: user
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-08-07

traitRate <- function(phy, seq, x, mainType = "Optimize_Model", n,
                      charModelParam1 = 0.5, charModelParam2 = 1, 
                      gammaParam = 0.5, seqModelParam1 = 2,
                      exec = "/Applications/traitRate-1.1/programs/traitRate"){
  
  ## check input data
  ## ----------------
  if ( !inherits(phy, "phylo") )
    stop("'phy' is not of class 'phylo'")
  if ( !is.ultrametric(phy) )
    stop("'phy' must be ultrametric")
  phy$node.label <- NULL # traitRate does not parse node labels
  if ( !inherits(seq, "DNAbin") )
    stop("'seq' is not of class 'DNAbin'")
  if ( !is.matrix(seq) )
    stop("'seq' must be a matrix")
  if ( ncol(seq) > 15000 )
    stop("traitRate cannot handle > 15000 bp")
  
  ## write input data
  ## ----------------
  fn <- c("model.tree", "in_msa.fasta", "in_chars.fasta")
  write.tree(phy, fn[1])
  write.fas(seq, fn[2])
  write.fas(x, fn[3])
  
  outDir <- "RESULTS"
  outFile <- "traitRate.res"
  
  ## what type of analysis?
  ## ----------------------
  mainType <- match.arg(mainType, 
                        c("Optimize_Model", 
                          "runTraitBootstrap"))
  
  ## parametric bootstrapping
  ## ------------------------
  if ( mainType == "runTraitBootstrap" ){
    res.fn <- paste(outDir, outFile, sep = "/")
    if ( !file.exists(res.fn) ) stop("cannot find rate estimates")
    ## parse estimates of rates of trait evolution
    r.est <- traitRateOutput(res.fn, "rates")
    charModelParam1 <- r.est[1]
    charModelParam2 <- r.est[2]
    ll <- traitRateOutput(res.fn, "likelihoods")
    outDir <- "RESULTS_BOOTSTRAP"
  } else {
    out <- NULL
  }
  ## write parameters file
  ## ---------------------
  traitRateParams(mainType = mainType, n, fn, 
                  outDir, outFile,
                  charModelParam1 = charModelParam1, 
                  charModelParam2 = charModelParam2,
                  gammaParam, seqModelParam1)
  
  call <- paste(exec, "traitRate.doubleRep", sep = "/")
  call <- paste(call, "params.txt")
  system(call)
  
  ## run traitRate on simulated replicates
  ## -------------------------------------
  if ( mainType == "runTraitBootstrap" ){
    for ( i in seq(from = 0, to = n - 1) ){
      od <- paste(outDir, "/sim_", i, sep = "")
      fn[3] <- paste(od, "simRandomChars.fasta", sep = "/")
      traitRateParams(mainType = "Optimize_Model", n, fn, 
                      outDir = od, outFile = outFile,
                      charModelParam1, charModelParam2,
                      gammaParam, seqModelParam1)
      
      call <- paste(exec, "traitRate.doubleRep", sep = "/")
      call <- paste(call, "params.txt")
      system(call)
      
    }
    res.fn <- paste("sim", 0:(n - 1), sep = "_")
    res.fn <- paste(outDir, res.fn, outFile, 
                    sep = "/")
    out <- lapply(res.fn, traitRateOutput)
    out <- do.call(rbind, out)
    out <- cbind(out, 
                 diff = out[, "logL"] - out[, "logL0"])
  }
  out
}