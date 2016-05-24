"deletion.gvlma" <-
function(gvlmaobj)
{
  if (!inherits(gvlmaobj, "gvlma"))
    stop("First argument must be a gvlma object.")
  GTest <- gvlmaobj$GlobalTest
  n <- length(fitted(gvlmaobj))
  DeletionStatistics <- NULL
  processstat <- function(stat1, stat2)
    {
      deltastat <- (stat1$Value - stat2$Value)/stat2$Value * 100
      pvalue <- stat1$pvalue
      c(deltastat, pvalue)
    }
  whichstats <- grep("Stat", names(GTest))
  for(i in 1:n) {
    args <- list(as.name("update"))
    args$object <- gvlmaobj
    args$subset <- -i
    args$warn <- FALSE
    argscall <- as.call(args)
    GTestTemp <- eval(argscall)
    GTestTemp <- GTestTemp$GlobalTest
    ##
    stati <- mapply(processstat, GTestTemp[whichstats],
                    GTest[whichstats], SIMPLIFY = TRUE)
    DeletionStatistics <- rbind(DeletionStatistics, c(stati))
  }
  z <- as.data.frame(DeletionStatistics)
  names(z) <- c("DeltaGlobalStat", "GStatpvalue",
                "DeltaStat1", "Stat1pvalue",
                "DeltaStat2", "Stat2pvalue",
                "DeltaStat3", "Stat3pvalue",
                "DeltaStat4", "Stat4pvalue")
  class(z) <- c("gvlmaDel", "data.frame")
  attr(z, "gvlmacall") <- gvlmaobj$GlobalTest$call
  attr(z, "lmcall") <- gvlmaobj$call
  attr(z, "timeseq") <- gvlmaobj$GlobalTest$timeseq
  z
}

