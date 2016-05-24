statshelpers.cv <-
function(diff) {
  # Descr:    Calculates the coefficient of variance
  # Deps:     (none)
  # I/p:      inD
  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.cv", fg="red"), sep="")
  }

  # Structure of "diff" at this point:
  #        [gene1]   [gene2]  [gene3]
  # [1,]  -2.72395  -9.48870 45.79565
  # [2,]  14.97560  38.06380 88.60285
  # [3,]   8.76280  21.09895 50.51950

#  inD = as.data.frame(diff)
#  inD = ifelse(is.nan(diff), NA, inD)
  Stdv = apply(diff, MARGIN=1, sd, na.rm=T)
  Mean = rowMeans(diff)
  cv = Stdv/Mean

  # TFL is CRITICAL, because there are occasional "Inf" in the matrix for 
  # descriptive statistics "NDC"
  is.na(cv) <- do.call(cbind, lapply(cv, is.infinite))
  return(cv)
}
