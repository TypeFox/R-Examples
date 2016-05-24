save_log <-
function(phaseL, checkL, typeL, SNPL = allSNPs, allSNPs = 1L, actionL, noteL = "", fileL) {
  write.table(t(c(phaseL, checkL, typeL, SNPL, round(100*SNPL/allSNPs, digits = 2), actionL, noteL)), paste0(fileL, "_log.txt"), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  return(invisible())
}
