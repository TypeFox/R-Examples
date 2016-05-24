loghelpers.prntmngr <-
function(inText, uprFlg=F, nwlFlg=T) {
  ##############################
  # Function "loghelpers.prnt" #
  ##############################
  # Descr:  print text to screen
  # Deps:   -
  # I/p:    inText

  #metadataBool = get("P2C2M_flg_metadataBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  # Logging to screen
  if (verboseBool) {
      if (uprFlg & nwlFlg) {cat("\n", toupper(inText), "\n", sep="")}
      if (!uprFlg & nwlFlg) {cat("\n", inText, "\n", sep="")}
      if (!uprFlg & !nwlFlg) {cat(inText)}
  }

  ## Logging to parameter file
  #if (metadataBool) {
  #  if (tofileFlg) {
  #    writeLines(paste("\n## ", toupper(inText), "\n", sep=""),
  #               prmFHndl)
  #  }
  #}

}
