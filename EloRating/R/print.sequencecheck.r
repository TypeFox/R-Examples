# print.sequencecheck 14_11_19

print.sequencecheck <- function(x, ...) {
  cat("presence data supplied, see below for details\n")
  # purely interaction sequence
  if(x$checksum["IDcheck"]==1) cat("IDs occur in the data with inconsistent capitalization (ignore if on purpose)...WARNING\n")
  if(x$checksum["selfinteractions"]==1) cat(x$selfinteractions, "\n", sep="")
  if(x$checksum["length"]==1) cat("your data vectors do not match in length\n")
  if(x$checksum["singledayobs"]==1) cat("the following individuals were observed only on one day: ", paste(x$singledaycases, collapse=", "), " ...WARNING\n", sep="" )
  
  if(x$checksum["selfinteractions"]+x$checksum["IDcheck"]+x$checksum["length"]+x$checksum["singledayobs"]==0) cat("Everything seems to be fine with the interaction sequence...OK\n")
  
  # presence related
  
  if(x$checksum["presence"]==0) cat("\n#####################################\n\n")
  
  if(x$checksum["datecol"]==1) cat("no 'Date' column found in presence...ERROR\n")
  
  if(x$checksum["startpresence1"]==1) cat("presence starts earlier than data...WARNING\n")
  if(x$checksum["startpresence2"]==1) cat("presence starts AFTER data...ERROR\n")
  if(x$checksum["endpresence1"]==1) cat("presence stops BEFORE data...ERROR\n")
  if(x$checksum["endpresence2"]==1) cat("presence continues beyond data...WARNING\n")
  
  if(x$checksum["IDmatch"]==1) { 
    cat("IDs in datasequence and presence do not match!\n")
    if(x$IDmatch1[1]!="none") { 
      cat("the following IDs occur in the presence data but NOT in the data sequence:...WARNING\n")
      cat("   ", paste(x$IDmatch1, collapse=", "), "\n")
    }
    if(x$IDmatch2[1]!="none") { 
      cat("the following IDs occur in the data sequence but NOT in the presence data:...ERROR\n")
      cat("   ", paste(x$IDmatch2, collapse=", "), "\n")
    }
  }
  
  if(x$checksum["IA_presencematch"]==1) {  
    cat("during the following interactions, IDs were absent according to the presence data:...ERROR\n")
    cat("   ", paste(x$IA_presencematchN, collapse=", "), "\n")
  }
  
  if(x$checksum["IA_presencematch"]==2) {  
    cat("your presence data does not match the interaction data and therefore no check could be performed as to whether IDs were actually present on interaction dates...WARNING\n")
  }
  
  
  if(x$checksum["presenceentries"]==1) {  cat("at least one presence entry is not 1 or 0...ERROR\n")
  }
  
  if(x$checksum["continouspres"]==1) {  cat("there appear to be gaps in your presence (days missing?)...ERROR\n")
  }
  
  
  if(sum(x$checksum[c("startpresence1", "startpresence2", "endpresence1", "endpresence2", "IDmatch", "IA_presencematch", "presenceentries", "datecol", "continouspres")])==0) cat("presence data seems to be fine and matches interaction sequence...OK\n")
  
  if(x$checksum["presence"]==0) cat("\n#####################################\n")
}

