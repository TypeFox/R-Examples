# ==>  prep_getannots&nl=xx
# key_name{|subkey_name} \n       nl such lines sent to server
# ....     \n
# <== code=0 \n
# This command must be used before using the getannots or the modifylist&type=scan commands to specify 
# what sorts of annotation records will be returned by the getannots command or will be scanned.
# nl: announces the number of key names that follow.
# key_name: an annotation key name. 
# subkey_name: optionally, an annotation sub-item name (e.g., CDS when key_name = FT)
# For the EMBL/SWISSPROT format, keys are: ALL, AC, DT, KW, OS, OC, OG, OH, 
# RN, RC, RP, RX, RA, RG, RT, RL, DR, AH, AS, CC, FH, FT, SQ, SEQ.
# For GenBank: ALL, ACCESSION, VERSION, KEYWORDS, SOURCE, ORGANISM, REFERENCE, AUTHORS, CONSRTM, 
# TITLE, JOURNAL, PUBMED, REMARK, COMMENT, FEATURES, ORIGIN, SEQUENCE. Names of annotation subitems 
# (e.g.,  JOURNAL) must be entered with their 2 leading space characters.
# For FT(embl,swissprot) and FEATURES(GenBank), one or more specific feature keys can be specified
# using lines with only uppercase and such as
# FEATURES|CDS
# FT|TRNA
# Keys ALL and SEQ/SEQUENCE stand for all annotation and sequence lines, respectively.
# For the scan operation, key ALL stand for the DE/DEFINITION lines, 
# and SEQ/SEQUENCE cannot be used (annotations but not sequence are scanned).


prepgetannots <- function(what = "all",
                       setfor = c("scan", "getannots"),
                       socket = autosocket(), verbose = FALSE){
  #
  # Default is to set for scan :
  #
  setfor <- setfor[1]
  if(!(setfor %in% c("scan", "getannots"))) stop("Wrong setfor argument")
  if(verbose) cat(paste("setfor is", setfor, "\n"))
  #
  # Get annotation lines names for current database:
  #
  annotlines <- cfl(socket = socket)$annotlines
  if(verbose) cat("annotlines:\n", annotlines, "\n")
  annotlines <- toupper(annotlines)
  #
  # Building the list of annotation lines:
  #
  if(what == "all"){
  	  if(verbose) cat("what == all, turning all to true\n")
  	  if(setfor == "scan"){
  	  	  # For scan all but the first ("ALL") and the last ("SEQ*")
  	  	  todo <- annotlines[2:(length(annotlines) - 1)]
  	  	} else {
  	  	  # For getannots all but the first ("ALL")
  	  	  todo <- annotlines[2:(length(annotlines))]
  	  	}  
  	} else { 
  		if(verbose) cat("what != all, user supplied list of annotation lines\n")
  		todo <- toupper(what)
  	}  	  	
  	#
  	# Sending request to server:
  	#  
  	if(verbose) cat("todo:\n", todo, "\n")
  	request <- paste(paste("prep_getannots&nl=", length(todo), "\n", sep = ""),
  	  	               paste(todo, collapse = "\n"), sep = "")
  	if(verbose) cat("Sending:\n", request)
  	writeLines(request, socket)
  answerFromServer <- readLines(socket, n = 1)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    warning("Empty answer from server")
    return(NA)
  }
  if(verbose) cat("... answer from server is:", answerFromServer, "\n")
  resitem <- parser.socket(answerFromServer)
  if(resitem[1] != "0"){
    stop(paste("error code returned by server :", resitem[1]))
  }
  if(verbose) cat("... everything is OK up to now\n")
  invisible(todo)
}

pga <- prepgetannots

