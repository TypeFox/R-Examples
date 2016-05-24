

################################################################################
# summary method for objects of class "din"                                    #
################################################################################
summary.din <-
function(object, top.n.skill.classes = 6, overwrite = FALSE, 
    ...){

# Call: generic
# Input: 
#   object of class din
#   top.n.skill.classes: a numeric, specifying the number of skill classes, 
#  starting with the most frequent, to be returned. Default value is 6.
#   log.file: an optional character vector, specifying the directory and/or filename for an extensive log file.
#   overwrite: an optional boolean, specifying wether or not the method is supposed to 
#   overwrite an existing log.file. If the log.file exists and overwrite is FALSE, 
#    the user is asked to confirm the overwriting.
# Output: a named list, of the class summary.din (to be passed to print.summary.din), 
#    consisting of the following five components
# 	CALL: a character specifying the model rule, the number of items and the 
#      number of attributes underlying the items.
#	  IDI: a vector giving the item discrimination index. (see help file)
#	  SKILL.CLASSES: a vector giving the top.n.skill.classes most frequent skill 
#             classes and the corresponding class probability.
#	  AIC: a numeric giving the AIC of the specified model object.
#	  BIC: a numeric giving the BIC of the specified model object.

################################################################################
# extract output from din object                                               #
################################################################################

	log.file <- NULL
#	osink( file = file , suffix = paste0( "__SUMMARY.Rout") )

	CALL <- paste(object$display,"on", ncol(object$data), "items for", nrow(object$skill.patt),"attributes")
	
	AIC <- round(object$AIC, 3)
	BIC <- round(object$BIC, 3)

	IDI <- t(matrix(round(object$IDI, 4)))
	rownames(IDI) <- ""
	colnames(IDI) <- rownames(object$item)
	# item parameters
	item <- data.frame( "item" = colnames(object$data) , "guess" = object$guess[,1] , 
				"slip" = object$slip[,1] , "IDI" = object$IDI , "rmsea" = object$itemfit.rmsea )
	for (vv in 2:5){ item[,vv] <- round( item[,vv] , 3 ) }
	
#	SKILL.CLASSES <- object$attribute.patt[order(object$attribute.patt[,1], decreasing = TRUE),][
#	1:min(top.n.skill.classes, 2^length(object$skill.patt)),]
#	SKILL.CLASSES <- round(t(SKILL.CLASSES)[1, ], 4)
	
#  	if(top.n.skill.classes > nrow(object$attribute.patt))
#  		warning("There are at most ", 2^length(object$skill.patt), 
#  		" different skill classes. Returning all skill classes.\n")
  	
################################################################################
# catch log file writing errors                                                #
################################################################################

  if(!is.null(log.file)){
    if(file.exists(log.file)){
      if(!overwrite){
        cat("Press 'y' to overwrite existing file: ")
        conf <- readLines(con = stdin(), n = 1)
        if(conf %in% c("y", "Y")){
          wrn <- getOption("warn"); options(warn = -1)
          err <- try({ff <- file(log.file); open(ff, "w"); close(ff)}, silent = TRUE)
          options("warn" = wrn)
          if(!is.null(err)){
            warning("'log.file' argument ignored due to ", err)
            log.file <- NULL
          }
        }
      }      
    }else{
      wrn <- getOption("warn"); options(warn = -1)
      err <- try({ff <- file(log.file); open(ff, "w"); close(ff)}, silent = TRUE)
      options("warn" = wrn)
      if(!is.null(err)){
        warning("'log.file' argument ignored due to ", err)
        log.file <- NULL
      }
    }
  }
 
################################################################################
# return list                                                                  #
################################################################################

# print(object)
	
  out <- list("CALL"=CALL,"IDI"=IDI, 
			"call" = deparse(object$call) , 
			"deviance" = -2*object$loglike,
#			"SKILL.CLASSES"=SKILL.CLASSES, 
			"AIC"=AIC, "BIC"=BIC, "item" = item , 
			"Npars" = object$Npars , 
			"log.file" = file, "din.object" = object,
			"start.analysis"=object$start.analysis ,
			"end.analysis"=object$end.analysis
			)			
	class(out) <- "summary.din"
	return(out)
}