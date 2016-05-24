#' toggle_Min_PMin transforms a set of parameters from Min, MinB or MinE to PMin, PminB or PminE, or reverse
#' @title Transform a set of parameters from Min, MinB or MinE to PMin, PminB or PminE, or reverse
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters Set of current parameters
#' @description This function is used to transform a set of parameters 
#' that uses Min, MinB or MinE to a set of parameters 
#' that uses PMin, PminB or PminE, or reverse.
#' @examples 
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", reference=refdate, format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot)
#' # Change the parameters to PMinB and PMinE
#' parg1<-toggle_Min_PMin(parameters=parg)
#' # And change back to MinB and MinE
#' parg2<-toggle_Min_PMin(parameters=parg1)
#' @export


toggle_Min_PMin <-
function(parameters=stop("A set of parameters must be indicated")) {
	
# si j'ai un PMin, je le transforme en Min=Max*PMin/100
# avec les suffixes de Min pris de Max
# si j'ai un Min, je le transforme en PMin=100*Min/Max
# J'ai un seul PMin, je prends la moyenne des PMin avec les Max de toutes les séries

# je fais la même chose avec les E et B

	p <- parameters
	n <- names(parameters)
	
	pM <- (substr(n, 1, 4)=="PMin") & (substr(n, 1, 5)!="PMinE") & (substr(n, 1, 5)!="PMinB")
	pMB <- (substr(n, 1, 5)=="PMinB")
	pME <- (substr(n, 1, 5)=="PMinE")
	M <- (substr(n, 1, 3)=="Min") & (substr(n, 1, 4)!="MinE") & (substr(n, 1, 4)!="MinB")
	MB <- (substr(n, 1, 4)=="MinB")
	ME <- (substr(n, 1, 4)=="MinE")
	
	mx <- (substr(n, 1, 3)=="Max")
	
	if (any(pM)) {
		min_p <- parameters[pM]/100*parameters[mx]
		names(min_p) <- paste0("Min", substring(names(parameters[mx]), 4))
		p <- p[!pM]
		p <- c(p, min_p)
	}
	
	if (any(pMB)) {
		min_pb <- parameters[pMB]/100*parameters[mx]
		names(min_pb) <- paste0("MinB", substring(names(parameters[mx]), 4))
		n <- names(p)
		p <- p[substr(n, 1, 5)!="PMinB"]
		p <- c(p, min_pb)
	}

	if (any(pME)) {
		min_pe <- parameters[pME]/100*parameters[mx]
		names(min_pe) <- paste0("MinE", substring(names(parameters[mx]), 4))
		n <- names(p)
		p <- p[substr(n, 1, 5)!="PMinE"]
		p <- c(p, min_pe)
	}

	if (any(M)) {
		pmin_p <- mean(100*parameters[M]/parameters[mx])
		names(pmin_p) <- "PMin"
		n <- names(p)
		p <- p[!((substr(n, 1, 3)=="Min") & (substr(n, 1, 4)!="MinE") & (substr(n, 1, 4)!="MinB"))]
		p <- c(p, pmin_p)
	}

	if (any(MB)) {
		pmin_pb <- mean(100*parameters[MB]/parameters[mx])
		names(pmin_pb) <- "PMinB"
		n <- names(p)
		p <- p[substr(n, 1, 4)!="MinB"]
		p <- c(p, pmin_pb)
	}

	if (any(ME)) {
		pmin_pe <- mean(100*parameters[ME]/parameters[mx])
		names(pmin_pe) <- "PMinE"
		n <- names(p)
		p <- p[substr(n, 1, 4)!="MinE"]
		p <- c(p, pmin_pe)
	}

		
	return(p)
}
