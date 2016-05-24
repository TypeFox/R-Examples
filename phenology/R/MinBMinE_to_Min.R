#' MinBMinE_to_Min transforms a set of parameters from MinB and MinE to Min
#' @title Transform a set of parameters from MinB and MinE to Min
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters Set of current parameters
#' @description This function is used to transform a set of parameters 
#' that uses MinB and MinE to a set of parameters 
#' that uses Min.
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
#' parg1<-MinBMinE_to_Min(parameters=parg)
#' @export


MinBMinE_to_Min <- function(parameters=stop("A set of parameters must be indicated")) {
  nm <- data.frame(nom=names(parameters[substr(names(parameters), 1, 3)=="Min"]), valeur=parameters[substr(names(parameters), 1, 3)=="Min"], stringsAsFactors=FALSE)
  nm <- cbind(nm, parametre=rep(NA, nrow(nm)), serie=rep(NA, nrow(nm)))
  l <- strsplit(nm[,1], "_")
  p <- unlist(lapply(l, function(y) y[1]))
  nm[,"parametre"] <- p
  p <- unlist(lapply(l, function(y) {
    g <- y[2]
    if (length(y)!=2) {
      for(j in 3:length(y)) g <- paste0(g, "_", y[j])
    }
    g
  }
  )
  )
  nm[,"serie"] <- p
  nm[,"serie"] <- as.factor(nm[,"serie"])
  newp <- NULL
  for(serie in levels(nm[,"serie"])) {
    v <- mean(nm[nm[,"serie"]==serie, "valeur"])
    names(v) <- paste0("Min_", serie)
    newp <- c(newp, v)
  }
  parameters <- parameters[substr(names(parameters), 1, 3)!="Min"]
  c(parameters, newp)
}
