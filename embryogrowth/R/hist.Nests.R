#' hist.Nests shows the histogram of temperatures with set of nests
#' @title Show the histogram of temperatures with set of nests
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return A list with an histogram object with information on histogram or 
#' NULL if no series was selected and the complete set of temperatures used.
#' @param x Data formated using formatdata.
#' @param series Series to be used, logical (TRUE ou FALSE), numbers or names. If "all", all series are used.
#' @param ... Parameters used by hist function
#' @description Show the histogram of temperatures with set of nests
#' hist(data)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' h <- hist(formated, series="all")
#' }
#' @method hist Nests
#' @export

hist.Nests <- function(x, series="all", ...) {

nids <- x

# Plus possible car fonction speciale
 if (inherits(nids, "NestsResult")) {
# j'ai un objet de resultat
# je prends les donnees
	nids <- nids$data
 }


if (!inherits(nids, "Nests")) {
#Je n'ai ni un objet result ni un objet formate. Je quitte
	print("'Nests' object obtained after FormatNests() or 'NestsResult' obtained after searchR() must be provided.")
	return(invisible())
}


NbTS <- nids[["IndiceT"]][3]
if (series[[1]]=="all") {
	series<-rep(TRUE, NbTS)
} else {
	if (any(!is.logical(series))) {
		if (is.numeric(series)) {
			seriesx <- rep(FALSE, NbTS)
			seriesx[series] <- TRUE
		} else {
			seriesx <- (names(nids[1:NbTS])==series)
		}
		series <- seriesx
	} else {
# c'est des valeurs logiques, je verifie si le bon nombre
		if (length(series)!=NbTS) {
			series <- rep(series, NbTS)
			series <- series[1:NbTS]
		}
	}
}

temptotal=NULL

for (j in 1:NbTS) {
	if (series[j]) {

		nidsx <- nids[[j]][, 1:2]
		colnames(nidsx) <- c("Time", "Temperature")

# Je dois recalculer toutes les informations heure par heure pour avoir la TSP en heure
# la duree d'incubation est= (0:(nids[,1][length(nids[,1])]%/%60)-1)*60
# et ensuite j'intercale les temperatures

# Je cree un tableau avec les donnees heures par heure et je rajoute la derniere donnee
		tl1 <- c((0:(nidsx[,1][length(nidsx[,1])]%/%60-1))*60, nidsx[,1][length(nidsx[,1])])
# je prends les vraies donnees
		tl2 <- nidsx[,1]

		it <- findInterval(tl1, tl2)

		temptotal<-append(temptotal, nidsx[it,"Temperature"])

	}
}


if (is.null(temptotal)) {
	print("No nest is selected !")
	a <- NULL
} else {
	x <- temptotal

	L <- modifyList(list(ylab="Frequencies among the nests", xlab=expression("Temperatures in " * degree * "C"), main="", freq=FALSE, las=1), modifyList(list(x=x), list(...))) 
  
	if (any(names(L)=="plot")) {
    if (!L$plot)
      L <- modifyList(list(x=x), list(...)) 
	}
  
	a <- do.call(hist, L) 

}

return(invisible(list(histogram=a, temperatures=x)))

}
