#' Allelic frequency histograms for an ecogen genetic data frame
#' 
#' @description This program computes the frequency of each allele and plots 
#' a histogram for the number of alleles with a given frequency. 
#' The distribution is expected to be L-shaped under mutation-drift equilibrium.
#' When a factor is given, the program plots a histogram for each group.
#' @param eco Object of class "ecogen".
#' @param grp Optional factor (column of the S slot) for plots by group.
#' @examples 
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' eco.alfreq(eco)
#' eco.alfreq(eco, "pop")
#' 
#' }
#' 
#' @references 
#' 
#' Luikart G., F. Allendorf, F. Cornuet, and W. Sherwin. 1998. 
#' Distortion of allele frequency distributions provides a test for recent 
#' population bottlenecks. Journal of Heredity, 89: 238-247.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export

setGeneric("eco.alfreq", function(eco, grp = NULL) {
	
	if(!is.null(grp)) {
		cual <- which(colnames(eco@S) == grp)
		grp.num <- as.numeric(levels(eco@S[, cual]))[eco@S[, cual]] 
		nfact <- max(grp.num)
	} else {
		dummy <- rep(1, nrow(eco@G))
		eco@S <- data.frame(dummy)
		cual <- which(colnames(eco@S) =="dummy")
		grp.num <- as.numeric(levels(as.factor(eco@S[, cual])))[eco@S[, cual]] 
		nfact <- 1
	}
	
  out <- list()
	for(i in 1:nfact) {
		eco2 <- eco[which(eco@S[, cual] == i)]
		clases<- as.numeric(eco2@INT@loc.fac)
		tabla <- eco2@A
		tabla <- 2 * tabla
		frecuencia <- apply(tabla, 2, sum)
		alelos.locus <- tapply(frecuencia, clases, sum)
		for( j in 1:length(clases)) {
			temp <- clases[j]
			clases[j] <- alelos.locus[temp]
		}
		frecuencia <- frecuencia / clases
		
		#tapply(frecuencia, clases, sum) #verification
		
		frecuencia <- as.data.frame(frecuencia)
		if(nfact >1) {
			tit <-paste("Pop",levels(eco@S[, cual])[i])
		} else {
			tit <-""
		}
		
		grafico<- ggplot2::ggplot(frecuencia, ggplot2::aes(frecuencia),
															fill = "black") + 
			ggplot2::geom_histogram(ggplot2::aes(y = ..density..)) + 
			ggplot2::geom_density(alpha = 0.5, fill = "red") +
			ggplot2::labs(title = tit) +
			ggplot2::xlab("Frequency class") +
			ggplot2::ylab("Density")
		out[[i]] <- grafico
	}
	out
})

