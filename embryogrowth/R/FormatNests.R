#' FormatNests creates a dataset of class "Nests" to be used with searchR
#' @title Create a dataset of class Nests to be used with searchR
#' @author Marc Girondot
#' @return A list with all the nests formated to be used with searchR.
#' @param data Data to be newly formated
#' @param previous Data already formated
#' @param simplify If TRUE, simply the time series by removing identical time series of temperatures
#' @param weight Named vector with weight for likelihhod
#' @description Will create a dataset of class Nests to be used with searchR\cr
#' FormatNests(nest, previous=x) with x being a previously formated data.\cr
#' The raw data must be organized being:\cr
#' First column is the time in minutes since the beginning of incubation\cr
#' Each column next is the trace of temperatures, one column for each nest.\cr
#' For example, for two nests:\cr
#' Time   Nest1    Nest2\cr
#' 0       29.8     27.6\cr
#' 90      30.2     28.8\cr
#' 120     30.4     30.7\cr
#' 180     31.2     32.6\cr
#' ...\cr
#' 65800   30.8     32.6\cr
#' 65890            30.2\cr
#' 65950            30.4\cr
#' \cr
#' The Nest1 ends incubation at 65800 minutes whereas Nest2 ends incubation at 65950 (last row\cr
#' with temperature for each).\cr
#' The parameter Weight is a vector: weight=c(Nest1=1, Nest2=1.2)\cr
#' It can be used to format database already formated with old format; in this case, just use data=xxx with xxx being the old format database.\cr
#' @examples 
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest, previous=NULL)
#' formated <- FormatNests(nest)
#' }
#' @export


FormatNests <-
function(data=stop("A dataset must be provided !"), previous=NULL, 
         simplify=TRUE, weight=NULL) {

# Je cree une une fonction qui prepare un fichier pour etre utilise
# Les differents nids seront des matrices dans une liste
# Dans chaque matrice on a le temps depuis le debut de l'incubation
# la temperature en C, la temperature en K, la valeur de r et la masse

# previous=NULL; simplify=TRUE; weight=NULL
  
if (class(data)=="Nests") {

  nidsEC <- data[1:data$IndiceT["NbTS"]]
  for(i in 1:length(nidsEC)) {
    ess <- nidsEC[[i]]
    if (dim(ess)[2]==5) {
      ess <- cbind(ess, IndiceK=NA)
      nidsEC[[i]] <- ess
    }
  }
  
} else {

nidsEC <- list()
# meme chose que nids=as.list(NULL)

for (j in 2:dim(data)[2]) {

# je veux le transformer en une liste avec a l'interieur un dataframe par enregistreur
# je commence par sortir la liste un par un qui n'ont pas de NA

# je stocke les donnees de l'enregistreur dans une matrice

newess2 <- as.numeric(subset(data[,j], !is.na(data[,j])))
newessT <- as.numeric(subset(data[,1], !is.na(data[,j])))

# 19/10/2012 je calcule les etats intermediaires en terme de temps

ess2 <- newess2[1]
essT <- newessT[1]

for(i in 1:(length(newess2)-1)) {
	ess2 <- c(ess2, newess2[i+1])
	minutes <- newessT[i]+(newessT[i+1]-newessT[i])/2
	essT <- c(essT, minutes)
	
}

ess2 <- c(ess2, newess2[length(newess2)])
essT <- c(essT, newessT[length(newess2)])

ess<-matrix(c(essT, ess2, ess2+273.15, rep(NA, 3*length(ess2))), ncol=6)

colnames(ess)<-c("Time", "Temperatures C", "Temperatures K", "r", "Mass", "IndiceK")


# ensuite je supprime les temps avec des valeurs de temperatures identiques
# sauf la derniere qui doit rester - 20/7/2012
# si que deux lignes, je ne fais rien - 27/7/2012
if ((dim(ess)[1]>2) & simplify) {
	for (i in 2:(dim(ess)[1]-1)) {
    	if (ess[i,2]==ess[i-1,2]) ess[i,1]=NA
	}
}

# nidsEC[[names(data[j])]]<-subset(ess, !is.na(ess[,1]))
# 23/4/2015
nidsEC[[colnames(data)[j]]]<-subset(ess, !is.na(ess[,1]))
}

}


if (!is.null(previous)) {
	nidsprevious<-previous[1:previous$IndiceT["NbTS"]]
	nidsEC<-c(nidsEC, nidsprevious)
}

# J'enregistre les tempmin et max car toujours pareil - Je n'en ai plus besoin; je les garde cependant

# j'enregistre les temp comme des facteurs
nbts <- length(nidsEC)
temp<-NULL
for (j in 1:nbts) {
	temp<- c(temp, nidsEC[[j]][, "Temperatures K"])
}
tempmin <- min(temp)
tempmax <- max(temp)

tempaf<- as.factor(temp)
templevels <- levels(tempaf)

nidsEC[["IndiceT"]] <- c(Tmin=tempmin,Tmax=tempmax, NbTS=nbts)
nidsEC[["Temperatures"]] <- templevels

# Je remplis la colonne IndiceK

for (j in 1:nbts) {
  temp<- nidsEC[[j]][, "Temperatures K"]
  nidsEC[[j]][, "IndiceK"] <- match(as.character(temp), nidsEC[["Temperatures"]])
}

nidsEC[["weight"]] <- weight

class(nidsEC) <- "Nests"

# 27/4/2015

if (any(duplicated(names(nidsEC)[1:(nidsEC[["IndiceT"]]["NbTS"])]))) {
  stop("Nests must have unique names")
} else {
  return(nidsEC)
}

}
