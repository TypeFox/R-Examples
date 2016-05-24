#' @exportMethod rings
setGeneric("rings",
    function(object, N = seq(100, 500, 50), R = Inf, progression = TRUE){
        standardGeneric("rings")
    }
)

#' Calculation of rings of equal number of observation and/or equal radius.
#' 
#' For each cluster, this function determines a ring of equal number of observations 
#' and/or equal radius and calculates several indicators from observations located inside that ring.
#' 
#' @param object object of class \code{\link[=prevR-class]{prevR}}.
#' @param N minimum number of observations.
#' @param R maximum rings radius (in kilometers if coordinates in decimal degrees, 
#'   in the unit of the projection otherwise).
#' @param progression show a progress bar?
#' 
#' @details For each ligne of the data frame \code{clusters} of \code{object}, \code{rings} determines 
#' a ring, centred on the cluster. It could be:\itemize{
#'   \item rings of eaqul number of observations if \code{N} is finite and \code{R=Inf};
#'   \item rings of equal radius if \code{N=Inf} and \code{R} is finite;
#'   \item a combination of both (see below) if \code{N} and \code{R} are finite.
#' }
#' For \emph{rings of equal number of observations}, \code{rings} selects the smallest 
#' ring containing at least \code{N} valid observations.\cr
#' For \emph{rings of equal radius}, \code{rings} selects all clusters located at a lower
#' distance than \code{R} from the central cluster.\cr For \emph{combination of both}, \code{rings} 
#' calculates firts the ring with the minimum number of observations and test if its radius is lower 
#' than \code{R} or not. If so, the ring is kept, otherwise the ring of maximum radius is calculated.
#' 
#' Different series of rings could be simultaneoulsy calculated by providing different values for \code{N} 
#' and \code{R}. \code{rings} will calculate rings corresponding to each couple (N,R).
#' 
#' @return Return \code{object} with the slot \code{rings} completed for each couple (N,R).
#' 
#' Each entry is composed of 3 elements: \code{N}, minimum number of observations per ring; \code{R}, 
#' maximum radius of rings and \code{estimates}, a data frame with the following variables:\itemize{
#'   \item "id" cluster ID.
#'   \item "r.pos" number of positive cases inside the ring.
#'   \item "r.n" number of valid observations inside the ring.
#'   \item "r.prev" observed prevalence (in \%) inside the ring (r.pos/r.n).
#'   \item "r.radius" ring radius (in kilometers if coordinates in decimal degrees, 
#'       in the unit of the projection otherwise).
#'   \item "r.clusters" number of clusters located inside the ring.
#'   \item "r.wpos" (optional) sum of weights of positive cases inside the ring.
#'   \item "r.wn" (optional) sum of weights of valid observations inside the ring.
#'   \item "r.wprev" (optional) weighted observed prevalence (in \%) inside the ring (r.wpos/r.wn).
#' }
#' Note: the list \code{rings} is named, the name of each element is N\emph{N_value}.R\emph{R_value}, 
#' for example \emph{N300.RInf}. 
#' 
#' Note 2: \emph{r.wpos}, \emph{r.wn} and \emph{r.wprev} are calculated only if the slot \code{clusters} 
#' of \code{object} contains weighted data.
#' 
#' @seealso \code{\link{prevR-class}}.
#' 
#' @references 
#' Larmarange Joseph, Vallo Roselyne, Yaro Seydou, Msellati Philippe and Meda Nicolas (2011) 
#' "Methods for mapping regional trends of HIV prevalence from Demographic and Health Surveys (DHS)", 
#' \emph{Cybergeo : European Journal of Geography}, no 558, \url{http://cybergeo.revues.org/24606}, 
#' DOI: 10.4000/cybergeo.24606.
#' 
#' Larmarange Joseph (2007) \emph{Prévalences du VIH en Afrique : validité d'une mesure}, 
#' PhD thesis in demography, directed by Benoît Ferry, université Paris Descartes, 
#' \url{http://tel.archives-ouvertes.fr/tel-00320283}.
#' 
#' @examples 
#' \dontrun{
#' print(fdhs)
#' dhs <- rings(fdhs,N=c(100,200,300,400,500))
#' print(dhs)
#' }
#' @aliases rings rings-methods rings,prevR-method
#' @keywords math spatial

setMethod("rings","prevR",
  function (object, N = seq(100, 500, 50), R = Inf, progression = TRUE)
  {
  ###############################################################################################
  # Cette fonction calcule la prevalence (et la prevalence ponderee) par la methode des cercles.
  # Cette fonction ajoute chaque ring defini par un couple N R a la list du slot rings de object
  # C'est a dire que si le slot rings de object n'est pas NULL alors les nouveaux resulats
  #    de calul viennent s'ajouter a la liste contenue dans l'element rings de object
  # Chaque ring est une liste contenant 3 elements
  #    N : contient la valeur de N pour laquelle ce ring a ete calcule
  #    R : contient la valeur de R pour laquelle ce ring a ete calcule
  #    estimates : contient un dataframe dont les colonnes sont
  #       id : l'identificateur du cercle (qui est egal a l'identificater du cluster centre du cercle)
  #       r.pos : le nombre de cas positifs dans le cercle
  #       r.n : l'effectif du cercle
  #       r.prev : la prevalence du cercle (100*r.pos/r.n)
  #       r.radius :  le rayon du cercle
  #       r.clusters : le nombre de clusters dans le cercle
  #       Si des donnees ponderes sont presentes 
  #       r.wpos : le nombre pondere de cas positifs dans le cercle
  #       r.wn : l'effectif pondere du cercle
  #       r.wprev : la prevalence ponderee du cercle (100*r.wpos/r.wn)
  # Remarque la list rings est nommee. Le nom de chacun de ses elements est de la forme  
  #  N(valeur de n ).R(valeur de R) : exemple N500.RInf 
  # Besoin du package fields pour rdist et rdist.earth   
  #
  ###############################################################################################

# On teste la coherence des donnees d'entree N et R Pour plus de precision regardez la fonction .isInputOk.prevR
    .isInputOk.prevR(N =N, R = R)
    clusters = slot(object,"clusters")
    boundary = slot(object,"boundary")
    coord.clust = clusters[,c("x","y")]

    projCRS = slot(object,"proj")
    proj    = slot(projCRS,"projargs")
# Si les donnees georeferencees sont exprimees en longitude latitude les distances entre clusters sont calculees en Km
# par la fonction rdist.earth du package fields
# Autrement les distances sont calculees par la fonction rdist (fields) en unite des donnees de depart
    if(regexpr("longlat",proj)==-1 && regexpr("latlong",proj)==-1){
      distances = fields::rdist(coord.clust)
    } else {
      distances = fields::rdist.earth(coord.clust,miles=F)
    }
    rings   = slot(object,"rings")
# On cree les couple N-R a partir desquels seront calclules les rayons des cercles
# et par consequent l'effectif et la prevalence dans chaque cercle 
    couples = cbind(N=N,R=R)
# Barre de progression
    if (progression) {
      message("Progress of calculations:", domain="R-prevR")
      barre = txtProgressBar(min=0, max=nrow(couples)*nrow(clusters), initial=0, style=3)
    }
# boucle sur les couples  N R 
    for (ic in 1:nrow(couples)) {
      one.N = couples[ic,"N"]
      one.R = couples[ic,"R"]

      result  = data.frame()
# boucle sur les clusters
      for (i in 1:nrow(clusters)) {
        one.clust = clusters[i, ]
        temp = clusters
        temp$dist = distances[i, ]
        temp = temp[order(temp$dist), ]

        temp$cum.n = cumsum(temp[["n"]])
        if (length(temp[temp$cum.n >= one.N, ]$cum.n) == 0){
          maxi = Inf
        } else {
          maxi = min(temp[temp$cum.n >= one.N, ]$cum.n)
        }
        temp = temp[temp$cum.n <= maxi, ]
#
# A ce niveau temp contient un sous dataframe de clusters qui est tel que
# Les clusters (c'est a dire les lignes de temp) sont classes en fonction de la distance au cluster d'origine (one.clust)
#    du plus pres au plus loin
# On ne garde que les clusters qui sont tels que l'effectif des clusters dans le cercle est de l'ordre de N sous la contrainte
#    que le rayon doit etre inferieur a R
# 

        temp2      = temp[temp$dist <= one.R, ]
        one.result = data.frame(id = one.clust[,"id"])
# r.pos : le nombre de cas positifs dans le cercle
# r.n : l'effectif du cercle
# r.prev : la prevalence du cercle (100*r.pos/r.n)
# r.radius :  le rayon du cercle
# r.clusters : le nombre de clusters dans le cercle 
        one.result$r.pos             = sum(temp2[["pos"]],na.rm=T)
        one.result$r.n               = sum(temp2[["n"]],na.rm=T)
        one.result$r.prev            = 100*one.result$r.pos/one.result$r.n
        if ( one.R != Inf && one.N == Inf) {
          one.result$r.radius          = R
        } else {
          one.result$r.radius          = max(temp2$dist)
        }
        one.result$r.clusters        = length(temp2$dist)

# Si des donnees ponderes sont presentes on calule
# r.wpos : le nombre pondere de cas positifs dans le cercle
# r.wn : l'effectif pondere du cercle
# r.wprev : la prevalence ponderee du cercle (100*r.wpos/r.wn)
        wn     = temp2[["wn"]]
        if(!is.null(wn)) {
          one.result$r.wpos  = sum(temp2[["wpos"]])
          one.result$r.wn    = sum(temp2[["wn"]])
          one.result$r.wprev = 100*one.result$r.wpos/one.result$r.wn
        }
        result = rbind(result, one.result)
        
        if (progression) {
          barre.cur = (ic-1)*nrow(clusters)+i
          setTxtProgressBar(barre,value=barre.cur)
        }
        
      }
      rings[[paste("N",one.N,".R",one.R,sep="")]] = list(N=one.N,R=one.R,estimates=result)
    }

    if (progression) {
      close(barre)
    }
    slot(object,"rings") = rings
    return(object)
  }
)