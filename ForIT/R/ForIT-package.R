

##' Functions for Italian Tree species Volume and Biomass
##' 
##' The ForIT package is the implementation of the biomass and volume models
##' carried out by Gasparini and Tabacchi (2011) and Tabacchi et al. (2011a)
##' during the 2nd Italian National Forest Inventory. An English description of
##' the methodology is provided by Tabacchi et al. (2011b). This package is
##' intended as the close translation in R of the literature cited above.
##' 
##' \tabular{ll}{ Package: \tab ForIT\cr Type: \tab Package\cr Version: \tab
##' 1.0\cr Date: \tab 2014-07-04\cr License: \tab What license is it under?\cr
##' }
##' 
##' @name ForIT-package
##' @aliases ForIT-package ForIT
##' @docType package
##' @author Nicola Puletti, Marco Mura, Cristiano Castaldi, Maurizio Marchi, Ugo Chiavetta, Roberto Scotti
##' 
##' Maintainer: Nicola Puletti \email{nicola.puletti@@gmail.com}
##' 
##' @references Gasparini, P., Tabacchi, G.(eds), 2011. \emph{L'Inventario
##' Nazionale delle Foreste e dei serbatoi forestali di Carbonio INFC 2005.
##' Secondo inventario forestale nazionale italiano. Metodi e risultati}.
##' Edagricole. 653 pp. [ITA, ita]
##' 
##' Tabacchi G., Di Cosmo L., Gasparini P., Morelli S., 2011a. \emph{Stima del
##' volume e della fitomassa delle principali specie forestali italiane.
##' Equazioni di previsione, tavole del volume e tavole della fitomassa arborea
##' epigea. Stima del volume e della fitomassa delle principali specie
##' forestali italiane. Equazioni di previsione, tavole del volume e tavole
##' della fitomassa arborea epigea}. 412 pp. [ITA, ita]
##' 
##' Tabacchi G., Di Cosmo L., Gasparini P., 2011b. \emph{Aboveground tree
##' volume and phytomass prediction equations for forest species in Italy}.
##' European Journal of Forest Research 130: 6 911-934 [ENG, eng]
##' @keywords package
##' @examples
##' 
##' # one single tree
##' INFCvpe('Acca', 22, 14, mod='v', freq=2, aggr=FALSE)
##' 
##' # a list with ten trees of the same specie
##' INFCvpe(rep('Acca',10),d=c(10,15,20,30,32,24,36,40,8,18),
##'     h=c(7,9,12,20,21,18,21,22,8,12), mod='v', aggr=TRUE)
##' 
##' # a list of different species
##' species <- rep(c('Abal','Piab'),2)
##' dbh <- c(10,41,20,30)
##' heigths <- c(12,14,13,15)
##' frequences <- c(2,6,5,4)
##' data.frame(species, dbh, heigths, frequences)
##' 
##' # single-tree estimates
##' INFCvpe(species, dbh, heigths, mod='v', frequences, aggr=FALSE)
##' 
##' # estimates aggregated at species level
##' INFCvpe(species, dbh, heigths, mod='v', frequences, aggr=TRUE)
##' 
NULL





##' The range of applicability of \code{INFCvpe()} function
##' 
##' A \code{data.frame} containing the "range of applicability" (or "domain")
##' of \code{INFCvpe()} function
##' 
##' @name INFCdomain
##' @docType data
##' @format A data frame with 18228 observations and 2 columns
##' @references Tabacchi G., Di Cosmo L., Gasparini P., Morelli S., 2011a.
##' \emph{Stima del volume e della fitomassa delle principali specie forestali
##' italiane. Equazioni di previsione, tavole del volume e tavole della
##' fitomassa arborea epigea. Stima del volume e della fitomassa delle
##' principali specie forestali italiane. Equazioni di previsione, tavole del
##' volume e tavole della fitomassa arborea epigea}. 412 pp. [ITA, ita]
##' @keywords datasets
##' 
##' 
NULL





##' Equation statistics
##' 
##' A \code{data.frame} containing values for variance and covariance matrices
##' 
##' The list of species and of their \code{spg}-codes is:\cr \tabular{ll}{Abies alba:\tab \code{Abal}\cr
##' Acer campestre:\tab \code{Acca}\cr Acer monspessolanum:\tab \code{Acmo}\cr
##' Acer opalus:\tab \code{Acop}\cr Acer pseudoplatanus:\tab \code{Acps}\cr Alnus cordata:\tab \code{Alco}\cr Alnus
##' glutinosa:\tab \code{Algl}\cr Carpinus orientalis:\tab \code{Caor}\cr
##' Cupressus spp:\tab \code{Cusp}\cr Eucalyptus occidentalis:\tab
##' \code{Euoc}\cr Fagus sylvatica:\tab \code{Fasy}\cr Fraxinus
##' angustifolia:\tab \code{Fran}\cr Fraxinus excelsior:\tab \code{Frex}\cr
##' Fraxinus ornus:\tab \code{Fror}\cr Laburnum alpinum:\tab \code{Laal}\cr
##' Larix decidua:\tab \code{Lade}\cr Ostrya carpinifolia:\tab \code{Ossp}\cr
##' Picea abies:\tab \code{Piab}\cr Pinus cembra:\tab \code{Pice}\cr Pinus
##' halepensis:\tab \code{Piha}\cr Pinus nigra var. laricio:\tab \code{Pila}\cr
##' Pinus nigra var. nigra:\tab \code{Pini}\cr Pinus pinaster:\tab
##' \code{Pips}\cr Pinus pinea:\tab \code{Pipi}\cr Pinus radiata:\tab
##' \code{Pira}\cr Pinus strobus:\tab \code{Pist}\cr Pinus sylvestris:\tab
##' \code{Pisy}\cr Populus canescens:\tab \code{Poca}\cr Populus nigra:\tab
##' \code{Poni}\cr Populus tremula:\tab \code{Potr}\cr Prunus avium:\tab
##' \code{Prav}\cr Pseudotsuga menziesii:\tab \code{Psme}\cr Quercus
##' cerris:\tab \code{Quce}\cr Quercus ilex:\tab \code{Quil}\cr Quercus
##' pubescens:\tab \code{Qupu}\cr Robinia pseudoacacia:\tab \code{Rops}\cr
##' Salix alba:\tab \code{Saal}\cr Salix caprea:\tab \code{Saca}\cr Sorbus
##' aria:\tab \code{Soar}\cr Tilia cordata:\tab \code{Tico}\cr Tilia
##' platyphyllos:\tab \code{Tipl}\cr Ulmus minor:\tab \code{Ulmi}\cr }
##' 
##' @name INFCstats
##' @docType data
##' @format A data frame with 220 observations on the following 21 variables.
##' @references Tabacchi G., Di Cosmo L., Gasparini P., Morelli S., 2011a.
##' \emph{Stima del volume e della fitomassa delle principali specie forestali
##' italiane. Equazioni di previsione, tavole del volume e tavole della
##' fitomassa arborea epigea. Stima del volume e della fitomassa delle
##' principali specie forestali italiane. Equazioni di previsione, tavole del
##' volume e tavole della fitomassa arborea epigea}. 412 pp. [ITA, ita]
##' @keywords datasets
##' 
##' 
NULL



