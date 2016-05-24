##' Volume and Phytomass Estimates
##' 
##' Estimate tree volume and phytomass per species from stem diameter at 1.3 m
##' height (DBH) and total height (HT)
##' 
##' 
##' @param spg a string or factor indicating the species. \cr Possible codes
##' are:\cr \tabular{ll}{Abies alba:\tab \code{Abal}\cr Acer campestre:\tab
##' \code{Acca}\cr Acer monspessolanum:\tab \code{Acmo}\cr Acer opalus:\tab
##' \code{Acop}\cr Acer pseudoplatanus:\tab \code{Acps}\cr Alnus cordata:\tab
##' \code{Alco}\cr Alnus glutinosa:\tab \code{Algl}\cr Carpinus orientalis:\tab
##' \code{Caor}\cr Cupressus spp:\tab \code{Cusp}\cr Eucalyptus
##' occidentalis:\tab \code{Euoc}\cr Fagus sylvatica:\tab \code{Fasy}\cr
##' Fraxinus angustifolia:\tab \code{Fran}\cr Fraxinus excelsior:\tab
##' \code{Frex}\cr Fraxinus ornus:\tab \code{Fror}\cr Laburnum alpinum:\tab
##' \code{Laal}\cr Larix decidua:\tab \code{Lade}\cr Ostrya carpinifolia:\tab
##' \code{Ossp}\cr Picea abies:\tab \code{Piab}\cr Pinus cembra:\tab
##' \code{Pice}\cr Pinus halepensis:\tab \code{Piha}\cr Pinus nigra var.
##' laricio:\tab \code{Pila}\cr Pinus nigra var. nigra:\tab \code{Pini}\cr
##' Pinus pinaster:\tab \code{Pips}\cr Pinus pinea:\tab \code{Pipi}\cr Pinus
##' radiata:\tab \code{Pira}\cr Pinus strobus:\tab \code{Pist}\cr Pinus
##' sylvestris:\tab \code{Pisy}\cr Populus canescens:\tab \code{Poca}\cr
##' Populus nigra:\tab \code{Poni}\cr Populus tremula:\tab \code{Potr}\cr
##' Prunus avium:\tab \code{Prav}\cr Pseudotsuga menziesii:\tab \code{Psme}\cr
##' Quercus cerris:\tab \code{Quce}\cr Quercus ilex:\tab \code{Quil}\cr Quercus
##' pubescens:\tab \code{Qupu}\cr Robinia pseudoacacia:\tab \code{Rops}\cr
##' Salix alba:\tab \code{Saal}\cr Salix caprea:\tab \code{Saca}\cr Sorbus
##' aria:\tab \code{Soar}\cr Tilia cordata:\tab \code{Tico}\cr Tilia
##' platyphyllos:\tab \code{Tipl}\cr Ulmus minor:\tab \code{Ulmi}\cr }
##' @param d a value or vector indicating the stem diameter at 1.3 m height
##' (DBH) [cm]
##' @param h a value or vector indicating the total stem height (HT) [m]
##' @param mod a character: \emph{v} for volume of the stem and large branches,
##' \emph{dw1} for phytomass of the stem and large braches, \emph{dw2} for
##' phytomass of the small branches, \emph{dw3} for phytomass of the stump,
##' \emph{dw4} for phytomass of the whole tree
##' @param freq the number of trees of the same spg with equal DBH and HT
##' @param aggr a flag allowing estimates aggregated at \code{spg} level (when
##' aggr is \code{TRUE}) or at tree level (when aggr is \code{FALSE})
##' @return Returns a \code{list} with the following objects: \cr
##' \item{$mainData }{a \code{data.frame} with the following columns
##' \tabular{ll}{ \code{spg} \tab a string with the species group code \cr \code{d130}\tab
##' a value indicating the stem diameter at 1.3 m height (DBH) [cm] \cr
##' \code{h_tot}\tab a value indicating the total stem height (HT) [m]\cr \code{freq}\tab
##' the number of trees of the same \code{spg} with equal DBH and HT"\cr \code{mod}\tab a
##' character, the same as \code{mod} in arguments. \cr \code{T_0}\tab a value of
##' the estimates for \code{mod} value.  \code{v} is expressed in \emph{dm^3}
##' while \code{dw1}, \code{dw2}, \code{dw3}, \code{dw4} are expressed in
##' \emph{kg}; \cr \code{SEE}\tab a value of Standard Error of the Estimates \cr
##' \code{dof}\tab the degree of freedom \cr \code{in.range}\tab tree inside (\code{y})
##' or out of the range (\code{n}) of the sampled trees in Gasparini and
##' Tabacchi (2011)\cr } } \item{out.of.range }{a \code{data.frame} listing the
##' trees out of the range of application (domain) }
##' @author Nicola Puletti \email{nicola.puletti@@gmail.com}, Marco Mura,
##' Cristiano Castaldi, Roberto Scotti
##' @references Gasparini, P., Tabacchi, G.(eds), 2011. \emph{L'Inventario
##' Nazionale delle Foreste e dei serbatoi forestali di Carbonio INFC 2005.
##' Secondo inventario forestale nazionale italiano. Metodi e risultati}.
##' Edagricole. 653 pp. [ITA, ita]
##' 
##' Tabacchi G., Di Cosmo L., Gasparini P., Morelli S., 2011a. \emph{Stima
##' del volume e della fitomassa delle principali specie forestali italiane.
##' Equazioni di previsione, tavole del volume e tavole della fitomassa arborea
##' epigea. Stima del volume e della fitomassa delle principali specie
##' forestali italiane. Equazioni di previsione, tavole del volume e tavole
##' della fitomassa arborea epigea}. 412 pp. [ITA, ita]
##' 
##' Tabacchi G., Di Cosmo L., Gasparini P., 2011b. \emph{Aboveground tree
##' volume and phytomass prediction equations for forest species in Italy}.
##' European Journal of Forest Research 130: 6 911-934 [ENG, eng]
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
##' 
##' @export INFCvpe
INFCvpe <-
function(spg, d, h, mod, freq, aggr=F) {
  
  aggr <- ifelse(sum(table(spg)>1)==0, F, aggr)
  
  if(missing(freq)) {
    cav <- data.frame(spg, d130=d, h_tot=h, frequency=rep(1,length(d)) )
  } else {
    cav <- data.frame(spg, d130=d, h_tot=h, frequency=freq)
  }
  
  grandezze <- data.frame(grandezza= c(
    "v: volume of the stem and large branches [dm^3]",
    "dw1: phytomass of the stem and large braches [kg]",
    "dw2: phytomass of the small branches [kg]",
    "dw3: phytomass of the stump [kg]",
    "dw4: phytomass of the whole tree [kg]"
  ))
  
  grandezze$cod <- substr(grandezze$grandezza,1,regexpr(":", grandezze$grandezza)-1)
  if (!(mod %in% grandezze$cod)) {cat("ERROR - Unrecognized esitimation quantity")}
  
  spl <- as.character(unique(cav$spg))
  
  # s <- spl[1]
  
  result <- NULL
  for (s in spl) {
    stats_s <- INFCstats[INFCstats$spg==s & INFCstats$mod==mod,]
    
    n_par <- stats_s$n_par
    b <- c(stats_s$b.1., stats_s$b.2.)
    b <- t(t(b))
    
    cav2 <- cav[cav$spg==s,]
    d2h <- with(cav2, d130^2 * h_tot * frequency)
    one <- rep(1, length(d2h))
    D_0 <- cbind(one, d2h)
    D_0_u <- one %*% D_0
    T_0 <- D_0_u %*% b
    if(aggr==F) T_0 <- D_0 %*% b
    
    if (n_par == 3) b[3] <- stats_s$b.3.
    b <- t(t(b))
    d130 <- cav2$d130
    d2h <- with(cav2, d130^2 * h_tot * frequency)
    one <- rep(1, length(d2h))
    if (n_par == 3) D_0 <- cbind(one, d2h, d130) # ho aggiunto un "if (n_par == 3)" per evitare l'incompatibilità
    D_0_u <- one %*% D_0
    T_0 <- D_0_u %*% b
    if(aggr==F) T_0 <- D_0 %*% b
    
    ##
    mvc <- NA
    if (n_par == 2) mvc <- matrix(c(stats_s$mvc.1., stats_s$mvc.3.,
                                    stats_s$mvc.3., stats_s$mvc.4.), 2, 2) # [i] row = [i] col
    
    if (n_par == 3) mvc <- matrix(c(stats_s$mvc.1., stats_s$mvc.4., stats_s$mvc.7.,
                                    stats_s$mvc.4., stats_s$mvc.5., stats_s$mvc.8.,
                                    stats_s$mvc.7., stats_s$mvc.8., stats_s$mvc.9.), 3, 3)
    
    sa2 <- stats_s$saq
    
    if (aggr==F) {Var_tot <- (D_0[1,] %*% mvc %*% t(t(D_0[1,]))) + (sa2*d2h^2) #rel. [22]
    } else {
      W_0 <- diag(x=d2h^2)
      Var_stima <- D_0_u %*% mvc %*% t(D_0_u)
      Var_resid <- one %*% W_0 %*% one * sa2
      Var_tot <- Var_stima + Var_resid
    }
    
    SEE <- sqrt(Var_tot)
    dof <- with(stats_s, n_oss-n_par)
    
    ##
    if(aggr==T){
      out <- data.frame(spg=s, mod=mod, T_0=round(T_0,1), N=sum(cav2$frequency), SEE, dof)
    } else {
      out <- data.frame(spg=cav2$spg, d130=cav2$d130, h_tot=cav2$h_tot, freq=cav2$frequency, mod=mod, T_0=round(T_0,1), SEE, dof)
    }
    
    result <- rbind(result, out)
    
  }
  result <- result[!duplicated(result), ]
  
  if(aggr==F) {
    result$key <- paste(result$spg, round(result$d130,0), round(result$h_tot,0),sep='-')
    
    result <- merge(result, INFCdomain, by = c("key", "key"), all.x=T)
    result$in.range <- as.character(result$in.range)
    result$in.range[is.na(result$in.range)] <- 'n'
    
    out.of.range <- result[result$in.range=='n',]
  }
  
  if(aggr==T) {
    cav$key <- paste(cav$spg, round(cav$d130,0), round(cav$h_tot,0),sep='-')
    
    cav <- merge(cav, INFCdomain, by = c("key", "key"), all.x=T)
    cav$in.range <- as.character(cav$in.range)
    cav$in.range[is.na(cav$in.range)] <- 'n'
    
    out.of.range <- cav[cav$in.range=='n',]
  }
  list(mainData = result, out.of.range=out.of.range)
}
