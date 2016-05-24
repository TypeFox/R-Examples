#' Classification of climate according to Koeppen - Geiger, of aridity indices, of continentality 
#' indices, of water balance after Thornthwaite, of viticultural bioclimatic indices. Drawing 
#' climographs: Thornthwaite, Peguy, Bagnouls-Gaussen.
#' 
#' @section ClimClass functions:
#' The \code{\link{koeppen_geiger}} function performs the Koeppen - Geiger's classification, as 
#' described in Trewartha (1980).  Almost all sub-classes have been considered, with the only 
#' exception of those whose attribution is based on qualitative assessment of climatic features.
#'
#' The package collects several criteria for climate classification. The most general is , implemented in function 
#' A classic graphical visualization of temperature and precipitation, according to Bagnouls and 
#' Gaussen (1953), is provided by function \code{\link{bagn_gau}}. A similar, but more sophisticated
#' representation of the same variable, is that of Walter - Lieth (Lieth et al., CD). This 
#' function is implemented in library \code{climatol} (http://www.climatol.eu/).
#'
#' The Function \code{\link{arid}} calculates a set of six annual aridity indices (Emberger, 1955; 
#' Lang, R., 1920; Rivas - Martinez, (website); and UNEP, 1997; De Martonne, 1925; Thornthwaite, 
#' 1948). For the latter two also a monthly index is calculated.
#'
#' A set of four continentality indices is proposed by function \code{\link{contin}} (Gorczynski, 
#' L., 1920; Conrad, 1946; Gams, 1932; Rivas - Martinez, web page).
#'
#' Thornthwaite's method for the assessment of soil water balance (Thornthwaite, 1948; 
#' Thornthwaite and Mather, 1955; Thornthwaite and Mather, 1957) makes use of monthly series to 
#' calculate the main quantities in water balance: evapotranspiration, soil water deficit, soil 
#' water surplus. From these series, quantiles are calculated for every month, to infer climatic 
#' features concerning soil water. Function \code{\link{thornthwaite}} provides such analysis, and 
#' function \code{\link{plot}} manages the plot of the quantiles of the relevant quantities.
#' The assessment of potential evapotranspiration by Thornthwaite and Mather's algorithm requires 
#' the estimation of extra-atmospheric radiation, which is calculated by function 
#' \code{\link{ExAtRa}}, based on the algorithm of Allen et al., 2005.
#'
#' The Function \code{\link{as.datcli}} tranforms a data frame as in example dataset 
#' \code{\link{Trent_climate}} into a data frame format like \code{datcli} in \code{climatol} 
#' package. It can be used to plot Walter - Lieth's climographs (see examples documentation).
#'
#' The Function \code{\link{oiv_ind}} calculates several bioclimatic indices for viticulture 
#' proposed by the International Organization of Viticuture, OIV (Resolution OIV-VITI 423-2012), 
#' plus one index (Branas). 
#'
#' One index of OIV's list, Riou's drought index, needing daily series, is calculated by another 
#' function, \code{\link{RDI}}.
#'
#' The data set included in the library is formed by monthly and daily time series of temperature 
#' and precipitation from Trentino, Italy (courtesy of Autonomous Province of Trento - 
#' Meteotrentino, and of Fondazione Edmund Mach, San Michele all'Adige). Climatic normals are 
#' calculated, too (output of function \code{\link{climate}}). The output of function 
#' \code{\link{thornthwaite}} is present in the data set \code{\link{Trent_climate}}, as input 
#' for function \code{\link{plot}}. Reference tables for aridity and continentality indices are 
#' provided as lists, to rank the classifications on standard scales (\code{\link{arid_ind_tables}} 
#' and \code{\link{continental_ind_tables}}, respectively). See a first application in Eccel et al.,
#'  2015.
#' @docType package
#' @name ClimClass
#' @references Allen, R.G., Walter, I.A., Elliott, R.L., Howell, T.A., Itenfisu, D., Jensen, M.E.,and Snyder, R.L. (eds.), 2005: ASCE Standardized Reference Evapotranspiration Equation. 216 pp.
#' @references Amerine, M.A., and Winkler, A.J., 1944: Composition and quality of musts and wines of California grapes. Hilgardia. 15(6): 493-673.
#' @references  Bagnouls, F., and Gaussen, H., 1953: Saison seche et indice xerothermique. Docum. pour les Cartes des Prod. Veget. Serie: Generalite, 1 (1953), pp. 1-49.
#' @references  Conrad, V. 1946: Usual formulas of continentality and their limits of validity. Transactions, American Geophysical Union, Volume 27, Issue 5, p. 663-664
#' @references De Martonne E., 1925: Traite de Geographie Physique: 3 tomes, Paris.
#' @references Eccel, E., Cordano, E., Zottele, F., Toller, GB., 2015: ClimClass and ClimClassMap: two R- packages for climatic and agro-bioclimatic indices. An application to Trentino. XVIII National Congress of Agrometeorology, 9-11 June 2015, San Michele all?Adige, Book of Extended Abstract (available from Autors).
#' @references Emberger, L., 1955: Une classification biogeographique des climats. Receuil des travaux des laboratoires de botanique, geologie et zoologie de la faculte des sciences de l'universite de Montpellier (Serie Botanique), Fascicule 7, 3-43.
#' @references Eynard, I. e Dal Masso, G., 1990: Viticoltura moderna. Manuale pratico. Hoepli Milano. 778 pp.
#' @references Fregoni, C., et  Pezzutto, S., 2000: Principes et premieres approches de l'indice bioclimatique de qualite Fregoni, Progr.Agric.Vitic. 117: 390-396.
#' @references Gams, H., 1932: Die klimatische Begrenzung von Pflanzenarealen und die Verteilung der hygrischen Kontinentalitaet in den Alpen. Zeitschr. Ges. Erdkunde, Berlin. 
#' @references  Gladstones, J.S., 2004: Climate and Australian Viticulture. In 'Viticulture. Volume 1-Resources'. (Eds Dry PR, Coombe BG) pp. 90-118.
#' @references  Huglin, M.P., 1978: Nouveau mode d'evaluation des possibilites heliothermiques d'un milieu viticole. Comptes Rendus de l'Academie de l'Agriculture de France. 64: 1117-1126.
#' @references Gorczynski, L., 1920: Sur le calcul du degre de continentalisme et son application dans la climatologie. Geografiska Annaler 2, 324-331.
#' @references Hargreaves, G.H., and Samani, Z.A., 1985: Reference crop evapotranspiratin from temperature. Applied Engineering in Agriculture, 1(2):96-99
#' @references Lang, R., 1920: Verwitterung und Bodenbildung als Einfuehrung in die Bodenkunde. Schweizerbart Science Publishers, Stuttgart
#' @references Lebourgeoise, F., 2010: Cours de bioclimatologie a l'usage des forestiers. Departement SIAFEE, UFR Forets, Arbres et Milieux Naturels. ENGREF, Nancy Cedex.
#' @references Lieth, H., Berlekamp, J., Fuest, S., and Riediger, S.: Walter-Lieth: Climate Diagram World Atlas, CD-Series I of Climate and Biosphere, 1st edit.
#' @references Michalet, R., and Souchier, B., 1991: Une approche synthetique biopedoclimatique des montagnes mediterraneennes: l'exemple du Maroc septemptrional. Thesis, Univ. J. Fourier, Grenoble, 273 pp
#' @references Rivas-Martinez: \url{http://www.globalbioclimatics.org/}
#' @references Thornthwaite, C. W., 1948: An Approach toward a Rational Classification of Climate. Geographical Review, Vol. 38, No. 1(Jan.):55-94.
#' @references Thornthwaite, C. W., and Mather, J.R., 1955: The water balance.  Publications in Climatology, Volume 8(1), Laboratory of Climatology
#' @references Thornthwaite, C. W., and Mather, J.R., 1957: Instructions and tables for computing potential evapotranspiration and the water balance.  Publications in climatology, Volume 10(3), Laboratory of Climatology.
#' @references Tonietto, J., and Carbonneau, A., 2004: A multicriteria climatic classification system for grape-growing regions worldwide. Agricultural and Forest Meteorology. 124(1/2): 81-97.
#' @references Trewartha, G.T. and Lyle, H.H., 1980: An Introduction to Climate. MacGraw - Hill, 5th Ed. Appendix: Koeppen's Classification of Climates.
#' @references UNEP (United Nations Environment Programme), 1997: World atlas of desertification 2ED. UNEP, London
NULL