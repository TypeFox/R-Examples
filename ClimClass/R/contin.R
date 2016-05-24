NULL
#' @description Calculates climate continentality / oceanicity according to several indices.
#' 
#' @param clim_norm climatic normals
#' @param latitude station latitude in degrees. Used in Gorczynski's and Conrad's classifications (indices 1 and 2). Default is \code{NULL}.
#' @param elevation station elevation in m. Used in Gams' classification (index 3). Default is \code{NULL}.
#' @param Michalet_correction logic: if \code{TRUE}, Michalet's correction is applied to index 3 (Gams). Default is \code{FALSE}.
#' @param indices set of aridity indices to be listed. Default is all indices (1 to 5).
#'
#' @title Continentality indices
#' 
#' @author Emanuele Eccel
#'
#' @return A single-line data frame with the desired continentality index(es).
#'
#' @details clim_norm is a monthly data frame of climate normals, with column names: "P", "Tn", "Tx", "Tm" (precipitation, minimum, maximum and mean temperature, respectively). It can be the output of function \code{\link{climate}}.
#' 
#' \code{indices}' values are the following:
#'
#' 1: Gorczynski - K.G. (Gorczynski, L., 1920).
#' 
#' 2: Conrad - K.C. (Conrad, 1946).
#' 
#' 3: Gams - alpha. (Gams, H., 1932). For Michalet's correction: Michalet and Souchier, 1991.
#' 
#' 4: Rivas-Martinez - Ic. (Rivas - Martinez, web page).
#' 
#' 5: Amann - H. (Amann, 1929)
#'
#' A reference for the continentality / oceanicity degree is given in the list object \code{continental_ind_tables} of data set \code{\link{Trent_climate}}.
#'
#' If Michalet's correction is applied to Gams' hygric continentality index, the value of precipitation is proportionally diminished for elevations below 900 m a.s.l. See also Lebourgeoise, 2010.
#'
#' @export
#'
#' @references 
#' Amann, J., 1929: L'hygrothermie du climat, facteur determinant la repartition des especes atlantiques. Revue Bryol., 56:126-133.
#'
#' Conrad, V., 1946: Usual formulas of continentality and their limits of validity. Transactions, American Geophysical Union, Volume 27, Issue 5, p. 663-664.
#' 
#' Gams, H., 1932: Die klimatische Begrenzung von Pflanzenarealen und die Verteilung der hygrischen Kontinentalitaet in den Alpen. Zeitschr. Ges. Erdkunde, Berlin.
#' 
#' Gorczynski, L., 1920: Sur le calcul du degre de continentalisme et son application dans la climatologie. Geografiska Annaler 2, 324-331.
#' 
#' Lebourgeoise, F., 2010: Cours de bioclimatologie a l'usage des forestiers. Departement SIAFEE, UFR Forets, Arbres et Milieux Naturels. ENGREF, Nancy Cedex.
#' 
#' Michalet, R., and Souchier, B., 1991: Une approche synthetique biopedoclimatique del montagnes mediterraneennes: l'exemple du Maroc septemptrional. Thesis, Univ. J. Fourier, Grenoble, 273 pp.
#' 
#' Rivas-Martinez: http://www.globalbioclimatics.org/.
#'
#'
#' @examples
#' 
#' data(Trent_climate)
#' 
#' 
#' # clima_81_10 is a list of data frames having climatic means of temperature and precipitation as 
#' # required by the aridity indices algorithms, each one referring to one station. 
#' # It can be the output of function climate.
#' 
#' # creates a data frame with all the continentality indices for all stations in clima_81_10
#' 
#' latit<-coord_elev$North
#' elev<-coord_elev$Elevation
#' 
#' contin_I<-NULL
#' for(i in 1:length(clima_81_10)) {
#'   contin_I[[i]]<-contin(clima_81_10[[i]], 
#'    latitude=latit[i], 
#'    elevation=elev[i], 
#'    Michalet_correction=TRUE)
#' }
#' names(contin_I)<-names(clima_81_10)
#'
#' @seealso \code{\link{climate}}


contin<- function(clim_norm,  latitude=NULL, elevation=NULL, Michalet_correction=FALSE, indices= 1:5)
{
  
  K.G <- NA
  K.C <- NA
  alpha <- NA
  Ic <- NA
  H <- NA
  
  # Gorczynski & Conrad
  if(!is.null(latitude)) 
  {
    K.G <- round((1.7 * (max(clim_norm$Tm)  - min(clim_norm$Tm)) / sin(abs(latitude)*2*pi/360) - 20.4),2)
    K.C <- round((1.7 * (max(clim_norm$Tm)  - min(clim_norm$Tm)) / sin(abs(latitude + 10)*2*pi/360) - 14),2)
  }
  # Gams
  if(!is.null(elevation)) 
    if(Michalet_correction==TRUE & elevation < 900)
      alpha <- round( atan(elevation / (sum(clim_norm$P) - (900 - elevation)/100 *sum(clim_norm$P)/10)) * 360 / (2*pi),2) else
        alpha <- round(    atan(elevation / sum(clim_norm$P))* 360 / (2*pi)  ,2)
  
  # Rivas-Martinez
  Ic<-round((max(clim_norm$Tm) - min(clim_norm$Tm)),2)
  
  # Amann
  
  H <- round(sum(clim_norm$P) * mean(clim_norm$Tm) / (max(clim_norm$Tm)  - min(clim_norm$Tm)), 1)
  
  continentality<-data.frame(K.G = K.G, K.C = K.C, alpha=alpha, Ic = Ic, H=H)
  if(length(indices) != 5)
    continentality<-continentality[indices]
  
  
  return(continentality)  
}
