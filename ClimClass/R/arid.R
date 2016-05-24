NULL
#' @description Calculates aridity according to several indices.
#' 
#' @param clim_norm climatic normals
#' @param coeff_rad mean monthly solar radiation; used only for Thornthwaite's annual index Im. Default is \code{NULL}
#' @param coeff_Hargr (vector of monthly) correction coefficient(s) for Hargreaves' equation
#' @param monthly logic. Sets calculation to the monthly mode if \code{TRUE}. Default is \code{FALSE}.
#' @param indices set of aridity indices to be listed. Default is all indices (1 to 6 for annual, 1 to 2 for monthly).
#'
#' @title Aridity indices
#' @author Emanuele Eccel
#'
#' @return Either a single-line data frame (when \code{monthly = FALSE}) with the desired aridity index(es), or a data frame (\code{monthly = TRUE}), with monthly values of the desired index(es).
#'
#' @details \code{clim_norm} is a monthly data frame of climate normals, with column names: "P", "Tn", "Tx", "Tm" (precipitation, minimum, maximum and mean temperature, respectively). It can be the output of function \code{\link{climate}}.
#' 
#' Monthly potential evapotranspiration (PE) is calculated via the Hargreaves' formula (Hargreaves and Samani, 1985):
#' 
#' PE = (0.0023*(clim_norm$Tx - clim_norm$Tn)^(0.5)*(clim_norm$Tm+17.8)*coeff_rad)* lmv * coeff_Hargr
#' 
#' where Tn, Tx, Tm are min, max, and mean temperatures, respectively, and lmv is the number of days in any month.
#' 
#' \code{coeff_rad} and \code{coeff_Hargr} are needed only by Thornthwaite's annual index \code{Im} and UNEP's \code{Ai} index, whose PE term is calculated via Hargreaves' equation.
#' 
#' \code{coeff_rad} corresponds to the mean monthly extra-atmospheric radiation (see function \code{\link{ExAtRa}}). 
#' 
#' \code{coeff_Hargr} is either a single value or a vector of 12 coefficients to adjust Hargreaves' estimation of potential evapotranspiration (implemented in \code{Im} and \code{Ai} indices). From calibration in 6 stations from the same network of \code{\link{Trent_climate}}, its average value is 0.75.
#' 
#' When \code{monthly} is \code{TRUE}, a data frame with monthly detail is generated for one station, instead of a synthetic single-line data frame.
#' 
#' \code{indices}' values are the following:
#'
#' 1 De Martonne - Ia (annual or monthly). De Martonne, 1925.
#'
#' 2 Thornthwaite - Im (annual or monthly). Thornthwaite, 1948.
#'
#' 3 Emberger - Q (annual only). Emberger, 1955.
#'
#' 4 Lang - R (annual only). Lang, R., 1920.
#'
#' 5 Rivas-Martinez - Io (annual only). Rivas - Martinez, website http://www.iao.florence.it/training/geomatics/BenSlimane/Marocco21_3_1_2.htm
#'
#' 6 UNEP - Ai (annual only). UNEP, 1997.
#'
#' A reference for the aridity degree for any index is given in the list object \code{arid_ind_tables} (see \code{\link{Trent_climate}}.
#'
#' @export
#'
#' @references 
#' De Martonne E., 1925: Traite de Geographie Physique: 3 tomes, Paris.
#'
#' Emberger, L., 1955. Une classification biogeographique des climats. Receuil des travaux des laboratoires de botanique, geologie et zoologie de la faculte des sciences de l'universite de Montpellier (Serie Botanique), Fascicule 7, 3-43.
#' 
#' Hargreaves, G.H., and Samani, Z.A., 1985. Reference crop evapotranspiratin from temperature. Applied Engineering in Agriculture, 1(2):96-99
#'
#' Lang, R., 1920. Verwitterung und Bodenbildung als Einfuehrung in die Bodenkunde. Schweizerbart Science Publishers, Stuttgart
#'
#' Rivas-Martinez - http://www.iao.florence.it/training/geomatics/BenSlimane/Marocco21_3_1_2.htm
#'
#' Thornthwaite, C. W., 1948: An Approach toward a Rational Classification of Climate. Geographical Review, Vol. 38, No. 1(Jan.):55-94.
#'
#' UNEP (United Nations Environment Programme), 1997. World atlas of desertification 2ED. UNEP, London.
#'
#' @examples
#' 
#' data(Trent_climate)
#' # clima_81_10 is a list of data frames having climatic means of temperature and precipitation 
#' # as required by the aridity indices algorithms, each one referring to one station. 
#' # It can be the output of function climate.
#' # coeff_rad is a monthly vector of average daily extra-atmospheric solar radiation, 
#' # calculated e.g. by function ExAtRa.
#' 
#' aridity_Y<-lapply(clima_81_10, coeff_rad=coeff_rad, FUN=arid, monthly=FALSE, indices=c(1,2,5))
#'
#' @seealso \code{\link{climate}}, \code{\link{ExAtRa}}


arid<- function(clim_norm, coeff_rad=NULL, coeff_Hargr = rep(0.75,12), monthly=FALSE, indices= 1:6)
{
  # yearly indices
  if(monthly == FALSE)
  {
    
    # De Martonne
    Ia_y<- sum(clim_norm$P) / mean(clim_norm$Tm + 10)
    # Lang
    R<- sum(clim_norm$P)/mean(clim_norm$Tm)
    # Emberger
    Q<- 2000*sum(clim_norm$P)/(max(clim_norm$Tx + 273.15)^2 - min(clim_norm$Tn  + 273.15)^2) 
    # Rivas-Martinez
    Io <- 10 * sum(clim_norm$P[clim_norm$Tm >0]) / sum(clim_norm$Tm[clim_norm$Tm >0] * 10)
    # Thornthwaite (calculation of ET according to Hargreaves)
    lmv<-c(31,28.25,31,30,31,30,31,31,30,31,30,31)
    ET_hargr<- (0.0023*(clim_norm$Tx - clim_norm$Tn)^(0.5)*(clim_norm$Tm+17.8)*coeff_rad)* lmv * coeff_Hargr
    Im_y<-100*(sum(clim_norm$P) / sum(ET_hargr) -1)
    # UNEP 1997
    Ai <- sum(clim_norm$P) / sum(ET_hargr) 

    aridity<-round(data.frame(Ia = Ia_y, Im = Im_y, Q=Q, R=R, Io=Io, Ai=Ai), 2)[indices] 
    
  } else  # monthly indices
  {
    month_names<-c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    aridity<-NULL
    Ia_m<- round(12* clim_norm$P / (clim_norm$Tm + 10), 2); names(Ia_m)<-1:12
    Im_m<-round(1.65*(clim_norm$P / (clim_norm$Tm + 12.2))^(10/9), 2); names(Im_m)<-1:12 
    aridity$Ia<- Ia_m; names(aridity$Ia)<- month_names
    aridity$Im<- Im_m; names(aridity$Im)<- month_names
    if(length(indices) == 1 & (indices == 1 | indices == 2)[1])
      aridity<-aridity[indices]
  }
  
  return(aridity)
}
