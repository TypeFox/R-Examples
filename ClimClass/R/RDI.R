NULL
#' @description Calculation of Riou's drought index described in OIV bio-climatic indices for viticulture (see references)
#' 
#' @param series series of mean monthly weather values
#' @param clim_norm the reference climatic values for each month, used for gap filling. Default is \code{NULL} (no replacement of missing values)
#' @param first.yr of the period over which water balance is calculated. Default is  \code{NULL} (calculations start with the first year of the series)
#' @param last.yr of the period over which water balance is calculated. Default is  \code{NULL} (calculations start with the last year of the series)
#' @param TAW total available water content of soil
#' @param coeff_rad vector of solar radiation coefficients (12 values) for calculation of potential evapotranspiration
#' @param coeff_Hargr (vector of monthly) correction coefficient(s) for Hargreaves' equation
#' @param quant vector of quantiles for the statistical ranking of the year representative for balance (0..1)
#' 
#' @title Riou's drought index
#' @author Emanuele Eccel
#'
#' @return A two-column table reporting Riou's drought indices for each quantile chosen (one line each, minimum is 1). Both "harvest time" and minimum values
#' are calculated (see details).
#' 
#'  
#' @details 
#' For full description of algorithm see OIV standards at http://www.oiv.int/oiv/info/enresolution2012?lang=en and the references:  Riou, 1994; Tonietto, 1999.
#' Evapotranspiration is calculated by Hargreaves' equation (see \code{\link{arid}}).
#' 
#' \code{series} is a data frame of the monthly series (means) of: cumulated precipitation (mm), minimum temperature, maximum temperature, mean temperature (optional) 
#' - all in deg. C. Includes the following columns (and names): "year", "month", "P", "Tn", "Tx", "Tm" (optional), for precipitation, minimum, maximum and 
#' mean temperature, respectively. If \code{Tm} is missing it is calculated as (Tn + Tx)/2. Format is the same of \code{lista_cli}.
#'
#' \code{clim_norm} is a monthly data frame of 12 climate normals, with the same column names of \code{series}, except "year".
#' It can be the output of function \code{\link{climate}}. If \code{clim_norm} is not  \code{NULL}, any missing value in 
#' the monthly series is substituted by the corresponding climatic value in \code{clim_norm}.
#' 
#' A default value of 200 mm for \code{TAW} is suggested by the authors of the index. It can be changed according to the known
#' pedological features of soil.
#' 
#' \code{coeff_rad} corresponds to the mean monthly extra-atmospheric radiation (see function \code{\link{ExAtRa}}). It is required in Hargreaves' equation.
#' 
#' \code{coeff_Hargr} is either a single value or a vector of 12 coefficients to adjust Hargreaves' estimation of potential evapotranspiration. 
#' From calibration in 6 stations from the same network of \code{Trent_climate}, its average value is 0.75.
#' 
#' \code{quant_vector} a vector of minimum one element. 0 yields minimum absolute case, 0.5 the median. Values range from 0 to 1 (inappropropriate if > 0.5).
#' 
#' The algorithm described in OIV assesses water balance at the last month of the ripenining period, early autumn. However, in
#' humid or sub-humid climates the driest period for soil generally falls in summer. For this reason, the output table reports 
#' both cases ("harvest" time value and monthly minimum over the season, "WB_harv" and "WB_min", respectfully). Harvest time is 
#' conventionally September (N emisphere) or March (S emisphere).
#' @importFrom stats quantile
#' @export
#' 
#' @references
#' Riou, C. 1994. Le determinisme climatique de la maturation du raisin: application au zonage de la teneur en sucre dans la Communaute Europeenne (E. Commission, ed.).
#' Office des Publications Officielles des Communautes Europeennes, Luxembourg, 322p.
#' 
#' Tonietto, J. 1999. Les Macroclimats Viticoles Mondiaux et l'Influence du Mesoclimat sur la Typicite de la Syrah et du Muscat de Hambourg dans le Sud de la France 
#' Methodologie de Caracterisation. These de doctorat, Ecole Nationale Superieure Agronomique de Montpellier, Montpellier (France), 216p.
#' 
#' @examples
#' data(Trent_climate)
#' RDI(lista_cli[[1]], clim_norm=clima_81_10[[1]], first.yr=1981, last.yr=2010, coeff_rad=coeff_rad)
#'
#' @seealso \code{\link{oiv_ind}},  \code{\link{arid}}

RDI <- function(series, clim_norm=NULL, first.yr = NULL,last.yr = NULL, TAW=200, coeff_rad, coeff_Hargr =  rep(0.75,12), quant = c(0, 0.1, 0.5))
  
{
  
  k_monthly<-c(0.1,0.3,0.5,0.5,0.5,0.5)
  lmv.12<-c(31, 28.25, 31,30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  # substitutes missing values with those of clim_norm, if passed to function
  if(!is.null(clim_norm))
    for(i in 1:nrow(series))
    {
      if(is.na(series$P[i])) series$P[i] <- clim_norm$P[clim_norm$month == series$month[i]]  
      if(is.na(series$Tn[i])) series$Tn[i] <- clim_norm$Tn[clim_norm$month == series$month[i]]  
      if(is.na(series$Tx[i])) series$Tx[i] <- clim_norm$Tx[clim_norm$month == series$month[i]]
      if(!is.null(series$Tm)) if(is.na(series$Tm[i])) series$Tm[i] <- clim_norm$Tm[clim_norm$month == series$month[i]]
    }
  
  
  # checks S or N emisphere (T in July <> T in January) and sets year's period
  S_emisph <- mean(series$Tn[series$month==1], na.rm=T) > mean(series$Tn[series$month==7], na.rm=T) 
  if(S_emisph) index_months<-c(10:12, 1:3) else index_months<-4:9  
  
  lmv <- lmv.12[index_months]
  coeff_rad <- coeff_rad[4:9] # for both S and N emisph.
  coeff_Hargr <- coeff_Hargr[index_months]
  
  if(S_emisph & !is.null(first.yr)) first.yr<-first.yr-1  
  if(!is.null(first.yr)) series<-series[series$year>= first.yr,]
  if(!is.null(last.yr)) series<-series[series$year<= last.yr,]
  
  series<-series[series$month %in% index_months,]
  
  if(S_emisph) 
  {
    series<- series[!(series$year==min(series$year) & series$month %in% (1:3)),]
    series<- series[!(series$year==max(series$year) & series$month %in% (10:12)),]
    series$year[series$month %in% (10:12)]<-series$year[series$month %in% (10:12)]+1   # ficticious attribution for water balance
  }
  
  if("Tm" %in% names(series) == FALSE) series<-data.frame(series[c("year","month","P","Tn","Tx")],  Tm=(series$Tn + series$Tx)/2)
  Tn_m<-series[,c("year", "month", "Tn")]
  Tx_m<-series[,c("year", "month", "Tx")]
  Tm_m<-series[,c("year", "month", "Tm")]
  P_m<-series[,c("year", "month", "P")]
  
  if(Tm_m$month[1] !=index_months[1] | P_m$month[1] !=index_months[1]) 
  {
    print("Warning: series does not start with first month (Apr. N-emisph., Oct. S-emisph.). First year will be ignored", quote=FALSE)
    Tm_m<-Tm_m[Tm_m$year > min(Tm_m$year),]
    Tn_m<-Tn_m[Tn_m$year > min(Tn_m$year),]
    Tx_m<-Tx_m[Tx_m$year > min(Tx_m$year),]
    P_m<-P_m[P_m$year > min(P_m$year),]
  }
  
  if(Tm_m$month[length(Tm_m$month)] !=index_months[6] | P_m$month[length(P_m$month)] !=index_months[6]) 
  {
    print("Warning: series does not end with last month (Sept. N-emisph., March S-emisph.). Last year will be ignored", quote=FALSE)
    Tm_m<-Tm_m[Tm_m$year < max(Tm_m$year),]
    Tn_m<-Tn_m[Tn_m$year < max(Tn_m$year),]
    Tx_m<-Tx_m[Tx_m$year < max(Tx_m$year),]
    P_m<-P_m[P_m$year < max(P_m$year),]
  }
  
  
  # Hargreaves' PET
  ET_hargr<- cbind(year=Tm_m$year, month=Tm_m$month, (0.0023 * (Tx_m - Tn_m)[-(1:2)]^(0.5) * (Tm_m + 17.8)[-(1:2)] * coeff_rad) * lmv * coeff_Hargr); names(ET_hargr)[3]<-"ET"
  
  # vine traspiration
  Tv<-cbind(year=ET_hargr$year, month=ET_hargr$month, ET_hargr[-(1:2)]*k_monthly)
  
  # soil traspiration
  # max  between P/5 and nr. days in a month
  Max<-data.frame(P_m[1:2], P_m[-(1:2)]/5)   
  for(r in 1:nrow(Max))
    Max$P[r]<-max(Max$P[r], lmv.12[Max$month[r]]) 
  Es<-data.frame(Tv[1:2],round(ET_hargr[-(1:2)] / lmv.12[index_months] * (1-k_monthly) * Max[-(1:2)], 1))
  
  # water balance
  WB<-data.frame(Es[1:2], WB=rep(NA, nrow(Es)))
  WB_y.harv<-data.frame(year=min(Es$year):max(Es$year), matrix(nrow=max(Es$year)-min(Es$year)+1, ncol=length(WB)-2)); WB_y.min<-WB_y.harv
  for(yr in min(WB$year):max(WB$year))
  {
    WB$WB[WB$year==yr & WB$month==index_months[1]]<-TAW
    for(m in 2:length(index_months))
      WB$WB[WB$year==yr & WB$month==index_months[m]]<-min(TAW, max(0,WB$WB[WB$year==yr & WB$month==index_months[m-1]])+P_m$P[P_m$year==yr &  P_m$month==index_months[m]] - Tv$ET[Tv$year==yr &  Tv$month==index_months[m]] - Es$ET[Es$year==yr &  Es$month==index_months[m]])
    WB_y.harv[WB_y.harv$year==yr,]<-WB[WB$year==yr & WB$month==index_months[length(index_months)],][-2]  
    WB_y.min[WB_y.min$year==yr,]<-apply(WB[WB$year==yr,], FUN=min, MARGIN=2)[-2]
  }
  WB_y.harv<-round(WB_y.harv, c(0,1)); names(WB_y.harv)<-c("year", "WB")
  WB_y.min<-round(WB_y.min, c(0,1)); names(WB_y.min)<-c("year", "WB")
  
  # quantile table
  RDI_harv<-round(apply(WB_y.harv[-1], FUN=quantile, na.rm=TRUE, probs=quant, MARGIN=2), 1)
  row.names(RDI_harv)[quant==0]<- "min"; row.names(RDI_harv)[quant==0.5]<- "median"; row.names(RDI_harv)[quant==1]<- "max"
  RDI_min<-round(apply(WB_y.min[-1], FUN=quantile, na.rm=TRUE, probs=quant,  MARGIN=2), 1)
  row.names(RDI_min)[quant==0]<- "min"; row.names(RDI_min)[quant==0.5]<- "median"; row.names(RDI_min)[quant==1]<- "max"
  
  RDI <- data.frame(RDI_harv,  RDI_min); names(RDI)<-c("WB_harv", "WB_min")
  
  return(RDI)
  
}
