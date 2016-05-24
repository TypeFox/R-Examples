NULL
#' @description Calculation of OIV bio-climatic indices for viticulture (ref.: http://www.oiv.int/oiv/info/enresolution2012?lang=en).
#' 
#' @param daily_Tn series of daily minimum temperature (data frame). Must include the following columns (and names): "year", "month", "day" and one or more station id(s), each column one Tn series.
#' @param daily_Tx series of daily maximum temperature (data frame). Must include the following columns (and names): "year", "month", "day" and one or more station id(s), each column one Tx series.
#' @param daily_Tm series of daily mean temperature (data frame). Must include the following columns (and names): "year", "month", "day" and one or more station id(s), each column one Tm series.
#' @param daily_P series of daily precipitation (data frame). Must include the following columns (and names): "year", "month", "day" and one or more station id(s), each column one P series.
#' @param first.yr first year of the period over which indices are calculated
#' @param last.yr last year of the period over which indices are calculated
#' @param subs_missing logical: if \code{TRUE} (default), missing values in input series are replaced by corresponding average values for each day in the series.
#' @param coeff_HI Huglin's daylength correction, as long as the number of stations, or a single coeffient; in this case it is recycled on all stations. See details.
#' @param quant_Tn_rest quantile (0..1) for the choice of the position of the representative year in the series, referred to the minimum temperature during "rest" period. Default is 0.5 (median)
#' @param quant_Tn_veg quantile (0..1) for the choice of the position of the representative year in the series, referred to the minimum temperature during the "vegetative" period. Default is 0.5 (median)
#' @param quant_Tx_veg quantile (0..1) for the choice of the position of the representative year in the series, referred to the maximum temperature during the "vegetative" period. Default is 0.5 (median)
#' @param indices set of OIV indices to be listed. Default is all indices (1 to 10). See details for numbered list of indices.
#' 
#' @title OIV bioclimatic indices for viticulture
#' @author Emanuele Eccel
#'
#' @return A table (one line per station) reporting OIV bioclimatic indices. A Branas' index is added.
#'  
#' @details General info about OIV indices can be sought at http://www.oiv.int/oiv/info/enresolution2012?lang=en.
#' In general, if \code{first.yr} or \code{last.yr} are NULL (default), the lowest and highest values in series are taken as begin and end of calculation period.
#' A coherence check is done on input of start / end years.
#' 
#' If any input is missing, corresponding indices are not be calculated. If \code{daily_Tm} is missing and both \code{daily_Tn} and  \code{daily_Tx} are passed to the function,
#'  \code{daily_Tm} is calculated as the mean of  \code{daily_Tn} and  \code{daily_Tx}.
#'
#' OIV indices are the following:
#' \code{Tm_veg}: 1 - mean temperature during vegetation period. The latter is April - October (N emisphere) or October - April (S emisphere). The case (N or S) is automatically recognised by temperature data.
#' 
#' \code{WI}: 2 - Winkler index (degree days with a 10 C base, summation over vegetative period), see note on \code{Tm_veg}. Ref: Amerine and Winkler, 1944.
#' 
#' \code{BEDD}: 3 - biologically active degree days. Degree days with a lower threhold of 10 C and an upper threshold at 19 C. Ref: Gladstone, 2004.
#' 
#' \code{HI}: 4 - heliothermic Huglin index. A degree day summation of the average between  mean and maximum temperature above 10 C, corrected by a coefficient of daylength duration.
#' The coefficient is given by the author in a table according to latitude. If one value is given, it is used for all stations. Default is 1.04 (lat. 44-46). Ref.: Huglin, P., 1978
#' \code{CNI}: 5 - cool night index. Mean of September (N emisphere) or March (S emisphere) minimum temperatures. Ref.: Tonietto and Carbonneau, 2004.
#' 
#' \code{FSI}: 6 - Fregoni's simplified index. Given by the product between the summation of thermal range (Tx - Tn) and the number of days with Tm > 10 C, for 30 days before ripening.
#' The period before ripening is taken as September (N emisphere) or March (S emisphere). Ref.: Fregoni et Pezzutto, 2000.
#' 
#' \code{BI}: 7 - Branas' hydrothermic index. The only index not included in OIV list, used for fungine infection proneness. 
#' It is given by the product of precipitation (mm) by the mean temperature (C) during the period April - October (N emisphere) or October - April (S emisphere). Ref: Eynard and Dal Masso, 1990.
#' 
#' \code{Tn_rest}: 8 - mean of minimum temperature during rest period. Useful for assessing winter severity. The rest period is November - March (N emisphere) or May - September (S emisphere). 
#' The case (N or S) is automatically recognised by temperature data.
#' 
#' \code{Tn_veg}: 9 - mean of minimum temperature during vegetative period. Useful for assessing spring frosts. See note on \code{Tm_veg} for periods.
#' 
#' \code{Tx_veg}: 10 - mean of maximum temperature during vegetative period. Useful for assessing summer hot spells. See note on \code{Tm_veg} for periods.
#' 
#' \code{quant_Tn_rest}, \code{quant_Tn_veg}, \code{quant_Tx_veg} define the statistical rank of the year to be chosen as representative for assessing \code{Tn_rest},
#' \code{Tn_veg}, and \code{Tx_veg}, respectfully. 0.5 (default) is the median year, 0 is the minimum (lowest temperature), 1 is the maximum (highest temperature).
#'  
#' The only missing index among those selected by OIV is Riou's Drought Index, which is calculated by function \code{RDI} on monthly series.
#' @importFrom stats aggregate
#' @export
#' 
#' @references Amerine, M.A., and Winkler, A.J., 1944. Composition and quality of musts and wines of California grapes. Hilgardia. 15(6): 493-673.
#' @references Eynard, I. e Dal Masso, G., 1990: Viticoltura moderna. Manuale pratico. Hoepli Milano. 778 pp.
#' @references Fregoni, C., et  Pezzutto, S., 2000 : Principes et premieres approches de l'indice bioclimatique de qualite Fregoni, Progr.Agric.Vitic. 117: 390-396.
#' @references Gladstones, J.S., 2004: Climate and Australian Viticulture. In 'Viticulture. Volume 1-Resources'. (Eds Dry PR, Coombe BG) pp. 90-118.
#' @references Huglin, M.P., 1978: Nouveau mode d'evaluation des possibilites heliothermiques d'un milieu viticole. Comptes Rendus de l'Academie de l'Agriculture de France. 64: 1117-1126.
#' @references Tonietto, J., and Carbonneau, A., 2004. A multicriteria climatic classification system for grape-growing regions worldwide. Agricultural and Forest Meteorology. 124(1/2): 81-97.

#' @examples
#' data(Trent_climate)
#' oiv_ind(daily_Tn=Tn,daily_Tx=Tx, daily_P=P, first.yr=1981, last.yr=2010, subs_missing=FALSE)
#'
#' @seealso \code{\link{RDI}}


oiv_ind<-function (daily_Tn=NULL,daily_Tx=NULL,daily_Tm=NULL, daily_P=NULL, first.yr=NULL, last.yr=NULL, subs_missing=TRUE, coeff_HI=1.04, quant_Tn_rest=0.5,quant_Tn_veg=0.5,quant_Tx_veg=0.5, indices=1:10)
  
{
  
  
  if(is.null(daily_Tm)) daily_Tm <- (daily_Tn + daily_Tx)/2
  
  # checks and sets series' first / last year
  
  first.yr.series<-vector(length=4)
  if(!is.null(daily_Tn)) first.yr.series[1]<-min(daily_Tn$year)
  if(!is.null(daily_Tx)) first.yr.series[2]<-min(daily_Tx$year)
  if(!is.null(daily_Tm)) first.yr.series[3]<-min(daily_Tm$year)
  if(!is.null(daily_P)) first.yr.series[4]<-min(daily_P$year)
  first.yr.series<-first.yr.series[first.yr.series!=0]
  short_series<-FALSE; if(!is.null(first.yr)) short_series<- first.yr < max(first.yr.series)
  if(is.null(first.yr)) {
    if(length(unique(first.yr.series))>1) {
      print("Warning! Non-univoque start year: please pass correct input to function", quote=F)
      first.yr<-max(first.yr.series) } else first.yr<-first.yr.series[1] }      else   if(short_series)
      {
        print("Warning! First year least than series' start year", quote=F)
        first.yr<-max(first.yr.series)
      }
  
  last.yr.series<-vector(length=4)
  if(!is.null(daily_Tn)) last.yr.series[1]<-max(daily_Tn$year)
  if(!is.null(daily_Tx)) last.yr.series[2]<-max(daily_Tx$year)
  if(!is.null(daily_Tm)) last.yr.series[3]<-max(daily_Tm$year)
  if(!is.null(daily_P)) last.yr.series[4]<-max(daily_P$year)
  last.yr.series<-last.yr.series[last.yr.series!=0]
  short_series<-FALSE; if(!is.null(last.yr)) short_series<- last.yr > min(last.yr.series)
  if(is.null(last.yr)) {
    if(length(unique(last.yr.series))>1) {
      print("Warning! Non-univoque end year: please pass correct input to function", quote=F)
      last.yr<-min(last.yr.series) } else last.yr<-last.yr.series[1] }         else if(short_series)
      {
        print("Warning! Last year greater than series' end year", quote=F)
        last.yr<-min(last.yr.series)
      }
  
  if(first.yr > last.yr) print("First year greater than last year!", quote=FALSE) else {
    
    
    if(!is.null(daily_Tn)) daily_Tn <- daily_Tn[daily_Tn$year >= first.yr & daily_Tn$year <= last.yr,]
    if(!is.null(daily_Tx)) daily_Tx <- daily_Tx[daily_Tx$year >= first.yr & daily_Tx$year <= last.yr,]
    if(!is.null(daily_Tm)) daily_Tm <- daily_Tm[daily_Tm$year >= first.yr & daily_Tm$year <= last.yr,]
    if(!is.null(daily_P)) daily_P <- daily_P[daily_P$year >= first.yr & daily_P$year <= last.yr,]
    
    
    # checks S or N emisphere (T in July <> T in Jaunary) and sets periods
    S_emisph <- mean(colMeans(daily_Tm[daily_Tm$month==1,-(1:3)], na.rm=T), na.rm=T) > mean(colMeans(daily_Tm[daily_Tm$month==7,-(1:3)], na.rm=T), na.rm=T)
    if(S_emisph) index_months<-c(10:12, 1:4) else index_months<-4:10
    
    # filters series by index months
    Tm_period<-daily_Tm[daily_Tm$month %in% index_months,]
    if(!is.null(daily_Tn)) Tn_period<-daily_Tn[daily_Tn$month %in% index_months,]
    if(!is.null(daily_Tx)) Tx_period<-daily_Tx[daily_Tx$month %in% index_months,]
    if(!is.null(daily_P)) P_period<-daily_P[daily_P$month %in% index_months,]
    
    # checks first and last years' completeness  
    if(Tn_period$month[1] != index_months[1])
    {
      first.yr<-first.yr+1
      print("Warning: first year incomplete! First year data will be ignored.", quote=FALSE)
      Tm_period<-Tm_period[Tm_period$year >= first.yr,]
      Tn_period<-Tn_period[Tn_period$year >= first.yr,]
      Tx_period<-Tx_period[Tx_period$year >= first.yr,]
      P_period<-P_period[P_period$year >= first.yr,]
    }  
    if(Tn_period$month[length(Tn_period$month)] != index_months[length(index_months)]) 
    {
      last.yr<-last.yr-1
      print("Warning: last year incomplete! Last year data will be ignored.", quote=FALSE)
      Tm_period<-Tm_period[Tm_period$year <= last.yr,]
      Tn_period<-Tn_period[Tn_period$year <= last.yr,]
      Tx_period<-Tn_period[Tx_period$year <= last.yr,]
      P_period<-P_period[P_period$year <= last.yr,]
    }  
    
    print(paste("Indices calculated over", first.yr,"-",last.yr), quote=FALSE)
    
    # gap filling
    if(subs_missing)
    {
      Tn_daily_climate<-aggregate(Tn_period, by=list(Tn_period$day, Tn_period$month), FUN=mean, na.rm=T)[-(1:3)]
      Tx_daily_climate<-aggregate(Tx_period, by=list(Tx_period$day, Tx_period$month), FUN=mean, na.rm=T)[-(1:3)]
      Tm_daily_climate<-aggregate(Tm_period, by=list(Tm_period$day, Tm_period$month), FUN=mean, na.rm=T)[-(1:3)]
      P_daily_climate<-aggregate(P_period, by=list(P_period$day, P_period$month), FUN=mean, na.rm=T)[-(1:3)]
      
      for(i in 1:nrow(Tn_period))
        for(sta in names(Tn_period)[-(1:3)])
        {
          if(is.na(Tn_period[i,sta])) Tn_period[i,sta] <- Tn_daily_climate[Tn_daily_climate$month== Tn_period$month[i] & Tn_daily_climate$day== Tn_period$day[i], sta]
          if(is.na(Tx_period[i,sta])) Tx_period[i,sta] <- Tx_daily_climate[Tx_daily_climate$month== Tx_period$month[i] & Tx_daily_climate$day== Tx_period$day[i], sta]
          if(is.na(Tm_period[i,sta])) Tm_period[i,sta] <- Tm_daily_climate[Tm_daily_climate$month== Tm_period$month[i] & Tm_daily_climate$day== Tm_period$day[i], sta]
        }
      for(i in 1:nrow(P_period))
        for(sta in names(P_period)[-(1:3)])
          if(is.na(P_period[i,sta])) P_period[i,sta] <- P_daily_climate[P_daily_climate$month== P_period$month[i] & P_daily_climate$day== P_period$day[i], sta]
    }
    
    
    
    
    # mean T in the vegetative period
    Tm_veg<-round(colMeans(Tm_period, na.rm=TRUE)[-(1:3)],1)
    
    # Winkler's Growing Degree-Days
    GDD<-function(x)
    {
      diff<- x-10
      diff[diff<0]<-0
      return(diff)
    }
    
    GDD_series<-data.frame(Tm_period[1:3], apply(Tm_period[-(1:3)], FUN=GDD, MARGIN=2))
    winkler<-round(aggregate(GDD_series, by=list(GDD_series$year),FUN=sum, na.rm=T)[-(2:4)],0); names(winkler)[1]<-"year"
    winkler_mean <- round(colMeans(winkler[-1], na.rm=T),0)  
    
    
    # Biologically Effective Degree Days
    BEDD<-function(x)
    {
      diff<- x-10
      diff[diff<0]<-0
      diff[diff>9]<-9
      return(diff)
    }
    BEDD_series<-data.frame(Tm_period[1:3], apply(Tm_period[-(1:3)], FUN=BEDD, MARGIN=2))
    BEDD_table<-round(aggregate(BEDD_series, by=list(BEDD_series$year),FUN=sum, na.rm=T)[-(2:4)],0); names(BEDD_table)[1]<-"year"
    BEDD_mean <- round(colMeans(BEDD_table[-1], na.rm=T),0)
    
    # Huglin's Heliothermic Index
    if(!is.null(daily_Tn) & !is.null(daily_Tx))  
    {
      HDD<-function(Tn,Tx)
      {
        diff_m<- (Tn + Tx)/2 -10
        diff_x<-  Tx -10
        diff_weighted<- (diff_m + diff_x)/2
        diff_weighted[diff_weighted<0]<-0
        return(diff_weighted)
      }
      if(S_emisph) 
      {
        Tn_period_HI<-Tn_period[Tn_period$month != 4,]
        Tx_period_HI<-Tx_period[Tx_period$month != 4,]
      } else
      {
        Tn_period_HI<-Tn_period[Tn_period$month != 10,]
        Tx_period_HI<-Tx_period[Tx_period$month != 10,]
      } 
      if(length(coeff_HI) ==1) correct.coeff<-rep(coeff_HI, (length(Tm_period)-3))
      HI_series_non.correct<-HDD(Tn=Tn_period_HI[-(1:3)], Tx=Tx_period_HI[-(1:3)])
      HI_series<-data.frame(Tn_period_HI[1:3], t(t(HI_series_non.correct)* correct.coeff))
      HI_table<-round(aggregate(HI_series, by=list(HI_series$year),FUN=sum, na.rm=T)[-(2:4)],0); names(HI_table)[1]<-"year"
      HI_mean <- round(colMeans(HI_table[-1], na.rm=T),0)
    } else HI_mean<-rep(NA,length(daily_Tm)-3)
    
    # cool night index
    if(!is.null(daily_Tn))
    {
      if(S_emisph) Tn_period_ripen <- Tn_period[Tn_period$month == 3,] else  Tn_period_ripen <- Tn_period[Tn_period$month == 9,]
      CNI_table<-round(aggregate(Tn_period_ripen, by=list(Tn_period_ripen$year),FUN=mean, na.rm=T)[-(2:4)],1); names(CNI_table)[1]<-"year"
      CNI_mean <- round(colMeans(CNI_table[-1], na.rm=T),1)
    } else CNI_mean<-rep(NA,length(daily_Tm)-3)
    
    
    # Fregoni's simplified index
    if(!is.null(daily_Tn) & !is.null(daily_Tx))
    {
      if(S_emisph) 
      {
        Tn_period_ripen <- Tn_period[Tn_period$month == 3,] 
        Tx_period_ripen <- Tx_period[Tx_period$month == 3,] 
        Tm_period_ripen <- Tm_period[Tm_period$month == 3,] 
      }  else  
      {
        Tn_period_ripen <- Tn_period[Tn_period$month == 9,]
        Tx_period_ripen <- Tx_period[Tx_period$month == 9,]
        Tm_period_ripen <- Tm_period[Tn_period$month == 9,]
      }
      delta.T<-data.frame(Tn_period_ripen[1:3], Tx_period_ripen[-(1:3)] - Tn_period_ripen[-(1:3)])
      delta.T_table <- round(aggregate(delta.T, by=list(delta.T$year),FUN=sum, na.rm=T)[-(2:4)],1); names(delta.T_table)[1]<-"year"
      days.less.than.10_table <- aggregate(Tm_period_ripen, by=list(Tm_period_ripen$year),FUN= function(x) {length(which(x<10))})[-(2:4)]
      FSI_period<-data.frame(delta.T_table[1], delta.T_table[-1] * days.less.than.10_table[-1])
      FSI_mean <- round(colMeans(FSI_period[-1], na.rm=T),1)
    } else FSI_mean<-rep(NA,length(daily_Tm)-3)
    
    
    # Branas' hydrothermic index
    if(!is.null(daily_P))
    {
      if(S_emisph) 
      {
        Tm_period_bran <- Tm_period[Tm_period$month >= 10 | Tm_period$month <= 2,] 
        P_period_bran <- P_period[P_period$month >= 10 | P_period$month <= 2,] 
      }  else  
      {
        Tm_period_bran <- Tm_period[Tm_period$month >= 4 & Tm_period$month <= 8,] 
        P_period_bran <- P_period[P_period$month >= 4 & P_period$month <= 8,] 
      }
      P_period_bran<-P_period_bran[names(P_period_bran) %in% names(Tm_period_bran)]
      Tm_table<-round(aggregate(Tm_period_bran, by=list(Tm_period_bran$year),FUN=mean, na.rm=T)[-(2:4)],1); names(Tm_table)[1]<-"year"
      P_table<-round(aggregate(P_period_bran, by=list(P_period_bran$year),FUN=mean, na.rm=T)[-(2:4)],1); names(P_table)[1]<-"year"
      P_table[-1]<-P_table[-1]*153 # nr. days during sensitivity period 
      # in case P and Tm table have different stations ids:
      P_table<-P_table[names(P_table) %in% names(Tm_table)]
      Tm_table<-Tm_table[names(Tm_table) %in% names(P_table)]
      BI_table<-data.frame(P_table[1], P_table[-1] * Tm_table[-1])
      BI_mean <- round(colMeans(BI_table[-1], na.rm=T),0)
    } else BI_mean<-rep(NA,length(daily_Tm)-3)
    
    # minimum temperature during vegetative rest
    if(!is.null(daily_Tn))
    {
      rest_months<- setdiff(1:12, index_months)
      Tn_period_rest<-daily_Tn[daily_Tn$month %in% rest_months,]
      Tn_rest_table <- round(aggregate(Tn_period_rest, by=list(Tn_period_rest$year),FUN=min, na.rm=T)[-(2:4)],1); names(Tn_rest_table)[1]<-"year"
      Tn_rest_quantile<-round(apply(Tn_rest_table, FUN=quantile, probs=quant_Tn_rest, MARGIN=2)[-1],1)
    } else Tn_rest_quantile<-rep(NA,length(daily_Tm)-3)
    
    
    # minimum temperature during vegetative period
    if(!is.null(daily_Tn))
    {
      Tn_period_veg<-daily_Tn[daily_Tn$month %in% index_months,]
      Tn_veg_table <- round(aggregate(Tn_period_veg, by=list(Tn_period_veg$year),FUN=min, na.rm=T)[-(2:4)],1); names(Tn_veg_table)[1]<-"year"
      Tn_veg_quantile<-round(apply(Tn_veg_table, FUN=quantile, probs=quant_Tn_veg, MARGIN=2)[-1],1)
    } else Tn_veg_quantile<-rep(NA,length(daily_Tm)-3)
    
    # maximum temperature during vegetative period
    if(!is.null(daily_Tx))
    {
      Tx_period_veg<-daily_Tx[daily_Tx$month %in% index_months,]
      Tx_veg_table <- round(aggregate(Tx_period_veg, by=list(Tx_period_veg$year),FUN=max, na.rm=T)[-(2:4)],1); names(Tx_veg_table)[1]<-"year"
      Tx_veg_quantile<-round(apply(Tx_veg_table, FUN=quantile, probs=quant_Tx_veg, MARGIN=2)[-1],1)
    } else Tx_veg_quantile<-rep(NA,length(daily_Tm)-3)
    
    # creates a general table
    oiv_indic_table<-cbind(Tm_veg,winkler_mean,BEDD_mean,HI_mean,CNI_mean,FSI_mean,BI_mean,Tn_rest_quantile,Tn_veg_quantile,Tx_veg_quantile)
    colnames(oiv_indic_table)<-c("Tm_veg", "WI", "BEDD", "HI", "CNI", "FSI", "BI", "Tn_rest", "Tn_veg", "Tx_veg")
    
    if(length(indices) != 10)
      oiv_indic_table<-oiv_indic_table[,indices]
    
    return(oiv_indic_table)
  }
}

