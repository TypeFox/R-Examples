NULL
#' @description Calculates Thornthwaite and Mather's water balance from monthly series of precipitation and temperature. Aimed at a classification of a site's climate according to its water balance features.
#' 
#' @param series the monthly series of temperature and precipitation. 
#' @param latitude latitude of the station in degrees.
#' @param clim_norm climatic normals.
#' @param first.yr first year of the period over which water balance is calculated. Default is \code{NULL} (calculations start with the first year of the series).
#' @param last.yr last year of the period over which water balance is calculated. Default is \code{NULL} (calculations stop with the last year of the series).
#' @param quant vector of quantiles for which water balance has to be assessed. Default is: min, 10th, 25th 50th, 75th, 90th, max.
#' @param snow.init initial water equivalent for snowpack (mm). Default is 20.
#' @param Tsnow maximum temperature (monthly mean) for precipitation to be treated as snowfall. Default is -1 degree C.
#' @param TAW maximum (field capacity) for soil water retention, and initial soil water content (mm). Default is 100.
#' @param fr.sn.acc  fraction of snow that contributes to snowpack (0-1). 1 - fr.sn.acc is treated as liquid monthly precipitation Default is 0.95.
#' @param snow_melt_coeff monthly coefficient(s) for snowmelt. Default is 1.
#' 
#'@import geosphere
#' 
#' 
#' @title Thornthwaite and Mather's water balance
#' 
#' @author Giambattista Toller and Emanuele Eccel
#'
#' @return A \code{thornthwaite} S3 object, consisting on a list of two lists. The first (name: W_balance) is a list of data frames containing the monthly series of all indices, the second (name: quantiles) the relevant quantiles. See details for meanings of single variables.
#'
#' @details The algorithm for the calculation of water balance is adapted from Thornthwaite, 1948; Thornthwaite and Mather, 1955; Thornthwaite and Mather, 1957.
#' 
#' \code{series} is a data frame with years, months, temperature and precipitation values. Names in series columns must include: year, month, Tn and Tx (minimum and maximum temperatures, respectively) or, as an alternative, Tm (mean temperatures), and P (mandatory).
#' 
#' \code{clim_norm} is a monthly data frame of climate normals, with column names: "P", "Tn", "Tx", "Tm" (precipitation, minimum, maximum and mean temperature, respectively). It can be the output of function \code{\link{climate}}. If \code{clim_norm} is not NULL, any missing value in the monthly series is substituted by the corresponding climatic value in \code{clim_norm}.
#' 
#' At any winter season, the maximum monthly snowpack height is attained in the last month before "spring" conditions (\code{Tm} >= \code{Tsnow}), even if a month with Tm < Tsnow may occur later.
#' 
#' \code{snow_melt_coeff} is (are) the coefficient(s) for snow melt fraction(s) at any month where the condition for melting exists. If \code{snow_melt_coeff} = 1 (default), all the melting occurs in the first month when \code{Tm >= Tsnow}; if it is a vector, melting is spread over more than one month. If the sum of coefficients is less than 1, the residual melting occurs in one further month.
#' 
#' The output function is a list of two lists of data frames (balance and quantile). In both lists, data frame (and names) are the following (all variables in mm):
#' 
#' \code{Precipitation} (repeats input values);
#' 
#' \code{Et0} (potential evapotranspiration);
#' 
#' \code{Storage} (water stored in soil);
#' 
#' \code{Prec. - Evap.} (difference between precipitation and potential evapotranspiration);
#' 
#' \code{Deficit} (difference between potential and real evapotranspiration, due to water unavailability in soil);
#' 
#' \code{Surplus} (water surplus in soil, routed to runoff).
#' 
#' Please, refer to the quoted references for details.
#' 
#' This function requires the function \code{\link{daylength}} (libr. \code{\link{geosphere}}).
#' @export
#'
#' @references 
#' Thornthwaite, C. W., 1948: An Approach toward a Rational Classification of Climate. Geographical Review, Vol. 38, No. 1(Jan.):55-94.
#' 
#' Thornthwaite, C. W., and Mather, J.R., 1955: The water balance.  Publications in Climatology, Volume 8(1), Laboratory of Climatology
#' 
#' Thornthwaite, C. W., and Mather, J.R., 1957: Instructions and tables for computing potential evapotranspiration and the water balance.  Publications in climatology, Volume 10(3), Laboratory of Climatology
#' 
#'
#' @examples
#' 
#' data(Trent_climate)
#' 
#' 
#' # lista_cli is a list of data frames of the type "series", 
#' # each one referring to one station - see function "climate".
#' # clima_81_10 is a list of data frames having climatic means 
#' # of temperature and precipitation, each one referring to one station. 
#' # It can be the output of function "climate".
#' library(geosphere) # required for function daylength
#' thornt_lst<-NULL
#' lista_cli <- lista_cli[1:3] ## lista_cli is reduced to diminish elapsed time of execution!
#' for(k in 1 : length(lista_cli[1:3])) {
#'   thornt_lst[[k]]<-thornthwaite(series=lista_cli[[k]], 
#'   clim_norm=clima_81_10[[k]],
#'   latitude = 46, first.yr=1981, 
#'   last.yr=2010, snow_melt_coeff=c(0.5,0.5 )  )
#' }
#' names(thornt_lst)<-names(lista_cli)
#'   
#' # splits list into two lists
#' W_balance<-NULL; quantiles<-NULL
#' for(k in 1 : length(lista_cli))
#' {
#'   W_balance[[k]]<-thornt_lst[[k]]$W_balance
#'   quantiles[[k]]<-thornt_lst[[k]]$quantiles
#'  }
#'  names(W_balance)<-names(thornt_lst); names(quantiles)<-names(thornt_lst)
#'  
#' @seealso \code{\link{climate}}, \code{\link{ExAtRa}}, \code{\link{plot.thornthwaite}}

thornthwaite <- function(series,  latitude, clim_norm=NULL, first.yr=NULL, last.yr=NULL, quant=c(0,0.10,0.25,0.50,0.75,0.90,1.00),
                         snow.init=20, Tsnow=-1,  TAW=100, fr.sn.acc=0.95,  snow_melt_coeff=1)
  
{
  month_names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  month_lengths<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  m<-c("min",paste(quant[- c(1, length(quant))]*100, "%", sep=""),"max")  
  # creates a set of DOYs centered in each month
  yDay  <- strptime(paste(15,"/",1:12,"/",2014,sep=""), format="%d/%m/%Y")$yday + 1  
  rel_daylength<-daylength(latitude,yDay)/12
  # calculates extra-atmospheric radiation for each day in the set
  chg<- ExAtRa(DOY=yDay,latitude=latitude, unit="mm")
  # calculates total available water
  fTAW<-(-0.000000105*TAW^2+0.000054225*TAW-0.00878)
  # first storage value = TAW
  ultStorPrA<-TAW
  # determines calculation period
  if(is.null(first.yr))  first.yr<-min(series$year) else first.yr<-max(first.yr, min(series$year))
  if(is.null(last.yr))  last.yr<-max(series$year) else last.yr<-min(last.yr, max(series$year))
  years<-first.yr:last.yr
  
  tPrec<-seq(1:12)
  tET0<-seq(1:12)
  tStor<-seq(1:12)
  tPmE<-seq(1:12)
  tDef<-seq(1:12)
  tSur<-seq(1:12)
  
  # creates a data frame
  
  snowpack_jan<-snow.init
  
  for(ii in 1:length(years))
  {
    zzz<-data.frame(series[series$year==years[ii],])
    zzz<-data.frame(zzz[3], (zzz$Tn + zzz$Tx)/2, zzz$Tn, zzz$Tx); names(zzz)<-c("Prec1", "Tmean1","Tn1", "Tx1")
    # substitutes missing values with those of clim_norm, if passed to function
    if(!is.null(clim_norm))
    {
      zzz$Prec1[is.na(zzz$Prec1)]<-clim_norm$P[is.na(zzz$Prec1)]
      zzz$Tmean1[is.na(zzz$Tmean1)]<-clim_norm$Tm[is.na(zzz$Tmean1)]
      zzz$Tn1[is.na(zzz$Tn1)]<-clim_norm$Tn[is.na(zzz$Tn1)]
      zzz$Tx1[is.na(zzz$Tx1)]<-clim_norm$Tx[is.na(zzz$Tx1)]
    }
    # check on residual missing values
    if(sum(is.na(zzz$Prec1))!=0 | sum(is.na(zzz$Tmean1))!=0) print(paste("Monthly values missing in year", first.yr+ii-1), quote=FALSE) else
    {
      Prec  <-zzz[,1]
      Tmean <-zzz[,2]
      Tn    <-zzz[,3]
      Tx    <-zzz[,4]
      
      Storage<-rep(TAW,length(Tmean))
      Percola<-rep(0,length(Tmean))
      
      # calculation of thermal index
      TmeanCor<-Tmean
      TmeanCor[TmeanCor<0]<-0
      ITM<-((TmeanCor/5)^1.514) # monthly
      ITA<-sum(ITM) # annual
      
      # calculation of PET    
      exp<- 0.000000675*ITA^3 - 0.0000771 *ITA^2 + 0.01792*ITA + 0.492
      PET<-(16*(10*TmeanCor/ITA)^exp)*rel_daylength
      # snow accumulation
      SnowAcc<-rep(0,12)
      SnowPrec<-Prec; SnowPrec[Tmean>=Tsnow]<-0 
      
      
      # calculates maximum snowpack before first thaw month
      month_max<-max(min(which(Tmean>=Tsnow)-1), 1) # month before first with Tmean >= Tsnow
      SnowAcc_wint<-NULL
      for(i in 1:month_max) 
      {
        if(i==1) SnowAcc_wint[i] <- snowpack_jan + SnowPrec[i]*fr.sn.acc  else
          SnowAcc_wint[i]<-SnowAcc_wint[i-1] + SnowPrec[i]*fr.sn.acc
      }
      snowpack_ref<-SnowAcc_wint[month_max]
      
      
      snow_depl<-NULL
      SnowAcc[1] <- SnowAcc_wint[1]
      snow_depl[1] <- SnowPrec[1]*(1-fr.sn.acc) # 0 if Tmean>=Tsnow
      if(Tmean[1]>=Tsnow)
      {
        snow_depl[1] <- snow_melt_coeff[1]*SnowAcc[1]
        SnowAcc[1]<- SnowAcc[1]-snow_depl[1]
      }
      
      
      # following months
      count_thaw<-0
      for(i in 2:12)
      {
        
        snow_depl[i]<-  SnowPrec[i]*(1-fr.sn.acc) # 0 if Tmean>=Tsnow
        SnowAcc[i]<- SnowAcc[i-1]+SnowPrec[i]*fr.sn.acc
        if(Tmean[i]>=Tsnow) 
        {
          count_thaw<-count_thaw+1
          if(count_thaw>length(snow_melt_coeff)) 
          {
            snow_depl[i]<- SnowAcc[i]
          } else
            snow_depl[i]<- snow_depl[i] + SnowAcc[i] *snow_melt_coeff[count_thaw]
        }
        SnowAcc[i]<- SnowAcc[i]-snow_depl[i]
        
      } 
      
      snowpack_jan<-SnowAcc[12]
      
      
      # calculation of difference between Prec and PET        
      Liq_Prec<- Prec; Liq_Prec[Tmean<Tsnow] <- Prec[Tmean<Tsnow]*(1-fr.sn.acc)
      P.minus.ET<- Liq_Prec + snow_depl - PET
      
      # water balance
      if(ii==1)   # first year, field capacity
        last_Storage<-TAW
      # first month
      if(!is.na(P.minus.ET[1]) & P.minus.ET[1]>0)
      {
        Storage[1]<-last_Storage+P.minus.ET[1]
        if(Storage[1]>TAW)
        {
          Percola[1]<-Storage[1]-TAW
          Storage[1]<-TAW
        }
      }  else
      {
        PETvir<-(log10(TAW)-log10(last_Storage))/fTAW # virtual PET
        Storage[1]<-TAW*10^(-(PETvir + P.minus.ET[1])*fTAW)
      }
      # following months
      for(i in 2:length(Storage))
      {
        if(!is.na(P.minus.ET[i]) & P.minus.ET[i]>0)
        {
          Storage[i]<-Storage[i-1]+P.minus.ET[i]
          if(Storage[i]>TAW)
          {
            Percola[i]<-Storage[i]-TAW
            Storage[i]<-TAW
          }
        } else
        {
          PETvir<-(log10(TAW)-log10(Storage[i-1]))/fTAW
          Storage[i]<-TAW*10^(-(PETvir+P.minus.ET[i])*fTAW)
        }
      } 
      # saves the last storage
      last_Storage<-Storage[12]
      # calculates real ET as difference between prec. and storage
      delta.sto<-c(0,diff(Storage))
      ETr<-(Liq_Prec + snow_depl - delta.sto)
      for(i in 1:length(ETr))
      {
        if(P.minus.ET[i]>0)
          ETr[i]<-PET[i]
      }  
      
      # deficit: difference between real and potential ET (non-negative only)
      Def <- PET - ETr; Def[Def<0]<-0
      
      # water surplus
      Sur <- Liq_Prec + snow_depl - ETr - (TAW - Storage)
      Sur[Sur<0]<-0
      
      
      tPrec<-data.frame(tPrec, Prec)
      tET0<-data.frame(tET0, PET)
      tStor<-data.frame(tStor, Storage)
      tPmE<-data.frame(tPmE, P.minus.ET)
      tDef<-data.frame(tDef, Def)
      tSur<-data.frame(tSur, Sur)
      
    } # end else (no missing monthly data)
  } # end cycle on years
  
  names(tPrec)<-c("Month", years); row.names(tPrec)<-month_names
  names(tET0)<-c("Month", years); row.names(tET0)<-month_names
  names(tStor)<-c("Month", years); row.names(tStor)<-month_names
  names(tPmE)<-c("Month", years); row.names(tPmE)<-month_names
  names(tDef)<-c("Month", years); row.names(tDef)<-month_names
  names(tSur)<-c("Month", years); row.names(tSur)<-month_names
  
  list_thornt<-list(round(tPrec[-1],1), round(tET0[-1],1), round(tStor[-1],1), round(tPmE[-1],1), round(tDef[-1],1), round(tSur[-1],1))
  names(list_thornt)<- c("Precipitation", "Et0", "Storage", "Prec. - Evap.", "Deficit", "Surplus")
  
  # quantile tables
  
  qPrec<-data.frame(quant)
  qET0<-data.frame(quant)
  qStor<-data.frame(quant)
  qPmE<-data.frame(quant)
  qDef<-data.frame(quant)
  qSur<-data.frame(quant)
  
  
  ttPrec<-t(tPrec)
  ttET0<-t(tET0)
  ttStor<-t(tStor)
  ttPmE<-t(tPmE)
  ttDef<-t(tDef)
  ttSur<-t(tSur)
  
  qPrec<-as.data.frame(apply(ttPrec[-1,],FUN=quantile,probs=quant,na.rm=T, MARGIN=2)); names(qPrec)<-month_names
  qET0<-as.data.frame(apply(ttET0[-1,],FUN=quantile,probs=quant,na.rm=T, MARGIN=2)); names(qET0)<-month_names
  qStor<-as.data.frame(apply(ttStor[-1,],FUN=quantile,probs=quant,na.rm=T, MARGIN=2)); names(qStor)<-month_names
  qPmE<-as.data.frame(apply(ttPmE[-1,],FUN=quantile,probs=quant,na.rm=T, MARGIN=2)); names(qPmE)<-month_names
  qDef<-as.data.frame(apply(ttDef[-1,],FUN=quantile,probs=quant,na.rm=T, MARGIN=2)); names(qDef)<-month_names
  qSur<-as.data.frame(apply(ttSur[-1,],FUN=quantile,probs=quant,na.rm=T, MARGIN=2)); names(qSur)<-month_names
  
  list_quant<-list(round(qPrec,1),round(qET0,1),round(qStor,1),round(qPmE,1),round(qDef,1),round(qSur,1))
  names(list_quant)<- names(list_thornt)
  
  list_2.lists<-list(W_balance=list_thornt, quantiles=list_quant)
  class(list_2.lists) <- "thornthwaite"

  return(list_2.lists)
  
}
