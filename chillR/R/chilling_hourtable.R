#' Add chilling and heat accumulation to table of hourly temperatures
#' 
#' This function calculates cumulative values for three chill metrics and one
#' heat metric for every hour of an hourly temperature record. The count is
#' restarted on a specified date each year.
#' 
#' Chill metrics are calculated as given in the references below. Chilling
#' Hours are all hours with temperatures between 0 and 7.2 degrees C. Units of
#' the Utah Model are calculated as suggested by Richardson et al. (1974)
#' (different weights for different temperature ranges, and negation of
#' chilling by warm temperatures). Chill Portions are calculated according to
#' Fishman et al. (1987a,b). More honestly, they are calculated according to an
#' Excel sheet produced by Amnon Erez and colleagues, which converts the
#' complex equations in the Fishman papers into relatively simple Excel
#' functions. These were translated into R. References to papers that include
#' the full functions are given below. Growing Degree Hours are calculated
#' according to Anderson et al. (1986), using the default values they suggest.
#' 
#' @param hourtemps a dataframe of stacked hourly temperatures (e.g. produced
#' by stack_hourly_temps). This data frame must have a column for Year, a
#' column for JDay (Julian date, or day of the year), a column for Hour and a
#' column for Temp (hourly temperature).
#' @param Start_JDay the start date (in Julian date, or day of the year) of the
#' calculation for the four metrics. The count is restarted on this date every
#' year.
#' @return data frame consisting of all the columns of the THourly input data
#' frame, plus the following additional columns: Chilling_Hours (cumulative
#' number of Chilling Hours since the last Start_JDay), Chill_Portions (same
#' for units of the Dynamic Models), Chill_Units (same for units of the Utah
#' Model) and GDH (same for Growing Degree Hours).
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Chilling Hours or Utah
#' Models, especially in warm climates! The Dynamic Model (Chill Portions),
#' though far from perfect, seems much more reliable.
#' @author Eike Luedeling
#' @references Model references:
#' 
#' Chilling Hours:
#' 
#' Weinberger JH (1950) Chilling requirements of peach varieties. Proc Am Soc
#' Hortic Sci 56, 122-128
#' 
#' Bennett JP (1949) Temperature and bud rest period. Calif Agric 3 (11), 9+12
#' 
#' Utah Model:
#' 
#' Richardson EA, Seeley SD, Walker DR (1974) A model for estimating the
#' completion of rest for Redhaven and Elberta peach trees. HortScience 9(4),
#' 331-332
#' 
#' Dynamic Model:
#' 
#' Erez A, Fishman S, Linsley-Noakes GC, Allan P (1990) The dynamic model for
#' rest completion in peach buds. Acta Hortic 276, 165-174
#' 
#' Fishman S, Erez A, Couvillon GA (1987a) The temperature dependence of
#' dormancy breaking in plants - computer simulation of processes studied under
#' controlled temperatures. J Theor Biol 126(3), 309-321
#' 
#' Fishman S, Erez A, Couvillon GA (1987b) The temperature dependence of
#' dormancy breaking in plants - mathematical analysis of a two-step model
#' involving a cooperative transition. J Theor Biol 124(4), 473-483
#' 
#' Growing Degree Hours:
#' 
#' Anderson JL, Richardson EA, Kesner CD (1986) Validation of chill unit and
#' flower bud phenology models for 'Montmorency' sour cherry. Acta Hortic 184,
#' 71-78
#' 
#' Model comparisons and model equations:
#' 
#' Luedeling E, Zhang M, Luedeling V and Girvetz EH, 2009. Sensitivity of
#' winter chill models for fruit and nut trees to climatic changes expected in
#' California's Central Valley. Agriculture, Ecosystems and Environment 133,
#' 23-31
#' 
#' Luedeling E, Zhang M, McGranahan G and Leslie C, 2009. Validation of winter
#' chill models using historic records of walnut phenology. Agricultural and
#' Forest Meteorology 149, 1854-1864
#' 
#' Luedeling E and Brown PH, 2011. A global analysis of the comparability of
#' winter chill models for fruit and nut trees. International Journal of
#' Biometeorology 55, 411-421
#' 
#' Luedeling E, Kunz A and Blanke M, 2011. Mehr Chilling fuer Obstbaeume in
#' waermeren Wintern? (More winter chill for fruit trees in warmer winters?).
#' Erwerbs-Obstbau 53, 145-155
#' 
#' Review on chilling models in a climate change context:
#' 
#' Luedeling E, 2012. Climate change impacts on winter chill for temperate
#' fruit and nut production: a review. Scientia Horticulturae 144, 218-229
#' 
#' The PLS method is described here:
#' 
#' Luedeling E and Gassner A, 2012. Partial Least Squares Regression for
#' analyzing walnut phenology in California. Agricultural and Forest
#' Meteorology 158, 43-52.
#' 
#' Wold S (1995) PLS for multivariate linear modeling. In: van der Waterbeemd H
#' (ed) Chemometric methods in molecular design: methods and principles in
#' medicinal chemistry, vol 2. Chemie, Weinheim, pp 195-218.
#' 
#' Wold S, Sjostrom M, Eriksson L (2001) PLS-regression: a basic tool of
#' chemometrics. Chemometr Intell Lab 58(2), 109-130.
#' 
#' Mevik B-H, Wehrens R, Liland KH (2011) PLS: Partial Least Squares and
#' Principal Component Regression. R package version 2.3-0.
#' http://CRAN.R-project.org/package0pls.
#' 
#' Some applications of the PLS procedure:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' 
#' Yu H, Luedeling E and Xu J, 2010. Stronger winter than spring warming delays
#' spring phenology on the Tibetan Plateau. Proceedings of the National Academy
#' of Sciences (PNAS) 107 (51), 22151-22156.
#' 
#' Yu H, Xu J, Okuto E and Luedeling E, 2012. Seasonal Response of Grasslands
#' to Climate Change on the Tibetan Plateau. PLoS ONE 7(11), e49230.
#' 
#' The exact procedure was used here:
#' 
#' Luedeling E, Guo L, Dai J, Leslie C, Blanke M, 2013. Differential responses
#' of trees to temperature variation during the chilling and forcing phases.
#' Agricultural and Forest Meteorology 181, 33-42.
#' 
#' The chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords chill and heat calculation
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2008),])
#' 
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' cht<-chilling_hourtable(hourtemps,20)
#' 
#' @export chilling_hourtable
chilling_hourtable <-
function (hourtemps,Start_JDay)             #hourtemps is a data frame with columns Year, JDay, Hour and Temp
{
  
  if((length(names(hourtemps))==2) & ("hourtemps" %in% names(hourtemps)) & ("QC" %in% names(hourtemps))) 
  {hourtemps<-hourtemps$hourtemps
  QC<-hourtemps$QC}
  
  cols<-colnames(hourtemps)
  
  hourtemps<-hourtemps[which(!is.na(hourtemps[,"Temp"])),]
  
  #Chilling Hours
  CH_range<-which(hourtemps$Temp<=7.2&hourtemps$Temp>=0)
  hourtemps[,"CH_weights"]<-0
  hourtemps[CH_range,"CH_weights"]<-1
  hourtemps[,"CH"]<-0
  #      for (nl in normal_lines) {hourtemps[nl,"CH"]<-hourtemps[nl-1,"CH"]+hourtemps[nl,"CH_weights"]}
  
  #Utah Model
  Utah_range_0.5<-which(hourtemps$Temp<=2.4&hourtemps$Temp>1.4|
    hourtemps$Temp<=12.4&hourtemps$Temp>9.1)
  Utah_range_1.0<-which(hourtemps$Temp<=9.1&hourtemps$Temp>2.4)
  Utah_range_min0.5<-which(hourtemps$Temp<=18.0&hourtemps$Temp>15.9)
  Utah_range_min1.0<-which(hourtemps$Temp>18.0)
  hourtemps[,"Utah_weights"]<-0
  hourtemps[Utah_range_0.5,"Utah_weights"]<-0.5
  hourtemps[Utah_range_1.0,"Utah_weights"]<-1
  hourtemps[Utah_range_min0.5,"Utah_weights"]<-(-0.5)
  hourtemps[Utah_range_min1.0,"Utah_weights"]<-(-1)
  
  
  #Dynamic Model
  e0<-4153.5
  e1<-12888.8
  a0<-139500
  a1<-2567000000000000000
  slp<-1.6
  tetmlt<-277
  aa<-a0/a1
  ee<-e1-e0
  
  
  hourtemps[,"TK"]<-hourtemps$Temp+273
  hourtemps[,"ftmprt"]<-slp*tetmlt*(hourtemps[,"TK"]-tetmlt)/hourtemps[,"TK"]
  hourtemps[,"sr"]<-exp(hourtemps[,"ftmprt"])
  hourtemps[,"xi"]<-hourtemps[,"sr"]/(1+hourtemps[,"sr"])
  hourtemps[,"xs"]<-aa*exp(ee/hourtemps[,"TK"])
  hourtemps[,"ak1"]<-a1*exp(-e1/hourtemps[,"TK"])
  hourtemps[1,"interE"]<-0
  
  memo<-new.env(hash=TRUE)
  
  
  
  posi<-1
  assign(x=paste(1),value=0,envir=memo)
  E=0
  
  xs<-hourtemps[,"xs"]
  xi<-hourtemps[,"xi"]
  ak1<-hourtemps[,"ak1"]
  S<-ak1
  S[1]<-0
  E<-S
  options(scipen=30)
  
  for (l in 2:nrow(hourtemps))  {if(E[l-1]<1)
  {S[l]<-E[l-1]
   E[l]<-xs[l]-(xs[l]-S[l])*exp(-ak1[l])} else
   {S[l]<-E[l-1]-E[l-1]*xi[l-1]
    E[l]<-xs[l]-(xs[l]-S[l])*exp(-ak1[l])}
  }
  hourtemps[,"interE"]<-E
  
  
  
  
  hourtemps[which(hourtemps$interE<1),"delt"]<-0
  hourtemps[which(hourtemps$interE>=1),"delt"]<-hourtemps[which(hourtemps$interE>=1),"interE"]*hourtemps[which(hourtemps$interE>=1),"xi"]
  
  Stress<-1
  Tb<-4
  Tu<-25
  Tc<-36
  
  hourtemps[,"GDH_weight"]<-0
  hourtemps[which(hourtemps$Temp>=Tb&hourtemps$Temp<=Tu),"GDH_weight"]<-Stress*(Tu-Tb)/2*
    (1+cos(pi+pi*(hourtemps[which(hourtemps$Temp>=Tb&hourtemps$Temp<=Tu),"Temp"]-Tb)/(Tu-Tb)))
  hourtemps[which(hourtemps$Temp>Tu&hourtemps$Temp<=Tc),"GDH_weight"]<-Stress*(Tu-Tb)*
    (1+cos(pi/2+pi/2*(hourtemps[which(hourtemps$Temp>Tu&hourtemps$Temp<=Tc),"Temp"]-Tu)/(Tc-Tu)))
  
  
  add_up_weights<-function(hourtemps,outcol,weightcol,SDay)
  {weights<-hourtemps[,weightcol]
   SD<-hourtemps$JDay==SDay
   temp<-weights
   temp[1]<-0
   nn<-nrow(hourtemps)
   for (l in 2:nn)
     if (SD[l]) {temp[l]<-0} else
     {temp[l]<-temp[l-1]+weights[l]}
   hourtemps[,outcol]<-temp  
   return(hourtemps)
  }
  
  hourtemps<-add_up_weights(hourtemps,"Chilling_Hours","CH_weights",Start_JDay)
  hourtemps<-add_up_weights(hourtemps,"Chill_Portions","delt",Start_JDay)
  hourtemps<-add_up_weights(hourtemps,"Chill_Units","Utah_weights",Start_JDay)
  hourtemps<-add_up_weights(hourtemps,"GDH","GDH_weight",Start_JDay)
  
  return(hourtemps[,c(cols,"Chilling_Hours","Chill_Portions","Chill_Units","GDH")])
  
}
