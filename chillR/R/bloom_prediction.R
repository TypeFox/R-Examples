#' Bloom prediction from chilling and forcing requirements, assumed to be
#' fulfilled strictly in sequence
#' 
#' This is a pretty rudimentary function to predict phenological dates from
#' chilling and forcing requirements and hourly chilling and forcing data. Note
#' that there are enormous uncertainties in these predictions, which are hardly
#' ever acknowledged. So please use this function with caution.
#' 
#' This function is a bit preliminary at the moment. It will hopefully be
#' refined later.
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
#' @param HourChillTable a data frame resulting from the chilling_hourtable
#' function.
#' @param Chill_model character string describing the chill model to use.
#' Options are "Chilling_Hours", "Chill_Portions" and "Chill_Units".
#' @param Chill_req numeric parameter indicating the chilling requirement of
#' the particular growth stage (in the unit specified by "Chill_model")
#' @param Heat_req numeric parameter indicating the heat requirement of the
#' particular growth stage (in Growing Degree Hours)
#' @param Start_JDay numeric parameter indicating the day when chill
#' accumulation is supposed to start
#' @return data frame containing the predicted dates of chilling requirement
#' fulfillment and timing of the phenological stage. Columns are Creqfull,
#' Creq_year, Crey_month, Creq_day and Creq_JDay (the row number, date and
#' Julian date of chilling requirement fulfillement), Hreqfull, Hreq_year,
#' Hreq_month, Hreq_day and Hreq_JDay (the row number, date and Julian date of
#' heat requirement fulfillment - this corresponds to the timing of the
#' phenological event.
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
#' @keywords bloom prediction
#' @examples
#' 
#' 
#' hourtemps<-stack_hourly_temps(fix_weather(KA_weather[which(KA_weather$Year>2008),]),latitude=50.4)
#' 
#' CT<-chilling_hourtable(hourtemps,Start_JDay=305)
#' 
#' bloom_prediction(CT,Chill_model="Chill_Portions",Chill_req=60,Heat_req=5000,Start_JDay=305)
#' 
#' 
#' @export bloom_prediction
bloom_prediction <-
function (HourChillTable, Chill_model, Chill_req, Heat_req,Start_JDay=305) 
{
  
  cchh<-HourChillTable$Chilling_Hours
  ccpp<-HourChillTable$Chill_Portions
  ccuu<-HourChillTable$Chill_Units
  sea<-HourChillTable$Season
  stdd<-HourChillTable$JDay
  for (s in unique(HourChillTable$Season))
      {cchh[which(sea==s)]<-cchh[which(sea==s)]-cchh[which(sea==s&stdd==round(Start_JDay))][1]
       ccpp[which(sea==s)]<-ccpp[which(sea==s)]-ccpp[which(sea==s&stdd==round(Start_JDay))][1]
       ccuu[which(sea==s)]<-ccuu[which(sea==s)]-ccuu[which(sea==s&stdd==round(Start_JDay))][1]
       }
  HourChillTable$Chilling_Hours<-cchh
  HourChillTable$Chill_Portions<-ccpp
  HourChillTable$Chill_Units<-ccuu
  
  results <- data.frame()
  HCT <- HourChillTable
  chill <- HCT[, Chill_model]
  chill2 <- chill[c(2:length(chill), 1)]
  Creqfull <- which((chill<=Chill_req)&(chill2>Chill_req)) + 1
  Creqfull <- Creqfull[!chill[Creqfull] == 0]
  Creqfull <- Creqfull[which(!is.na(Creqfull))]
  results <- data.frame(Creqfull = Creqfull)
  results[, c("Creq_year", "Creq_month", "Creq_day", "Creq_JDay")] <- HCT[Creqfull, 
                                                                          c("Year", "Month", "Day", "JDay")]
  results[, "Hreqfull"]<-rep(NA,nrow(results))
  results[, "Hreq_year"]<-rep(NA,nrow(results))
  results[, "Hreq_month"]<-rep(NA,nrow(results))
  results[, "Hreq_day"]<-rep(NA,nrow(results))
  results[, "Hreq_JDay"]<-rep(NA,nrow(results))
  
  if(nrow(results)>0)
  {for (i in 1:nrow(results)) {
    if (!i == nrow(results)) 
      tabend <- results$Creqfull[i + 1]
    else tabend <- nrow(HCT)
    temp <- HCT[results$Creqfull[i]:tabend, ]
    temp[, "GDH"] <- temp[, "GDH"] - temp[1, "GDH"]
    GDH <- temp[, "GDH"]
    GDH2 <- GDH[c(2:length(GDH), 1)]
    Hreqfull <- which((GDH - Heat_req) * (GDH2 - Heat_req) < 
                        0) + 1
    Hreqfull <- Hreqfull[!GDH[Hreqfull] <= 0]
    Hreqfull <- Hreqfull[which(!is.na(Hreqfull))][1]
    results[i, c("Hreqfull", "Hreq_year", "Hreq_month", "Hreq_day", 
                 "Hreq_JDay")] <- c(Hreqfull, temp[Hreqfull, c("Year", 
                                                               "Month", "Day", "JDay")])
  }
   results <- results[which(results[, "Hreq_year"] - results[, 
                                                             "Creq_year"] < 2), ]
  }
  return(results)
}
