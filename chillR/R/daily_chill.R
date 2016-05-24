#' Calculation of daily chill and heat accumulation
#' 
#' This function calculates daily chill (with three models) and heat
#' accumulation for every day of an hourly temperature record (best generated
#' with stack_hourly_temps). It includes the option to include calculation of a
#' running mean, which smoothes accumulation curves. Especially for the Dynamic
#' Model, this may be advisable, because it does not accumulate chill smoothly,
#' but rather in steps.
#' 
#' Temperature metrics are calculated according to the specified models. They
#' are computed based on hourly temperature records and then summed to produce
#' daily chill accumulation rates.
#' 
#' @param hourtemps a dataframe of stacked hourly temperatures (e.g. produced
#' by stack_hourly_temps). This data frame must have a column for Year, a
#' column for JDay (Julian date, or day of the year), a column for Hour and a
#' column for Temp (hourly temperature).
#' @param running_mean what running mean should be applied to smooth the chill
#' and heat accumulation curves? This should be an odd integer. Use 1 (default)
#' for no smoothing.
#' @param models named list of models that should be applied to the hourly
#' temperature data. These should be functions that take as input a vector of
#' hourly temperatures. This defaults to the set of models provided by the
#' chilling function.
#' @param THourly hourtemps was called THourly in an earlier version of this
#' package. So in order to allow function calls written before the 0.57 update
#' to still work, this is included here.
#' @return a daily chill object consisting of the following elements
#' \item{object_type}{a character string "daily_chill" indicating that this is
#' a daily_chill object} \item{daily_chill}{data frame consisting of the
#' columns YYMMDD, Year, Month, Day and Tmean, plus one column for each model
#' that is evaluated. The latter columns have the name given to the model in
#' the models list and they contain daily total accumulations of the computed
#' metrics.}
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Chilling Hours or Utah
#' Models, especially in warm climates! The Dynamic Model (Chill Portions),
#' though far from perfect, seems much more reliable.
#' @author Eike Luedeling
#' @references Model references for the default models:
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
#' models<-list(CP=Dynamic_Model,CU=Utah_Model,GDH=GDH)
#' 
#' dc<-daily_chill(stack_hourly_temps(fix_weather(KA_weather[which(KA_weather$Year>2009),]),
#'  latitude=50.4),11,models)
#' 
#' 
#' @export daily_chill
daily_chill <-
function(hourtemps=NULL,running_mean=1,models=list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,
                                          Chill_Portions=Dynamic_Model,GDH=GDH),THourly=NULL)
{  if(is.null(hourtemps) & !is.null(THourly)) hourtemps<-THourly
if((length(names(hourtemps))==2) & ("hourtemps" %in% names(hourtemps))) {QC<-hourtemps$QC; hourtemps<-hourtemps$hourtemps} else QC<-NULL


  for(m in 1:length(models))
    hourtemps[,names(models)[m]]<-do.call(models[[m]], list(HourTemp=hourtemps[,"Temp"],summ=F), quote = FALSE, envir = parent.frame())
  
hourtemps$YYMMDD<-hourtemps$Year*10000+hourtemps$Month*100+hourtemps$Day
  
dc<-aggregate(x=hourtemps[,names(models)],by=list(hourtemps$YYMMDD),FUN=sum)
colnames(dc)<-c("YYMMDD",names(models))
  dc$Year<-trunc(dc[,1]/10000)
  dc$Month<-trunc((dc[,1]-dc$Year*10000)/100)
  dc$Day<-trunc((dc[,1]-dc$Year*10000-dc$Month*100))
  dc<-dc[,c("YYMMDD","Year","Month","Day",names(models))]
  dc[,"Tmean"]<-aggregate(x=hourtemps[,"Temp"],by=list(hourtemps$YYMMDD),FUN=mean)[,2]

  dc[,names(models)]<-sapply(names(models),function (x) runn_mean(dc[,x],running_mean))
  
  if("no_Tmin" %in% colnames(hourtemps)) dc[,"no_Tmin"]<-aggregate(x=hourtemps[,"no_Tmin"],by=list(hourtemps$YYMMDD),FUN=sum)[,2]>0
  if("no_Tmax" %in% colnames(hourtemps)) dc[,"no_Tmax"]<-aggregate(x=hourtemps[,"no_Tmax"],by=list(hourtemps$YYMMDD),FUN=sum)[,2]>0
  
  if(!is.null(QC)) if(is.data.frame(QC)) QC<-QC else QC<-NA else QC<-NA
    #chillout<-cbind(chillout,QC[which(QC$End_year %in%chillout$End_year),(which(colnames(QC)=="Data_days")+1):ncol(QC)])
  
  
 return(list(object_type="daily_chill",daily_chill=dc,QC=QC))
}
