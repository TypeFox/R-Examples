#' Calculation of cumulative temperature metric according to a user-defined
#' stepwise weight function
#' 
#' This function calculates heat for temperate trees according to a stepwise
#' model provided by the user.
#' 
#' Temperature-based metric calculated according to the user-defined model.
#' 
#' @param HourTemp Vector of hourly temperatures.
#' @param df data.frame with three columns: lower, upper and weight. lower
#' should contain the lower boundary of a chilling weight interval and upper
#' should contain the upper boundary. weight indicates the weighting to be
#' applied to the respective temperature interval.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative
#' temperature metric over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @keywords chill and heat calculation
#' @examples
#' 
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' 
#' stack<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' df=data.frame(
#'   lower=c(-1000,1,2,3,4,5,6),
#'   upper=c(1,2,3,4,5,6,1000),
#'   weight=c(0,1,2,3,2,1,0))
#' 
#' custom<-function(x) step_model(x,df)
#' 
#' custom(stack$Temp)
#' 
#' models<-list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,
#' Chill_Portions=Dynamic_Model,GDH=GDH,custom=custom)
#' 
#' tempResponse(stack,Start_JDay = 305,End_JDay = 60,models)
#' 
#' @export step_model
step_model<-function(HourTemp,
                     df=data.frame(
                       lower=c(-1000,1.4,2.4,9.1,12.4,15.9,18),
                       upper=c(1.4,2.4,9.1,12.4,15.9,18,1000),
                       weight=c(0,0.5,1,0.5,0,-0.5,-1)),summ=TRUE)
    {lower<-df$lower;upper<-df$upper;weight<-df$weight
      if (summ==TRUE) return(cumsum(sapply(HourTemp,function(x) weight[which(x>lower&x<=upper)]))) else
                      return(sapply(HourTemp,function(x) weight[which(x>lower&x<=upper)]))
    }

#df=data.frame(
#  lower=c(-1000,1,2,3,4,5,6),
#  upper=c(1,2,3,4,5,6,1000),
#  weight=c(0,1,2,3,2,1,0),
#  lower_include_equal=c(F,F,T,T,T,T,T),
#  upper_include_equal=c(F,F,T,T,T,T,T))




#' Calculation of cumulative chill according to the Utah Model
#' 
#' This function calculates winter chill for temperate trees according to the
#' Utah Model.
#' 
#' Units of the Utah Model are calculated as suggested by Richardson et al.
#' (1974) (different weights for different temperature ranges, and negation of
#' chilling by warm temperatures).
#' 
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Utah
#' Chill Units over the entire duration of HourTemp.
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Utah Model, especially in
#' warm climates! The Dynamic Model (Chill Portions), though far from perfect,
#' seems much more reliable.
#' @author Eike Luedeling
#' @references Utah Model reference:
#' 
#' Richardson EA, Seeley SD, Walker DR (1974) A model for estimating the
#' completion of rest for Redhaven and Elberta peach trees. HortScience 9(4),
#' 331-332
#' @keywords chill and heat calculation
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' 
#' stack<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' Utah_Model(stack$hourtemps$Temp)
#' 
#' @export Utah_Model
Utah_Model<-function(HourTemp,summ=TRUE)
  return(step_model(HourTemp,df=data.frame(lower=c(-1000,1.4,2.4,9.1,12.4,15.9,18),upper=c(1.4,2.4,9.1,12.4,15.9,18,1000),weight=c(0,0.5,1,0.5,0,-0.5,-1)),summ=summ))
  


#Utah_Model<-function(HourTemp)
#  {#Utah Model
#  Utah_range_0.5<-which(HourTemp<=2.4&HourTemp>1.4|
#                          HourTemp<=12.4&HourTemp>9.1)
#  Utah_range_1.0<-which(HourTemp<=9.1&HourTemp>2.4)
#  Utah_range_min0.5<-which(HourTemp<=18.0&HourTemp>15.9)
#  Utah_range_min1.0<-which(HourTemp>18.0)
#  Utah_weights<-rep(0,length(HourTemp))
#  Utah_weights[Utah_range_0.5]<-0.5
#  Utah_weights[Utah_range_1.0]<-1
#  Utah_weights[Utah_range_min0.5]<-(-0.5)
#  Utah_weights[Utah_range_min1.0]<-(-1)
#  return(cumsum(Utah_weights))
#}
  


#' Calculation of cumulative chill according to the Chilling Hours Model
#' 
#' This function calculates winter chill for temperate trees according to the
#' Chilling Hours Model.
#' 
#' Chilling Hours are calculated as suggested by Bennett (1949) (all hours with
#' temperatures between 0 and 7.2 degrees C are considered as one Chilling
#' Hour.
#' 
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Chilling
#' Hours over the entire duration of HourTemp.
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Chilling Hours, especially
#' in warm climates! The Dynamic Model (Chill Portions), though far from
#' perfect, seems much more reliable.
#' @author Eike Luedeling
#' @references Chilling Hours references:
#' 
#' Weinberger JH (1950) Chilling requirements of peach varieties. Proc Am Soc
#' Hortic Sci 56, 122-128
#' 
#' Bennett JP (1949) Temperature and bud rest period. Calif Agric 3 (11), 9+12
#' @keywords chill and heat calculation
#' @examples
#' 
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' 
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' Chilling_Hours(hourtemps$hourtemps$Temp)
#' 
#' @export Chilling_Hours
Chilling_Hours<-function(HourTemp,summ=TRUE)
{
  CH_range<-which(HourTemp<=7.2&HourTemp>=0)
  CH_weights<-rep(0,length(HourTemp))
  CH_weights[CH_range]<-1
  if(summ==TRUE) return(cumsum(CH_weights)) else
    return(CH_weights)
}



#' Calculation of cumulative chill according to the Dynamic Model
#' 
#' This function calculates winter chill for temperate trees according to the
#' Dynamic Model.
#' 
#' Chill Portions are calculated as suggested by Erez et al. (1990).
#' 
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Chill
#' Portions over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @references Dynamic Model references:
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
#' @keywords chill and heat calculation
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' 
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' Dynamic_Model(hourtemps$hourtemps$Temp)
#' 
#' @export Dynamic_Model
Dynamic_Model<-function(HourTemp,summ=TRUE)
  {#Dynamic Model
  e0<-4153.5
  e1<-12888.8
  a0<-139500
  a1<-2567000000000000000
  slp<-1.6
  tetmlt<-277
  aa<-a0/a1
  ee<-e1-e0
  
  TK<-HourTemp+273
  ftmprt<-slp*tetmlt*(TK-tetmlt)/TK
  sr<-exp(ftmprt)
  xi<-sr/(1+sr)
  xs<-aa*exp(ee/TK)
  ak1<-a1*exp(-e1/TK)
  interE<-0
  
  memo<-new.env(hash=TRUE)
  
  posi<-1
  assign(x=paste(1),value=0,envir=memo)
  E=0
  
  S<-ak1
  S[1]<-0
  E<-S
  options(scipen=30)
  
  for (l in 2:length(HourTemp))  {if(E[l-1]<1)
  {S[l]<-E[l-1]
  E[l]<-xs[l]-(xs[l]-S[l])*exp(-ak1[l])} else
  {S[l]<-E[l-1]-E[l-1]*xi[l-1]
  E[l]<-xs[l]-(xs[l]-S[l])*exp(-ak1[l])}
  }
  interE<-E
  delt<-rep(0,length(HourTemp))
  delt[which(interE>=1)]<-interE[which(interE>=1)]*xi[which(interE>=1)]
  if (summ==TRUE) return(cumsum(delt)) else
    return(delt)
}



#' Calculation of cumulative heat according to the Growing Degree Hours Model
#' 
#' This function calculates heat for temperate trees according to the Growing
#' Degree Hours Model.
#' 
#' Growing Degree Hours are calculated as suggested by Anderson et al. (1986).
#' 
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Growing
#' Degree Hours over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @references Growing Degree Hours reference:
#' 
#' Anderson JL, Richardson EA, Kesner CD (1986) Validation of chill unit and
#' flower bud phenology models for 'Montmorency' sour cherry. Acta Hortic 184,
#' 71-78
#' @keywords chill and heat calculation
#' @examples
#' 
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' 
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' GDH(hourtemps$hourtemps$Temp)
#' 
#' @export GDH
GDH<-function(HourTemp,summ=TRUE)
  {Stress<-1
  Tb<-4
  Tu<-25
  Tc<-36
  
  GDH_weight<-rep(0,length(HourTemp))
  GDH_weight[which(HourTemp>=Tb&HourTemp<=Tu)]<-Stress*(Tu-Tb)/2*
    (1+cos(pi+pi*(HourTemp[which(HourTemp>=Tb&HourTemp<=Tu)]-Tb)/(Tu-Tb)))
  GDH_weight[which(HourTemp>Tu&HourTemp<=Tc)]<-Stress*(Tu-Tb)*
    (1+cos(pi/2+pi/2*(HourTemp[which(HourTemp>Tu&HourTemp<=Tc)]-Tu)/(Tc-Tu)))
  if (summ) return(cumsum(GDH_weight)) else
    return(GDH_weight)
}


#Johann_model<-function(HourTemp,summ=TRUE)
#  if(summ) return(cumsum(0.6702*(exp(-0.148*HourTemp)))) else
#    return(0.6702*(exp(-0.148*HourTemp)))


#models<-list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,Chill_Portions=Dynamic_Model,GDH=GDH,experimental=step_model,Johann=Johann_model)

#tempResponse(stack,Start_JDay = 305,End_JDay = 60,models)


