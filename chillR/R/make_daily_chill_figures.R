#' Produce image of daily chill and heat accumulation
#' 
#' Function to make figures showing the mean rate of chill and heat
#' accumulation for each day of the year, as well as as the standard deviation.
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
#' This function uses the Kendall package.
#' 
#' @param daily_chill a daily chill object. This should be generated with the
#' daily_chill function.
#' @param file_path the path where data should be saved. Can either end with
#' '/' or include a prefix for all images that are produced.
#' @param models column names of the data.frame stored in daily_chill's
#' daily_chill object that contain the metrics to be plotted. Defaults to four
#' standard metrics of interest in fruit tree phenology analysis.
#' @param labels labels to be used in the plots for the metrics listed under
#' models. This defaults to NA, which means that the character strings given in
#' models are used for the figures. If alternative labels are to be used, these
#' should be given as a vector of length length(models).
#' @return data frame containing all information used to make the figures that
#' are saved. For each Julian Date, means and standard deviations of all chill
#' and heat metrics are saved. In addition, Mann-Kendall tests are performed
#' for daily accumulations of all metrics. p and tau values from this test
#' indicate the level of statistical significance. This non-parametric test is
#' reliable for time series data.
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
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2005),])
#' 
#' dc<-daily_chill(stack_hourly_temps(weather,50.4), 11,models=list(Chill_Portions=Dynamic_Model))
#' 
#' md<-make_daily_chill_figures(dc, paste(getwd(),"/daily_chill_",sep=""),models="Chill_Portions",
#'  labels="Chill Portions")
#' 
#' 
#' @export make_daily_chill_figures
make_daily_chill_figures <-
function(daily_chill,file_path,models=c("Chilling_Hours","Utah_Chill_Units","Chill_Portions","GDH"),
         labels=NA)
{
  if(daily_chill[1]=="daily_chill")
    {dc<-daily_chill$daily_chill
      
     dc[,"JDay"]<-strptime(paste(dc$Month,"/",dc$Day,"/",dc$Year,sep=""),"%m/%d/%Y")$yday+1
     df<-data.frame(JDay=1:365)
     
     Kendall_0<-function(x,y)
     {if(sum(y)==0) out<-list(sl=NA,tau=NA) else out<-Kendall(x,y)
      return(out)}
     
     for(m in models)
       for(i in 1:365) {df[i,paste(m,"_mean",sep="")]<-mean(dc[which(dc$JDay == i),m])
                        df[i,paste(m,"_sd",sep="")]<-sd(dc[which(dc$JDay == i),m])
                        df[i,paste(m,"_Kendall_p",sep="")]<-Kendall_0(dc[which(dc$JDay == i),"Year"],dc[which(dc$JDay == i), m])$sl[1]
                        df[i,paste(m,"_Kendall_tau",sep="")]<-Kendall_0(dc[which(dc$JDay == i),"Year"],dc[which(dc$JDay == i), m])$tau[1]
       }
     

     make_plot<-function(JDay,means,sds,name,label,file_name)
     {
       months<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
       leg<-strptime(strptime(paste("01/",1,"/2001",sep=""),format="%d/%m/%Y")-86400+JDay*86400,"%Y-%m-%d")
       tick_marks<-leg[sapply(strsplit(as.character(leg),"-"),"[",3)=="01"]
       tick_label_pos<-leg[sapply(strsplit(as.character(leg),"-"),"[",3)=="15"]
       tick_labels<-as.Date(tick_marks,origin="1999-1-2")
       tick_labels<-as.POSIXlt(tick_labels)$mon
       tick_labels<-months[tick_labels+1]
       tick_marks<-as.numeric(format(strptime(tick_marks,"%Y-%m-%d"),format="%j"))
       ticks_labels<-as.numeric(format(strptime(tick_label_pos,"%Y-%m-%d"),format="%j"))
       
       png(file_name,width=2000,height=1400,pointsize = 20)
       par(mar=c(6.1,9.1,6.1,2.1))
       plot(JDay,means,main=paste(label,"accumulation"),ylab=NA,xlab=NA,xaxs="i",yaxs="i",xaxt="n",cex.lab=4,cex.axis=4,cex.main=5,type="l",lwd=3,col="BLACK",ylim=c(min(means-sds),max(means+sds)))
       arrows(JDay,means+sds,JDay,means-sds, angle=90, code=3,lwd=6,length=0,col="GRAY")
       lines(JDay,means,lwd=5)
       axis(1,lwd.ticks=3,labels=FALSE,at=tick_marks,cex.axis=4,padj=1);axis(2,lwd.ticks=3,labels=FALSE);box(which="plot",lwd=3);
       axis(1,lwd.ticks=0,at=ticks_labels,labels=tick_labels,cex.axis=4,padj=1);axis(2,lwd.ticks=3,labels=FALSE);box(which="plot",lwd=3);
       mtext(side=2,text=paste(label,"per day"),line=5,cex=4)
       dev.off()
     }
  
     if(is.na(labels[1])) labels<-models
     
     for(m in models)
          make_plot(JDay=df$JDay,means=df[,paste(m,"_mean",sep="")],
               sds=df[,paste(m,"_sd",sep="")],name=m,label=labels[which(models==m)],file_name=paste(file_path,"unit_accumulation_",m,".png",sep=""))
  
     return(list(daily_chill_figure_summary=df))
      
    } else {"Error: not a daily chill object; use function daily_chill to make one from hourly temperature data"}
    
}
