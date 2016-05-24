#' Output of Partial Least Squares analysis results of phenology vs. daily mean
#' temperatures
#' 
#' This function produces figures that illustrate statistical correlations
#' between temperature variation during certain phases and the timing of
#' phenological event, based on a PLS analysis conducted with the PLS_pheno or
#' the PLS_chill_force function.
#' 
#' Ths figure illustrates results from the PLS_pheno function, which uses
#' Partial Least Squares (or Projection to Latent Structures) regression to
#' examine the relationship between mean daily temperatures and the timing of
#' an annual biological event. It produces a plot (as a bmp image) with three
#' panels: the top panel shows the value of the VIP score for each day of the
#' year; the middle panel shows the model coefficients and the bottom panel
#' shows the mean temperature and its standard deviation. In the top plot, all
#' days with VIP scores above VIP_threshold are shown in blue. In the other two
#' panels, values for the same days are shown in red, which high VIP scores
#' coincide with negative model coefficients, and in green for positive
#' coefficients. This function does not produce an output, but as side effects
#' it produces a bmp image and a table that summarizes all data used for making
#' the figure in the specified folder.
#' 
#' @param PLS_output a PLS_output object - the output of the PLS_pheno
#' function. This object is a list with a list element called PLS_summary (and
#' an optional object called PLS_output). This element is a data.frame with the
#' following columns: Date, JDay, Coef, VIP, Tmean, Tstdev. Date is the day of
#' the year in MDD format. JDay is the Julian day (day of the year) of the year
#' in which the biological event is observed; since the analysis will often
#' start in the year before the event, this column often starts with negative
#' numbers. Coef is the coefficient of the PLS regression output. VIP is the
#' Variable Importance in the Projection, another output of the PLS regression.
#' Tmean is the mean observed temperature of the respective day of the year,
#' for the duration of the phenology record. Tstdev is the standard deviation
#' of temperature on a given day of the year over the length of the phenology
#' record.
#' @param PLS_results_path the path where analysis outputs should be saved.
#' Should include the file name, but without suffix.
#' @param VIP_threshold the VIP threshold, above which a variable is considered
#' important. Defaults to 0.8.
#' @param colorscheme color scheme used for plotting. For grayscale image, this
#' should be set to "bw". Otherwise a color plot is produced.
#' @param plot_bloom boolean variable specifying whether the range of bloom
#' dates should be shown in the plots. If set to TRUE, this range is shown by a
#' semi-transparent gray rectangle. The median bloom date is shown as a dashed
#' line. This only works if the full range of bloom dates is visible in the
#' plot, and it should be set to FALSE if anything other than Julian dates are
#' used as dependent variables.
#' @param fonttype font style to be used for the figure. Can be 'serif'
#' (default) or 'sans'.
#' @param add_chill option for indicating the chilling period in the plot. This
#' should be a numeric vector: c(start_chill,end_chill).
#' @param add_heat option for indicating the forcing period in the plot. This
#' should be a numeric vector: c(start_heat,end_heat).
#' @author Eike Luedeling
#' @references The method is described here:
#' 
#' Luedeling E and Gassner A, 2012. Partial Least Squares Regression for
#' analyzing walnut phenology in California. Agricultural and Forest
#' Meteorology 158, 43-52.
#' 
#' Wold S, 1995. PLS for multivariate linear modeling. In: van der Waterbeemd H
#' (ed) Chemometric methods in molecular design: methods and principles in
#' medicinal chemistry, vol 2. Chemie, Weinheim, pp 195-218.
#' 
#' Wold S, Sjostrom M, Eriksson L, 2001. PLS-regression: a basic tool of
#' chemometrics. Chemometr Intell Lab 58(2), 109-130.
#' 
#' Mevik B-H, Wehrens R, Liland KH, 2011. PLS: Partial Least Squares and
#' Principal Component Regression. R package version 2.3-0.
#' http://CRAN.R-project.org/package0pls.
#' 
#' Some applications:
#' 
#' Guo L, Dai J, Wang M, Xu J, Luedeling E, 2015. Responses of spring phenology
#' in temperate zone trees to climate warming: a case study of apricot
#' flowering in China. Agricultural and Forest Meteorology 201, 1-7.
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
#' @keywords phenology analysis
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
#' #Plots look much better with weather<-fix_weather(KA_weather)
#' #but that takes to long to run for passing CRAN checks
#' 
#' PLS_results<-PLS_pheno(
#'   weather_data=weather$weather,
#'   split_month=6,   #last month in same year
#'   bio_data=KA_bloom)
#'   
#' PLS_results_path<-paste(getwd(),"/PLS_output",sep="")
#'   
#' plot_PLS(PLS_results,PLS_results_path)
#' #plot_PLS(PLS_results,PLS_results_path,add_chill=c(307,19),add_heat=c(54,109))
#' 
#' dc<-daily_chill(stack_hourly_temps(weather,50.4), 11)
#' plscf<-PLS_chill_force(daily_chill_obj=dc, bio_data_frame=KA_bloom, split_month=6)
#' 
#' plot_PLS(plscf,PLS_results_path)
#' #plot_PLS(plscf,PLS_results_path,add_chill=c(307,19),add_heat=c(54,109))
#' 
#' 
#' 
#' @export plot_PLS
plot_PLS<-function (PLS_output, PLS_results_path, VIP_threshold = 0.8, 
                    colorscheme = "color",plot_bloom=TRUE,fonttype='serif',
                    add_chill=c(NA,NA),add_heat=c(NA,NA)) 
{
  get_leg<-function(x,leg)
  {leg_SC<-x-leg[1]$yday
  if(leg_SC<0) leg_SC<-leg_SC+365
  return(leg[leg_SC])}

  pheno<-suppressWarnings(as.numeric(as.character(PLS_output$pheno$pheno)))
  
  if (!"object_type" %in% names(PLS_output)) 
    return("Error: not a PLS results object produced with chillR.")
  if (!PLS_output$object_type %in% c("PLS_Temp_pheno", "PLS_chillforce_pheno")) 
    return("Error: not a PLS results object produced with chillR.")
  if (PLS_output$object_type == "PLS_Temp_pheno") {
    PLS_obj <- PLS_output["PLS_summary"]
    pnames <- colnames(PLS_obj$PLS_summary)
    PLS_obj <- as.data.frame(PLS_obj)
    colnames(PLS_obj) <- pnames
    VIPs <- PLS_obj[, "VIP"]
    Coefs <- PLS_obj[, "Coef"]
    lc <- nrow(PLS_obj)
    leg <- 1:nrow(PLS_obj)
    
    yearswitch<-which(trunc(PLS_obj$Date/100)>trunc(PLS_obj$Date/100)[c(2:nrow(PLS_obj),1)])
    if(length(yearswitch>0)) years<-c(rep("2001",yearswitch),
                                      rep("2002",nrow(PLS_obj)-(yearswitch+1)+1)) else years=rep("2001",nrow(PLS_obj))
    leg<-strptime(paste(PLS_obj$Date%%100,"/",trunc(PLS_obj$Date/100),"/",years,sep=""),format = "%d/%m/%Y")
    
    tick_marks <- leg[sapply(strsplit(as.character(leg), 
                                      "-"), "[", 3) == "01"]
    tick_labels <- as.Date(tick_marks, origin = "1999-1-2")
    tick_labels <- as.POSIXlt(tick_labels)$mon
    tick_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[tick_labels + 
                                                                 1]
    mean_clim <- PLS_obj$Tmean
    dev_clim <- PLS_obj$Tstdev
    tick_label_pos <- leg[sapply(strsplit(as.character(leg), 
                                          "-"), "[", 3) == "15"]
    png(paste(PLS_results_path, ".png", sep = ""), width = 2000, 
        height = 2000, pointsize = 20)
    par(family=fonttype)  
    par(mfcol = c(3, 1))
    par(mar = c(6.1, 9.4, 6.1, 2.1))
    par(mgp = c(4, 1.5, 0))
    if (colorscheme == "bw") 
      color_bars <- color_bar_maker(VIPs[1:lc], VIPs[1:lc], 
                                    threshold = VIP_threshold, col1 = "BLACK", col2 = "BLACK", 
                                    col3 = "GREY")  else color_bars <- color_bar_maker(VIPs[1:lc], VIPs[1:lc], 
                                       threshold = VIP_threshold, col1 = "DARK BLUE", col2 = "DARK BLUE", 
                                       col3 = "DARK GREY")
    plot(leg, VIPs[1:lc], main = "VIP", xlab = NA, ylab = NA, 
         xaxs = "i", xaxt = "n", yaxs = "i", cex.lab = 4, 
         cex.axis = 4, cex.main = 5, type = "h", col = color_bars, 
         lwd = 6)
    tick_marks <- grconvertX(tick_marks, from = "user", to = "user")
    X_coords <- grconvertX(leg, from = "user", to = "user")
    ticks_labels <- grconvertX(tick_label_pos, from = "user", 
                               to = "user")
    tick_labels<-tick_labels[1:length(ticks_labels)]
    tick_label_pos<-tick_label_pos[1:length(ticks_labels)]    
    
    axis(1, lwd.ticks = 3, at = tick_marks, labels = FALSE, 
         cex.axis = 4, padj = 1)
    axis(2, lwd.ticks = 3, labels = FALSE)
    box(which = "plot", lwd = 3)
    axis(1, lwd.ticks = 0, at = ticks_labels, labels = tick_labels, 
         cex.axis = 4, padj = 1)
    mtext(side = 2, text = "VIP", line = 6, cex = 3)
    if (colorscheme == "bw") 
      color_bars <- color_bar_maker(VIPs[1:lc], Coefs[1:lc], 
                                    threshold = VIP_threshold, col1 = "BLACK", col2 = "#CECECE", 
                                    col3 = "#7B7B7B") else color_bars <- color_bar_maker(VIPs[1:lc], Coefs[1:lc], 
                                       threshold = VIP_threshold, col1 = "RED", col2 = "DARK GREEN", 
                                       col3 = "DARK GREY")
    if(!is.na(add_chill[1])&!is.na(add_chill[2]))
      rect(get_leg(add_chill[1],leg),-10000,get_leg(add_chill[2],leg),10000,col = rgb(204/255,229/255,1,1/2),border=NA)
    if(!is.na(add_heat[1])&!is.na(add_chill[1]))
      rect(get_leg(add_heat[1],leg),-10000,get_leg(add_heat[2],leg),10000,col = rgb(1,204/255,204/255,1/2),border=NA)
    
    if(plot_bloom) {
      rect(get_leg(get_last_date(pheno,first=T),leg),-10000,
           get_leg(get_last_date(pheno),leg),10000,col = rgb(0.5,0.5,0.5,1/3),border=NA)
      lines(list(x=rep(get_leg(median(pheno,na.rm=TRUE),leg),2),y=c(-10000,10000)),lwd=3,lty="dashed")}
    par(new=TRUE)    
    plot(leg, VIPs[1:lc], main = "VIP", xlab = NA, ylab = NA, 
         xaxs = "i", xaxt = "n", yaxs = "i", cex.lab = 4, 
         cex.axis = 4, cex.main = 5, type = "h", col = color_bars, 
         lwd = 6)

    plot(leg, Coefs[1:lc], main = "Model coefficients", ylab = NA, 
         xlab = NA, xaxs = "i", xaxt = "n", yaxs = "i", cex.lab = 4, 
         cex.axis = 4, cex.main = 5, type = "h", col = color_bars, 
         lwd = 6)
    axis(1, lwd.ticks = 3, at = tick_marks, labels = FALSE, 
         cex.axis = 4, padj = 1)
    axis(2, lwd.ticks = 3, labels = FALSE)
    box(which = "plot", lwd = 3)
    axis(1, lwd.ticks = 0, at = ticks_labels, labels = tick_labels, 
         cex.axis = 4, padj = 1)
    mtext(side = 2, text = "Model coefficient", line = 6, 
          cex = 3)
    if (colorscheme == "bw") 
      color_bars <- color_bar_maker(VIPs[1:lc], Coefs[1:lc], 
                                    threshold = VIP_threshold, col1 = "BLACK", col2 = "#CECECE", 
                                    col3 = "#7B7B7B")
    else color_bars <- color_bar_maker(VIPs[1:lc], Coefs[1:lc], 
                                       threshold = VIP_threshold, col1 = "RED", col2 = "DARK GREEN", 
                                       col3 = "DARK GREY")
    if(!is.na(add_chill[1])&!is.na(add_chill[2]))
      rect(get_leg(add_chill[1],leg),-10000,get_leg(add_chill[2],leg),10000,col = rgb(204/255,229/255,1,1/2),border=NA)
    if(!is.na(add_heat[1])&!is.na(add_chill[1]))
      rect(get_leg(add_heat[1],leg),-10000,get_leg(add_heat[2],leg),10000,col = rgb(1,204/255,204/255,1/2),border=NA)
    if(plot_bloom) {
      rect(get_leg(get_last_date(pheno,first=T),leg),-10000,
           get_leg(get_last_date(pheno),leg),10000,col = rgb(0.5,0.5,0.5,1/3),border=NA)
      lines(list(x=rep(get_leg(median(pheno,na.rm=TRUE),leg),2),y=c(-10000,10000)),lwd=3,lty="dashed")}
    par(new=TRUE)
    plot(leg, Coefs[1:lc], main = "Model coefficients", ylab = NA, 
         xlab = NA, xaxs = "i", xaxt = "n", yaxs = "i", cex.lab = 4, 
         cex.axis = 4, cex.main = 5, type = "h", col = color_bars, 
         lwd = 6)
    
    plot(leg, mean_clim[1:lc], main = "Mean temperature", 
         ylab = NA, xlab = NA, xaxs = "i", yaxs = "i", xaxt = "n", 
         cex.lab = 4, cex.axis = 4, cex.main = 5, type = "l", 
         lwd = 3, col = "BLACK", ylim = c(min(mean_clim[1:lc] - 
                                                dev_clim[1:lc]), max(mean_clim[1:lc] + dev_clim[1:lc])))
    arrows(X_coords, mean_clim[1:lc] + dev_clim[1:lc], X_coords, 
           mean_clim[1:lc] - dev_clim[1:lc], angle = 90, code = 3, 
           lwd = 6, length = 0, col = color_bars)
    lines(leg, mean_clim[1:lc], lwd = 3)
    axis(1, lwd.ticks = 3, at = tick_marks, labels = FALSE, 
         cex.axis = 4, padj = 1)
    axis(2, lwd.ticks = 3, labels = FALSE)
    box(which = "plot", lwd = 3)
    axis(1, lwd.ticks = 0, at = ticks_labels, labels = tick_labels, 
         cex.axis = 4, padj = 1)
    mtext(side = 2, text = expression("Mean temperature ("^"o" * 
                                        "C)"), line = 5, cex = 3)
    if(!is.na(add_chill[1])&!is.na(add_chill[2]))
      rect(get_leg(add_chill[1],leg),-10000,get_leg(add_chill[2],leg),10000,col = rgb(204/255,229/255,1,1/2),border=NA)
    if(!is.na(add_heat[1])&!is.na(add_chill[1]))
      rect(get_leg(add_heat[1],leg),-10000,get_leg(add_heat[2],leg),10000,col = rgb(1,204/255,204/255,1/2),border=NA)
    
    if(plot_bloom) {
      rect(get_leg(get_last_date(pheno,first=T),leg),-10000,
           get_leg(get_last_date(pheno),leg),10000,col = rgb(0.5,0.5,0.5,1/3),border=NA)
      lines(list(x=rep(get_leg(median(pheno,na.rm=TRUE),leg),2),y=c(-10000,10000)),lwd=3,lty="dashed")}
    par(new=TRUE)    
    arrows(X_coords, mean_clim[1:lc] + dev_clim[1:lc], X_coords, 
           mean_clim[1:lc] - dev_clim[1:lc], angle = 90, code = 3, 
           lwd = 6, length = 0, col = color_bars)
    lines(leg, mean_clim[1:lc], lwd = 3)
    
    dev.off()
    write.csv(PLS_obj, paste(PLS_results_path, ".csv", sep = ""), 
              row.names = FALSE)
  }
  if (PLS_output$object_type == "PLS_chillforce_pheno") {
    chill_models<-names(PLS_output)[3:length(PLS_output)]
    heat_models<-names(PLS_output[[3]])
    HM <- heat_models
    for (CM in chill_models)
      for(HM in heat_models)
        { PLS_obj <- PLS_output[[CM]][[HM]]["PLS_summary"]
          pnames <- colnames(PLS_obj$PLS_summary)
          PLS_obj <- as.data.frame(PLS_obj)
          colnames(PLS_obj) <- pnames
          lc <- nrow(PLS_obj)/2
          leg <- 1:(nrow(PLS_obj)/2)
          yearswitch<-which(trunc(PLS_obj$Date/100)>trunc(PLS_obj$Date/100)[c(2:nrow(PLS_obj),1)])
          if(length(yearswitch>0)) years<-c(rep("2001",yearswitch[1]),
                                            rep("2002",lc-(yearswitch[1]+1)+1)) else years=rep("2001",lc)
          leg<-strptime(paste(PLS_obj$Date%%100,"/",trunc(PLS_obj$Date/100),"/",years,sep=""),format = "%d/%m/%Y")[1:lc]
          
          tick_marks <- leg[sapply(strsplit(as.character(leg), 
                                            "-"), "[", 3) == "01"]
          tick_labels <- as.Date(tick_marks, origin = "1999-1-2")
          tick_labels <- as.POSIXlt(tick_labels)$mon
          tick_labels <- c("Jan", "Feb", "Mar", "Apr", "May", 
                           "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[tick_labels + 
                                                                              1]
          
          mean_clim <- PLS_obj$CHmean
          dev_clim <- PLS_obj$CHstdev
          tick_label_pos <- leg[sapply(strsplit(as.character(leg), 
                                                "-"), "[", 3) == "15"]
          
          png(paste(PLS_results_path, "_", CM, "_", HM, ".png", 
                    sep = ""), width = 4000, height = 2000, pointsize = 20)
          par(family=fonttype)  
          par(mfcol = c(3, 2))
          par(mar = c(6.1, 9.1, 6.1, 2.1))
          par(mgp = c(4, 1.5, 0))
          for (graph in c("Chill", "Heat")) {
            if (graph == "Chill") {
              VIPs <- PLS_obj[1:lc, "VIP"]
              Coefs <- PLS_obj[1:lc, "Coef"]
              mean_climg <- PLS_obj[1:lc, "MetricMean"]
              dev_climg <- PLS_obj[1:lc, "MetricStdev"]
            }
            if (graph == "Heat") {
              VIPs <- PLS_obj[-(1:lc), "VIP"]
              Coefs <- PLS_obj[-(1:lc), "Coef"]
              mean_climg <- PLS_obj[-(1:lc), "MetricMean"]
              dev_climg <- PLS_obj[-(1:lc), "MetricStdev"]
            }
            if (colorscheme == "bw") 
              color_bars <- color_bar_maker(VIPs, VIPs, threshold = VIP_threshold, 
                                            col1 = "BLACK", col2 = "BLACK", col3 = "GREY") else
                                              color_bars <- color_bar_maker(VIPs, VIPs, 
                                                                            threshold = VIP_threshold, col1 = "DARK BLUE", 
                                                                            col2 = "DARK BLUE", col3 = "DARK GREY")
                                            plot(leg, VIPs, main = "VIP", xlab = NA, ylab = NA, 
                                                 xaxs = "i", xaxt = "n", yaxs = "i", cex.lab = 4, 
                                                 cex.axis = 4, cex.main = 5, type = "h", col = color_bars, 
                                                 lwd = 6)
                                            tick_marks <- grconvertX(tick_marks, from = "user", 
                                                                     to = "user")
                                            X_coords <- grconvertX(leg, from = "user", to = "user")
                                            ticks_labels <- grconvertX(tick_label_pos, from = "user", 
                                                                       to = "user")
                                            tick_labels<-tick_labels[1:length(ticks_labels)]
                                            tick_label_pos<-tick_label_pos[1:length(ticks_labels)]
                                            axis(1, lwd.ticks = 3, at = tick_marks, labels = FALSE, 
                                                 cex.axis = 4, padj = 1)
                                            axis(2, lwd.ticks = 3, labels = FALSE)
                                            box(which = "plot", lwd = 3)
                                            axis(1, lwd.ticks = 0, at = ticks_labels, labels = tick_labels, 
                                                 cex.axis = 4, padj = 1)
                                            mtext(side = 2, text = "VIP", line = 6, cex = 3)
                                            if (colorscheme == "bw") 
                                              color_bars <- color_bar_maker(VIPs, Coefs, 
                                                                            threshold = VIP_threshold, col1 = "BLACK", 
                                                                            col2 = "#CECECE", col3 = "#7B7B7B")
                                            else color_bars <- color_bar_maker(VIPs, Coefs, 
                                                                               threshold = VIP_threshold, col1 = "RED", col2 = "DARK GREEN", 
                                                                               col3 = "DARK GREY")
                                            if(!is.na(add_chill[1])&!is.na(add_chill[2]))
                                              rect(get_leg(add_chill[1],leg),-10000,get_leg(add_chill[2],leg),10000,col = rgb(204/255,229/255,1,1/2),border=NA)
                                            if(!is.na(add_heat[1])&!is.na(add_chill[1]))
                                              rect(get_leg(add_heat[1],leg),-10000,get_leg(add_heat[2],leg),10000,col = rgb(1,204/255,204/255,1/2),border=NA)
                                            
                                            if(plot_bloom) {
                                              rect(get_leg(get_last_date(pheno,first=T),leg),-10000,
                                                   get_leg(get_last_date(pheno),leg),10000,col = rgb(0.5,0.5,0.5,1/3),border=NA)
                                              lines(list(x=rep(get_leg(median(pheno,na.rm=TRUE),leg),2),y=c(-10000,10000)),lwd=3,lty="dashed")}
                                            
                                            par(new=TRUE)   
                                            plot(leg, VIPs, main = "VIP", xlab = NA, ylab = NA, 
                                                 xaxs = "i", xaxt = "n", yaxs = "i", cex.lab = 4, 
                                                 cex.axis = 4, cex.main = 5, type = "h", col = color_bars, 
                                                 lwd = 6)
                                            
                                            plot(leg, Coefs, main = "Model coefficients", 
                                                 ylab = NA, xlab = NA, xaxs = "i", xaxt = "n", 
                                                 yaxs = "i", cex.lab = 4, cex.axis = 4, cex.main = 5, 
                                                 type = "h", col = color_bars, lwd = 6)
                                            axis(1, lwd.ticks = 3, at = tick_marks, labels = FALSE, 
                                                 cex.axis = 4, padj = 1)
                                            axis(2, lwd.ticks = 3, labels = FALSE)
                                            box(which = "plot", lwd = 3)
                                            axis(1, lwd.ticks = 0, at = ticks_labels, labels = tick_labels, 
                                                 cex.axis = 4, padj = 1)
                                            mtext(side = 2, text = "Model coefficient", line = 6, 
                                                  cex = 3)
                                            if (colorscheme == "bw") 
                                              color_bars <- color_bar_maker(VIPs, Coefs, 
                                                                            threshold = VIP_threshold, col1 = "BLACK", 
                                                                            col2 = "#CECECE", col3 = "#7B7B7B")
                                            else color_bars <- color_bar_maker(VIPs, Coefs, 
                                                                               threshold = VIP_threshold, col1 = "RED", col2 = "DARK GREEN", 
                                                                               col3 = "DARK GREY")
                                            if(!is.na(add_chill[1])&!is.na(add_chill[2]))
                                              rect(get_leg(add_chill[1],leg),-10000,get_leg(add_chill[2],leg),10000,col = rgb(204/255,229/255,1,1/2),border=NA)
                                            if(!is.na(add_heat[1])&!is.na(add_chill[1]))
                                              rect(get_leg(add_heat[1],leg),-10000,get_leg(add_heat[2],leg),10000,col = rgb(1,204/255,204/255,1/2),border=NA)
                                            
                                            if(plot_bloom) {
                                              rect(get_leg(get_last_date(pheno,first=T),leg),-10000,
                                                   get_leg(get_last_date(pheno),leg),10000,col = rgb(0.5,0.5,0.5,1/3),border=NA)
                                              lines(list(x=rep(get_leg(median(pheno,na.rm=TRUE),leg),2),y=c(-10000,10000)),lwd=3,lty="dashed")   
                                            }
                                            par(new=TRUE)
                                            plot(leg, Coefs, main = "Model coefficients", 
                                                 ylab = NA, xlab = NA, xaxs = "i", xaxt = "n", 
                                                 yaxs = "i", cex.lab = 4, cex.axis = 4, cex.main = 5, 
                                                 type = "h", col = color_bars, lwd = 6)
                                            
                                            if (graph == "Chill") 
                                              plot(leg, mean_climg, main = "Chill accumulation", 
                                                   ylab = NA, xlab = NA, xaxs = "i", yaxs = "i", 
                                                   xaxt = "n", cex.lab = 4, cex.axis = 4, cex.main = 5, 
                                                   type = "l", lwd = 3, col = "BLACK", ylim = c(min(mean_climg - 
                                                                                                      dev_climg), max(mean_climg + dev_climg)))
                                            if (graph == "Heat") 
                                              plot(leg, mean_climg, main = "Heat accumulation", 
                                                   ylab = NA, xlab = NA, xaxs = "i", yaxs = "i", 
                                                   xaxt = "n", cex.lab = 4, cex.axis = 4, cex.main = 5, 
                                                   type = "l", lwd = 3, col = "BLACK", ylim = c(min(mean_climg - 
                                                                                                      dev_climg), max(mean_climg + dev_climg)))
                                            arrows(X_coords, mean_climg + dev_climg, X_coords, 
                                                   mean_climg - dev_climg, angle = 90, code = 3, 
                                                   lwd = 6, length = 0, col = color_bars)
                                            lines(leg, mean_climg[1:lc], lwd = 3)
                                            axis(1, lwd.ticks = 3, at = tick_marks, labels = FALSE, 
                                                 cex.axis = 4, padj = 1)
                                            axis(2, lwd.ticks = 3, labels = FALSE)
                                            box(which = "plot", lwd = 3)
                                            axis(1, lwd.ticks = 0, at = ticks_labels, labels = tick_labels, 
                                                 cex.axis = 4, padj = 1)
                                            if (graph == "Chill") {
                                              ccc <- unlist(strsplit(CM, "_"))
                                              cctemp <- ccc[1]
                                              if (length(ccc) > 1) {
                                                for (i in 2:length(ccc)) cctemp <- paste(cctemp, 
                                                                                         ccc[i])
                                              }
                                              mtext(side = 2, text = paste(cctemp, "per day"), 
                                                    line = 5, cex = 3)
                                            }
                                            if (graph == "Heat") {
                                              ccc <- unlist(strsplit(HM, "_"))
                                              cctemp <- ccc[1]
                                              if (length(ccc) > 1) {
                                                for (i in 2:length(ccc)) cctemp <- paste(cctemp, 
                                                                                         ccc[i])
                                              }
                                              mtext(side = 2, text = paste(cctemp, "per day"), 
                                                    line = 5, cex = 3)
                                            }
                                            if(!is.na(add_chill[1])&!is.na(add_chill[2]))
                                              rect(get_leg(add_chill[1],leg),-10000,get_leg(add_chill[2],leg),10000,col = rgb(204/255,229/255,1,1/2),border=NA)
                                            if(!is.na(add_heat[1])&!is.na(add_chill[1]))
                                              rect(get_leg(add_heat[1],leg),-10000,get_leg(add_heat[2],leg),10000,col = rgb(1,204/255,204/255,1/2),border=NA)
                                            
                                            if(plot_bloom) {
                                              rect(get_leg(get_last_date(pheno,first=T),leg),-10000,
                                                   get_leg(get_last_date(pheno),leg),10000,col = rgb(0.5,0.5,0.5,1/3),border=NA)
                                              lines(list(x=rep(get_leg(median(pheno,na.rm=TRUE),leg),2),y=c(-10000,10000)),lwd=3,lty="dashed")}
                                            par(new=TRUE)
                                            arrows(X_coords, mean_climg + dev_climg, X_coords, 
                                                   mean_climg - dev_climg, angle = 90, code = 3, 
                                                   lwd = 6, length = 0, col = color_bars)
                                            lines(leg, mean_climg[1:lc], lwd = 3)
                                            
          }
          dev.off()
          write.csv(PLS_obj, paste(PLS_results_path, "_", CM, 
                                   "_", HM, ".csv", sep = ""), row.names = FALSE)
        }
  }
}
