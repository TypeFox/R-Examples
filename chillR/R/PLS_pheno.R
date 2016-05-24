#' Partial Least Squares analysis of phenology vs. daily mean temperatures
#' 
#' This function conducts a Partial Least Squares (PLS) regression analysis
#' relating an annual biological phenomenon, e.g. fruit tree flowering or leaf
#' emergence, to mean daily temperatures of the preceding 12 months. It
#' produces figures that illustrate statistical correlations between
#' temperature variation during certain phases and the timing of phenological
#' event.
#' 
#' PLS regression is useful for exploring the relationship between daily
#' temperature data and biological phenomena that only occur once per year. The
#' statistical challenge is that a normally quite small number of observations
#' must be related to variation in a much larger number (365) of daily
#' temperatures, which are also highly autocorrelated. Most regression
#' approaches are not suitable for this, but PLS regression offers a potential
#' solution. The method is frequently used in chemometrics and hyperspectral
#' remote sensing, where similar statistical challenges are encountered. The
#' basic mechanism is that PLS first constructs latent factors (similar to
#' principal components) from the independent data (temperatures) and then uses
#' these components for the regression. The contribution of each individual
#' variable to the PLS model is then evaluated with two main metrics: the
#' Variable Importance in the Projection statistic (VIP) indicates how much
#' variation in a given independent variable is correlated with variation in
#' the dependent variable. A threshold of 0.8 is often used for determining
#' importance. The standardized model coefficients of the PLS model then give
#' an indication of the direction and strength of the effect, e.g. if
#' coefficients are positive and high, high values for the respective
#' independent variable are correlated with high values of the dependent
#' variable (e.g. late occurrence of a phenological stage). This procedure was
#' inspired by the challenge of explaining variation in bloom and leaf
#' emergence dates of temperate fruit trees in Mediterranean climates. These
#' are generally understood to result from (more of less) sequential
#' fulfillment of a chilling and a forcing requirement. During the chilling
#' phase, cool temperatures are needed; during the forcing phase, trees need
#' heat. There is no easily visible change in tree buds that would indicate the
#' transition between these two phases, making it difficult to develop a good
#' model of these processes. Where long-term phenology data are available and
#' can be couple with daily temperature records, PLS regression allows
#' detection of the chilling/forcing transition. This procedure has not often
#' been applied to biological phenomena at the time of writing this, and there
#' may be constraints to how generally applicable it is. Yet is has passed the
#' test of scientific peer review a few times, and it has produced plausible
#' results in a number of settings. This package draws heavily from the pls
#' package. It also incorporates very helpful comments from Sabine Guesewell of
#' ETH Zurich (Switzerland), who pointed out some errors in the PLS procedure
#' and made suggestions for improvement.
#' 
#' @param weather_data a dataframe containing daily minimum and maximum
#' temperature data (in columns called Tmin and Tmax, respectively), and/or
#' mean daily temperature (in a column called Tmean). There also has to be a
#' column for Year and one for JDay (the Julian date, or day of the year).
#' Alternatively, the date can also be given in three columns (Years, Month and
#' Day).
#' @param bio_data a data frame that contains information on the timing of
#' phenological events by year. It should consist of two columns called Year
#' and pheno. Data in the pheno column should be in Julian date (day of the
#' year).
#' @param split_month the procedure analyzes data by phenological year, which
#' can start and end in any month during the calendar year (currently only at
#' the beginning of a month). This variable indicates the last month (e.g. 5
#' for May) that should be included in the record for a given phenological
#' year. All subsequent months are assigned to the following phenological year.
#' @param runn_mean application of a running mean function to daily mean
#' temperatures before running the PLS procedure substantially enhances the
#' clarity of outputs. runn_mean requires an odd integer value specifying how
#' many days should be included in this running mean. runn_mean=11 has usually
#' produced good results.
#' @param expl.var percentage of the variation in the dependent variable that
#' the PLS model should explain. This is used as a threshold in finding the
#' appropriate number of components in the PLS regression procedure.
#' @param ncomp.fix fixed number of components for the PLS model. Defaults to
#' NULL, so that the number is automatically determined, but is can also be set
#' by the user.
#' @param use_Tmean boolean variable indicating whether or not the column Tmean
#' from the weather_data_frame should be used as input for the PLS analysis. If
#' this is set to FALSE, Tmean is calculated as the arithmetic mean of Tmin and
#' Tmax.
#' @param return.all boolean variable indicating whether or not the full set of
#' PLS results should be output from the function. If this is set to TRUE, the
#' function output is a list with two elements: PLS_summary and PLS_output; if
#' it is set to FALSE, only the PLS_summary is returned.
#' @param crossvalidate character variable indicating what kind of validation
#' should be performed by the PLS procedure. This defaults to "none", but the
#' plsr function (of the pls package) also takes "CV" and "LOO" as inputs. See
#' the documentation for the plsr function for details.
#' @param end_at_pheno_end boolean variable indicating whether the analysis
#' should disregard temperatures after the last date included in the
#' bio_data_frame dataset. If set to TRUE, only temperatures up this date are
#' considered. Phenology data is extracted from the PLS output files. If this
#' parameter is assigned a numeric value, only data up to the Julian date
#' specified by this number are considered.
#' @return \item{object_type}{ the character string "PLS_Temp_pheno". This is
#' only needed for choosing the correct method for the plot_PLS function.}
#' \item{pheno}{ a data frame containing the phenology data used for the PLS
#' regression, with columns Year and pheno.} \item{PLS_summary}{ a data frame
#' containing all important outputs of the PLS regression. Columns are Date (in
#' MMDD format), JDay (Julian date, or day of the year), Coefficient (the PLS
#' model coefficient for each daily temperature variable), and VIP (the
#' Variable Importance in the Projection score). The columns Tmean and Tstdev
#' contain the means and standard deviations of temperature for each day of the
#' year.} \item{PLS_output}{ this is the complete output of the plsr function
#' of the pls package. See the documentation for that package for further
#' details.}
#' @author Eike Luedeling, with contributions from Sabine Guesewell
#' @references The method is described here:
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
#' Some applications:
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
#' PLS_results<-PLS_pheno(
#'   weather_data=KA_weather,
#'   split_month=6,   #last month in same year
#'   bio_data=KA_bloom)
#'   
#' PLS_results_path<-paste(getwd(),"/PLS_output",sep="")
#'   
#' plot_PLS(PLS_results,PLS_results_path)
#' 
#' @export PLS_pheno
PLS_pheno <- function (weather_data,bio_data,split_month=7,runn_mean=11,expl.var=30,
                       ncomp.fix=NULL,use_Tmean=FALSE,return.all=FALSE,crossvalidate="none",end_at_pheno_end=TRUE) {

  # Interpolate missing temperatures and calculate Tmean from Tmin and Tmax - NEW: only if necessary
    if (!use_Tmean) { 
      # Remove missing temperature data at the beginning or end of the time series
      cl.dat <- which(!is.na(weather_data[,"Tmax"])|!is.na(weather_data[,"Tmin"]))
      weather_file <- weather_data[min(cl.dat):max(cl.dat),]
      weather_file [which(weather_file$Tmin > weather_file$Tmax), c("Tmin", "Tmax")] <- NA
      Tmin_gaps <- interpolate_gaps(weather_file$Tmin)
      weather_file$Tmin <- Tmin_gaps[[1]]
      Tmax_gaps <- interpolate_gaps(weather_file$Tmax)
      weather_file$Tmax <- Tmax_gaps[[1]]
      weather_file[, "Tmean"] <- (weather_file$Tmax + weather_file$Tmin)/2 }
    
    if (use_Tmean) {
      # Remove missing temperature data at the beginning or end of the time series
      cl.dat <- which(!is.na(weather_data[,"Tmean"]))
      weather_file <- weather_data[min(cl.dat):max(cl.dat),]
      Tmean_gaps <- interpolate_gaps(weather_file$Tmean)
      weather_file$Tmean <- Tmean_gaps[[1]] }
    
    # Adjust dates  - "Season" = current year before the cutpoint, following year after the cutpoint
    weather_file[weather_file$Month<=split_month,"Season"] <- weather_file[weather_file$Month<=split_month,"Year"]
    weather_file[weather_file$Month>split_month,"Season"] <- weather_file[weather_file$Month>split_month,"Year"] + 1
    weather_file[, "Date"] <- weather_file$Month*100 + weather_file$Day
    weather_file[, "JDay"] <- strptime(paste(weather_file$Month,"/", weather_file$Day,"/", weather_file$Year, sep = ""), 
        "%m/%d/%Y")$yday + 1
    ww <- weather_file[, "Tmean"]
    rr <- weather_file[, "Tmean"]
    
    # Calculate running means of temperatures
    for (dd in 1:length(ww)) {
        if (dd < ceiling(runn_mean/2)) { rr[dd] <- mean(ww[1:(dd + floor(runn_mean/2))]) }
        if ((dd >= ceiling(runn_mean/2)) & (dd <= length(ww)-ceiling(runn_mean/2))) {
            rr[dd] <- mean(ww[(dd-floor(runn_mean/2)):(dd+floor(runn_mean/2))])  }
        if (dd > (length(ww) - ceiling(runn_mean/2))) {rr[dd] <- mean(ww[(dd-floor(runn_mean/2)):length(ww)])} }
        weather_file[, "runn"] <- rr
    
    pls_ncomp <- function(indep,dep,threshold) {  ### NEW: Variable threshold of explained variation
        dat <- data.frame(dep)
        dat$runn <- indep
        if (length(dep) > 15) {
            suppressWarnings(pls_out <- plsr(dep~runn, data=dat, ncomp=10, validation="none"))  ### NEW: no cross-validation
            ncomp <- which(cumsum(explvar(pls_out)) > threshold)[1]
            if (is.na(ncomp)) ncomp <- 10 }
        else ncomp <- 2
        return(ncomp) }
        
    # Matrix with temperatures for each day of the year (columns) and each season (rows)    
    seasons <- unique(weather_file$Season)    # One-year periods between cutpoints
    for (yy in seasons) {
    	yearweather <- weather_file[weather_file$Season==yy, ]
        weathervector <- yearweather$runn[1:365]     # runn = running mean temperature
        if (yy == seasons[1]) 
            year_res <- weathervector
        else year_res <- rbind(year_res, weathervector)
        if (nrow(yearweather) == 365) {
            labdates <- yearweather$Date
            labJdates <- yearweather$JDay  }  }
    
    
    
    colnames(year_res) <- paste("runn_", 1:365, sep="")
    year_res <- cbind(Season = seasons, year_res)
    data_set <- year_res
    full_seasons <- which(!is.na(rowMeans(data_set)))
    data_set <- data_set[full_seasons, ]
    newseasons <- data_set[, "Season"]
    suppressWarnings(bio_data <- bio_data[which(bio_data[, 
                                                         "Year"] %in% newseasons), ])
    suppressWarnings(bio_data <- bio_data[which(!is.na(as.numeric(as.character(bio_data$pheno)))), 
                                          ])
    #all_outputs[["pheno"]]<-bio_data
    suppressWarnings(bio <- as.numeric(as.character(bio_data$pheno)))
    indep <- as.matrix(data_set[which(data_set[, "Season"] %in% 
                                        bio_data$Year), ])
    indep <- indep[, 2:ncol(indep)]
    
    pheno_end<-get_last_date(as.numeric(as.character(bio_data$pheno)))
    
    if(is.numeric(end_at_pheno_end)) pheno_end<-end_at_pheno_end
    
    if(end_at_pheno_end) if (pheno_end %in% labJdates) {dayskeep<-1:which(labJdates==pheno_end)
    labJdates<-labJdates[dayskeep]
    labdates<-labdates[dayskeep]
    
    indep<-indep[,dayskeep]}
    
    
    
   # colnames(year_res) <- paste("runn_", 1:365, sep="")
  #  data_set <- cbind(Season=seasons, year_res)
  #  data_set <- data_set[which(!is.na(rowMeans(data_set))), ]  # Only seasons with temperatures from first to last day
  #  newseasons <- data_set[, "Season"]  # vector with full seasons
    
  #  # Reduce phenological data to those of full seasons; Reduce both datasets to those with phenology data
  #  bio_data <- bio_data[!is.na(bio_data$pheno), ]
  #  bio_data <- subset(bio_data,bio_data$Year %in% newseasons)
  #  bio <- as.numeric(as.character(bio_data$pheno))
  #  indep <- data_set[data_set[,"Season"] %in% bio_data$Year, ]
  #  indep <- indep[,-1]
    
  #  pheno_end<-max(as.numeric(as.character(bio_data$pheno)),na.rm=TRUE)
    
  #  if(is.numeric(end_at_pheno_end)) pheno_end<-end_at_pheno_end
  #  
  #  if(end_at_pheno_end) if (pheno_end %in% labJdates) {dayskeep<-1:which(labJdates==pheno_end)
  #  labJdates<-labJdates[dayskeep]
  #  labdates<-labdates[dayskeep]
  #  
  #  indep<-indep[,dayskeep]}
    
    if (is.null(ncomp.fix)) {ncomp <- pls_ncomp(indep=indep,dep=bio,threshold=expl.var) } else {ncomp <- ncomp.fix}   ### NEW, works with fixed ncomp
   
    ### PLS regression
    # Temperatures of each day are scaled by their standard deviations, so that all regressors have same variation
    sdindep <- apply(indep,2,sd)  # SD of temperatures on each day of the season
    sdindep[which(sdindep == 0)] <- 1    # Avoids impossible cases when scaling by standard deviation
    PLS_output <- plsr(bio~indep,ncomp,validation=crossvalidate,method="oscorespls",scale=sdindep)
    d1 <- switch(split_month,31,59,89,120,151,181,212,243,274,304,335,365)
	labJdates[labJdates>d1] <- labJdates[labJdates>d1]-365
    out <- data.frame()
    tablength<-length(coef(PLS_output))
    out[1:tablength,"Date"] <- labdates
    out[1:tablength,"JDay"] <- labJdates
    out[1:tablength,"Coef"] <- coef(PLS_output)
    out[1:tablength,"VIP"] <- VIP(PLS_output)[ncomp, ]
    out[1:tablength,"Tmean"] <- colMeans(indep)
    out[1:tablength,"Tstdev"] <- apply(indep, 2, sd, na.rm = TRUE)
    if (return.all) return(list(object_type="PLS_Temp_pheno",pheno=bio_data,PLS_summary=out,PLS_output=PLS_output)) else return(list(object_type="PLS_Temp_pheno",pheno=bio_data,PLS_summary=out))  }

