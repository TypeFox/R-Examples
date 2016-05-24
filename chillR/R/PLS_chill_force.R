#' Partial Least Squares analysis of phenology vs. accumulated daily chill and
#' heat
#' 
#' This function conducts a Partial Least Squares (PLS) regression analysis
#' relating an annual biological phenomenon, e.g. fruit tree flowering or leaf
#' emergence, to mean daily rates of chill (with three models) and heat
#' accumulation of the preceding 12 months. It produces figures that illustrate
#' statistical correlations between temperature variation during certain phases
#' and the timing of phenological events.
#' 
#' PLS regression is useful for exploring the relationship between daily chill
#' and heat accumulation rates and biological phenomena that only occur once
#' per year. The statistical challenge is that a normally quite small number of
#' observations must be related to variation in a much larger number (730) of
#' daily chill and heat values, which are also highly autocorrelated. Most
#' regression approaches are not suitable for this, but PLS regression offers a
#' potential solution. The method is frequently used in chemometrics and
#' hyperspectral remote sensing, where similar statistical challenges are
#' encountered. The basic mechanism is that PLS first constructs latent factors
#' (similar to principal components) from the independent data (daily chill and
#' heat accumulation) and then uses these components for the regression. The
#' contribution of each individual variable to the PLS model is then evaluated
#' with two main metrics: the Variable Importance in the Projection statistic
#' (VIP) indicates how much variation in a given independent variable is
#' correlated with variation in the dependent variable. A threshold of 0.8 is
#' often used for determining importance. The standardized model coefficients
#' of the PLS model then give an indication of the direction and strength of
#' the effect, e.g. if coefficients are positive and high, high values for the
#' respective independent variable are correlated with high values of the
#' dependent variable (e.g. late occurrence of a phenological stage). This
#' procedure was inspired by the challenge of explaining variation in bloom and
#' leaf emergence dates of temperate fruit trees in Mediterranean climates.
#' These are generally understood to result from (more of less) sequential
#' fulfillment of a chilling and a forcing requirement. During the chilling
#' phase, cool temperatures are needed; during the forcing phase, trees need
#' heat. There is no easily visible change in tree buds that would indicate the
#' transition between these two phases, making it difficult to develop a good
#' model of these processes. Where long-term phenology data are available and
#' can be coupled with daily chill and heat records (derived from daily
#' temperature data), PLS regression allows detection of the chilling/forcing
#' transition. This procedure has not often been applied to biological
#' phenomena at the time of writing this, and there may be constraints to how
#' generally applicable it is. Yet is has passed the test of scientific peer
#' review a few times, and it has produced plausible results in a number of
#' settings. This package draws heavily from the pls package.
#' 
#' Per default, chill metrics used are the ones given in the references below.
#' Chilling Hours are all hours with temperatures between 0 and 7.2 degrees C.
#' Units of the Utah Model are calculated as suggested by Richardson et al.
#' (1974) (different weights for different temperature ranges, and negation of
#' chilling by warm temperatures). Chill Portions are calculated according to
#' Fishman et al. (1987a,b). More honestly, they are calculated according to an
#' Excel sheet produced by Amnon Erez and colleagues, which converts the
#' complex equations in the Fishman papers into relatively simple Excel
#' functions. These were translated into R. References to papers that include
#' the full functions are given below. Growing Degree Hours are calculated
#' according to Anderson et al. (1986), using the default values they suggest.
#' 
#' It is possible, however, for the user to specify other metrics to be
#' evaluated. These should be indicated by the chill_models and heat_models
#' parameters, which should contain the names of the respective columns of the
#' daily_chill_obj$daily_chill data frame.
#' 
#' @param daily_chill_obj a daily chill object. This should be generated with
#' the daily_chill function.
#' @param bio_data_frame a data frame that contains information on the timing
#' of phenological events by year. It should consist of two columns called Year
#' and pheno. Data in the pheno column should be in Julian date (day of the
#' year).
#' @param split_month the procedure analyzes data by phenological year, which
#' can start and end in any month during the calendar year (currently only at
#' the beginning of a month). This variable indicates the last month (e.g. 5
#' for May) that should be included in the record for a given phenological
#' year. All subsequent months are assigned to the following phenological year.
#' @param expl.var percentage of the variation in the dependent variable that
#' the PLS model should explain. This is used as a threshold in finding the
#' appropriate number of components in the PLS regression procedure.
#' @param ncomp.fix fixed number of components for the PLS model. Defaults to
#' NULL, so that the number is automatically determined, but it can also be set
#' by the user.
#' @param return.all boolean variable indicating whether or not the full set of
#' PLS results should be returned by the function. If this is set to TRUE, the
#' function output is a list with two elements (besides the object_type
#' string): PLS_summary and PLS_output; if it is set to FALSE, only the
#' PLS_summary is returned.
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
#' @param chill_models Character vector containing names of chill models that
#' should be considered in the PLS regression. These names should correspond to
#' column names of daily_chill. This defaults to c("Chilling_Hours",
#' "Utah_Chill_Units", "Chill_Portions").
#' @param heat_models Character vector containing names of heat models that
#' should be considered in the PLS regression. These names should correspond to
#' column names of daily_chill. This defaults to c("GDH").
#' @return \item{object_type}{ the character string "PLS_chillforce_pheno".
#' This is only needed for choosing the correct method for the plot_PLS
#' function.} \item{pheno}{ a data frame containing the phenology data used for
#' the PLS regression, with columns Year and pheno.}
#' \item{<chill_model>$<heat_model>}{ for each combination of elements from
#' chill_models and heat_models, a list element is generated, which contains a
#' list with elements PLS_summary and (if(return.all=TRUE) PLS_output. These
#' contain the results of the PLS analysis that used the respective chill and
#' heat metrics as independent variables.}
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Chilling Hours or Utah
#' Models, especially in warm climates! The Dynamic Model (Chill Portions),
#' though far from perfect, seems much more reliable.
#' @author Eike Luedeling, with contributions from Sabine Guesewell
#' @references Model references, for the default option:
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
#' @keywords phenology analysis chill and heat calculation
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
#' #Plots look much better with weather<-fix_weather(KA_weather)
#' #but that takes to long to run for passing CRAN checks
#' 
#' dc<-daily_chill(stack_hourly_temps(weather,50.4), 11)
#' plscf<-PLS_chill_force(daily_chill_obj=dc, bio_data_frame=KA_bloom, split_month=6)
#' 
#' #PLS_results_path<-paste(getwd(),"/PLS_output",sep="")
#' #plot_PLS(plscf,PLS_results_path)
#' #plot_PLS(plscf,PLS_results_path,add_chill=c(307,19),add_heat=c(54,109))
#' 
#' 
#' 
#' 
#' @export PLS_chill_force
PLS_chill_force<-function (daily_chill_obj, bio_data_frame, split_month, expl.var = 30, 
                           ncomp.fix = NULL, return.all = FALSE, crossvalidate = "none",end_at_pheno_end=TRUE,
                           chill_models=c("Chilling_Hours","Utah_Chill_Units","Chill_Portions"),
                           heat_models=c("GDH")) 
{
  if (daily_chill_obj[1] == "daily_chill") {
    dc <- daily_chill_obj$daily_chill
    weather_file <- daily_chill_obj$daily_chill
    bio_data <- bio_data_frame
    weather_file[which(weather_file$Month <= split_month), 
                 "Season"] <- weather_file[which(weather_file$Month <= 
                                                   split_month), "Year"]
    weather_file[which(weather_file$Month > split_month), 
                 "Season"] <- weather_file[which(weather_file$Month > 
                                                   split_month), "Year"] + 1
    weather_file[, "Date"] <- weather_file$Month * 100 + 
      weather_file$Day
    weather_file[, "JDay"] <- strptime(paste(weather_file$Month, 
                                             "/", weather_file$Day, "/", weather_file$Year, sep = ""), 
                                       "%m/%d/%Y")$yday + 1
    pls_ncomp_old <- function(indep, dep, threshold = 30) {
      dat <- data.frame(dep)
      dat$runn <- indep
      if (length(dep) > 15) {
        suppressWarnings(pls_out <- plsr(dep ~ runn, 
                                         data = dat, ncomp = 10, validation = "CV"))
        suppressWarnings(pls_cv <- crossval(pls_out, 
                                            segments = 10))
        res <- data.frame(ncomp = c(1:10), explvar = explvar(pls_out), 
                          cumul = NA)
        res$cumul[1] <- res$explvar[1]
        for (i in 2:nrow(res)) {
          res$cumul[i] <- res$cumul[i - 1] + res$explvar[i]
        }
        ncomp <- which(res$cumul > threshold)[1]
        if (is.na(ncomp)) 
          ncomp <- 10
      }
      else ncomp <- 2
      return(ncomp)
    }
    pls_ncomp <- function(indep, dep, threshold) {
      dat <- data.frame(dep)
      dat$runn <- indep
      if (length(dep) > 15) {
        suppressWarnings(pls_out <- plsr(dep ~ runn, 
                                         data = dat, ncomp = 10, validation = "none"))
        ncomp <- which(cumsum(explvar(pls_out)) > threshold)[1]
        if (is.na(ncomp)) 
          ncomp <- 10
      }
      else ncomp <- 2
      return(ncomp)
    }
    seasons <- unique(weather_file$Season)
    #chill_models <- c("Chilling_Hours", "Utah_Chill_Units", 
    #                  "Chill_Portions")
    #heat_models <- c("GDH")
    all_outputs <- list(object_type = "PLS_chillforce_pheno")
    for (CM in chill_models) for (HM in heat_models) {
      for (yy in seasons) {
        yearweather <- weather_file[which(weather_file$Season == 
                                            yy), ]
        weathervector <- c(yearweather[1:365, CM], yearweather[1:365, 
                                                               HM])
        if (yy == seasons[1]) 
          year_res <- weathervector else year_res <- rbind(year_res, weathervector)
        if (nrow(yearweather) == 365) {
          labdates <- yearweather$Date
          labJdates <- yearweather$JDay
        }
      }
      colnames(year_res) <- c(paste("Chill_", 1:365, sep = ""), 
                              paste("Heat_", 1:365, sep = ""))
      year_res <- cbind(Season = seasons, year_res)
      data_set <- year_res
      full_seasons <- which(!is.na(rowMeans(data_set)))
      data_set <- data_set[full_seasons, ]
      newseasons <- data_set[, "Season"]
      suppressWarnings(bio_data <- bio_data[which(bio_data[, 
                                                           "Year"] %in% newseasons), ])
      suppressWarnings(bio_data <- bio_data[which(!is.na(as.numeric(as.character(bio_data$pheno)))), 
                                            ])
      all_outputs[["pheno"]]<-bio_data
      suppressWarnings(bio <- as.numeric(as.character(bio_data$pheno)))
      indep <- as.matrix(data_set[which(data_set[, "Season"] %in% 
                                          bio_data$Year), ])
      indep <- indep[, 2:ncol(indep)]
      
      #pheno_end<-max(as.numeric(as.character(bio_data$pheno)))
      
      pheno_end<-get_last_date(as.numeric(as.character(bio_data$pheno)))
      
      if(is.numeric(end_at_pheno_end)) pheno_end<-end_at_pheno_end
      
      if(end_at_pheno_end) if (pheno_end %in% labJdates) {dayskeep<-1:which(labJdates==pheno_end)
      labJdates<-labJdates[dayskeep]
      labdates<-labdates[dayskeep]
      
      indep<-indep[,c(dayskeep,365+dayskeep)]}
      
      if (is.null(ncomp.fix)) {
        ncomp <- pls_ncomp(indep = indep, dep = bio, 
                           threshold = expl.var)
      }      else {        ncomp <- ncomp.fix      }
      sdindep <- apply(indep, 2, sd)
      sdindep[which(sdindep == 0)] <- 1
      PLS_output <- plsr(bio ~ indep, ncomp, validation = "none", 
                         method = "oscorespls", scale = sdindep)
      d1 <- switch(split_month, 31, 59, 89, 120, 151, 181, 
                   212, 243, 274, 304, 335, 365)
      labJdates[labJdates > d1] <- labJdates[labJdates > 
                                               d1] - 365
      out <- data.frame()
      tablength<-length(coef(PLS_output))
      out[1:tablength, "Date"] <- labdates
      out[1:tablength, "Type"] <- c(rep("Chill", tablength/2), rep("Heat", 
                                                                   tablength/2))
      out[1:tablength, "JDay"] <- labJdates
      out[1:tablength, "Coef"] <- coef(PLS_output)
      if(ncomp>1) out[1:tablength, "VIP"] <- VIP(PLS_output)[ncomp, ] else out[1:tablength, "VIP"] <- VIP(PLS_output)[ncomp]
      out[1:tablength, "MetricMean"] <- colMeans(indep)
      out[1:tablength, "MetricStdev"] <- apply(indep, 2, sd, na.rm = TRUE)
      if (return.all) 
        all_outputs[[CM]][[HM]] <- list(PLS_summary = out, 
            PLS_output = PLS_output)  else all_outputs[[CM]][[HM]] <- list(PLS_summary = out)
    }
    return(all_outputs)
  }
  else {
    "Error: not a daily chill object; use function daily_chill to make one from hourly temperature data"
  }
}

