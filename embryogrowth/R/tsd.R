#' tsd estimates the parameters that best describe temperature-dependent sex determination
#' @title Estimate the parameters that best describe temperature-dependent sex determination
#' @author Marc Girondot
#' @return A list the pivotal temperature, transitional range of temperatures and their SE
#' @param males A vector with male numbers
#' @param females A vector with female numbers
#' @param N A vector with total numbers
#' @param temperatures The constant incubation temperatures used to fit sex ratio
#' @param durations The duration of incubation or TSP used to fit sex ratio
#' @param df A dataframe with at least two columns named males, females or N and temperatures, Incubation.temperature or durations column
#' @param l The limit to define TRT (see Girondot, 1999)
#' @param parameters.initial Initial values for P, S or K search as a vector, ex. c(P=29, S=-0.3)
#' @param fixed.parameters Parameters that will not be changed
#' @param SE Standard errors for parameters
#' @param males.freq If TRUE data are shown as males frequency
#' @param equation Could be "logistic", "Hill", "Richards", "Hulin", "Double-Richards" or "GSD"
#' @param replicates Number of replicates to estimate confidence intervals
#' @param range.CI The range of confidence interval for estimation, default=0.95
#' @param limit.low.TRT.minimum Minimum lower limit for TRT
#' @param limit.high.TRT.maximum Maximum higher limit for TRT
#' @param temperatures.plot Sequences of temperatures that will be used for plotting. If NULL, does not estimate them
#' @param durations.plot Sequences of durations that will be used for plotting. If NULL, does not estimate them
#' @param print Do the results must be printed at screen? TRUE (default) or FALSE
#' @description Estimate the parameters that best describe temperature-dependent sex determination
#' @references Girondot, M. 1999. Statistical description of temperature-dependent sex determination using maximum likelihood. Evolutionary Ecology Research, 1, 479-486.
#' @references Godfrey, M.H., Delmas, V., Girondot, M., 2003. Assessment of patterns of temperature-dependent sex determination using maximum likelihood model selection. Ecoscience 10, 265-272.
#' @references Hulin, V., Delmas, V., Girondot, M., Godfrey, M.H., Guillon, J.-M., 2009. Temperature-dependent sex determination and global change: are some species at greater risk? Oecologia 160, 493-506.
#' @family Functions for temperature-dependent sex determination
#' @examples
#' \dontrun{
#' CC_AtlanticSW <- subset(STSRE_TSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & Sexed!=0)
#' tsdL <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="logistic"))
#' tsdH <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Hill"))
#' tsdR <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Richards"))
#' tsdDR <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Double-Richards"))
#' gsd <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="GSD"))
#' compare_AIC(Logistic_Model=tsdL, Hill_model=tsdH, Richards_model=tsdR, 
#'                DoubleRichards_model=tsdDR, GSD_model=gsd)
#' ##############
#' eo <- subset(STSRE_TSD, Species=="Emys orbicularis", c("Males", "Females", 
#'                                        "Incubation.temperature"))
#'                                        
#' eo_Hill <- with(eo, tsd(males=Males, females=Females, 
#'                                        temperatures=Incubation.temperature,
#'                                        equation="Hill"))
#' eo_Hill <- tsd(df=eo, equation="Hill")
#' eo_logistic <- tsd(eo)
#' eo_Richards <- with(eo, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature, 
#'                                  equation="Richards"))
#' ### The Hulin model is a modification of Richards (See Hulin et al. 2009)
#' ### limit.low.TRT and limit.high.TRT must be setup for Hulin equation
#' par <- eo_Richards$par
#' names(par)[which(names(par)=="K")] <- "K2"
#' par <- c(par, K1=0)
#' eo_Hulin <- with(eo, tsd(males=Males, females=Females, 
#'                                  parameters.initial=par, 
#'                                  temperatures=Incubation.temperature, 
#'                                  equation="Hulin", 
#'                                  limit.low.TRT.minimum=25, 
#'                                  limit.high.TRT.maximum=35))
#' ### The Double-Richards model is a Richards model with K1 and K2 using the two values
#' ### below and above P
#' par <- eo_Richards$par
#' names(par)[which(names(par)=="K")] <- "K2"
#' par <- c(par, K1=as.numeric(par["K2"])-0.1)
#' par["K1"] <- par["K1"]-0.1
#' eo_Double_Richards <- with(eo, tsd(males=Males, females=Females,
#'                                  parameters.initial=par,
#'                                  temperatures=Incubation.temperature,
#'                                  equation="Double-Richards"))
#' compare_AIC(Logistic=eo_logistic, Hill=eo_Hill, Richards=eo_Richards, 
#'              Hulin=eo_Hulin, Double_Richards=eo_Double_Richards)
#' ### Note the asymmetry of the Double-Richards model
#' predict(eo_Double_Richards, 
#'        temperatures=c(eo_Double_Richards$par["P"]-0.2, eo_Double_Richards$par["P"]+0.2))
#' predict(eo_Double_Richards)
#' ### It can be used also for incubation duration
#' CC_AtlanticSW <- subset(STSRE_TSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & Sexed!=0)
#' tsdL_IP <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  durations=IP.mean, 
#'                                  equation="logistic"))
#' plot(tsdL_IP, xlab="Incubation durations in days")
#' }
#' @export


tsd <- function(df=NULL, males=NULL, females=NULL, N=NULL, 
                temperatures=NULL, 
                durations=NULL,
                l=0.05, parameters.initial=c(P=NA, S=-0.5, K=0, K1=1, K2=0), 
                males.freq=TRUE, 
                fixed.parameters=NULL, SE=NULL,
                equation="logistic", replicates=1000, range.CI=0.95, 
                limit.low.TRT.minimum=5, limit.high.TRT.maximum=90, print=TRUE, 
                temperatures.plot=seq(from=20, to=40, by=0.1), 
                durations.plot=seq(from=15, to=100, by=0.1)) {
  
  # df=NULL; males=NULL; females=NULL; N=NULL; temperatures=NULL; durations=NULL; l=0.05; parameters.initial=c(P=NA, S=-0.5, K=0, K1=1, K2=0); males.freq=TRUE; fixed.parameters=NULL; SE=NULL; equation="logistic"; replicates=1000; range.CI=0.95; limit.low.TRT.minimum=5; limit.high.TRT.maximum=90; print=TRUE; temperatures.plot=seq(from=20, to=40, by=0.1); durations.plot=seq(from=15, to=100, by=0.1)
  # females=females; males=males; durations=durations; parameters.initial=c(P=54.0445435044487, S=0.0567693543096269); equation="Hill"
  # eo <- subset(STSRE_TSD, Species=="Emys orbicularis", c("Males", "Females","Incubation.temperature"))

  # males=eo$Males; females=eo$Females; N <- Males+Females; temperatures=eo$Incubation.temperature; equation<- "Double-Richards" 
  
  # males <- c(10, 14, 7, 4, 3, 0, 0) ; females <- c(0, 1, 2, 4, 15, 10, 13); temperatures <- c(25, 26, 27, 28, 29, 30, 31)
  
  # males <- c(15L, 23L, 16L, 0L); females <- c(0L, 1L, 3L, 19L)
  # temperatures <- c(70.4, 57.7, 53.9, 51.7)
  # parameters.initial <- c(P=53, S=0.05)
  # col.TRT="gray"; col.TRT.CI=rgb(0.8, 0.8, 0.8, 0.5); col.PT.CI=rgb(0.8, 0.8, 0.8, 0.5)
  # N <- NULL; males <- CC_AtlanticSW$Males;females <- CC_AtlanticSW$Females; temperatures <- CC_AtlanticSW$Incubation.temperature-CC_AtlanticSW$Correction.factor
  # equation <- "Richards"
  # equation <- "Hill"
  # equation <- "Double-Richards"
  # equation <- "logistic"
  
  range.CI.qnorm <- qnorm(1-((1-range.CI)/2))
  
  if (!is.null(df)) {
    if (class(df)!="data.frame" & class(df)!="matrix") {
      stop("df parameter must be a data.frame or a matrix")
    }
    
    namesdf <- tolower(colnames(df))
    males <- NULL
    females <- NULL
    N <- NULL
    
    if (any(namesdf=="males")) males <- df[, which(namesdf=="males")]
    if (any(namesdf=="male")) males <- df[, which(namesdf=="male")]
    if (any(namesdf=="females")) females <- df[, which(namesdf=="females")]
    if (any(namesdf=="female")) females <- df[, which(namesdf=="female")]
    if (any(namesdf=="n")) N <- df[, which(namesdf=="n")]
    if (any(namesdf=="sexed")) N <- df[, which(namesdf=="sexed")]
    if (any(namesdf=="temperatures")) temperatures <- df[, which(namesdf=="temperatures")]
    if (any(namesdf=="temperature")) temperatures <- df[, which(namesdf=="temperature")]
    if (any(namesdf=="incubation.temperature")) temperatures <- df[, which(namesdf=="incubation.temperature")]
    if (any(namesdf=="durations")) durations <- df[, which(namesdf=="durations")]
    if (any(namesdf=="duration")) durations <- df[, which(namesdf=="duration")]
    if (any(namesdf=="IP.mean")) durations <- df[, which(namesdf=="IP.mean")]
  }
  
  if (is.null(durations) & is.null(temperatures)) {
    stop("Or temperatures or durations must be provided")
  }
  
  if (!is.null(durations) & !is.null(temperatures)) {
    stop("Temperatures and durations cannot be both provided")
  }
  
  
  # si je retire des lignes, c'est decale
  #  if (!is.null(males)) males <- na.omit(males)
  #  if (!is.null(females)) females <- na.omit(females)
  #  if (!is.null(N)) N <- na.omit(N)
  if (is.null(temperatures)) {
    IP <- TRUE
    temperatures <- durations
    temperatures.plot <- durations.plot
  } else {
    IP <- FALSE
  }
  
  cpt1 <- ifelse(is.null(males), 0, 1)+ifelse(is.null(females), 0, 1)+ifelse(is.null(N), 0, 1)
  if (cpt1<2) {
    stop("At least two informations from males, females or N must be provided")
  }
  
  if (is.null(males)) males <- N-females
  if (is.null(females)) females <- N-males
  if (is.null(N)) N <- males+females
  
  if (length(temperatures)!=length(N)) {
    stop("A temperature or duration must be provided for each experiment")
  }
  
  if (equation!="GSD" & equation!="Hulin" & equation!="logistic" & equation!="Hill" & equation!="Richards" & equation!="Double-Richards") {
    stop("Equations supported are GSD, logistic, Hill, Hulin, Double_Richards and Richards")
  }
  
  
  if (equation=="GSD") {
    result <- list(par = NULL, SE=NULL, hessian=NULL, 
                   TRT=NULL,
                   SE_TRT=NULL,
                   message=NULL,
                   convergence=0,
                   counts=NULL,
                   value=-sum(dbinom(x=males, size=males+females, prob=0.5, log=TRUE)),
                   AIC=-2*sum(dbinom(x=males, size=males+females, prob=0.5, log=TRUE)))
#    limit.low.TRT <- min(temperatures)
#    limit.high.TRT <- max(temperatures)	
  } else {
    par <- parameters.initial
    if (!is.null(par)) {
    if (is.na(par["P"])) {
      par["P"] <- temperatures[which.min(abs((males/(males+females))-0.5))]
      if (IP) par["S"] <- abs(par["S"])
    }
      
    }
    if (equation!="Richards") par <- par[which(names(par)!="K")]
    if (equation!="Hulin" & equation!="Double-Richards") par <- par[which(names(par)!="K1")]
    if (equation!="Hulin" & equation!="Double-Richards") par <- par[which(names(par)!="K2")]
    
    repeat {
     # result  <- optim(par, embryogrowth:::.tsd_fit, fixed.parameters=fixed.parameters, males=males, N=N, temperatures=temperatures, equation=equation, method="BFGS", hessian=TRUE, control = list(maxit=1000))
     result  <- optim(par, getFromNamespace(".tsd_fit", ns="embryogrowth"), fixed.parameters=fixed.parameters, males=males, N=N, temperatures=temperatures, equation=equation, method="BFGS", hessian=TRUE, control = list(maxit=1000))

      if (result$convergence==0) break
      par<-result$par
      if (print) print("Convergence is not acheived. Optimization continues !")
    }

    par <- c(result$par, fixed.parameters)
    
    if (!is.null(SE)) {
      res <- SE
    } else {
      mathessian <- result$hessian
      
      inversemathessian <- try(solve(mathessian), silent=TRUE)
      
      if (inherits(inversemathessian, "try-error")) {
        res <- rep(NA, length(par))
        names(res) <- names(par)
      } else {
        res <- abs(diag(inversemathessian))
      }
    }
    
    result$SE <- res
    result$AIC <- 2*result$value+2*length(par)
  }
  
  result$range.CI <- range.CI
  result$l <- l
  
  if (equation!="GSD") {
    l.l.TRT.c <- NULL
    l.h.TRT.c <- NULL
    TRT.c <- NULL
    #   plist <- NULL
    
    rep <- replicates-1
    df_par <- data.frame(P=unname(c(par["P"], rnorm(rep, par["P"], res["P"]))), 
                         S=unname(c(par["S"], rnorm(rep, par["S"], res["S"]))),
                         K=unname(ifelse(is.na(par["K"]), rep(NA, replicates), c(par["K"], rnorm(rep, par["K"], ifelse(is.na(res["K"]), 0, res["K"]))))),
                         K1=unname(ifelse(is.na(par["K1"]), rep(NA, replicates), c(par["K1"], rnorm(rep, par["K1"], ifelse(is.na(res["K1"]), 0, res["K1"]))))),
                         K2=unname(ifelse(is.na(par["K2"]), rep(NA, replicates), c(par["K2"], rnorm(rep, par["K2"], ifelse(is.na(res["K2"]), 0, res["K2"]))))))

    out_tsd <- apply(df_par, 1, function(par) {
      # marche en temperatures mais par en IP - corrige 9/9/2014
      limit.low.TRT <- limit.low.TRT.minimum
      limit.high.TRT <- limit.high.TRT.maximum
      for (i in 1:4) {
        temperatures.se <- seq(from=limit.low.TRT, to=limit.high.TRT, length=20)
        p <- getFromNamespace(".modelTSD", ns="embryogrowth")(par, temperatures.se, equation)
        
        if (sign(par["S"])<0) {
          limit.low.TRT <- temperatures.se[tail(which(p>(1-l)), n=1)]
          limit.high.TRT <- temperatures.se[which(p < l)[1]]
        } else {
          limit.low.TRT <- temperatures.se[tail(which(p<l), n=1)]
          limit.high.TRT <- temperatures.se[which(p >(1-l))[1]]
        }
        limit.low.TRT <- ifelse(identical(limit.low.TRT, numeric(0)), NA, limit.low.TRT)
        limit.high.TRT <- ifelse(identical(limit.high.TRT, numeric(0)), NA, limit.high.TRT)
        if (is.na(limit.low.TRT) | is.na(limit.high.TRT)) break
      }
      return(data.frame(l.l.TRT.c=limit.low.TRT, l.h.TRT.c=limit.high.TRT, TRT.c=limit.high.TRT-limit.low.TRT))
    }
    )
    
    out_tsd2 <- t(sapply(out_tsd, c))
    
    if (!is.null(temperatures.plot)) {
      out_tsd_plot <- apply(df_par, 1, function(par) {
        p <- getFromNamespace(".modelTSD", ns="embryogrowth")(par, temperatures.plot, equation)
        return(list(sr=p))
      }
      )
      
      out_tsd2_plot <- sapply(out_tsd_plot, c, simplify = TRUE)
      # dans un df chaque colonne est une temperature
      out_tsd3_plot <- t(sapply(out_tsd2_plot, c, simplify = TRUE))
#      dim(out_tsd3_plot)
      outquant <- apply(out_tsd3_plot, 2, quantile, probs = c(0.025, 0.975))
      outquant <- rbind(outquant, mean=out_tsd3_plot[1,], temperatures=temperatures.plot)
      result$CI <- outquant
      
    }
    
    
    result$TRT <- as.numeric(out_tsd2[1, "TRT.c"])
    result$TRT_limits <- as.numeric(c(out_tsd2[1, "l.l.TRT.c"], out_tsd2[1, "l.h.TRT.c"]))
    result$SE_TRT <- sd(as.numeric(out_tsd2[, "TRT.c"]), na.rm = TRUE)
    result$SE_TRT_limits <- c(sd(as.numeric(out_tsd2[, "l.l.TRT.c"]), na.rm = TRUE), sd(as.numeric(out_tsd2[, "l.h.TRT.c"]), na.rm = TRUE))
    limit.low.TRT <- result$TRT_limits[1]-range.CI.qnorm*result$SE_TRT_limits[1]
    limit.high.TRT <- result$TRT_limits[2]+range.CI.qnorm*result$SE_TRT_limits[2]
  }
  
  saturated <- 0
  for (i in 1:length(N)) saturated <- saturated- dbinom(males[i], N[i], males[i]/N[i], log=TRUE)
  LRT <- 2*(result$value-saturated)
  result$GOF <- 1-pchisq(LRT, length(N)-length(par))
  
  result$males <- males
  result$females <- females
  result$N <- N
  result$temperatures <- temperatures
  result$males.freq <- males.freq
  result$equation <- equation
  result$l <- l
  result$fixed.parameters <- fixed.parameters
  
  if (print) print(paste("The goodness of fit test is", sprintf("%.5f",result$GOF)))
  
  if (equation!="GSD" & print) {
    print(paste("The pivotal is", sprintf("%.3f",par["P"]), "SE", sprintf("%.3f",res["P"])))
    print(paste("The transitional range of temperatures is", sprintf("%.3f",result$TRT), "SE", sprintf("%.3f",result$SE_TRT)))
    print(paste("The lower limit of transitional range of temperatures is", sprintf("%.3f",result$TRT_limits[1]), "SE", sprintf("%.3f",result$SE_TRT_limits[1])))
    print(paste("The higher limit of transitional range of temperatures is", sprintf("%.3f",result$TRT_limits[2]), "SE", sprintf("%.3f",result$SE_TRT_limits[2])))
  }
  
  class(result) <- "tsd"
  return(invisible(result))
}
