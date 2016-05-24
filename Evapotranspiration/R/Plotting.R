ETPlot <- function(results, type = "Aggregation", OBS, OBSplot = FALSE,  Sdate = time(results$ET.Daily)[1], Edate = time(results$ET.Daily)[length(results$ET.Daily)]) { # plot estimations and observations
  
  # Aggregation plot (default)
  if (type == "Aggregation") {
    if (is.null(Sdate)) {
      Sdate = time(results$ET.Daily)[1]
    }
    if (is.null(Edate)) {
      Edate = time(results$ET.Daily)[length(results$ET.Daily)]
    }
    # Plots - Aggregation
    par(ask=FALSE)
    plot.new()
    x.Date <- as.Date(time(results$ET.Daily[which(time(results$ET.Daily) == as.Date(Sdate)):which(time(results$ET.Daily) == as.Date(Edate))]))
    if (!is.null(results$ET.Daily)) {
      plot(results$ET.Daily[which(as.Date(time(results$ET.Daily)) == x.Date[1]) : which(as.Date(time(results$ET.Daily)) == x.Date[length(x.Date)])], main = paste("Daily", results$ET_formulation, results$ET_type), xlab = "Year", ylab = list(c(results$ET_type, "mm/day")))
      if (OBSplot == TRUE) {
        if (!is.null(OBS)) {
          if (!is.null(OBS$E_obs.Daily)) {
            lines(OBS$E_obs.Daily[which(as.Date(time(results$ET.Daily)) == x.Date[1]) : which(as.Date(time(results$ET.Daily)) == x.Date[length(x.Date)])], main = "Observed evaporation", type = "o", pch = ".", col = "RED", xlab = "Year", ylab = "Observed evaporation mm/day")
            legend("topright", inset = .05, c(paste(results$ET_formulation,results$ET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
          }
        } else {
          warning("No observed data available for plotting, only the estimated values are plotted")
        }
      }
    } else {
      warning("No daily results obtained from ", results$ET_formulation)
    }
    par(ask=TRUE)
    
    plot(results$ET.Monthly[which(as.yearmon(time(results$ET.Monthly)) == as.yearmon(x.Date[1])):which(as.yearmon(time(results$ET.Monthly)) == as.yearmon(x.Date[length(x.Date)]))], ylim = c(0,max(results$ET.Monthly)*1.5), main = paste("Monthly", results$ET_formulation, results$ET_type), type = "o", pch = ".", xlab = "Year", ylab = list(c(results$ET_type, "mm/month")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.Monthly[which(as.yearmon(time(results$ET.Monthly)) == as.yearmon(x.Date[1])):which(as.yearmon(time(results$ET.Monthly)) == as.yearmon(x.Date[length(x.Date)]))], main = "Observed evaporation", type = "o", pch = ".", ylim = c(0,400), col = "RED", xlab = "Year", ylab = "Observed evaporation mm/month")
        legend("topright", inset = .05, c(paste(results$ET_formulation,results$ET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    plot(results$ET.Annual[which(floor(as.numeric(as.yearmon(time(results$ET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[1])))):which(floor(as.numeric(as.yearmon(time(results$ET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[length(x.Date)]))))],  main = paste("Annual", results$ET_formulation, results$ET_type), type = "o", xlab = "Year", ylim = c(0,3000), ylab = list(c(results$ET_type, "mm/year")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.Annual[which(floor(as.numeric(as.yearmon(time(results$ET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[1])))):which(floor(as.numeric(as.yearmon(time(results$ET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[length(x.Date)]))))], main = "Observed evaporation", type = "o", col = "RED", ylim = c(0,3000), xlab = "Year", ylab = "Observed evaporation mm/year")
        legend("topright", inset = .05, c(paste(results$ET_formulation,results$ET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    
    par(ask=FALSE)
    paste("Completed plotting aggregations for results calculated by", results$ET_formulation, "formulation")
  }
  
  # Aggregation plot (default)
  if (type == "Average") {
    Sdate = time(results$ET.Daily)[1]
    Edate = time(results$ET.Daily)[length(results$ET.Daily)]
    
    par(ask=FALSE)
    plot.new()
    # Plots - Average
    plot(results$ET.MonthlyAve~unique(1:12), main = paste("Monthly average", results$ET_formulation, results$ET_type), type = "o", ylim = c(0, max(results$ET.MonthlyAve)*1.5), xlab = "Month", ylab = list(c(paste("Monthly average", results$ET_type), "mm/day")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.MonthlyAve~unique(1:12), main = "Observed evaporation", type = "o", ylim = c(0,10), col = "RED", xlab = "Month", ylab = "Monthly average observed evaporation mm")
        legend("topright", inset = .05, c(paste(results$ET_formulation,results$ET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    
    par(ask=TRUE)
    
    plot(results$ET.AnnualAve~unique(as.POSIXlt(data$Date.daily)$year + 1900), main = paste("Annual average", results$ET_formulation, results$ET_type), type = "o", ylim = c(0, max(results$ET.AnnualAve)*1.5), xlab = "Year", ylab = list(c(paste("Annual average", results$ET_type), "mm/day")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.AnnualAve~unique(as.POSIXlt(OBS$Date.OBS)$year + 1900), main = "Observed evaporation", type = "o", ylim = c(0,10), col = "RED", xlab = "Year", ylab = "Annual average observed evaporation mm")
        legend("topright", inset = .05, c(paste(results$ET_formulation,results$ET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    
    par(ask=FALSE)
    paste("Completed plotting averages for results calculated by", results$ET_formulation, "formulation")
  }
  par(ask=FALSE)
}

#-------------------------------------------------------------------------------------


ETComparison <- function(results1, results2, results3 = NULL, results4 = NULL, results5 = NULL, results6 = NULL, results7 = NULL, 
                         labs, Sdate=NULL, Edate=NULL, type = "Monthly", ylim=rep(NA,2)) { # plot up to 7 estimations
  # check number of sets of results to compare and compile list of legend text
  
  if (exists("results1")==F | exists("results2")==F) {
    stop("Please provide at least results1 and results2 for producing comparison plot")
  }
  if (is.null(results1) | is.null(results2)) {
    stop("Please provide at least results1 and results2 for producing comparison plot")
  }
  
  Ncomp <- 2
  legtext <- vector(mode = "character", length = Ncomp)
  legtext[1] <- paste(results1$ET_formulation, results1$ET_type, labs[1])
  legtext[2] <- paste(results2$ET_formulation, results2$ET_type, labs[2])
  
  if (!is.null(results3)) {
    Ncomp <- 3
    length(legtext) <- 3
    legtext[3] <- paste(results3$ET_formulation, results3$ET_type, labs[3])
  }
  if (!is.null(results4)) {
    Ncomp <- 4
    length(legtext) <- 4
    legtext[4] <- paste(results4$ET_formulation, results4$ET_type, labs[4])
  }
  if (!is.null(results5)) {
    Ncomp <- 5
    length(legtext) <- 5
    legtext[5] <- paste(results5$ET_formulation, results5$ET_type, labs[5])
  }
  if (!is.null(results6)) {
    Ncomp <- 6
    length(legtext) <- 6
    legtext[6] <- paste(results6$ET_formulation, results6$ET_type, labs[6])
  }
  if (!is.null(results7)) {
    Ncomp <- 7
    length(legtext) <- 7
    legtext[7] <- paste(results7$ET_formulation, results7$ET_type, labs[7])
  }
  
  if (length(labs)!= Ncomp) {
    stop("Please ensure that the number of labels is the same as the number of sets of results to compare")
  }
  
  
  # Time series
  if (type == "Daily") {
    # Daily
    par(ask=FALSE)
    if (is.null(results1$ET.Daily) | is.null(results2$ET.Daily)) {
      stop("Unable to compare results because estimations of daily ET are not available in both results1 and results2")
    }
    
    if (is.null(Sdate)) {
      Sdate = time(results1$ET.Daily)[1]
    }
    if (is.null(Edate)) {
      Edate = time(results1$ET.Daily)[length(results1$ET.Daily)]
    }
    x.Date <- as.Date(time(results1$ET.Daily[which(time(results1$ET.Daily) == as.Date(Sdate)):which(time(results1$ET.Daily) == as.Date(Edate))]))
    
    comp <- matrix(nrow=length(x.Date),ncol=Ncomp+1)
    comp[,1] <- x.Date
    comp[,2] <- results1$ET.Daily[which(as.Date(time(results1$ET.Daily)) == x.Date[1]) : which(as.Date(time(results1$ET.Daily)) == x.Date[length(x.Date)])]
    comp[,3] <- results2$ET.Daily[which(as.Date(time(results2$ET.Daily)) == x.Date[1]) : which(as.Date(time(results2$ET.Daily)) == x.Date[length(x.Date)])]
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$ET.Daily)) {
          comp[,4] <- results3$ET.Daily[which(as.Date(time(results3$ET.Daily)) == x.Date[1]) : which(as.Date(time(results3$ET.Daily)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated daily ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$ET.Daily)) {
          comp[,5] <- results4$ET.Daily[which(as.Date(time(results4$ET.Daily)) == x.Date[1]) : which(as.Date(time(results4$ET.Daily)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated daily ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$ET.Daily)) {
          comp[,6] <- results5$ET.Daily[which(as.Date(time(results5$ET.Daily)) == x.Date[1]) : which(as.Date(time(results5$ET.Daily)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated daily ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$ET.Daily)) {
          comp[,7] <- results6$ET.Daily[which(as.Date(time(results6$ET.Daily)) == x.Date[1]) : which(as.Date(time(results6$ET.Daily)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated daily ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$ET.Daily)) {
          comp[,8] <- results7$ET.Daily[which(as.Date(time(results7$ET.Daily)) == x.Date[1]) : which(as.Date(time(results7$ET.Daily)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated daily ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.na(ylim[1])) {
      ylim[1] <- 0
    }
    if (is.na(ylim[2])) {
      ylim[2] <- max(results1$ET.Daily)*1.5
    }
    
    plot(comp[,2]~as.Date(comp[,1]), type = "l", main = "Daily estimations of ET",  ylim = ylim, col = 2, xlab = "Year", ylab = "Estimated ET mm/day", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~as.Date(comp[,1]), type = "o", pch = ".", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(x.Date),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the daily ET", xlim = c(0,1), ylim = ylim, col = 2, xlab = "Non-exceedance probability", ylab = "Estimated daily ET mm/day", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the daily estimations of ET", xlab = "", ylab = "Distribution of daily PET in mm", show.names = FALSE, pars = list(ylim = ylim, cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), fill = c(seq(from=2,to=Ncomp+1)), bty="o")
  }
  
  if (type == "Monthly") {
    if (is.null(Sdate)) {
      Sdate = time(results1$ET.Monthly)[1]
    }
    if (is.null(Edate)) {
      Edate = time(results1$ET.Monthly)[length(results1$ET.Monthly)]
    }
    
    # Monthly
    par(ask=FALSE)
    if (is.null(results1$ET.Monthly) | is.null(results2$ET.Monthly)) {
      stop("Unable to compare results because estimations of monthly ET are not available in both results1 and results2")
    }
    x.Date <- as.yearmon(time(results1$ET.Monthly[which(time(results1$ET.Monthly) == as.yearmon(Sdate)):which(time(results1$ET.Monthly) == as.yearmon(Edate))]))
    
    comp <- matrix(nrow=length(x.Date),ncol=Ncomp+1)
    comp[,1] <- x.Date
    comp[,2] <- results1$ET.Monthly[which(as.yearmon(time(results1$ET.Monthly)) == x.Date[1]) : which(as.yearmon(time(results1$ET.Monthly)) == x.Date[length(x.Date)])]
    comp[,3] <- results2$ET.Monthly[which(as.yearmon(time(results2$ET.Monthly)) == x.Date[1]) : which(as.yearmon(time(results2$ET.Monthly)) == x.Date[length(x.Date)])]
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$ET.Monthly)) {
          comp[,4] <- results3$ET.Monthly[which(as.yearmon(time(results3$ET.Monthly)) == x.Date[1]) : which(as.yearmon(time(results3$ET.Monthly)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated monthly ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$ET.Monthly)) {
          comp[,5] <- results4$ET.Monthly[which(as.yearmon(time(results4$ET.Monthly)) == x.Date[1]) : which(as.yearmon(time(results4$ET.Monthly)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated monthly ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$ET.Monthly)) {
          comp[,6] <- results5$ET.Monthly[which(as.yearmon(time(results5$ET.Monthly)) == x.Date[1]) : which(as.yearmon(time(results5$ET.Monthly)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated monthly ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$ET.Monthly)) {
          comp[,7] <- results6$ET.Monthly[which(as.yearmon(time(results6$ET.Monthly)) == x.Date[1]) : which(as.yearmon(time(results6$ET.Monthly)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated monthly ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$ET.Monthly)) {
          comp[,8] <- results7$ET.Monthly[which(as.yearmon(time(results7$ET.Monthly)) == x.Date[1]) : which(as.yearmon(time(results7$ET.Monthly)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated monthly ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ylim[1])) {
      ylim[1] <- 0
    }
    if (is.null(ylim[2])) {
      ylim[2] <- max(results1$ET.Monthly)*1.5
    }
    
    plot(comp[,2]~as.yearmon(comp[,1]), type = "l", main = "Monthly estimations of ET",  ylim = ylim, col = 2, xlab = "Year", ylab = "Monthly aggregations of estimated ET mm/month", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~as.yearmon(comp[,1]), type = "o", pch = ".", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(x.Date),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the monthly ET", xlim = c(0,1), ylim = ylim, col = 2, xlab = "Non-exceedance probability", ylab = "Monthly aggregations of estimated ET mm/month", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), main = "Box plot of the monthly estimations of ET", boxwex = 0.5, range = 0, xlab = "", ylab = "Distribution of monthly PET in mm", show.names = FALSE, pars=list(ylim = ylim, cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), fill = c(seq(from=2,to=Ncomp+1)), bty="o")
    
  }
  
  if (type == "Annual") {
    if (is.null(Sdate)) {
      Sdate = time(results1$ET.Annual)[1]
    }
    if (is.null(Edate)) {
      Edate = time(results1$ET.Annual)[length(results1$ET.Annual)]
    }
    
    # Annual
    par(ask=FALSE)
    if (is.null(results1$ET.Annual) | is.null(results2$ET.Annual)) {
      stop("Unable to compare results because estimations of annual ET are not available in both results1 and results2")
    }
    x.Date <- floor(as.numeric(as.yearmon(time(results1$ET.Annual[which(time(results1$ET.Annual) == floor(as.numeric(as.yearmon(Sdate)))):which(time(results1$ET.Annual) == floor(as.numeric(as.yearmon(Edate))))]))))
    
    comp <- matrix(nrow=length(x.Date),ncol=Ncomp+1)
    comp[,1] <- x.Date
    comp[,2] <- results1$ET.Annual[which(floor(as.numeric(as.yearmon(time(results1$ET.Annual)))) == x.Date[1]) : which(floor(as.numeric(as.yearmon(time(results1$ET.Annual)))) == x.Date[length(x.Date)])]
    comp[,3] <- results2$ET.Annual[which(floor(as.numeric(as.yearmon(time(results2$ET.Annual)))) == x.Date[1]) : which(floor(as.numeric(as.yearmon(time(results2$ET.Annual)))) == x.Date[length(x.Date)])]
    
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$ET.Annual)) {
          comp[,4] <- results3$ET.Annual[which(as.Date(time(results3$ET.Annual)) == x.Date[1]) : which(as.Date(time(results3$ET.Annual)) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated annual ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$ET.Annual)) {
          comp[,5] <- results4$ET.Annual[which(floor(as.numeric(as.yearmon(time(results4$ET.Annual)))) == x.Date[1]) : which(floor(as.numeric(as.yearmon(time(results4$ET.Annual)))) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated annual ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$ET.Annual)) {
          comp[,6] <- results5$ET.Annual[which(floor(as.numeric(as.yearmon(time(results5$ET.Annual)))) == x.Date[1]) : which(floor(as.numeric(as.yearmon(time(results5$ET.Annual)))) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated annual ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$ET.Annual)) {
          comp[,7] <- results6$ET.Annual[which(floor(as.numeric(as.yearmon(time(results6$ET.Annual)))) == x.Date[1]) : which(floor(as.numeric(as.yearmon(time(results6$ET.Annual)))) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated annual ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$ET.Annual)) {
          comp[,8] <- results7$ET.Annual[which(floor(as.numeric(as.yearmon(time(results7$ET.Annual)))) == x.Date[1]) : which(floor(as.numeric(as.yearmon(time(results7$ET.Annual)))) == x.Date[length(x.Date)])]
        } else {
          warning("No estimated annual ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ylim[1])) {
      ylim[1] <- 0
    }
    if (is.null(ylim[2])) {
      ylim[2] <- max(results1$ET.Annual)*1.5
    }
    
    plot(comp[,2]~floor(as.numeric(as.yearmon(comp[,1]))), type = "l", main = "Annual estimations of ET",  ylim = ylim, col = 2, xlab = "Year", ylab = "Annual aggregations of estimated ET mm/year", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~floor(as.numeric(as.yearmon(comp[,1]))), type = "l", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), fill = c(seq(from=2,to=Ncomp+1)), bty="o")
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(x.Date),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the annual ET", xlim = c(0,1), ylim = ylim, col = 2, xlab = "Non-exceedance probability", ylab = "Annual aggregations of estimated ET mm/year", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the annual estimations of ET", xlab = "", ylab = "Distribution of annual PET in mm", show.names = FALSE, pars = list(ylim = ylim, cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), fill = c(seq(from=2,to=Ncomp+1)), bty="o")
  }
  
  if (type == "MonthlyAve") {
    if (is.null(Sdate)) {
      Sdate = time(results1$ET.MonthlyAve)[1]
    }
    if (is.null(Edate)) {
      Edate = time(results1$ET.MonthlyAve)[length(results1$ET.MonthlyAve)]
    }
    
    # Monthly average
    par(ask=FALSE)
    if (is.null(results1$ET.MonthlyAve) | is.null(results2$ET.MonthlyAve)) {
      stop("Unable to compare results because estimations of monthly averaged ET are not available in both results1 and results2")
    }
    
    x.Date <- floor(as.numeric(as.yearmon(time(results1$ET.MonthlyAve[which(time(results1$ET.MonthlyAve) == floor(as.numeric(as.yearmon(Sdate)))):which(time(results1$ET.MonthlyAve) == floor(as.numeric(as.yearmon(Edate))))]))))
    
    comp <- matrix(nrow=length(results1$ET.MonthlyAve),ncol=Ncomp+1)
    comp[,1] <- x.Date
    comp[,2] <- results1$ET.MonthlyAve
    comp[,3] <- results2$ET.MonthlyAve
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$ET.MonthlyAve)) {
          comp[,4] <- results3$ET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$ET.MonthlyAve)) {
          comp[,5] <- results4$ET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$ET.MonthlyAve)) {
          comp[,6] <- results5$ET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$ET.MonthlyAve)) {
          comp[,7] <- results6$ET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$ET.MonthlyAve)) {
          comp[,8] <- results7$ET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ylim[1])) {
      ylim[1] <- 0
    }
    if (is.null(ylim[2])) {
      ylim[2] <- max(results1$ET.MonthlyAve)*1.5
    }
    
    plot(comp[,2]~unique(1:12), type = "l", main = "Monthly average ET",  ylim = ylim, col = 2, xlab = "Month", ylab = "Monthly averages of estimated daily ET mm/day", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~unique(1:12), type = "o", pch = ".", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(results1$ET.MonthlyAve),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the monthly average ET", xlim = c(0,1), ylim = ylim, col = 2, xlab = "Non-exceedance probability", ylab = "Monthly averages of estimated daily ET mm/day", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the monthly average ET", xlab = "", ylab = "Distribution of monthly averaged PET in mm/day", show.names = FALSE, pars = list(ylim = ylim, cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), fill = c(seq(from=2,to=Ncomp+1)), bty="o")
  }
  
  if (type == "AnnualAve") {
    
    # Annual average
    par(ask=FALSE)
    if (is.null(results1$ET.AnnualAve) | is.null(results2$ET.AnnualAve)) {
      stop("Unable to compare results because estimations of annually averaged ET are not available in both results1 and results2")
    }
    
    x.Date <- floor(as.numeric(as.yearmon(time(results1$ET.AnnualAve[which(time(results1$ET.AnnualAve) == floor(as.numeric(as.yearmon(Sdate)))):which(time(results1$ET.AnnualAve) == floor(as.numeric(as.yearmon(Edate))))]))))
    
    comp <- matrix(nrow=length(results1$ET.AnnualAve),ncol=Ncomp+1)
    comp[,1] <- x.Date
    comp[,2] <- results1$ET.AnnualAve
    comp[,3] <- results2$ET.AnnualAve
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$ET.AnnualAve)) {
          comp[,4] <- results3$ET.AnnualAve
        } else {
          warning("No estimated annually averaged ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$ET.AnnualAve)) {
          comp[,5] <- results4$ET.AnnualAve
        } else {
          warning("No estimated annually averaged ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$ET.AnnualAve)) {
          comp[,6] <- results5$ET.AnnualAve
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$ET.AnnualAve)) {
          comp[,7] <- results6$ET.AnnualAve
        } else {
          warning("No estimated annually averaged ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$ET.AnnualAve)) {
          comp[,8] <- results7$ET.AnnualAve
        } else {
          warning("No estimated annually averaged ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ylim[1])) {
      ylim[1] <- 0
    }
    if (is.null(ylim[2])) {
      ylim[2] <- max(results1$ET.AnnualAve)*1.5
    }
    
    plot(comp[,2]~floor(as.numeric(as.yearmon(time(results1$ET.Annual)))), type = "l", main = "Annually average ET",  ylim = ylim, col = 2, xlab = "Month", ylab = "Annual averages of estimated daily ET mm/day", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~floor(as.numeric(as.yearmon(time(results1$ET.Annual)))), type = "l", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(results1$ET.AnnualAve),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the annual average ET", xlim = c(0,1), ylim = ylim, col = 2, xlab = "Non-exceedance probability", ylab = "Annual averages of estimated daily ET mm/day", cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the annual average ET", xlab = "", ylab = "Distribution ofannually averaged PET in mm/day", show.names = FALSE, pars = list(ylim = ylim, cex.lab = 1.3, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (1), fill = c(seq(from=2,to=Ncomp+1)), bty="o")
  }
  par(ask=FALSE)
}


#-------------------------------------------------------------------------------------

ETForcings <- function(data, results, forcing) {
  par(ask=FALSE)
  plot.new()
  # define forcing names
  Fnames <- list(v=c("Tmax","Tmin","u2","uz","Rs","n","Precip","Epan","RHmax","RHmin","Tdew"), 
                 n=c("maximum temperature","minimum temperature", "average wind speed at 2m", "average wind speed", "daily solar radiation", "daily sunshine hours", 
                     "precipitation", "Class-A pan evaporation", "maximum relative humidity", "minimum relative humidity", "average dew point temeprature"),
                 u=c("degree Celcius","degree Celcius", "m/s", "m/s", "MJ.per.sqm.per.day", "hours", "mm/day", "mm/day", "%", "%", "degree Celcius"))
  
  # check forcing value
  if (forcing %in% Fnames$v) {
    par(mfrow=c(2,2))
    
    # Aggregated PET vs forcing
    
    # Plot daily relationship
    if (!is.null(results$ET.Daily)) {
      interval <- "daily"
      if (!is.null(data[[forcing]])) {
        plot(results$ET.Daily~data[[forcing]], main = results$ET_formulation, xlab = paste(interval, Fnames$n[Fnames$v == forcing], Fnames$u[Fnames$v == forcing]), ylab = list(c(results$ET_type, "mm/day")))
      } else {
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(1, paste("no", forcing, "data to plot"))
      }
    } else {
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1, "no daily ET estimation to plot")
    }
    
    # Plot monthly relationship
    if (!is.null(results$ET.Monthly)) {
      interval <- "monthly"
      if (!is.null(data[[forcing]])) {
        if (forcing != "Precip") {
          Fmonthly <- aggregate(data[[forcing]], as.yearmon(data$Date.daily, "%m/%y"), mean)
        } else {
          Fmonthly <- aggregate(data[[forcing]], as.yearmon(data$Date.daily, "%m/%y"), sum)
        }
        plot(results$ET.Monthly~Fmonthly, main = results$ET_formulation, xlab = paste(interval, Fnames$n[Fnames$v == forcing], Fnames$u[Fnames$v == forcing]), ylab = list(c(results$ET_type, "mm/month")))
      } else {
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(1, paste("no", forcing, "data to plot"))
      }
    } else {
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1, "no monthly ET estimation to plot")
    }
    
    
    # Plot annual relationship
    if (!is.null(results$ET.Annual)) {
      interval <- "annual"
      if (!is.null(data[[forcing]])) {
        if (forcing != "Precip") {
          Fannual <- aggregate(Fmonthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), mean)
        } else {
          Fannual <- aggregate(Fmonthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), sum)
        }
        plot(results$ET.Annual~Fannual, main = results$ET_formulation, xlab = paste(interval, Fnames$n[Fnames$v == forcing], Fnames$u[Fnames$v == forcing]), ylab = list(c(results$ET_type, "mm/year")))
      } else {
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(1, paste("no", forcing, "data to plot"))
      }
    } else {
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1, "no annual ET estimation to plot")
    }
    
    
  } else {
    stop("Please select forcing from 'Tmax','Tmin','u2','uz','Rs','n','Precip','Epan','RHmax','RHmin','Tdew'") # forcing is not defined
  }
  
}

#-------------------------------------------------------------------------------------