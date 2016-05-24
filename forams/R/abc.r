`abc` <-
  function (df, Perm = 1000, confInt = 95) {
    Abundance <- df[order(df[, 1], decreasing=TRUE), ][1]
    Biomass <- df[order(df[, 2], decreasing=TRUE), ][2]
    PercAbun <- NA
    PercAbun[1] <- Abundance[1,] / sum(Abundance) * 100
    for (i in 2:nrow(Abundance)) {
      PercAbun[i] <- (Abundance[i,] / sum(Abundance)) * 100 + PercAbun[i-1]
    }
    PercAbun <- data.frame(PercAbun)
    colnames(PercAbun) <- 'Accum.Abund'
    rownames(PercAbun) <- rownames(Abundance)
    PercBio <- NA
    PercBio[1] <- Biomass[1,] / sum(Biomass) * 100
    for (i in 2:nrow(Biomass)) {
      PercBio[i] <- (Biomass[i,] / sum(Biomass)) * 100 + PercBio[i-1]
    }
    PercBio <- data.frame(PercBio)
    colnames(PercBio) <- 'Accum.Biomass'
    rownames(PercBio) <- rownames(Biomass)
    BiAi <- NA
    for(i in 2:nrow(df)) {
      BiAi[i] <- PercBio[i,] - PercAbun[i,]
    }
    abc <- cbind(PercAbun, PercBio[match(rownames(PercAbun),
                                         rownames(PercBio)),])
    abc[3] <- abc[2] - abc[1]
    colnames(abc) <- c('Accum.Abund', 'Accum.Biomass', 'BiAi')
    rownames(abc) <- rownames(PercAbun)
    
    wci <- function(x, permutations = Perm) {
      W <- function(x) { 
        return(round(sum(x)/ (50 * (length(x) - 1)), 4))
      }
      boots <- numeric(permutations)
      for (i in 1:length(boots)) {
        boots[i] <- W(sample(x, replace=T))
      }
      bias <- mean(boots) - W(x)
      
      CI <- new("numeric")
      
      if (confInt == 95) {
        CI[1] <- W(x) - bias - 1.95996*sqrt(var(boots))     
        CI[3] <- W(x) - bias + 1.95996*sqrt(var(boots))
      }
      else if (confInt == 90) {
        CI[1] <- W(x) - bias - 1.64485*sqrt(var(boots))     
        CI[3] <- W(x) - bias + 1.64485*sqrt(var(boots))
      }
      else if (confInt == 99) {
        CI[1] <- W(x) - bias - 2.57583*sqrt(var(boots))     
        CI[3] <- W(x) - bias + 2.57583*sqrt(var(boots))
      }
      else {
        stop("Invalid value, choose between 90, 95 or 99 percent 
             Confidence Interval.")
      }
      
      CI[2] <- W(x)
      
      if (confInt == 95) {
        names(CI) <- c("Left 95% CI", "W Statistic", "Right 95% CI")
      }
      else if (confInt == 90) {
        names(CI) <- c("Left 90% CI", "W Statistic", "Right 90% CI")
      }
      else if (confInt == 99) {
        names(CI) <- c("Left 99% CI", "W Statistic", "Right 99% CI")
      }
      else {
        stop("Invalid value, choose between 90, 95 or 99 percent 
             Confidence Interval.")
      }
      
      return(CI)
    }    
    Table <- new("abc")
    Table@abc <- as.data.frame(abc[1:3])
    rownames(Table@abc) <- attr(abc, "row.names")
    Table@abc <- Table@abc[sort(rownames(Table@abc)),]
    Table@W.Stat <- wci(Table@abc$BiAi)
    
    class(Table) <- "abc"
    return(Table)
  }
