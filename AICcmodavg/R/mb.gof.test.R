##MacKenzie and Bailey goodness of fit test for single season occupancy models
##generic function to compute chi-square
mb.chisq <- function(mod, print.table = TRUE, ...){
  UseMethod("mb.chisq", mod)
}


mb.chisq.default <- function(mod, print.table = TRUE, ...){
  stop("\nFunction not yet defined for this object class\n")
}


##for single-season occupancy models of class unmarkedFitOccu
mb.chisq.unmarkedFitOccu <- function(mod, print.table = TRUE, ...) {
  
##step 1:
##extract detection histories
y.raw <- mod@data@y

##if some rows are all NA and sites are discarded, adjust sample size accordingly
N.raw <- nrow(y.raw)
#if(all NA) {N - number of rows with all NA}
##identify sites without data
na.raw <- apply(X = y.raw, MARGIN = 1, FUN = function(i) all(is.na(i)))
  
##remove sites without data
y.data <- y.raw[!na.raw, ]

N <- N.raw - sum(na.raw)

#N is required for computations in the end
Ts <- ncol(y.data)

det.hist <- apply(X = y.data, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))

##compute predicted values of occupancy
preds.psi <- predict(mod, type = "state")$Predicted



##compute predicted values of p
preds.p <- matrix(data = predict(mod, type = "det")$Predicted,
                  ncol = Ts, byrow = TRUE)

##assemble in data.frame
out.hist <- data.frame(det.hist, preds.psi)

##identify unique histories
un.hist <- unique(det.hist)
n.un.hist <- length(un.hist)

##identify if missing values occur
na.vals <- length(grep(pattern = "NA", x = un.hist)) > 0

if(na.vals) {

  ##identify each history with NA
  id.na <- grep(pattern = "NA", x = un.hist)
  id.det.hist.na <- grep(pattern = "NA", x = det.hist)
  
  ##cohorts with NA
  cohort.na <- sort(un.hist[id.na])
  n.cohort.na <- length(cohort.na)
  ##determine cohorts that will be grouped together (same missing value)
  unique.na <- gsub(pattern = "NA", replacement = "N", x = cohort.na)
  ##determine which visit has missing value
  na.visits <- sapply(strsplit(x = unique.na, split = ""), FUN = function(i) paste(ifelse(i == "N", 1, 0), collapse = ""))
  ##add cohort labels for histories
  names(cohort.na) <- na.visits  
  ##number of histories in each cohort
  n.hist.missing.cohorts <- table(na.visits)
  ##number of missing cohorts
  n.missing.cohorts <- length(n.hist.missing.cohorts)
  out.hist.na <- out.hist[id.det.hist.na, ]
  out.hist.na$det.hist <- droplevels(out.hist.na$det.hist)

  ##groupings in out.hist.na
  just.na <- sapply(X = out.hist.na$det.hist, FUN = function(i) gsub(pattern = "1", replacement = "0", x = i))
  out.hist.na$coh <- sapply(X = just.na, FUN = function(i) gsub(pattern = "NA", replacement = "1", x = i))                                           
  ##number of sites in each missing cohort
  freqs.missing.cohorts <- table(out.hist.na$coh)
    
  ##number of sites with each history
  na.freqs <- table(det.hist[id.det.hist.na]) 
  preds.p.na <- preds.p[id.det.hist.na, ]
  
  ##cohorts without NA
  cohort.not.na <- sort(un.hist[-id.na])
  out.hist.not.na <- out.hist[-id.det.hist.na, ]
  out.hist.not.na$det.hist <- droplevels(out.hist.not.na$det.hist)
  n.cohort.not.na <- length(cohort.not.na)
  n.sites.not.na <- length(det.hist) - length(id.det.hist.na)
  preds.p.not.na <- preds.p[-id.det.hist.na, ]

  
} else {
  cohort.not.na <- sort(un.hist)
  out.hist.not.na <- out.hist
  preds.p.not.na <- preds.p
  n.cohort.not.na <- length(cohort.not.na)
  n.sites.not.na <- length(det.hist)
}
    
##for each missing data cohort, determine number of sites for each
  
  ##iterate over each site for each unique history
if(n.cohort.not.na > 0) {  ##expected frequencies for non-missing data
  exp.freqs <- rep(NA, n.cohort.not.na)
  names(exp.freqs) <- cohort.not.na ########################################SORT ENCOUNTER HISTORIES CHECK THAT ORDER IS IDENTICAL TO OBSERVED FREQS
  
  ##iterate over detection histories
  for (i in 1:n.cohort.not.na) {
    eq.solved <- rep(NA, n.sites.not.na)
    select.hist <- cohort.not.na[i]
    ##strip all values
    strip.hist <- unlist(strsplit(select.hist, split = ""))
    
    ##translate each visit in probability statement
    hist.mat <- matrix(NA, nrow = n.sites.not.na, ncol = Ts)

    ##iterate over sites
    for(j in 1:n.sites.not.na) {

      ##in extreme cases where only a single cohort occurs without missing values
      if(n.sites.not.na == 1) {
        hist.mat[j, ] <- ifelse(strip.hist == "1", preds.p.not.na[j],
                                ifelse(strip.hist == "0", 1 - preds.p.not.na[j],
                                       0))
      } else {
        hist.mat[j, ] <- ifelse(strip.hist == "1", preds.p.not.na[j, ],
                                ifelse(strip.hist == "0", 1 - preds.p.not.na[j, ],
                                       0))
      }
        
      
      ##combine into equation
      combo.p <- paste(hist.mat[j, ], collapse = "*")
      ##for history without detection
      if(sum(as.numeric(strip.hist)) == 0) {
        combo.first <- paste(c(out.hist.not.na[j, "preds.psi"], combo.p), collapse = "*")
        combo.psi.p <- paste((1 - out.hist.not.na[j, "preds.psi"]), "+", combo.first)
      } else {
        combo.psi.p <- paste(c(out.hist.not.na[j, "preds.psi"], combo.p), collapse = "*")       
      }
      eq.solved[j] <- eval(parse(text = as.expression(combo.psi.p)))
    }
    
    exp.freqs[i] <- sum(eq.solved, na.rm = TRUE)
  }


  ##for each detection history, compute observed frequencies
  freqs <- table(out.hist.not.na$det.hist)
  
  out.freqs <- matrix(NA, nrow = n.cohort.not.na, ncol = 4)
  colnames(out.freqs) <- c("Cohort", "Observed", "Expected", "Chi-square")
  rownames(out.freqs) <- names(freqs)
  ##cohort
  out.freqs[, 1] <- 0
  ##observed
  out.freqs[, 2] <- freqs
  ##expected
  out.freqs[, 3] <- exp.freqs
  ##chi-square
  out.freqs[, 4] <- ((out.freqs[, "Observed"] - out.freqs[, "Expected"])^2)/out.freqs[, "Expected"]
}

##if missing values
if(na.vals) {
##create list to store the chisquare for each cohort
  missing.cohorts <- list( )
  ##check if preds.p.na has only 1 row and change to matrix
  if(!is.matrix(preds.p.na)) {preds.p.na <- matrix(data = preds.p.na, nrow = 1)}
  
  for(m in 1:n.missing.cohorts) {
    ##select cohort
    select.cohort <- out.hist.na[which(out.hist.na$coh == names(freqs.missing.cohorts)[m]), ]
    select.preds.p.na <- preds.p.na[which(out.hist.na$coh == names(freqs.missing.cohorts)[m]), ]
    ##replace NA's with 1 to remove from likelihood
    if(!is.matrix(select.preds.p.na)) {select.preds.p.na <- matrix(data = select.preds.p.na, nrow = 1)}
    select.preds.p.na[, gregexpr(pattern = "N", text = gsub(pattern = "NA", replacement = "N", x = select.cohort$det.hist[1]))[[1]]] <- 1
    n.total.sites <- nrow(select.cohort)
    freqs.na <- table(droplevels(select.cohort$det.hist))
    cohort.na.un <- sort(unique(select.cohort$det.hist))
    n.hist.na <- length(freqs.na)
    exp.na <- rep(NA, n.hist.na)
    names(exp.na) <- cohort.na.un  

          for(i in 1:n.hist.na) {
            ##number of sites in given history
            n.sites.hist <- freqs.na[i] ##this should be number of sites for each history
            eq.solved <- rep(NA, n.total.sites)
            ##replace NA's with N
            select.hist <- gsub(pattern = "NA", replacement = "N", x = cohort.na.un[i])
            ##strip all values
            strip.hist <- unlist(strsplit(select.hist, split = ""))
            ##translate each visit in probability statement
            hist.mat <- matrix(NA, nrow = n.total.sites, ncol = Ts)
            ##iterate over sites
            for(j in 1:n.total.sites) {
##################
              ##modified
              if(Ts == 1) {
                hist.mat[j, ] <- ifelse(strip.hist == "1", select.preds.p.na[j],
                                        ifelse(strip.hist == "0", 1 - select.preds.p.na[j], 1))
              } else {
                hist.mat[j, ] <- ifelse(strip.hist == "1", select.preds.p.na[j, ],
                                        ifelse(strip.hist == "0", 1 - select.preds.p.na[j, ], 1))
              }
              ##modified
##################  
              ##replace NA by 1 (missing visit is removed from likelihood)
###################################################
###for missing value, remove occasion
###################################################      
              ##combine into equation
              combo.p <- paste(hist.mat[j, ], collapse = "*")
              ##for history without detection
              if(sum(as.numeric(gsub(pattern = "N", replacement = "0", x = strip.hist))) == 0) {
                combo.first <- paste(c(select.cohort[j, "preds.psi"], combo.p), collapse = "*")
                combo.psi.p <- paste((1 - select.cohort[j, "preds.psi"]), "+", combo.first)
              } else {
                combo.psi.p <- paste(c(select.cohort[j, "preds.psi"], combo.p), collapse = "*")
              }  

              eq.solved[j] <- eval(parse(text = as.expression(combo.psi.p)))
            }
            
            exp.na[i] <- sum(eq.solved, na.rm = TRUE)
          }

    ##compute chisq for missing data cohorts
    
    ##for each detection history, compute observed frequencies
    out.freqs.na <- matrix(NA, nrow = n.hist.na, ncol = 4)
    colnames(out.freqs.na) <- c("Cohort", "Observed", "Expected", "Chi-square")
    rownames(out.freqs.na) <- cohort.na.un
    ##cohort
    out.freqs.na[, 1] <- m
    ##observed
    out.freqs.na[, 2] <- freqs.na
    ##expected
    out.freqs.na[, 3] <- exp.na
    ##chi-square
    out.freqs.na[, 4] <- ((out.freqs.na[, "Observed"] - out.freqs.na[, "Expected"])^2)/out.freqs.na[, "Expected"]
    
    missing.cohorts[[m]] <- list(out.freqs.na = out.freqs.na)
  }
  
  
}
##test statistic is chi-square for all possible detection histories
##for observed detection histories, chisq = sum((obs - exp)^2)/exp = Y

##to avoid computing all possible detection histories, it is possible to obtain it by subtraction:
#for unobserved detection histories, chisq = sum((obs - exp)^2)/exp = sum(0 - exp)^2/exp = sum(exp) = X
##X = N.sites - sum(exp values from observed detection histories)
##Thus, test statistic = Y + X = Y + N - sum(exp values from observed detection histories)

##compute partial chi-square for observed detection histories (Y)
#chisq.obs.det <- sum(((out.freqs[, "observed"] - out.freqs[, "expected"])^2)/out.freqs[, "expected"])

##compute partial chi-square for unobserved detection histories (X)
if(na.vals) {
  chisq.missing <- do.call("rbind", lapply(missing.cohorts, FUN = function(i) i$out.freqs.na))
  if(n.cohort.not.na > 0) {
    chisq.unobs.det <- N - sum(out.freqs[, "Expected"]) - sum(chisq.missing[, "Expected"])
    chisq.table <- rbind(out.freqs, chisq.missing)
  } else {
    chisq.unobs.det <- N - sum(chisq.missing[, "Expected"])
    chisq.table <- chisq.missing
  }  
} else {
  chisq.unobs.det <- N - sum(out.freqs[, "Expected"])
  chisq.na <- 0
  chisq.table <- out.freqs
}

##test statistic (Y + X = Y - sum(exp observed detection histories)
chisq <- sum(chisq.table[, "Chi-square"]) + chisq.unobs.det

if(print.table) {
  out <- list(chisq.table = chisq.table, chi.square = chisq,
              model.type = "single-season")
} else {
  out <- list(chi.square = chisq, model.type = "single-season")
}

class(out) <- "mb.chisq"
return(out)
}



##dynamic occupancy models of class unmarkedFitColExt
mb.chisq.unmarkedFitColExt <- function(mod, print.table = TRUE, ...) {
  
  ##information on data set
  ##extract number of seasons
  orig.data <- mod@data
  detections <- orig.data@y
  n.seasons <- orig.data@numPrimary
  n.sites <- nrow(detections)

  ##total visits
  total.visits <- ncol(detections)
  ##number of visits per season
  n.visits <- total.visits/n.seasons

  ##determine if certain sites have no data for certain visits
  
  
  
  ##split encounter history for each season
  starts <- seq(1, total.visits, by = n.visits)

  ##create list to hold seasons
  y.seasons <- list( )
  
  for(k in 1:n.seasons) {
    first.col <- starts[k]
    y.seasons[[k]] <- detections[, first.col:(first.col+n.visits-1)]
  }


#################################
  ##from model
  ##predict init psi
  psi.init.pred <- predict(mod, type = "psi")$Predicted

  ##predict gamma
  gam.pred <- matrix(data = predict(mod, type = "col")$Predicted,
                     ncol = n.seasons, byrow = TRUE)[, 1:(n.seasons - 1)]

  ##predict epsilon
  eps.pred <- matrix(data = predict(mod, type = "ext")$Predicted,
                     ncol = n.seasons, byrow = TRUE)[, 1:(n.seasons - 1)]

  ##compute predicted values of psi for each season
  ##predicted init psi
  psi.init.pred <- predict(mod, type = "psi")$Predicted

  ##predicted p
  p.pred <- matrix(data = predict(mod, type = "det")$Predicted,
                   ncol = total.visits,
                   nrow = n.sites, byrow = TRUE)

  ##divide p's for each season
  p.seasons <- list( )
  for(k in 1:n.seasons) {
    first.col <- starts[k]
    p.seasons[[k]] <- p.pred[, first.col:(first.col+n.visits-1)]
  }


  ##compute psi for each season
  psi.seasons <- list( )
  ##add first year in list
  psi.seasons[[1]] <- psi.init.pred

  ##compute psi recursively for each year given psi(t - 1), epsilon, and gamma
  if(n.seasons == 2) {
      psi.seasons[[2]] <- psi.seasons[[1]] * (1 - eps.pred) + (1 - psi.seasons[[1]]) * gam.pred
    } else {
      for(m in 2:(n.seasons)){
        psi.seasons[[m]] <- psi.seasons[[m-1]] * (1 - eps.pred[, m-1]) + (1 - psi.seasons[[m-1]]) * gam.pred[, m-1]
      }
    }
  

################################################
###the following is adapted from mb.chisq

  ##create list to hold table for each year
  out <- vector(mode = "list", length = n.seasons)

  ##add label for seasons
  season.labels <- paste("Season", 1:n.seasons, sep = "")
  names(out) <- season.labels
  
  #names(all.chisq) <- season.labels

  ##iterate over seasons
  for (season in 1:n.seasons) {
    
    ##step 1:
    ##extract detection histories
    y.raw <- y.seasons[[season]]

    ##if some rows are all NA and sites are discarded, adjust sample size accordingly
    N.raw <- nrow(y.raw)
    ##if(all NA) {N - number of rows with all NA}
    ##identify sites without data
    na.raw <- apply(X = y.raw, MARGIN = 1, FUN = function(i) all(is.na(i)))
  
    ##remove sites without data
    #y.data <- y.raw[!na.raw, ] #in mb.chisq( ) for single season data, sites without data are removed
    ##with multiseason model, missing data in some years create an imbalance - better to maintain number of sites constant across years
    ##this creates a new cohort with only missing values
    y.data <- y.raw

    ##number of observed detection histories (excludes cases with all NA's)
    N <- N.raw - sum(na.raw)

    ##check if any columns are empty
    nodata.raw <- apply(X = y.data, MARGIN = 2, FUN = function(i) all(is.na(i)))
    with.data <- which(nodata.raw == FALSE)
    ##select only columns with data
    if(any(nodata.raw)) {
      y.data <- y.data[, with.data]
    }

    ##T is required for computations in the end
    ##if only a single visit, ncol( ) returns NULL
###########################
    ##modified
    if(is.vector(y.data)){
      Ts <- 1
    } else { 
      Ts <- ncol(y.data)
    }

    ##if single visit, returns error
    if(Ts == 1) {
      det.hist <- paste(y.data, sep = "")
    } else {
      det.hist <- apply(X = y.data, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
    }
    ##modified
###########################
    
    ##compute predicted values of occupancy for season i
    preds.psi <- psi.seasons[[season]] ##MODIFIED FROM mb.chisq - change iteration number

    ##extract matrix of p's for season i
    if(any(nodata.raw)) {
      preds.p <- p.seasons[[season]][, with.data]
    } else {
      preds.p <- p.seasons[[season]]
    }

    ##assemble in data.frame
    out.hist <- data.frame(det.hist, preds.psi)

    ##identify unique histories
    un.hist <- unique(det.hist)
    n.un.hist <- length(un.hist)
    
    ##identify if missing values occur
    na.vals <- length(grep(pattern = "NA", x = un.hist)) > 0

    if(na.vals) {
    
      ##identify each history with NA
      id.na <- grep(pattern = "NA", x = un.hist)
      id.det.hist.na <- grep(pattern = "NA", x = det.hist)
  
      ##cohorts with NA
      cohort.na <- sort(un.hist[id.na])
      n.cohort.na <- length(cohort.na)
      ##determine cohorts that will be grouped together (same missing value)
      unique.na <- gsub(pattern = "NA", replacement = "N", x = cohort.na)
      ##determine which visit has missing value
      na.visits <- sapply(strsplit(x = unique.na, split = ""), FUN = function(i) paste(ifelse(i == "N", 1, 0), collapse = ""))
      ##add cohort labels for histories
      names(cohort.na) <- na.visits  
      ##number of histories in each cohort
      n.hist.missing.cohorts <- table(na.visits)
      ##number of missing cohorts
      n.missing.cohorts <- length(n.hist.missing.cohorts)
      out.hist.na <- out.hist[id.det.hist.na, ]
      out.hist.na$det.hist <- droplevels(out.hist.na$det.hist)

      ##groupings in out.hist.na
      just.na <- sapply(X = out.hist.na$det.hist, FUN = function(i) gsub(pattern = "1", replacement = "0", x = i))
      out.hist.na$coh <- sapply(X = just.na, FUN = function(i) gsub(pattern = "NA", replacement = "1", x = i))                                           
      ##number of sites in each missing cohort
      freqs.missing.cohorts <- table(out.hist.na$coh)
      
      ##number of sites with each history
      na.freqs <- table(det.hist[id.det.hist.na]) 
#####################
      ##modified
      if(Ts == 1) {
        preds.p.na <- preds.p[id.det.hist.na]
      } else {
        preds.p.na <- preds.p[id.det.hist.na, ]
      }
      ##modified
#####################
  
      ##cohorts without NA
      cohort.not.na <- sort(un.hist[-id.na])
      out.hist.not.na <- out.hist[-id.det.hist.na, ]
      out.hist.not.na$det.hist <- droplevels(out.hist.not.na$det.hist)
      n.cohort.not.na <- length(cohort.not.na)
      n.sites.not.na <- length(det.hist) - length(id.det.hist.na)
#####################
      ##modified
      if(Ts == 1) {
        preds.p.not.na <- preds.p[-id.det.hist.na]
      } else {
        preds.p.not.na <- preds.p[-id.det.hist.na, ]
      }
      ##modified
#####################

      
    } else {
      cohort.not.na <- sort(un.hist)
      out.hist.not.na <- out.hist
      preds.p.not.na <- preds.p
      n.cohort.not.na <- length(cohort.not.na)
      n.sites.not.na <- length(det.hist)
    }
    
    ##for each missing data cohort, determine number of sites for each
    
    ##iterate over each site for each unique history
    if(n.cohort.not.na > 0) {  ##expected frequencies for non-missing data
      exp.freqs <- rep(NA, n.cohort.not.na)
      names(exp.freqs) <- cohort.not.na ########################################SORT ENCOUNTER HISTORIES CHECK THAT ORDER IS IDENTICAL TO OBSERVED FREQS
  
      ##iterate over detection histories
      for (i in 1:n.cohort.not.na) {
        eq.solved <- rep(NA, n.sites.not.na)
        select.hist <- cohort.not.na[i]
        ##strip all values
        strip.hist <- unlist(strsplit(select.hist, split = ""))
    
        ##translate each visit in probability statement
        hist.mat <- matrix(NA, nrow = n.sites.not.na, ncol = Ts)
        
        ##iterate over sites
        for(j in 1:n.sites.not.na) {

          ##in extreme cases where only a single cohort occurs without missing values
############
          ##modified
          if(n.sites.not.na == 1 || Ts == 1) {
            ##modified
#############            
            hist.mat[j, ] <- ifelse(strip.hist == "1", preds.p.not.na[j],
                                    ifelse(strip.hist == "0", 1 - preds.p.not.na[j],
                                           0))
          } else {
            hist.mat[j, ] <- ifelse(strip.hist == "1", preds.p.not.na[j, ],
                                    ifelse(strip.hist == "0", 1 - preds.p.not.na[j, ],
                                           0))
          }
          
      
          ##combine into equation
          combo.p <- paste(hist.mat[j, ], collapse = "*")
          ##for history without detection
          if(sum(as.numeric(strip.hist)) == 0) {
            combo.first <- paste(c(out.hist.not.na[j, "preds.psi"], combo.p), collapse = "*")
            combo.psi.p <- paste((1 - out.hist.not.na[j, "preds.psi"]), "+", combo.first)
          } else {
            combo.psi.p <- paste(c(out.hist.not.na[j, "preds.psi"], combo.p), collapse = "*")       
          }
          eq.solved[j] <- eval(parse(text = as.expression(combo.psi.p)))
        }
    
        exp.freqs[i] <- sum(eq.solved, na.rm = TRUE)
      }

      
      ##for each detection history, compute observed frequencies
      freqs <- table(out.hist.not.na$det.hist)
  
      out.freqs <- matrix(NA, nrow = n.cohort.not.na, ncol = 4)
      colnames(out.freqs) <- c("Cohort", "Observed", "Expected", "Chi-square")
      rownames(out.freqs) <- names(freqs)
      ##cohort
      out.freqs[, 1] <- 0
      ##observed
      out.freqs[, 2] <- freqs
      ##expected
      out.freqs[, 3] <- exp.freqs
      ##chi-square
      out.freqs[, 4] <- ((out.freqs[, "Observed"] - out.freqs[, "Expected"])^2)/out.freqs[, "Expected"]
    }

    ##if missing values
    if(na.vals) {
      ##create list to store the chisquare for each cohort
      missing.cohorts <- list( )
      ##check if preds.p.na has only 1 row and change to matrix
      if(!is.matrix(preds.p.na)) {preds.p.na <- matrix(data = preds.p.na, nrow = 1)}
    
      for(m in 1:n.missing.cohorts) {
        ##select cohort
        select.cohort <- out.hist.na[which(out.hist.na$coh == names(freqs.missing.cohorts)[m]), ]
#######################
        ##modified
        if(Ts == 1) {
          select.preds.p.na <- preds.p.na[which(out.hist.na$coh == names(freqs.missing.cohorts)[m])]
        } else {
          select.preds.p.na <- preds.p.na[which(out.hist.na$coh == names(freqs.missing.cohorts)[m]), ]
        }
        ##modified
#######################
        ##replace NA's with 1 to remove from likelihood
        if(!is.matrix(select.preds.p.na)) {select.preds.p.na <- matrix(data = select.preds.p.na, nrow = 1)}
        select.preds.p.na[, gregexpr(pattern = "N", text = gsub(pattern = "NA", replacement = "N", x = select.cohort$det[1]))[[1]]] <- 1
        n.total.sites <- nrow(select.cohort)
        freqs.na <- table(droplevels(select.cohort$det.hist))
        cohort.na.un <- sort(unique(select.cohort$det.hist))
        n.hist.na <- length(freqs.na)
        exp.na <- rep(NA, n.hist.na)
        names(exp.na) <- cohort.na.un  

        for(i in 1:n.hist.na) {
          ##number of sites in given history
          n.sites.hist <- freqs.na[i] ##this should be number of sites for each history
          eq.solved <- rep(NA, n.total.sites)
          ##replace NA's with N
          select.hist <- gsub(pattern = "NA", replacement = "N", x = cohort.na.un[i])
          ##strip all values
          strip.hist <- unlist(strsplit(select.hist, split = ""))
          ##translate each visit in probability statement
          hist.mat <- matrix(NA, nrow = n.total.sites, ncol = Ts)
          ##iterate over sites
          for(j in 1:n.total.sites) {
            hist.mat[j, ] <- ifelse(strip.hist == "1", select.preds.p.na[j, ],
                                    ifelse(strip.hist == "0", 1 - select.preds.p.na[j, ], 1))
            ##replace NA by 1 (missing visit is removed from likelihood)
###################################################
###for missing value, remove occasion
###################################################      
            ##combine into equation
            combo.p <- paste(hist.mat[j, ], collapse = "*")
            ##for history without detection
            if(sum(as.numeric(gsub(pattern = "N", replacement = "0", x = strip.hist))) == 0) {
              combo.first <- paste(c(select.cohort[j, "preds.psi"], combo.p), collapse = "*")
              combo.psi.p <- paste((1 - select.cohort[j, "preds.psi"]), "+", combo.first)
            } else {
              combo.psi.p <- paste(c(select.cohort[j, "preds.psi"], combo.p), collapse = "*")
            }  
            
            eq.solved[j] <- eval(parse(text = as.expression(combo.psi.p)))
          }
            
          exp.na[i] <- sum(eq.solved, na.rm = TRUE)
        }

        ##compute chisq for missing data cohorts
        
        ##for each detection history, compute observed frequencies
        out.freqs.na <- matrix(NA, nrow = n.hist.na, ncol = 4)
        colnames(out.freqs.na) <- c("Cohort", "Observed", "Expected", "Chi-square")
        rownames(out.freqs.na) <- cohort.na.un
        ##cohort
        out.freqs.na[, 1] <- m
        ##observed
        out.freqs.na[, 2] <- freqs.na
        ##expected
        out.freqs.na[, 3] <- exp.na
        ##chi-square
        out.freqs.na[, 4] <- ((out.freqs.na[, "Observed"] - out.freqs.na[, "Expected"])^2)/out.freqs.na[, "Expected"]
        
        missing.cohorts[[m]] <- list(out.freqs.na = out.freqs.na)
      }
      
  
    }
    ##test statistic is chi-square for all possible detection histories
    ##for observed detection histories, chisq = sum((obs - exp)^2)/exp = Y

    ##to avoid computing all possible detection histories, it is possible to obtain it by subtraction:
    ##for unobserved detection histories, chisq = sum((obs - exp)^2)/exp = sum(0 - exp)^2/exp = sum(exp) = X
    ##X = N.sites - sum(exp values from observed detection histories)
    ##Thus, test statistic = Y + X = Y + N - sum(exp values from observed detection histories)

    ##compute partial chi-square for observed detection histories (Y)
    ##chisq.obs.det <- sum(((out.freqs[, "observed"] - out.freqs[, "expected"])^2)/out.freqs[, "expected"])

    ##compute partial chi-square for unobserved detection histories (X)
    if(na.vals) {
      chisq.missing <- do.call("rbind", lapply(missing.cohorts, FUN = function(i) i$out.freqs.na))

      ##check for sites never sampled in table 
      sites.never <- rownames(chisq.missing)
      never.sampled <- grep(pattern = paste(rep(NA, Ts), collapse = ""), x = sites.never)

      if(length(never.sampled) > 0) {
        ##remove row for site never sampled
        chisq.missing <- chisq.missing[-never.sampled, ]
      }

      if(n.cohort.not.na > 0) {
        chisq.unobs.det <- N - sum(out.freqs[, "Expected"]) - sum(chisq.missing[, "Expected"])
        chisq.table <- rbind(out.freqs, chisq.missing)
      } else {
        chisq.unobs.det <- N - sum(chisq.missing[, "Expected"])
        chisq.table <- chisq.missing
      }  
    } else {
      chisq.unobs.det <- N - sum(out.freqs[, "Expected"])
      chisq.na <- 0
      chisq.table <- out.freqs
    }

    ##test statistic (Y + X = Y - sum(exp observed detection histories)
    chisq <- sum(chisq.table[, "Chi-square"]) + chisq.unobs.det

    
    if(print.table) {
      out[[season]] <- list(chisq.table = chisq.table, chi.square = chisq) #change to iteration number
    } else {
      out[[season]] <- list(chi.square = chisq)
    }
  }

  ##add a final named element combining the chi-square
  all.chisq <- unlist(lapply(out, FUN = function(i) i$chi.square))

  ##add label for seasons
  new.season.labels <- paste("Season", 1:n.seasons, sep = " ")
  names(all.chisq) <- new.season.labels
  
  ##print table or only test statistic
  if(print.table) {
    out.nice <- list(tables = out, all.chisq = all.chisq,
                     n.seasons = n.seasons, model.type = "dynamic")
  } else {
    out.nice <- list(all.chisq = all.chisq, n.seasons = n.seasons, model.type = "dynamic")
  }
  
  class(out.nice) <- "mb.chisq" 
  return(out.nice)
}




##simulating data from model to compute P-value of test statistic
##create generic mb.gof.test 
mb.gof.test <- function(mod, nsim = 5, plot.hist = TRUE, ...){
  UseMethod("mb.gof.test", mod)
}

mb.gof.test.default <- function(mod, nsim = 5, plot.hist = TRUE, ...){
  stop("\nFunction not yet defined for this object class\n")
}


##for single-season occupancy models of class unmarkedFitOccu
mb.gof.test.unmarkedFitOccu <- function(mod, nsim = 5, plot.hist = TRUE, ...){#more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract table from fitted model
  mod.table <- mb.chisq(mod)

  ##compute GOF P-value
  out <- parboot(mod, statistic = function(i) mb.chisq(i)$chi.square,
                 nsim = nsim)
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
    p.display <- paste("<", 1/nsim)
  } else {
    p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
  if(plot.hist) {
  hist(out@t.star, main = paste("Bootstrapped MacKenzie and Bailey fit statistic (", nsim, " samples)", sep = ""),
       xlim = range(c(out@t.star, out@t0)), xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""))
  title(main = bquote(paste(italic(P), " ", .(p.display))), line = 0.5)
  abline(v = out@t0, lty = "dashed", col = "red")
}
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)


##assemble result
gof.out <- list(model.type = mod.table$model.type, chisq.table = mod.table$chisq.table, chi.square = mod.table$chi.square,
                t.star = out@t.star, p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "mb.chisq"
  return(gof.out)  
}



##dynamic occupancy models of class unmarkedFitColExt
mb.gof.test.unmarkedFitColExt <- function(mod, nsim = 5, plot.hist = TRUE, plot.seasons = FALSE, ...){#more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract table from fitted model
  mod.table <- mb.chisq(mod)
  n.seasons <- mod.table$n.seasons
  n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    
  ##compute GOF P-value
  out <- parboot(mod, statistic = function(i) mb.chisq(i)$all.chisq, #extract chi-square for each year
                 nsim = nsim)

  ##list to hold results
  p.vals <- list( )

  

  ##if only season-specific plots are requested
  if(!plot.hist && plot.seasons) {
    ##determine arrangement of plots in matrix
    if(plot.seasons && n.seasons > 12) {
      n.seasons.adj <- 12
      warning("\nOnly first 12 seasons are plotted\n")
    }
    
    if(plot.seasons && n.seasons.adj <= 12) {

      ##if n.seasons < 12
      ##if 12, 11, 10 <- 4 x 3
      ##if 9, 8, 7 <- 3 x 3
      ##if 6, 5 <- 3 x 2
      ##if 4 <- 2 x 2
      ##if 3 <- 3 x 1
      ##if 2 <- 2 x 1
    
      if(n.seasons.adj >= 10) {
        par(mfrow = c(4, 3))
      } else {

        if(n.seasons.adj >= 7) {
          par(mfrow = c(3, 3))
        } else {

          if(n.seasons.adj >= 5) {
            par(mfrow = c(3, 2))
          } else {
            if(n.seasons.adj == 4) {
              par(mfrow = c(2, 2))
            } else {
              if(n.seasons.adj == 3) {
                par(mfrow = c(3, 1))
              } else {
                par(mfrow = c(2, 1))
              }
            }
          }
        }
      }
    }
  }

  
  ##if both plots for seasons and summary are requested
  if(plot.hist && plot.seasons){
    ##determine arrangement of plots in matrix
    if(plot.seasons && n.seasons > 12) {
      n.seasons.adj <- 11
      warning("\nOnly first 11 seasons are plotted\n")
    }

    if(plot.seasons && n.seasons.adj <= 11) {

      if(n.seasons.adj >= 9) {
        par(mfrow = c(4, 3))
      } else {

        if(n.seasons.adj >= 6) {
          par(mfrow = c(3, 3))
        } else {

          if(n.seasons.adj >= 4) {
            par(mfrow = c(3, 2))
          } else {
            if(n.seasons.adj == 3) {
              par(mfrow = c(2, 2))
            } else {
              if(n.seasons.adj == 2) {
                par(mfrow = c(3, 1))
              }
            }
          }
        }
      }
    }
  }

    
  ##determine significance for each season
  for(k in 1:n.seasons) {

    p.value <- sum(out@t.star[, k] >= out@t0[k])/nsim
    if(p.value == 0) {
      p.display <- paste("<", 1/nsim)
    } else {
      p.display = paste("=", round(p.value, digits = 4))
    }
    
    p.vals[[k]] <- list("p.value" = p.value, "p.display" = p.display)
  }  
      
  ##create plot for first 12 plots
  if(plot.seasons) {
    for(k in 1:n.seasons.adj) {
      hist(out@t.star[, k], main = paste("Bootstrapped MacKenzie and Bailey fit statistic (", nsim, " samples) - season ", k, sep = ""),
           xlim = range(c(out@t.star[, k], out@t0[k])),
           xlab = paste("Simulated statistic ", "(observed = ", round(out@t0[k], digits = 2), ")", sep = ""))
      title(main = bquote(paste(italic(P), " ", .(p.vals[[k]]$p.display))), line = 0.5)
      abline(v = out@t0[k], lty = "dashed", col = "red")
      
    }
  }


  ##estimate c-hat
  obs.chisq <- sum(mod.table$all.chisq)
  boot.chisq <- sum(colMeans(out@t.star))
  c.hat.est <- obs.chisq/boot.chisq
  
  all.p.vals <- lapply(p.vals, FUN = function(i) i$p.value)
  ##lapply(mod.table, FUN = function(i) i$chisq.table)

  ##compute P-value for obs.chisq
  sum.chisq <- rowSums(out@t.star)
  p.global <- sum(sum.chisq >= obs.chisq)/nsim 

  if(p.global == 0) {
    p.global.display <- paste("< ", 1/nsim)
  } else {
    p.global.display <- round(p.global, digits = 4)
  }

  ##optionally show sum of chi-squares
  ##create plot
  if(plot.hist) {
    hist(sum.chisq, main = paste("Bootstrapped sum of chi-square statistic (", nsim, " samples)", sep = ""),
         xlim = range(c(sum.chisq, obs.chisq)),
         xlab = paste("Simulated statistic ", "(observed = ", round(obs.chisq, digits = 2), ")", sep = ""))
    title(main = bquote(paste(italic(P), " ", .(p.vals[[k]]$p.display))), line = 0.5)
    abline(v = out@t0[k], lty = "dashed", col = "red")
  }

  
  ##assemble result
  gof.out <- list(model.type = mod.table$model.type, chisq.table = mod.table,
                  chi.square = obs.chisq, t.star = sum.chisq,
                  p.value = all.p.vals, p.global = p.global, c.hat.est = c.hat.est,
                  nsim = nsim, n.seasons = n.seasons)
  class(gof.out) <- "mb.chisq"
  return(gof.out)  
}



##print function
print.mb.chisq <- function(x, digits.vals = 2, digits.chisq = 4, ...) {
  ##single-season occupancy models
  if(identical(x$model.type, "single-season")) {
    cat("\nMacKenzie and Bailey goodness-of-fit for single-season occupancy model\n")
    if(any(names(x) == "chisq.table")) {
      cat("\nPearson chi-square table:\n\n")
      print(round(x$chisq.table, digits = digits.vals))
    }
    cat("\nChi-square statistic =", round(x$chi.square, digits = digits.chisq), "\n")
    if(any(names(x) == "c.hat.est")) {
      cat("Number of bootstrap samples =", x$nsim)
      cat("\nP-value =", x$p.value)
      cat("\n\nQuantiles of bootstrapped statistics:\n")
      print(quantile(x$t.star), digits = digits.vals)
      cat("\nEstimate of c-hat =", round(x$c.hat.est, digits = digits.vals), "\n")
    }
    cat("\n")
  }

    ##dynamic occupancy models
    if(identical(x$model.type, "dynamic")) {
      cat("\nGoodness-of-fit for dynamic occupancy model\n")
      cat("\nNumber of seasons: ",  x$n.seasons, "\n")
      ##cat("\nPearson chi-square table:\n\n")
      ##print(round(x$chisq.table, digits = digits.vals))
      ##x$chisq.table
      cat("\nChi-square statistic:\n")
      if(any(names(x) == "all.chisq")) {
      print(round(x$all.chisq, digits = digits.chisq))
      cat("\nTotal chi-square =", round(sum(x$all.chisq), digits = digits.chisq), "\n")
    } else {
      print(round(x$chisq.table$all.chisq, digits = digits.chisq))
      cat("\nTotal chi-square =", round(sum(x$chi.square), digits = digits.chisq), "\n")
      cat("Number of bootstrap samples =", x$nsim)
      cat("\nP-value =", x$p.global)
      cat("\n\nQuantiles of bootstrapped statistics:\n")
      print(quantile(x$t.star), digits = digits.vals)
      cat("\nEstimate of c-hat =", round(x$c.hat.est, digits = digits.vals), "\n")
    }
      cat("\n")
    }
}

