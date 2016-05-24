  # drops electrodes that are not specified
.cluster.seg <- function(data, electrodes, window, ...) {
    elec.cluster <- subset(data, select = c("Subject", "Stimulus", "Time", electrodes))
}

  # --------------------------------------------------------------------------------

  # averages selected electrodes together (useful for dense array)
.avg.subs <- function(data, electrodes, window, toAvg, Time.range, trial.types) {
    if (ncol(toAvg) <= 4) {
        colnames(toAvg)[4] <- "Means"
        avg.sub = toAvg
    } else {
        avg.sub <- data.frame(toAvg[, 1:3], Means = rowMeans(toAvg[, 4:ncol(toAvg)]))
    }
}

# --------------------------------------------------------------------------------

  # creates a list of means for each subject for each loaded condition
.ind.by.cond <- function(data, electrodes, window, Time.range,
                          avgsub, trial.types, Stimulus, num.subs, num.conditions) {
  with(avgsub, { # see .get.ga.mamps comment
  means.cond <- vector("list")  # individual means by condition
  means.cond.sub <- vector("list")  # individual means for each subject, by condition
    for (i in 1:num.conditions) {
      means.cond[[i]] <- subset(avgsub, avgsub$Stimulus == trial.types[i], select = Means)
    }
    for (j in 1:num.conditions) {
      means.cond.sub[[j]] <- matrix(unlist(means.cond[j]), ncol = num.subs)
    }
  return(means.cond.sub)
  }
  )
}  # close grand.average

# --------------------------------------------------------------------------------

# creates a data frame of grand average voltages by condition
.grand.average <- function(data, electrodes, window, Time.range,
                        avgsub, trial.types, Stimulus, num.subs, num.conditions) {
  Means <- NULL # to appease R CMD check no visible binding NOTE
  means.cond <- vector("list")  # individual means by condition
  means.cond.sub <- vector("list")  # individual means for each subject, by condition
  ga.means.cond <- vector("list")  # list of grand averaged means for each condition
  for (i in 1:num.conditions) {
    means.cond[[i]] <- subset(avgsub, Stimulus == trial.types[i], select = Means)
  }
  for (j in 1:num.conditions) {
    means.cond.sub[[j]] <- matrix(unlist(means.cond[j]), ncol = num.subs)
  }
  for (h in 1:num.conditions) {
    ga.means.cond[[h]] <- rowMeans(as.data.frame(means.cond.sub[h]))
  }
  grand.avg <- data.frame(Time.range, Stimulus, unlist(ga.means.cond))
  new.names <- c("Time", "Stimulus", "Means")
  colnames(grand.avg) = new.names
  return(grand.avg)
}  # close grand.average

# --------------------------------------------------------------------------------

########################################
# group of functions for peak measures #
########################################

# function that gets grand average peak amplitudes for each condition
.get.ga.pamps = function (z, grand.avg, Stimulus, Time.range, win1, win2, num.pts) {
  with(grand.avg, { # see .get.ga.mamps comment
    values = subset(grand.avg, Stimulus==z & Time.range >= win1 &
                          Time.range <= win2, select = Means)
  values = as.vector(unlist(values))
  times = subset(grand.avg, Stimulus == z, select = Time)
  rowinfo1 <- max(which(abs(times-win1) == min(abs(times-win1))))
  rowinfo2 <- max(which(abs(times-win2) == min(abs(times-win2))))
  winnew1 = rownames(times)[rowinfo1]
  winnew2 = rownames(times)[rowinfo2]
  vector = NULL
    for (i in winnew1:winnew2) {
        start <- i - num.pts
        end <- i + num.pts
        less <- i - 1
        more <- i + 1
        test1 <- grand.avg[i, 3] > grand.avg[start:less, 3]
        test2 <- grand.avg[i, 3] > grand.avg[more:end, 3]
        test3 <- grand.avg[i, 3] < grand.avg[start:less, 3]
        test4 <- grand.avg[i, 3] < grand.avg[more:end, 3]
      if (all(test1) == TRUE & all(test2) == TRUE | all(test3) == TRUE & all(test4) == TRUE) {
          vector <- c(vector, grand.avg[i, 3])
        }
    }
      if (is.null(vector)) {
        peaks.ga <- values[which.max(abs(values))]
      } else {
        vector <- unlist(vector)
        peaks.ga <- vector[which.max(abs(vector))]
      }
  }
  )
} # close main function

# function that gets the condition associated with each grand average peak measure
.get.peak.ga.cond = function(GAcondy, grand.avg, Stimulus, Time.range, win1, win2) {
  peak.ga.cond = grand.avg[which(grand.avg[,3]==GAcondy), 2]
}

# function that gets peak latency measures
.get.peak.ga.latency = function(GAlate, grand.avg, Stimulus, Time.range, win1, win2) {
  peak.ga.lat = grand.avg[which(grand.avg[,3]==GAlate), 1]
}

# function that gets peak amplitudes for each subject, for each condition
.get.peak.amps = function(x, y, avgsub, Stimulus, Time.range, win1, win2, num.pts) {
  with(avgsub, { # see .get.ga.mamps comment
    values = subset(avgsub, Subject == x & Stimulus == y & Time.range >= win1 &
                         Time.range <= win2, select=Means)
    all.values <- subset(avgsub, Subject == x & Stimulus == y)
    rownames(all.values) <- NULL # this resets the row indexes each time the program loops, else
      # the which command refers to unexpected rows
    values = as.vector(unlist(values))
    times = subset(all.values, Stimulus == y, select = Time)
    rowinfo1 <- max(which(abs(times-win1) == min(abs(times-win1))))
    rowinfo2 <- max(which(abs(times-win2) == min(abs(times-win2))))
    winnew1 = rownames(times)[rowinfo1]
    winnew2 = rownames(times)[rowinfo2]
    vector = NULL
    for (i in winnew1:winnew2) {
      start <- i - num.pts
      end <- i + num.pts
      less <- i - 1
      more <- i + 1
      test1 <- all.values[i, 4] > all.values[start:less, 4]
      test2 <- all.values[i, 4] > all.values[more:end, 4]
      test3 <- all.values[i, 4] < all.values[start:less, 4]
      test4 <- all.values[i, 4] < all.values[more:end, 4]
      if (all(test1) == TRUE & all(test2) == TRUE | all(test3) == TRUE & all(test4) == TRUE) {
        vector <- c(vector, all.values[i, 4])
      }
    }
    if (is.null(vector)) {
      peaks <- values[which.max(abs(values))]
    } else {
      vector <- unlist(vector)
      peaks <- vector[which.max(abs(vector))]
    }
  }
  )
} # close main function

# function that gets subject IDs for peak measures
.get.peak.sub = function(subs, avgsub, Stimulus, Time.range, win1, win2) {
  peak.sub = avgsub[which(avgsub[,4]==subs), 1]
}

# function that gets the condition associated with each peak measure
.get.peak.cond = function(condy, avgsub, Stimulus, Time.range, win1, win2) {
  peak.cond = avgsub[which(avgsub[,4]==condy), 2]
}

# function that gets peak latency measures
.get.peak.latency = function(late, avgsub, Stimulus, Time.range, win1, win2) {
  peak.lat = avgsub[which(avgsub[,4]==late), 3]
}

# --------------------------------------------------------------------------------

##################################################
# group of functions for mean amplitude measures #
##################################################

# function that gets grand average mean amplitudes for each condition
.get.ga.mamps <- function(z, grand.avg, Stimulus, Time.range,
                          win1, win2) {
  with(grand.avg, { # added to appease R CMD Note: no visibile binding

  means.ga <- colMeans(subset(grand.avg,
                              Stimulus == z & # Stimulus
                                Time.range >= win1 &
                                Time.range <= win2,
                              select = Means)) # Means
  } # close with
  ) # close with
} # close function

# function that gets grand average standard deviations for each condition
.get.ga.msds <- function(z, grand.avg, Stimulus, Time.range,
                         win1, win2) {
  with(grand.avg, { # see .get.ga.mamps comment
  sd.l <- subset(grand.avg,
                 Stimulus == z &
                   Time.range >= win1 &
                   Time.range <= win2,
                 select = Means)
  sd.ga <- sd(unlist(sd.l))
  }
  )
}

# function that gets mean amplitudes for each subject for each condition
.get.mean.amps <- function(x, y, avgsub, Stimulus, Time.range, win1, win2) {
  with(avgsub, { # see .get.ga.mamps comment
  means.ind <- colMeans(subset(avgsub,
                               Subject == x &
                                 Stimulus == y &
                                 Time.range >= win1 &
                                 Time.range <= win2,
                               select = Means))
  }
  )
}

# function that gets individual standard deviations for each condition
.get.mean.msds <- function(x, y, avgsub, Stimulus, Time.range, win1, win2) {
  with(avgsub, { # see .get.ga.mamps comment
  sd.mean <- sapply(subset(avgsub,
                           Subject == x &
                             Stimulus == y &
                             Time.range >= win1 &
                             Time.range <= win2,
                           select = Means),
                    sd)
  }
  )
}
