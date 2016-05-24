# auxilliary functions

#fast p-value computation for simple marginal lm fit test
pValComp <- function(x,y,n,suma){
  a <- lm.fit(cbind(x,1),y)
  b <- sum(a$residuals^2)
  1-pf((suma-b)/b*n,1, n)
}


# Function to replace missing values with mean for that col
replace_na_with_mean <- function(x) {
  x_bar <- mean(x, na.rm = TRUE)
  ifelse(is.na(x), x_bar, x)
}

#estimate noise in lm
estimate_noise <- function (X, y, intercept = TRUE) {
  n = nrow(X)
  if (intercept)
    X = cbind(rep(1, n), X)
  p = ncol(X)
  fit = lm.fit(X, y)
  sqrt(sum(fit$residuals^2)/(n - p))
}

#create clumping plot.data
create_clumping_plot_data <- function(x){
  plot.data <- NULL
  for(i in 1L:length(x$SNPclumps)){
    plot.data <- rbind(plot.data,
                       cbind(as.numeric(x$X_info[x$selectedSnpsNumbersScreening[x$SNPclumps[[i]]],1]),
                             as.numeric(x$X_info[x$selectedSnpsNumbersScreening[x$SNPclumps[[i]]],3]),
                             i, -log(x$pVals[x$selectedSnpsNumbersScreening[x$SNPclumps[[i]]]])))
  }
  rownames(plot.data) <- NULL
  plot.data <- data.frame(plot.data)
  colnames(plot.data) <- c("chromosome", "snp", "clump", "val")
  plot.data <- cbind(plot.data,
                     representatives = unlist(x$SNPclumps) %in% unlist(x$SNPnumber))
  granice <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
  granice_max <- cumsum(granice$x)
  granice$x <- c(0,head(cumsum(granice$x),-1))
  for(i in unique(plot.data$chromosome)){
    plot.data$snp[plot.data$chromosome==i] <- granice$x[i] +
      plot.data$snp[plot.data$chromosome==i]
  }
  plot.data$val[is.infinite(plot.data$val)] <- -log(2e-16) #R precision
  return(plot.data)
}

#create slopeResult plot.data
create_slope_plot_data <- function(x){
  plot.data <- NULL
  for(i in 1L:length(x$selectedClumps)){
    plot.data <- rbind(plot.data,
                       cbind(as.numeric(x$X_info[x$screenedSNPsNumbers[x$selectedClumps[[i]]],1]),
                             as.numeric(x$X_info[x$screenedSNPsNumbers[x$selectedClumps[[i]]],3]),
                             i, x$effects[i]^2/var(as.vector(x$y))))
  }
  rownames(plot.data) <- NULL
  plot.data <- data.frame(plot.data)
  colnames(plot.data) <- c("chromosome", "snp", "clump", "val")
  plot.data <- cbind(plot.data,
                     representatives = unlist(x$selectedClumps) %in% unlist(x$selectedSNPs))
  granice <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
  granice_max <- cumsum(granice$x)
  granice$x <- c(0,head(cumsum(granice$x),-1))
  for(i in unique(plot.data$chromosome)){
    plot.data$snp[plot.data$chromosome==i] <- granice$x[i] +
      plot.data$snp[plot.data$chromosome==i]
  }
  plot.data$val[plot.data$representatives] <- (x$effects^2/var(x$y))
  plot.data
}
