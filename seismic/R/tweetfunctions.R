#' Predicting information cascade by self-exciting point process model
#'
#' This package implements a self-exciting point process model for information cascades. 
#' An information cascade occurs when many people engage in the same acts after observing 
#' the actions of others. Typical examples are post/photo resharings on Facebook and retweets 
#' on Twitter. The package provides functions to estimate the infectiousness of an 
#' information cascade and predict its popularity given the observed history. 
#' For more information, see 
#' \url{http://snap.stanford.edu/seismic/}.
#'
#' @docType package
#' @name seismic
#' @references SEISMIC: A Self-Exciting Point Process Model for Predicting Tweet Popularity by Q. Zhao, M. Erdogdu, H. He, A. Rajaraman, J. Leskovec, ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD), 2015.
NULL

#' Memory kernel
#' 
#' Probability density function and complementary cumulative distribution function
#' for the human reaction time.
#' @keywords internal
#' 
#' @param t time
#' @param theta exponent of the power law
#' @param cutoff the cutoff value where the density changes from constant to power law
#' @param c the constant density when t is less than the cutoff
#' @return the density at t
#' @details default values are measured from a real Twitter data set.
#' @return \code{memory.pdf} returns the density function at t. 
#' \code{memory.ccdf} returns the ccdf (probabilty of greater than t).
#' 
memory.pdf <- function(t, theta=0.2314843, cutoff=300, c=0.0006265725) {
  if (t < cutoff)
    return(c)
  else
    return(c*exp((log(t) - log(cutoff))*(-(1+theta))))
}

#' @describeIn memory.pdf
memory.ccdf <- function(t, theta=0.2314843, cutoff=300, c=0.0006265725) {
  t[t<0] <- 0
  index1 <- which(t <= cutoff)
  index2 <- which(t > cutoff)
  ccdf <- rep(0, length(t))
  ccdf[index1] <- 1 - c*t[index1]
  ccdf[index2] <-  c*cutoff^(1+theta)/theta*(t[index2]^(-theta))
  ccdf
}

#' Integration with respect to locally weighted kernel
#' 
#' @keywords internal
#' 
#' @param t1 a vector of integral lower limit
#' @param t2 a vector of integral upper limit
#' @param ptime the time (a scalar) to estimate infectiousness and predict for popularity
#' @param slope slope of the linear kernel
#' @param window size of the linear kernel
#' @inheritParams memory.pdf
#' @inheritParams get.infectiousness
#' @return \code{linear.kernel} returns the integral from vector t1 to vector t2 of 
#' c*[slope(t-ptime) + 1]; 
#' \code{power.kernel} returns the integral from vector t1 to vector 2 of c*((t-share.time)/cutoff)^(-(1+theta))[slope(t-ptime) + 1]; 
#' \code{integral.memory.kernel} returns the vector with ith entry being integral_-inf^inf phi_share.time[i]*kernel(t-p.time)
#' @seealso \code{\link{memory.pdf}}
linear.kernel <- function(t1, t2, ptime, slope, c=0.0006265725){
  ## indefinite integral is c*(t-ptime*slope*t+(slope*t^2)/2)
  return(c*(t2-ptime*slope*t2+slope*t2^2/2) - c*(t1-ptime*slope*t1+slope*t1^2/2))
}

#' @describeIn linear.kernel
power.kernel <- function(t1, t2, ptime, share.time, slope, theta=0.2314843, cutoff=300, c=0.0006265725){
  return(c*cutoff^(1+theta)*(t2-share.time)^(-theta)*(share.time*slope-theta+(theta-1)*ptime*slope-theta*slope*t2+1)/((theta-1)*theta) - c*cutoff^(1+theta)*(t1-share.time)^(-theta)*(share.time*slope-theta+(theta-1)*ptime*slope-theta*slope*t1+1)/((theta-1)*theta))
}

#' @describeIn linear.kernel
integral.memory.kernel <- function(p.time, share.time, slope, window, theta=0.2314843, cutoff=300, c=0.0006265725){
  index1 <- which(p.time <= share.time)
  index2 <- which(p.time > share.time & p.time <= share.time + cutoff)
  index3 <- which(p.time > share.time + cutoff & p.time <= share.time + window)
  index4 <- which(p.time > share.time + window & p.time <= share.time + window + cutoff)
  index5 <- which(p.time > share.time + window + cutoff)
  integral <- rep(NA, length(share.time))
  integral[index1] <- 0
  integral[index2] <- linear.kernel(share.time[index2], p.time, p.time, slope)
  integral[index3] <- linear.kernel(share.time[index3], share.time[index3] + cutoff, p.time, slope) + power.kernel(share.time[index3]+cutoff, p.time, p.time, share.time[index3], slope)
  integral[index4] <- linear.kernel(p.time-window, share.time[index4]+cutoff, p.time, slope) + power.kernel(share.time[index4]+cutoff, p.time, p.time, share.time[index4], slope)
  integral[index5] <- power.kernel(p.time-window, p.time, p.time, share.time[index5], slope)
  return(integral)
}

#' Estimate the infectiousness of an information cascade
#' 
#' @param share.time observed resharing times, sorted, share.time[1] =0
#' @param degree observed node degrees
#' @param p.time equally spaced vector of time to estimate the infectiousness, p.time[1]=0
#' @param max.window maximum span of the locally weight kernel 
#' @param min.window minimum span of the locally weight kernel
#' @param min.count the minimum number of resharings included in the window
#' @details Use a triangular kernel with shape changing over time. At time p.time, use a triangluer kernel with slope = min(max(1/(\code{p.time}/2), 1/\code{min.window}), \code{max.window}).
#' @return a list of three vectors: \itemize{
#' \item infectiousness. the estimated infectiousness 
#' \item p.up. the upper 95 percent approximate confidence interval
#' \item p.low. the lower 95 percent approximate confidence interval
#' }
#' @export
#' @examples
#' data(tweet)
#' pred.time <- seq(0, 6 * 60 * 60, by = 60)
#' infectiousness <- get.infectiousness(tweet[, 1], tweet[, 2], pred.time)
#' plot(pred.time, infectiousness$infectiousness)
get.infectiousness <- function(share.time,
                      degree,
                      p.time, 
                      max.window = 2 * 60 * 60, 
                      min.window = 300, 
                      min.count = 5) {
  
  ix <- sort(share.time, index.return=TRUE)$ix
  share.time <- share.time[ix]
  
  slopes <- 1/(p.time/2)
  slopes[slopes < 1/max.window] <- 1/max.window
  slopes[slopes > 1/min.window] <- 1/min.window
  
  windows <- p.time/2
  windows[windows > max.window] <- max.window
  windows[windows < min.window] <- min.window
  
  for(j in c(1:length(p.time))) {
    ind <- which(share.time >= p.time[j] - windows[j] & share.time < p.time[j])
    if(length(ind) < min.count) {
      ind2 <- which(share.time < p.time[j])
      lcv <- length(ind2)
      ind <- ind2[max((lcv-min.count),1):lcv]
      slopes[j] <- 1/(p.time[j] - share.time[ind[1]])
      windows[j] <- p.time[j] - share.time[ind[1]]
    }
  }
  
  M.I <- matrix(0,nrow=length(share.time),ncol=length(p.time))
  for(j in 1:length(p.time)){
    M.I[,j] <- degree*integral.memory.kernel(p.time[j], share.time, slopes[j], windows[j])
  }
  infectiousness.seq <- rep(0, length(p.time))
  p.low.seq <- rep(0, length(p.time))
  p.up.seq <- rep(0, length(p.time))
  share.time <- share.time[-1]          #removes the original tweet from retweet
  for(j in c(1:length(p.time))) {
    share.time.tri <- share.time[which(share.time >= p.time[j] - windows[j] & share.time < p.time[j])]
    rt.count.weighted <- sum(slopes[j]*(share.time.tri - p.time[j]) + 1)
    #print(paste("p.time[i]", p.time[j], "rt.num", length(share.time.tri)))
    I <- sum(M.I[,j])
    rt.num <- length(share.time.tri)
    if (rt.count.weighted==0)
      next
    else {
      infectiousness.seq[j] <- (rt.count.weighted)/I
      p.low.seq[j] <- infectiousness.seq[j] * qchisq(0.05, 2*rt.num) / (2*rt.num)
      p.up.seq[j] <- infectiousness.seq[j] * qchisq(0.95, 2*rt.num) / (2*rt.num)
    }
  }
  ## p.low.seq[is.nan(p.low.seq)] <- 0
  ## p.up.seq[is.nan(p.up.seq)] <- 0
  list(infectiousness = infectiousness.seq, p.up = p.up.seq, p.low = p.low.seq)
}

#' Predict the popularity of information cascade
#' 
#' @param infectiousness a vector of estimated infectiousness, returned by \code{\link{get.infectiousness}}
#' @param n.star the average node degree in the social network
#' @param features.return if TRUE, returns a matrix of features to be used to further calibrate the prediction
#' @inheritParams get.infectiousness
#' @return a vector of predicted populatiry at each time in \code{p.time}.
#' @export
#' @examples
#' data(tweet)
#' pred.time <- seq(0, 6 * 60 * 60, by = 60)
#' infectiousness <- get.infectiousness(tweet[, 1], tweet[, 2], pred.time)
#' pred <- pred.cascade(pred.time, infectiousness$infectiousness, tweet[, 1], tweet[, 2], n.star = 100)
#' plot(pred.time, pred)
pred.cascade <- function(p.time, infectiousness, share.time, degree, n.star=100, features.return = FALSE){
  
  # n.star should a vector of the same length as p.time
  if (length(n.star) == 1) {
    n.star <- rep(n.star, length(p.time))
  }

  # to train for best n.star, we get feature matrices
  features <- matrix(0, length(p.time), 3)
  
  prediction <- matrix(0, length(p.time), 1)
  for (i in 1:length(p.time)) {
    share.time.now <- share.time[share.time <= p.time[i]]
    nf.now <- degree[share.time <= p.time[i]]
    rt0 <- sum(share.time <= p.time[i]) - 1
    rt1 <- sum(nf.now * infectiousness[i] * memory.ccdf(p.time[i] - share.time.now))
    prediction[i] <- rt0 + rt1 / (1 - infectiousness[i]*n.star[i])
    features[i, ] <- c(rt0, rt1, infectiousness[i])
    if (infectiousness[i] > 1/n.star[i]) {
        prediction[i] <- Inf
    }
  }
  
  colnames(features) <- c("current.rt", "numerator", "infectiousness")
  
  if (!features.return) {
    prediction
  } else {
    list(prediction = prediction, features = features)
  }
}

#' An example information cascade
#' 
#' A dataset containing all the (relative) resharing time and node degree of a tweet. The original Twitter ID is 127001313513967616.
#' 
#' \itemize{
#'   \item relative_time_second. resharing time in seconds
#'   \item number_of_followers. number of followers
#' }
#' 
#' @format A data frame with 15563 rows and 2 columns
#' @source \url{http://board.muse.mu/archive/index.php/t-85075.html}
#' @name tweet
NULL