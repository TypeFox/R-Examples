
## generate single behavior stream ####

r_behavior_stream_single <- function(mu, lambda, F_event, F_interim, stream_length,
                                     equilibrium, p0, tuning) {
  # initial condition
  start_state <- rbinom(1, 1, p0)
  
  # draw initial time
  if (equilibrium) {
    if (start_state) {
      b_stream <- F_event$r_eq(1, mu)
    } else {
      b_stream <- F_interim$r_eq(1, lambda)
    } 
  } else {
    if (start_state) {
      b_stream <- F_event$r_gen(1, mu)
    } else {
      b_stream <- F_interim$r_gen(1, lambda)
    }    
    
  }
  
  cum_length <- b_stream 
  cum_size <- 1
  
  # add event durations and interim times until total length exceeds stream length
  while (cum_length < stream_length) {
    
    # generate random event durations and interim times
    extend_size <- ceiling(tuning * (stream_length - cum_length) / (mu + lambda))
    event_times <- F_event$r_gen(n=extend_size, mean = mu)
    interim_times <- F_interim$r_gen(n=extend_size, mean = lambda)
    
    # lengthen behavior stream vector
    b_stream <- append(b_stream, cum_length + cumsum(
      if (start_state) c(rbind(interim_times, event_times)) 
      else c(rbind(event_times, interim_times))))
    
    # update totals
    cum_size <- cum_size + 2 * extend_size
    cum_length <- b_stream[cum_size]
  }
  
  list(start_state=start_state, b_stream = b_stream[b_stream < stream_length])
}



## generate behavior streams ####

#' @title Generates random behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution.
#' 
#' @param n number of behavior streams to generate
#' @param mu vector of mean event durations
#' @param lambda vector of mean interim time
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param stream_length length of behavior stream
#' @param equilibrium logical; if \code{TRUE}, then equilibrium initial conditions are used; 
#' if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
#' distributions are used for event durations and interim times.
#' @param p0 vector of initial state probabilities. Only used if \code{equilibrium = FALSE}, in which case
#' default is zero (i.e., behavior stream always starts with an interim time).
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time .
#' 
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. The vectors \code{mu}, \code{lambda}, and \code{p0} are 
#' recycled to length \code{n}.
#' @export
#' 
#' @return An object of class \code{behavior_stream} containing two elements.
#' 
#' @examples
#' # default equilibrium initial conditions
#' r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                   F_event = F_exp(), F_interim = F_exp(), 
#'                   stream_length = 100)
#'                   
#' # non-equilibrium initial conditions
#' r_behavior_stream(n = 5, mu = 3, lambda = 10,
#'                   F_event = F_gam(3), F_interim = F_gam(3),
#'                   stream_length = 100, 
#'                   equilibrium = FALSE, p0 = 0.5)

r_behavior_stream <- function(n, mu, lambda, F_event, F_interim, stream_length, 
                              equilibrium = TRUE, p0 = 0, tuning = 2) {
  
  mu_vec <- rep(mu, length.out = n)
  lambda_vec <- rep(lambda, length.out = n)
  p0_vec <- if (equilibrium) mu_vec / (mu_vec + lambda_vec) else rep(p0, length.out = n)
  
  BS <- list(stream_length = stream_length,
       b_streams = mapply(r_behavior_stream_single, 
                          mu = mu_vec,
                          lambda = lambda_vec,
                          p0 = p0_vec,
                          MoreArgs = list(F_event = F_event, F_interim = F_interim, 
                                      stream_length = stream_length, 
                                      equilibrium = equilibrium,
                                      tuning = tuning),
                          SIMPLIFY = FALSE))
  class(BS) <- "behavior_stream"
  return(BS)
}


#' @title Generates random partial interval recording behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' which are then coded as partial interval recording data with given interval length
#' and rest length.
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param stream_length length of behavior stream
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param interval_length total interval length
#' @param rest_length length of any recording time in each interval
#' @param summarize logical value indicating whether the behavior streams should by summarized by taking their mean
#' @param equilibrium logical; if \code{TRUE}, then equilibrium initial conditions are used; 
#' if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
#' distributions are used for event durations and interim times.
#' @param p0 Initial state probability. Only used if \code{equilibrium = FALSE}, in which case
#' default is zero (i.e., behavior stream always starts with an interim time).
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time .
#' 
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. Then applies a partial interval recording filter to the generated behavior streams.
#' 
#' @return If \code{summarize = FALSE}, a matrix with rows equal to \code{n} and a number of columns equal to the number intervals per session. If \code{summarize = TRUE} a vector of means of length \code{n}.
#' @export
#' 
#' @examples
#' 
#' # An unsummarized set of PIR observations
#' r_PIR(n = 5, mu = 2, lambda = 4, stream_length = 20, 
#'        F_event = F_exp(), F_interim = F_exp(), 
#'        interval_length = 1, rest_length = 0)
#'       
#' # A summarized set of of PIR observations
#' r_PIR(n = 5, mu = 2, lambda = 4, stream_length = 20, 
#'        F_event = F_exp(), F_interim = F_exp(), 
#'        interval_length = 1, rest_length = 0,
#'        summarize = TRUE)
#'        
#' @author Daniel Swan <dswan@@utexas.edu>

r_PIR <- function(n, mu, lambda, stream_length, F_event, F_interim, 
                  interval_length, rest_length = 0, summarize = FALSE, 
                  equilibrium = TRUE, p0 = 0, tuning = 2){
  
  if (equilibrium) p0 <- mu / (mu + lambda)
  
  n_intervals <- floor(stream_length / interval_length)
  start_time <- interval_length * (0:(n_intervals - 1))
  end_time <- start_time + interval_length - rest_length
  
  samples <- replicate(n, {
    BS <- r_behavior_stream_single(mu = mu, lambda = lambda, 
                            F_event = F_event, F_interim = F_interim, 
                            stream_length = stream_length, 
                            equilibrium = equilibrium,
                            p0 = p0, tuning = tuning)    
    IntRec_single(b_stream = BS, start_time = start_time, end_time = end_time)
  })
  
  if (summarize) colMeans(samples) else t(samples)
}
  
#' @title Generates random whole interval recording behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' which are then coded as whole interval recording data with given interval length
#' and rest length.
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param stream_length length of behavior stream
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param interval_length total interval length
#' @param rest_length length of any recording time in each interval
#' @param summarize logical value indicating whether the behavior streams should by summarized by taking their mean
#' @param equilibrium logical; if \code{TRUE}, then equilibrium initial conditions are used; 
#' if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
#' distributions are used for event durations and interim times.
#' @param p0 Initial state probability. Only used if \code{equilibrium = FALSE}, in which case
#' default is zero (i.e., behavior stream always starts with an interim time).
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time .
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. Then applies a whole interval recording filter to the generated behavior streams.
#' @return If \code{summarize = FALSE}, a matrix with rows equal to \code{n} and a number of columns equal to the number intervals per session. If \code{summarize = TRUE} a vector of means of length \code{n}.
#' @export
#' 
#' @examples
#' 
#' # An unsummarized set of WIR observations
#' r_WIR(n = 5, mu = 2, lambda = 4, stream_length = 20, 
#'        F_event = F_exp(), F_interim = F_exp(), 
#'        interval_length = 1, rest_length = 0)
#'       
#' # A summarized set of of WIR observations
#' r_WIR(n = 5, mu = 2, lambda = 4, stream_length = 20, 
#'        F_event = F_exp(), F_interim = F_exp(), 
#'        interval_length = 1, rest_length = 0,
#'        summarize = TRUE)
#'        
#' @author Daniel Swan <dswan@@utexas.edu>

r_WIR <- function(n, mu, lambda, stream_length, F_event, F_interim, 
                  interval_length, rest_length = 0, summarize = FALSE, 
                  equilibrium = TRUE, p0 = 0, tuning = 2){
  
  if (equilibrium) p0 <- mu / (mu + lambda)
  
  n_intervals <- floor(stream_length / interval_length)
  start_time <- interval_length * (0:(n_intervals - 1))
  end_time <- start_time + interval_length - rest_length
  
  samples <- replicate(n, {
    BS <- r_behavior_stream_single(mu = mu, lambda = lambda, 
                                   F_event = F_event, F_interim = F_interim, 
                                   stream_length = stream_length, 
                                   equilibrium = equilibrium,
                                   p0 = p0, tuning = tuning)
    IntRec_single(b_stream = BS, start_time = start_time, end_time = end_time, partial = FALSE)
  })
  
  if (summarize) colMeans(samples) else t(samples)
}

#' @title Generates random momentary time sampling behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' which are then coded as momentary time sampling data with given interval length
#' between moments.
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param stream_length length of behavior stream
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param interval_length length of time between moments
#' @param summarize logical value indicating whether the vector of moments should be summarized by taking their mean, excluding the first moment in each row.
#' @param equilibrium logical; if \code{TRUE}, then equilibrium initial conditions are used; 
#' if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
#' distributions are used for event durations and interim times.
#' @param p0 Initial state probability. Only used if \code{equilibrium = FALSE}, in which case
#' default is zero (i.e., behavior stream always starts with an interim time).
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time.
#' 
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. Then applies a momentary time sampling filter to the generated behavior streams.
#' 
#' @return If \code{summarize = FALSE}, a matrix of logicals with rows equal to \code{n} and length equal to \code{(stream_length/interval_length) + 1}. If \code{summarize = TRUE}, a vector of means of length \code{n}.
#' @export
#' 
#' @examples
#' 
#' # A set of unsummarized MTS observations
#' r_MTS(n = 5, mu = 2, lambda = 4, stream_length = 20, 
#'        F_event = F_exp(), F_interim = F_exp(), interval_length = 1)
#'       
#' # A set of summarized MTS observations
#' r_MTS(n = 5, mu = 2, lambda = 4, stream_length = 20, 
#'        F_event = F_exp(), F_interim = F_exp(), 
#'        interval_length = 1, summarize = TRUE)
#'        
#' @author Daniel Swan <dswan@@utexas.edu>

r_MTS <- function(n, mu, lambda, stream_length, F_event, F_interim, 
                  interval_length, summarize = FALSE, equilibrium = TRUE, 
                  p0 = 0, tuning = 2) {
  
  if (equilibrium) p0 <- mu / (mu + lambda)
  
  moments <- seq(interval_length * summarize, stream_length, interval_length)
  samples <- replicate(n, {
    BS <- r_behavior_stream_single(mu = mu, lambda = lambda, 
                             F_event = F_event, F_interim = F_interim, 
                             stream_length = stream_length, 
                             equilibrium = equilibrium,
                             p0 = p0, tuning = tuning)
    MTS_single(b_stream = BS, moments = moments)
    })
  
  if(summarize) colMeans(samples) else t(samples)
}

#' @title Generates random samples of continuously recorded behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' summarized as the total proportion of time the behavior of interest occurred.
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param stream_length length of behavior stream
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param equilibrium logical; if \code{TRUE}, then equilibrium initial conditions are used; 
#' if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
#' distributions are used for event durations and interim times.
#' @param p0 Initial state probability. Only used if \code{equilibrium = FALSE}, in which case
#' default is zero (i.e., behavior stream always starts with an interim time).
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time .
#' 
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. Then applies a continuous recording filter to the generated behavior streams.
#' 
#' @return A vector of proportions of length \code{n}.
#' @export
#' 
#' @examples
#' 
#' r_continuous_recording(n = 5, mu = 2, lambda = 4, stream_length = 20,
#'                        F_event = F_exp(), F_interim = F_exp())
#' 
#' @author Daniel Swan <dswan@@utexas.edu>
                        
r_continuous_recording <- function(n, mu, lambda, stream_length, F_event, F_interim, 
                                   equilibrium = TRUE, p0 = 0, tuning = 2) {
  
  if (equilibrium) p0 <- mu / (mu + lambda)
  
  samples <- replicate(n, {
    BS <- r_behavior_stream_single(mu = mu, lambda = lambda, 
                                   F_event = F_event, F_interim = F_interim, 
                                   stream_length = stream_length, 
                                   equilibrium = equilibrium,
                                   p0 = p0, tuning = tuning)
    CDR_single(b_stream = BS, stream_length = stream_length)
  })
  
  samples
}

#' @title Generates random samples of event counts
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' summarized as the the total number of behaviors that began during the recording
#' session
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param stream_length length of behavior stream
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param equilibrium logical; if \code{TRUE}, then equilibrium initial conditions are used; 
#' if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
#' distributions are used for event durations and interim times.
#' @param p0 Initial state probability. Only used if \code{equilibrium = FALSE}, in which case
#' default is zero (i.e., behavior stream always starts with an interim time).
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time .
#' 
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. Then applies an event counting filter to the generated behavior streams.
#' 
#' @return A vector of behavior counts of length \code{n}.
#' @export
#' 
#' @examples
#' 
#' r_event_counting(n = 5, mu = 2, lambda = 4, stream_length = 20,
#'                      F_event = F_exp(), F_interim = F_exp())
#'                      
#' @author Daniel Swan <dswan@@utexas.edu>

r_event_counting <- function(n, mu, lambda, stream_length, F_event, F_interim, 
                                 equilibrium = TRUE, p0 = 0, tuning = 2) {
  
  if (equilibrium) p0 <- mu / (mu + lambda)
  
  samples <- replicate(n,{
    BS <- r_behavior_stream_single(mu = mu, lambda = lambda, 
                                   F_event = F_event, F_interim = F_interim, 
                                   stream_length = stream_length, 
                                   equilibrium = equilibrium,
                                   p0 = p0, tuning = tuning)
    floor((length(BS$b_stream) + 1 - BS$start_state)/2)
  })
  
  samples
}
