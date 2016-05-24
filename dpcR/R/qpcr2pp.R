#' qPCR to Poisson Process
#' 
#' Describes qPCR as Poisson process.
#' 
#' @details 
#' Selected platforms (e.g., Open Array) are real-time platforms. dPCR can be
#' described by Poisson statistics. The function \code{qpcr2pp} takes a step
#' further and interprets the dPCR as a Poisson process if it is analyzed as a
#' "time" based process.
#' 
#' The dPCR Technology breaks fundamentally with the previous concept of
#' nucleic acid quantification. dPCR can be seen as a next generation nucleic
#' acid quantification method based on PCR. The key difference between dPCR and
#' traditional PCR lies in the method of measuring (absolute) nucleic acids
#' amounts. This is possible after ``clonal DNA amplification'' in thousands of
#' small separated partitions (e.g., droplets, nano chambers).  Partitions with
#' no nucleic acid remain negative and the others turn positive. Selected
#' technologies (e.g., OpenArray(R) Real-Time PCR System) monitor amplification
#' reactions in the chambers in real-time. Cq values are calculated from the
#' amplification curves and converted into discrete events by means of positive
#' and negative partitions and the absolute quantification of nucleic acids is
#' done by Poisson statistics.
#' 
#' PCR data derived from a qPCR experiment can be seen as a series of events
#' over time. We define t_i as the time between the first (i - 1)^st and the
#' i^th event. Therefore, the time \eqn{S_n}{S_n} is the sum of all
#' \eqn{t_i}{t_i} from \eqn{i = 1}{i = 1} to \eqn{i = n}{i = n}. This is the
#' time to the n^th event. \eqn{S(t)}{S(t)} is the number of events in \eqn{[0,
#' t]}{[0, t]}. This can be seen as a Poisson process. The Poisson statistics
#' is the central theorem to random processes in digital PCR.
#' 
#' The function \code{qpcr2pp} is used to model random point events in time
#' units (PCR cycles), such as the increase of signal during a qPCR reaction in
#' a single compartment. A Poisson process can be used to model times at which
#' an event occurs in a "system". The \code{qpcr2pp} (quantitative Real-Time
#' PCR to Poisson process) function transforms the qPCR amplification curve
#' data to quantification points (Cq), which are visualized as Poisson process.
#' This functions helps to spot differences between replicate runs of digital
#' PCR experiments. In ideal scenarios the \code{qpcr2pp} plots are highly
#' similar.
#' 
#' This tool might help to spot differences between experiments (e.g.,
#' inhibition of amplification reactions, influence of the chip arrays). The
#' qPCR is unique because the amplification of conventional qPCRs takes place
#' in discrete steps (cycles: 1, 2 ... 45), but the specific Cq values are
#' calculated with continuous outcomes (Cq: 18.2, 25.7, ...). Other
#' amplification methods such as isothermal amplifications are time based and
#' thus better suited for Poisson process.
#' 
#' @inheritParams limit_cq
#' @param NuEvents "number of expected events" within a time frame
#' (interval).
#' @param delta difference "time (cycles) points" e.g., Cycle 18 and 25.
#' @param exper The id of experiments.
#' @param replicate The id of technical replicates.
#' @param assay The name or id of assays.
#' @param type object of class \code{"character"} defining type of data. Could
#' be \code{"np"} (status (positive (1) or negative(0)) of each droplet) or\code{"ct"} 
#' (threshold cycle).
#' @return An object of \code{\linkS4class{qdpcr}} class.
#' @author Stefan Roediger, Michal Burdukiewicz.
#' @keywords Poisson Process qPCR
#' @export qpcr2pp
#' @examples
#' 
#' library(qpcR)
#' test <- cbind(reps[1L:45, ], reps2[1L:45, 2L:ncol(reps2)], 
#' 	      reps3[1L:45, 2L:ncol(reps3)])
#' 
#' # before interpolation qPCR experiment must be converted into dPCR
#' qpcrpp <- qpcr2pp(data = test, cyc = 1, fluo = NULL, Cq_range = c(20, 30), 
#'                   model = l5, delta = 5)
#' summary(qpcrpp)


qpcr2pp <- function(data, cyc = 1, fluo = NULL,
                    Cq_range = c(min(data[cyc]) + 6, max(data[cyc]) - 6), model = l5, 
                    SDM = TRUE, NuEvents = 1, delta = 1, exper = "qPCR1", replicate = 1, 
                    assay = "Unknown", type = "np") {
  
  if(!(type %in% c("np", "ct")))
    stop("'type' must have value 'ct' or 'np'.")
  
  cq_dat <- limit_cq(data = data, cyc = cyc, fluo = fluo, 
                     Cq_range = Cq_range, model = model, SDM = SDM, pb = FALSE)
  
  res_qPCR <- cq_dat[order(cq_dat[[1]]), ]
  res_qPCR <- cbind(res_qPCR, cumsum(res_qPCR[, 2]))
  colnames(res_qPCR) <- c("Cycles", "result", "lambda") 
  
  # do not know if this is correct, WIP
  # cycle_time should give the "average time" between the occurrence of a 
  # positive reaction and another positive reaction
  dens_tmp <- dpcr_density(sum(res_qPCR[, 2]), nrow(res_qPCR), plot = FALSE)
  cycle_time <- exp(1/dens_tmp[["k"]])
  # Determine probaility how often a events occur according to Poisson process 
  # with a certia rate per time frame (interval).
  # NuEvents is "number of expected events" within a time frame (interval)
  # dens_tmp$k gives the rate of the process according to dpcr_density()
  # delta is the difference "time (cycles) points" e.g., Cycle 18 and 25
  # mu is the expected number of events in defined interval (but this is somewhat
  # stupid since the intervals are discrete ... and so on)
  # cyc_occ gives the occurrence in an interval
  mu <- dens_tmp[["k"]] * delta
  fact_NuEvents <- factorial(NuEvents)
  if (fact_NuEvents != "Inf") {
    cyc_occ <- (exp(-mu) * mu^NuEvents)/fact_NuEvents
  } else cyc_occ <- "too large"
  # END WIP
  
  res <- construct_dpcr(data = cq_dat[, ifelse(type == "np", 2, 1)], n = nrow(cq_dat), exper = exper, 
                        replicate = replicate, assay = assay, type = type)
  class(res) <- "qdpcr"
  slot(res, "qpcr") <- data.matrix(res_qPCR)
  slot(res, "mu") <- mu
  slot(res, "CT") <- cycle_time
  slot(res, "CO") <- cyc_occ
  
  res
}