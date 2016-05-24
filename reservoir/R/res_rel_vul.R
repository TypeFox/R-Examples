#' @title Reliability, resilience, and vulnerability analysis for water supply reservoirs
#' @description Computes time-based, annual, and volumetric reliability, as well as resilience and dimensionless vulnerability for a single reservoir.
#' @param Q                  vector or time series object. Net inflow totals to the reservoir. Recommended units: Mm^3 (Million cubic meters).
#' @param target             numerical. The target release constant. Recommended units: Mm^3 (Million cubic meters).
#' @param capacity           numerical. The reservoir capacity. Should be same volumetric unit as Q and R.
#' @param surface_area       numerical. The reservoir water surface area at maximum capacity. Recommended units: km^2 (square kilometers).
#' @param max_depth          numerical. The maximum water depth of the reservoir at maximum capacity. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only. Recommended units: meters.
#' @param evap               vector or time series object of length Q, or a numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param double_cycle       logical. If TRUE the input series will be replicated and placed end-to-end to double the simulation. (Recommended if the critical period occurs at the end of the recorded inflow time series) 
#' @param plot               logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param S_initial          numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param policy             list. The output of the SDP function. If omitted, Standard Operating Policy is assumed.
#' @return Returns reliability, resilience and vulnerability metrics based on supply deficits.
#' @examples # Compare reliability, resilience and vulnerability for two operating policies (SOP and SDP).
#' rrv(resX$Q_Mm3, capacity = 20*resX$cap_Mm3, target = 0.95 * mean(resX$Q_Mm3))
#' pol_Markov <- sdp_supply(resX$Q_Mm3, capacity = 20 * resX$cap_Mm3,
#' target = 0.95 * mean(resX$Q_Mm3), Markov = TRUE)
#' rrv(resX$Q_Mm3, capacity = 20*resX$cap_Mm3, target = 0.95 * mean(resX$Q_Mm3), policy = pol_Markov)
#' @import stats
#' @importFrom graphics abline lines
#' @export
rrv <- function(Q, target, capacity, double_cycle = FALSE,
                surface_area, max_depth, evap,
                plot = TRUE, S_initial = 1, policy) {
    
    frq <- frequency(Q)
    if(length(target)==1){
      target <- rep(target, length(Q))
    }
    if(length(target) != length(Q))
      stop("target must be numerical constant or time series of length Q")
    
    target <- ts(target, start = start(Q), frequency = frequency(Q))
    
    
    x <- simRes(Q = Q, target = target, capacity = capacity, double_cycle = double_cycle, policy = policy,
                surface_area = surface_area, max_depth = max_depth, evap = evap, S_initial = S_initial, plot = FALSE)
    R <- x$releases
    S <- x$storage
    Spill <- x$spill

    #===============================================================================
    
    # COMPUTE METRICS FROM SIMULATION RESULTS---------------------------------------
    
    deficit <- ts(round(1 - (R / target),5), start = start(Q), frequency = frq)
    rel_ann <- sum(aggregate(deficit, FUN = mean) == 0) /
      length(aggregate(deficit, FUN = mean))
    rel_time <- sum(deficit == 0) / length(deficit)
    rel_vol <- sum(R) / sum(target)
    fail.periods <- which(deficit > 0)
    if (length(fail.periods) == 0) {
        resilience <- NA
        vulnerability <- NA
    } else {
        if (length(fail.periods) == 1) {
            resilience <- 1
            vulnerability <- max(deficit)
        } else {
            resilience <- (sum(diff(which(deficit > 0)) > 1) + 1) / (length(which(deficit > 0)))
            fail.refs <- vector("numeric", length = length(fail.periods))
            fail.refs[1] <- 1
            for (j in 2:length(fail.periods)) {
                if (fail.periods[j] > (fail.periods[j - 1] + 1)) {
                  fail.refs[j] <- fail.refs[j - 1] + 1
                } else {
                  fail.refs[j] <- fail.refs[j - 1]
                }
            }
            n.events <- max(fail.refs)
            event.starts <- by(fail.periods, fail.refs, FUN = min)
            event.ends <- by(fail.periods, fail.refs, FUN = max)
            max.deficits <- vector("numeric", length = n.events)
            for (k in 1:n.events) {
                max.deficits[k] <- max(deficit[event.starts[k]:event.ends[k]])
            }
            vulnerability <- mean(max.deficits)
        }
    }

    #===============================================================================
  
    results <- list(rel_ann, rel_time, rel_vol, resilience, vulnerability,
                    S, R, Spill)
    names(results) <- c("annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability", "storage", "releases", "spill")
    if (plot) {
        plot(R, ylab = "Release", ylim = c(0, max(target)), main = paste0("Ann rel. = ", round(rel_ann, 2),
                                                                         "; Time rel. = ", round(rel_time, 2),
                                                                         "; Vol rel. = ", round(rel_vol, 2),
                                                                         "; Resil. = ", round(resilience, 2),
                                                                         "; Vuln. = ", round(vulnerability, 2))); lines(target, lty=3)
        plot(S, ylab = "Storage", ylim = c(0, capacity))
        plot(Spill, ylab = "Uncontrolled Spill")
    }
    
    return(results)
}
