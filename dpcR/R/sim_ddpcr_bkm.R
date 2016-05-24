# #' Simulate Droplet Digital PCR
# #'
# #' Simulates results of a droplet digital PCR.
# #'
# #' @param m is either the concentration of templates molecules in the raw sample 
# #' (copies/microliter)
# #' or the expected number of template molecules per droplet.
# #' See \code{mexp} parameter for details
# #' Must be a (vector of) positive integers.
# #' 
# #' @param n the expected number of droplets per experiment. Must be a positive integer.
# #' Default 20000 based on the Bio-Rad ddPCR QX100 theoretical expected values
# #' 
# #' @param mexp If \code{TRUE}, m is the expected number of template molecules per droplet
# #' If \code{FALSE}, m is the concentration of the raw sample
# #' Default \code{TRUE} as in Jacobs et al.
# #' 
# #' @param n_exp the number of experiments that are simulated by the function for each given 
# #' \code{m}.
# #' Default 8 for eight replicates for each given \code{m} as in Jacobs et al. 2014
# #' 
# #' @param type Object of class \code{"character"} defining type of data. Could
# #' be \code{"tnp"} (total number of positive wells in the panel), \code{"fluo"} 
# #' (fluorescence) or \code{"np"} (status (positive (1) or negative(0)) of each droplet).
# #' 
# #' @param fluo_range defines expected space between two 
# #' consecutive measured droplets. Used only when parameter \code{type} has value 
# #' \code{"fluo"}. Values between 10-20 give nice results.
# #' 
# #' @param sddropc standard deviation of the number of droplets generated.
# #' Must be a real number between 0 and \code{n} divided by 10.
# #' Default 0 for constant number of droplets.
# 
# #' @param mudropr average proportion (between 0 and 1) of retained partitions.
# #' Must be a real number between 0 and 1.
# #' Default 1 for no loss.
# 
# #' @param sddropr relative standard deviation of the proportion of retained partitions.
# #' Must be a real number, preferably close to 0.
# #' Default 0 for a constant loss.
# 
# #' @param Pvar If \code{TRUE}, number of copies in constant volume follows P(c) distribution.
# #' If \code{FALSE},  number of copies in constant volume is constant.
# #' Default \code{TRUE} for the realistic Poisson model.
# 
# #' @param piperr coefficient of variation of the actual pipetted volume from the raw material.
# #' Must be a positive real number, preferably close to 0 (0.1 = 10% is very large).
# #' Default 0 for constant volume equal to the expected volume
# 
# #' @param dropsd relative variability of the droplet volume.
# #' parameter sigma of a lognormal distribution with mu = 0.
# #' Must be a positive real number, preferably close to 0.
# #' Default 0 for constant droplet size.
# 
# #' @param falpos probability that a partition containing no copy gives a positive result.
# #' Must be a real number between 0 and 1.
# #' Default 0 for no false positives.
# #' Only used if \code{fluo} is \code{NULL}.
# #' 
# #' @param falneg probability that a partition containing at least one copy gives a negative 
# #' result
# #' Must be a real number between 0 and 1
# #' Default 0 for no false negatives
# #' Only used with \code{fluo} is \code{NULL}
# #' 
# #' @param rain parameter that defines how much inhibition is enforced on positive droplets.
# #' Must be a real number between 0 and 1 with 0 no rain, 1 positive droplets follow same 
# #' distribution as negative droplets
# #' Default 0 for no rain. Used only when \code{type} is \code{fluo}.
# #' 
# #' @details sim_ddpcr_bkm is based on the R code from Jacobs et al. (2014) (see references).
# #' @references
# #' Jacobs B, Goetghebeur E, Clement L \emph{Impact of variance components on reliability of 
# #' absolute quantification using digital PCR} BMC Bioinformatics, 2014.
# #' @export
# #' @examples
# #' #two concentration, each 3 repetitions
# #' dat2_3 <- sim_ddpcr_bkm(c(0.5, 0.6), n_exp = 3, type = "tnp")
# 
# 
# sim_ddpcr_bkm <- function(m, n = 20000L, mexp = TRUE, n_exp = 8L, type = "np", 
#                           fluo_range = NULL, sddropc = 0, mudropr = 1, 
#                           sddropr = 0, Pvar = TRUE,
#                           piperr = 0, dropsd = 0, falpos = 0, falneg = 0, 
#                           rain = 0) {
#   
#   ##############
#   ### checks ###
#   ##############
#   
#   if(!is.logical(mexp)) stop("mexp must be a logical argument (TRUE or FALSE).", call. = TRUE, domain = NA)
#   if(max(!is.finite(m))) stop("Concentrations should all be numeric.", call. = TRUE, domain = NA)
#   if(min(m) < 0) stop("Concentrations cannot be negative.", call. = TRUE, domain = NA)
#   lambda <- if(mexp)
#     m else m * 0.89 / 1000
#   
#   if(!is.numeric(n)) stop("Number of droplets must have a numeric argument.", call. = TRUE, domain = NA)
#   if(n < 10) stop("Number of droplets must be larger than 10.", call. = TRUE, domain = NA)
#   if(!is.integer(n)) {
#     warning("Number of droplets will be rounded up to the next integer.", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
#     n <- ceiling(n)
#   }
#   
#   if(!(type %in% c("tnp", "fluo", "np")))
#     stop("Invalid value of 'type' parameter.", call. = TRUE, domain = NA)
#   
#   if(!is.numeric(n_exp)) stop("number of replicates must have a numeric argument.", call. = TRUE, domain = NA)
#   if(n_exp < 1) stop("number of replicates must be at least 1.", call. = TRUE, domain = NA)
#   if(n_exp != as.integer(n_exp)) {
#     warning("number of replicates 'n_exp' will be rounded up to the next integer.", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
#     n_exp <- as.integer(n_exp)
#   }
#   
#   if(!is.numeric(sddropc)) stop("sddropc must have a numeric argument.", call. = TRUE, domain = NA)
#   if(sddropc < 0) {warning("sddropc will be set to 0.", call. = TRUE, domain = NA)
#                    sddropc <- 0}
#   if(sddropc > n / 5) {warning("sddropc will be set to n/5.", call. = TRUE, domain = NA)
#                        sddropc <- n / 5}
#   
#   if(!is.numeric(mudropr)) stop("mudropr must have a numeric argument.", call. = TRUE, domain = NA)
#   if(mudropr * n < 10) stop("mudropr too small, too few droplets will be returned.", call. = TRUE, domain = NA)
#   if(mudropr > 1) warning("mudropr will be set to 1.", call. = TRUE, domain = NA) # happens in code
#   
#   if(mudropr < 1){if(!is.numeric(sddropr)) stop("sddropr must have a numeric argument.", call. = TRUE, domain = NA)
#                   if(sddropr < 0) {warning("sddropr will be set to 0.", call. = TRUE, domain = NA)
#                                    sddropr <- 0}
#                   if((sddropr >= mudropr) | ((sddropr + mudropr) >= 1)) stop("sddropr too large.", call. = TRUE, domain = NA)
#   }
#   
#   if(!is.numeric(piperr)) stop("pipette error must have a numeric argument.", call. = TRUE, domain = NA)
#   if(piperr < 0) stop("pipette error should be positive or 0.", call. = TRUE, domain = NA)
#   
#   if(!is.numeric(dropsd)) stop("dropsd must have a numeric argument.", call. = TRUE, domain = NA)
#   if(dropsd < 0) {
#     warning("dropsd will be set to 0.", call. = TRUE, domain = NA)
#     dropsd <- 0
#   }
#   
#   if(!is.logical(Pvar)) stop("Pvar must be a logical argument (TRUE or FALSE).", call. = TRUE, domain = NA)
#   
#   
#   if(!is.numeric(falpos)) stop("falpos must have a numeric argument.", call. = TRUE, domain = NA)
#   
#   if(falpos < 0) {
#     warning("falpos will be set to 0.", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
#     falpos <- 0
#   }
#   
#   if(!is.numeric(falneg)) stop("falneg must have a numeric argument.", call. = TRUE, domain = NA)
#   if(falneg < 0) {
#     warning("falneg will be set to 0.", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
#     falneg <- 0
#   }
#   if(falpos >= 1) stop("falpos too large, set a number between 0 and 1.", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
# 
#   if(!is.numeric(rain)) stop("rain must have a numeric argument.", call. = TRUE, domain = NA)
#   
#   if(rain < 0) {
#     warning("rain will be set to 0.", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
#     rain <- 0
#   }
#   
#   if(rain >= 1) {
#     warning("rain will be set to 1.", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
#     rain <- 1
#   }
#   
#   #add check for fluo
#   
#   ###############
#   ### repfunc ###
#   ###############
#   
#   # Same procedure for all replicates
#   # repfunc is called internally in samfunc
#   repfunc <- function(repdat) {
#     dropmem <- sample(repdat[1], repdat[2], replace = TRUE, 
#                       prob = rlnorm(repdat[1], 0, dropsd))
#     # droplet membership, probability proportional to size, size following a lognormal distribution
#     dropn <- ifelse(mudropr >= 1, repdat[1], 
#                     round(repdat[1]*plogis(rnorm(1, log(mudropr/(1 - mudropr)), log((mudropr + sddropr)/(mudropr - sddropr)*(1 - mudropr + sddropr)/(1 - mudropr - sddropr))/2))))
#     # number of droplets retained
#     dropmem <- dropmem[dropmem <= dropn]
#     # only retain copies of which the droplet is retained (lower rank)
#     dropvec <- 1L:dropn %in% dropmem
#     # vector with TRUE for positive droplets and FALSE for negative
#     
#     if(type != "fluo") {
#       return_drops <- 1L:length(dropvec)
#       return_drops[dropvec] <- rbinom(sum(dropvec), 1, 1 - falneg)
#       return_drops[!dropvec] <- rbinom(sum(!dropvec), 1, falpos)
#       # vector TRUE for positive signal and FALSE for negative signal
# 
#       return_fluo <- NULL
#     } else {
#       fluopeaks <- rnorm(dropn, 1000, 100) + 
#         8000*dropvec*(1 - runif(dropn)^(1/rain - 1))*(1 - rain^2) + 
#         2000*(1 - dropvec)*(1 - runif(dropn)^rain)*(1 - rain^2)
#       # random variation+downward rain+upward rain
#       dropfin <- (fluopeaks > 2500)
#       # hard threshold as in most software these days
#       # vector TRUE for positive signal and FALSE for negative signal
#       
#       return_drops <- list(dropfin = dropfin, dropn = dropn)
#       
#       return_fluo <- if (is.null(fluo_range)) {
#         fluopeaks
#       } else {
#         fluopos <- (1L:dropn)*fluo_range + 9 + runif(dropn)*2 + rnorm(dropn)
#         # vector of positions where the peak was found
#         # peaks on average 10 positions away from each other.
#         fluox <- 1L:(round(dropn*fluo_range + 90, -2))
#         fluoy <- rnorm(length(fluox), 50, 10)
#         # define fluo vectors with random background
#         j <- 1L
#         
#         #adjusted fluo_range
#         adj_fluo_range <- 30/fluo_range
#         for (i in fluox) {
#           if (j <= dropn) {
#             if (fluopos[j] < (i - 20)) {
#               j <- j + 1L
#             }
#             for(k in j:round(j + adj_fluo_range)){
#               # move j such that only the influence of peaks close by are counted
#               # influence of peaks further away would be marginally small anyway
#               if(k <= dropn) {
#                 dist <- fluopos[k] - i
#                 # distance between peak and current location
#                 fluoy[i] <- fluoy[i] + dnorm(dist/2 + rnorm(1, 0, 0.1))/dnorm(0)*fluopeaks[k]*0.95
#                 # add fluorescence signal stemming from this specific droplet
#               }
#             }
#           }
#         }
#         fluoy
#       }
#     }
#     
#     list(drop = return_drops, fluo = return_fluo)
#     # returns list with first element vector of droplets 1/0
#     # or pair of number of pos droplets and total number of droplets
#     # and second element either NULL, peak fluorescence of droplets
#     # or vector with continuous fluorescence output 
#   }
#   
#   
#   ###############
#   ### samfunc ###
#   ###############
#   
#   # Same procedure for all simulations
#   # samfunc is called in sim_ddpcr
#   samfunc <- function(lambdan){
#     dropstart <- round(rnorm(n_exp,n,sddropc))
#     # number of droplets
#     copyvar <- lambdan*rnorm(n_exp, 1, piperr)*dropstart
#     copyvar[copyvar < 0] <- 0
#     # expected number of copies after pipette variation
#     copyn <- ifelse(rep(Pvar, n_exp), rpois(n_exp, copyvar), round(copyvar))
#     # number of copies
#     lamdummy <- rep(lambdan, n_exp)
#     repdat <- data.frame(dropstart, copyn, lamdummy)
#     # number of droplets and copies in a list with n_exp elements, all pairs
#     # droplets is the first element, copies the second
#     apply(repdat, 1, repfunc)
#     # returns a list with n_exp*2 elements with
#     # first element vector of droplets 1/0 or
#     #   pair of number of pos droplets and total number of droplets
#     # second element either NULL, peak fluorescence of droplets or
#     #   vector with continuous fluorescence output 
#     # and so on for each replicate
#   }
#   
#   
#   ###############
#   ### execute ###
#   ###############
#   
#   out <- unlist(lapply(lambda, samfunc), recursive = FALSE)
#     
#   res <- suppressMessages(bind_dpcr(lapply(out, function(single_run) {
#     switch(type,
#            tnp = create_dpcr(sum(single_run[["drop"]]), length(single_run[["drop"]]), 
#                              NULL, type = "tnp", threshold = 2000, adpcr = FALSE),
#            np = create_dpcr(single_run[["drop"]], length(single_run[["drop"]]), NULL, 
#                             type = "np", threshold = 2000, adpcr = FALSE),
#            fluo = create_dpcr(single_run[["fluo"]], length(single_run[["drop"]]), NULL, 
#                               type = "fluo", threshold = 2000, adpcr = FALSE))
#   })))
# 
#   #crude solution to naming of the experiments
#   colnames(res) <- sapply(1L:length(lambda), function(single_lambda) 
#     paste0(single_lambda, ".", 1L:n_exp))
#   
#   res
# }
# 
