#' Simulate SONE values for scenario.
#'
#' @description A wrapper that takes a scenario, and produces the Significance Of iNdicator Exclusion (SONE) values for each exclusion and calculates efficacy. Used by  \code{\link{optimal_p}}.
#' @inheritParams scale_sim
#' @param n_sim number of simulations. 1000 is a start, 10000 was used in paper, but takes a long time
#' @param sizes An array of sample sizes to be simulated. Can be single value.
#' @param ... further tweaking of the scale simulator, see  \code{\link{scale_sim}} for details.
#' @return Returns a list of SONE values and related efficacy. See example for details
#' \enumerate{
#'   \item SONE results. Feed this to  \code{\link{scenario_plot}} or \code{\link{scenario_plot80}} (see examples)
#'   \item Summary efficacy. Such data comprises Table 1 in Vainik, Mottus et al., 2015 EJP
#'   \item Full efficacy data.
#' }
#' @encoding utf-8

#' @examples
#'
#' #A scenario with 8 items relating to outcome, testing 2 different samples
#'sizes=c(250,500)
#'n_sim=100
#'to_n=8
#'scen1=scenario_sim(sizes,n_sim,to_n)  # takes a few seconds..
#'scenario_plot80(scen1[[1]],sizes,n_sim)
#'
#'# A scenario with 2 out of 8 items relating to outcome, 2 different samples
#'to_n=2
#'scen2=scenario_sim(sizes,n_sim,to_n)  # takes a few seconds..
#'scenario_plot(scen2[[1]],sizes,n_sim,to_n)
#'



#' @export

scenario_sim <- function(sizes, n_sim, to_n, tn_n = 8 - to_n, ...) {
    
    # A wrapper that takes a scenario, and produces the SONE values for each exclusion and calculates efficacy
    # calculate means and CI-s for data
    meanci <- function(dat) {
        # calculate means and CI-s for data
        
        out_dat <- array(0, dim = c(nrow(dat), 3))
        out_dat[, 1] <- rowMeans(dat)
        out_dat[, 2:3] <- t(apply(dat, 1, quantile, probs = c(0.025, 0.975)))
        # if (boot==T) out_dat[,2:3]= t(apply(dat, 1, bsci, B=bootn)) else
        colnames(out_dat) <- c("mean", "2.5%", "97.5%")
        return(out_dat)
        
    }
    
    # calculate efficacy
    
    efficacy_correct <- function(eff) {
        # calculate efficacy. display percentages of correctly detected indicators for Study1 for the condition where no
        # more indicators should be discarded, the script displays the average percentage. For details see description of
        # efficacy_correct_single
        
        out <- array("", dim = nrow(eff))
        for (i in nrow(eff):1) {
            lowest <- summary(ordered(eff[i, ]))  # summary of indicators excluded across all simulations, for a given iteration
            lowest.exc <- round(((lowest/sum(lowest)) * 100), 2)  # convertthe summary into percentages
            if (i == nrow(eff)) {
                out[i] <- paste0(round(mean(lowest.exc), 2), "(", round(sd(lowest.exc), 2), ")")  # assess the efficacy of 8+0 scenario
            } else {
                # out[i]=as.character(sum(lowest.exc[1:(nrow(eff)-i)]))
                out[i] <- as.character(sum(lowest.exc[length(lowest.exc) - ((nrow(eff) - i - 1):0)]))  # get the percentage of largest numbers that belong to to/t2
                
            }
        }
        return(out)
    }
    
    
    efficacy_correct_single <- function(eff) {
        # calculate efficacy for the 8+0 Scale. Displays the mean and SD of percent of indicators discarded in 8+0 scale.
        # Ideally, the percentage should equal 1/nindicators, and the SD be as low as possible We ended up not using it
        # for the study 1
        
        
        lowest <- summary(ordered(eff))  # summary of indicators excluded across all simulations, for a given iteration
        lowest.exc <- round(((lowest/sum(lowest)) * 100), 2)  # convertthe summary into percentages
        out <- paste0(round(mean(lowest.exc), 2), "(", round(sd(lowest.exc), 2), ")")
        
        return(out)
    }
    
    
    
    
    
    
    if (length(sizes) == 1) 
        sizes <- c(sizes, 20)  # temporary fix to deal with a single 'sizes' value. Currently the script assumes that 'sizes' have at least two values, and the script breaks down otherwise. Therefore, a small secondary sample size is added, and later removed from the results. This fix will be addressed in a future release.
    
    # possibly could be fixed by a predefined matrix.
    if (tn_n == 0) {
        # for the case where all indicators relate to an outcome
        out <- array(0, dim = c(length(sizes), n_sim))
        excl_list <- array(0, c(length(sizes), n_sim))
        for (h in 1:length(sizes)) {
            for (i in 1:n_sim) {
                it_out <- scale_sim(n = sizes[h], to_n = to_n, tn_n = tn_n, ...)
                i8 <- ind_excl_step(it_out[[1]], it_out[[2]])
                excluded <- as.double(rownames(i8)[1])
                excl_list[h, i] <- excluded
                out[h, i] <- i8[1, 2]
            }
        }
        out_dat <- meanci(out)
        
        out_excl <- array("", dim = length(sizes))
        names(out_excl) <- sizes
        
        for (i in 1:length(sizes)) {
            out_excl[i] <- efficacy_correct_single(excl_list[i, ])
            
        }
        
        # cut the dummy results for the 6-person sample size, if needed
        if (sizes[2] == 20) {
            output <- list(out_dat[1, ], out_excl[1], excl_list[1, ])
        } else {
            output <- list(out_dat, out_excl, excl_list)
        }
        
        
        return(output)
    } else {
        excl <- to_n + 1
        out <- array(0, dim = c(excl, length(sizes), n_sim))
        excl_list <- array(0, dim = c(excl, length(sizes), n_sim))
        
        
        for (h in 1:length(sizes)) {
            for (i in 1:n_sim) {
                it_out <- scale_sim(n = sizes[h], to_n = to_n, tn_n = tn_n, ...)  # simulate trait
                i8 <- ind_excl_step(it_out[[1]], it_out[[2]])  # run 1-indicator-out
                out[1, h, i] <- i8[1, 2]  # export the highest p-value
                excluded <- as.double(rownames(i8)[1])  # find the indicator with highest pvalue
                excl_list[1, h, i] <- excluded
                indicators_excl <- it_out[[1]][, -excluded]  # exclude it from next iterations
                for (j in 2:excl) {
                  iexc <- ind_excl_step(indicators_excl, it_out[[2]])  # run 1-indicator-out again on smaller indicator sets
                  out[j, h, i] <- iexc[1, 2]
                  excluded <- as.double(rownames(iexc)[1])
                  excl_list[j, h, i] <- excluded
                  indicators_excl <- indicators_excl[, -excluded]
                }
            }
        }
        out_dat <- array(0, dim = c(excl, length(sizes), 3))
        
        for (i in 1:excl) {
            out_dat[i, , ] <- meanci(out[i, , ])
        }
        
        out_excl <- array("", dim = c(excl, length(sizes)))
        colnames(out_excl) <- sizes
        temprow <- c(paste0("+", to_n:1), " mean(SD)")
        rownames(out_excl) <- paste0(tn_n, temprow)  # before 16.01.19     rownames(out_excl)=paste0(nindicators-to,temprow)
        
        for (i in 1:length(sizes)) {
            out_excl[, i] <- efficacy_correct(excl_list[, i, ])
            
        }
        
        # cut the dummy results for the 6-person sample size, if needed
        if (sizes[2] == 20) {
            output <- list(out_dat[, 1, ], out_excl[, 1], excl_list[, 1, ])
        } else {
            output <- list(out_dat, out_excl, excl_list)
        }
        
        return(output)
    }
    
} 
