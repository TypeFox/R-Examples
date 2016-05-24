################################################################################
# fstOnly: a memory efficient function to calculate WC Fst and Fit
################################################################################
#' @export
fstOnly <- function(infile = NULL, outfile = NULL, gp = 3, 
                    bs_locus = FALSE, bs_pairwise = FALSE, 
                    bootstraps = 0, parallel = FALSE){
  .Deprecated(new = "diffCalc", msg = "This function is no longer in use. Please use 'diffCalc' instead, \nSee ?diffCalc for usage details.", 
              old = "fstOnly")
}
################################################################################
# END
################################################################################