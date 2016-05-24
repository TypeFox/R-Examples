setup.ncore <- function(ncore, bigmem = FALSE) {
  if(is.null(ncore) || ncore > 1) {
    os1 <- .Platform$OS.type
    if(os1 == "windows") {
       if(is.null(ncore)) 
          ncore = 1
       else
          stop("Multicore is NOT supported in Windows (Set ncore = 1 or NULL)")
    } else {
       if(bigmem) {
         oops <- requireNamespace("bigmemory", quietly = TRUE)
         if(!oops) {
           if(is.null(ncore))
             ncore <- 1
           else
             stop("Please install the bigmemory package from CRAN for running with multicore")
         }
       }
       if(is.null(ncore))
         ncore = parallel::detectCores()

       # Following lines check R internal varible for potential limit on multicore usage
       # Normally it does nothing, but will be helpful in running `R CMD check --as-cran`
       if(ncore > 1) {
          chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
          if (nzchar(chk) && (chk != "false")) ncore = 1
       }

    }
  }
  options(mc.cores = ncore)
  return(ncore)
}
