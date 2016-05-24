
# calling message from .onLoad gives a NOTE on R CMD CHECK, so we circumvent that here.
mymsg <- message

.onLoad <- function(libname, pkgname){

  nthread = parallel::detectCores()

  if ( is.na(nthread) || !is.numeric(nthread) ){
    nthread <- 1L
    mymsg("Could not detect number of cores, defaulting to 1.")
  }

  omp_thread_limit = as.numeric(Sys.getenv("OMP_THREAD_LIMIT"))
  if ( is.na(omp_thread_limit) ) omp_thread_limit <- nthread
  
  nthread = min(omp_thread_limit,nthread)
  if (nthread >= 4) nthread <- nthread - 1
 
  options(hashr_num_thread=nthread)
}


