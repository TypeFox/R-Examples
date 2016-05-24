## ----eval = F, engine='Rcpp'---------------------------------------------
## #include <sitmo.h> // SITMO PPRNG

## ----mm, engine='Rcpp',eval = F------------------------------------------
## // Generate engine called e.
## sitmo::prng_engine e;

## ----engine='Rcpp', eval = F---------------------------------------------
## #include <sitmo.h> // SITMO PPRNG
## 
## //' Example RNG Draws with sitmo
## //'
## //' Shows a basic setup and use case for sitmo.
## //'
## //' @param n A \code{unsigned int} is a .
## //' @return A \code{vec} with random sequences.
## //' @examples
## //' n = 10
## //' sitmo_draws(n)
## // [[Rcpp::export]]
## Rcpp::NumericVector sitmo_draws(unsigned int n) {
## 
##   Rcpp::NumericVector o(n);
## 
##   // Create a prng engine
##   sitmo::prng_engine eng;
## 
##   // Draw from base engine
##   for (unsigned int i=0; i< n ; ++i){
##     o(i) = eng();
##   }
## 
##   return o;
## }

## ----engine='Rcpp', eval = F---------------------------------------------
## #include <sitmo.h> // SITMO PPRNG
## 
## //' Example Seed Set and RNG Draws with sitmo
## //'
## //' Shows how to set a seed in sitmo.
## //'
## //' @param n    An \code{unsigned int} that dictates how many realizations occur.
## //' @param seed An \code{unsigned int} that controls the rng seed.
## //' @return A \code{vector} with random sequences.
## //' @examples
## //' n = 10
## //' a = sitmo_engine_seed(n, 1337)
## //' b = sitmo_engine_seed(n, 1337)
## //' c = sitmo_engine_seed(n, 1338)
## //'
## //' isTRUE(all.equal(a,b))
## //' isTRUE(all.equal(a,c))
## // [[Rcpp::export]]
## Rcpp::NumericVector sitmo_engine_seed(unsigned int n, unsigned int seed) {
## 
##   // Create Rcpp Matrix
##   Rcpp::NumericVector o(n);
## 
##   // Create a prng engine with a specific seed
##   sitmo::prng_engine eng(static_cast<uint32_t>(seed));
## 
##   // Draw from base engine
##   for (unsigned int i=0; i < n; ++i){
##     o(i) = eng();
##   }
## 
##   return o;
## }

## ----engine='Rcpp', eval = F---------------------------------------------
## #include <sitmo.h> // SITMO PPRNG
## 
## //' Example Seed Set and RNG Draws with sitmo
## //'
## //' Shows how to set a seed in sitmo.
## //'
## //' @param n    An \code{unsigned int} that dictates how many realizations occur.
## //' @param seed An \code{unsigned int} that controls the rng seed.
## //' @return A \code{matrix} with random sequences.
## //' @examples
## //' n = 10
## //' a = sitmo_engine_seed(n, 1337)
## //'
## //' isTRUE(all.equal(a[,1],a[,2]))
## // [[Rcpp::export]]
## Rcpp::NumericMatrix sitmo_engine_reset(unsigned int n, unsigned int seed) {
## 
##   // Create Rcpp Vector
##   Rcpp::NumericMatrix o(n,2);
## 
##   // Create a prng engine with a specific seed
##   sitmo::prng_engine eng(static_cast<uint32_t>(seed));
## 
##   // Draw from base engine
##   for (unsigned int i=0; i < n ; ++i){
##     o(i,0) = eng();
##   }
## 
##   // Reset seed
##   eng.seed();
## 
##   // Draw from base engine
##   for (unsigned int i=0; i< n ; ++i){
##     o(i,1) = eng();
##   }
## 
##   return o;
## }

## ----engine='Rcpp',eval = F----------------------------------------------
## #include <sitmo.h> // SITMO PPRNG
## 
## //' Two RNG engines running side-by-side
## //'
## //' Shows how to create two separate RNGs and increase them together.
## //'
## //' @param n     An \code{unsigned int} that dictates how many realizations occur.
## //' @param seeds A \code{vec} containing two integers greater than 0.
## //' @return A \code{matrix} with random sequences.
## //' @examples
## //' n = 10
## //' a = sitmo_engine_seed(n, c(1337,1338))
## //'
## //' b = sitmo_engine_seed(n, c(1337,1337))
## //'
## //' isTRUE(all.equal(a[,1],a[,2]))
## //'
## //' isTRUE(all.equal(b[,1],b[,2]))
## //'
## //' isTRUE(all.equal(a[,1],b[,1]))
## // [[Rcpp::export]]
## Rcpp::NumericMatrix sitmo_two_seeds(unsigned int n, Rcpp::NumericVector seeds) {
## 
##   if(seeds.size() != 2) Rcpp::stop("Need exactly two seeds for this example.");
## 
##   // Create Rcpp Matrix
##   Rcpp::NumericMatrix o(n,3);
## 
##   // Create a prng engine with a specific seed
##   sitmo::prng_engine eng1;
##   eng1.seed(seeds(0));
## 
##   sitmo::prng_engine eng2;
##   eng1.seed(seeds(1));
## 
##   // Draw from base engine
##   for (unsigned int i=0; i< n ; ++i){
##     o(i,0) = eng1();
##     o(i,1) = eng2();
##   }
## 
##   return o;
## }

## ----engine='Rcpp', eval = F---------------------------------------------
## //' Random Uniform Number Generator with sitmo
## //'
## //' The function provides an implementation of sampling from a random uniform distribution
## //'
## //' @param n    An \code{unsigned integer} denoting the number of realizations to generate.
## //' @param min  A \code{double} indicating the minimum \eqn{a} value
## //'               in the uniform's interval \eqn{\left[a,b\right]}
## //' @param max  A \code{double} indicating the maximum \eqn{b} value
## //'               in the uniform's interval \eqn{\left[a,b\right]}
## //' @param seed A special \code{unsigned integer} containing a single seed.
## //' @return A \code{vec} containing the realizations.
## //' @export
## //' @examples
## //' a = runif_sitmo(10)
## // [[Rcpp::export]]
## Rcpp::NumericVector runif_sitmo(unsigned int n, double min = 0.0, double max = 1.0, uint32_t seed = 1) {
##   Rcpp::NumericVector o(n);
## 
##   // Create a prng engine
##   sitmo::prng_engine eng(seed);
##   // Obtain the range between max and min
##   double dis = max - min;
## 
##   for(int i = 0; i < n; ++i) {
##     // Sample from the RNG and divide it by the maximum value possible (can also use SITMO_RAND_MAX, which is 4294967295)
##     // Apply appropriate scale (MAX-MIN)
##     o[i] = min + ((double) eng() / (sitmo::prng_engine::max())) * (dis);
##   }
## 
##   return o;
## }

## ----engine='Rcpp', eval = F---------------------------------------------
## #ifdef _OPENMP
## #include <omp.h>
## #endif

## ----engine='Rcpp', eval = F---------------------------------------------
## #ifdef _OPENMP
## // multithreaded OpenMP version of code
## #else
## // single-threaded version of code
## #endif

## ----engine='asis', eval = F---------------------------------------------
## PKG_LIBS =  $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CFLAGS)
## PKG_CFLAGS = $(SHLIB_OPENMP_CFLAGS)

## ----engine='Rcpp', eval = F---------------------------------------------
## #include <Rcpp.h>
## #include <sitmo.h> // SITMO PPRNG
## 
## #ifdef _OPENMP
## #include <omp.h>
## #endif
## 
## //' Test Generation using sitmo and C++11
## //'
## //' The function provides an implementation of creating realizations from the default engine.
## //'
## //' @param n An \code{unsigned integer} denoting the number of realizations to generate.
## //' @param seeds A \code{vec} containing a list of seeds. Each seed is run on its own core.
## //' @return A \code{vec} containing the realizations.
## //' @details
## //' The following function's true power is only accessible on platforms that support OpenMP (e.g. Windows and Linux).
## //' However, it does provide a very good example as to how to make ones code applicable across multiple platforms.
## //'
## //' With this being said, how we determine how many cores to split the generation to is governed by the number of seeds supplied.
## //' In the event that one is using OS X, only the first seed supplied is used.
## //'
## //' @export
## //' @examples
## //' a = sitmo_parallel(10, 5.0, c(1))
## //'
## //' b = sitmo_parallel(10, 5.0, c(1))
## //'
## //' c = sitmo_parallel(10, 5.0, c(2))
## //'
## //' isTRUE(all.equal(a,b))
## //'
## //' isTRUE(all.equal(a,c))
## // [[Rcpp::export]]
## Rcpp::NumericVector sitmo_parallel(unsigned int n, Rcpp::NumericVector& seeds){
## 
##   unsigned int ncores = seeds.size();
## 
##   Rcpp::NumericVector q(n);
## 
##   #ifdef _OPENMP
##   #pragma omp parallel num_threads(ncores) if(ncores > 1)
##   {
##   #endif
## 
##     // Engine requires uint32_t inplace of unsigned int
##     uint32_t active_seed;
## 
##     // Write the active seed per core or just write one of the seeds.
##     #ifdef _OPENMP
##       active_seed = static_cast<uint32_t>(seeds[omp_get_thread_num()]);
##     #else
##       active_seed = static_cast<uint32_t>(seeds[0]);
##     #endif
## 
##     sitmo::prng_engine eng( active_seed );
## 
##     // Parallelize the Loop
##     #ifdef _OPENMP
##     #pragma omp for schedule(static)
##     #endif
##     for (unsigned int i = 0; i < n; i++){
##       q[i] = eng().
##     }
## 
##   #ifdef _OPENMP
##   }
##   #endif
## 
##   return q;
## }

