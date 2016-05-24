
library(inline)

fgsl <- cxxfunction(signature(ns="integer", ws="numeric"),
                    plugin="RcppGSL",
                    include="#include <gsl/gsl_randist.h>",  body='

   unsigned int n = Rcpp::as<int>(ns);
   Rcpp::NumericVector w = Rcpp::NumericVector(ws);
   unsigned int k = w.size();

   //void rng::Multinomial(unsigned n, unsigned k, const double* w, unsigned* X)
   //{
   //	gsl_ran_multinomial(pWorkspace, k, n, w, X);
   //}

   const gsl_rng_type* type;
   gsl_rng* pWorkspace;

   gsl_rng_env_setup();
   type = gsl_rng_default;
   pWorkspace = gsl_rng_alloc(gsl_rng_default);

   unsigned int *z = new unsigned int[k];

   gsl_ran_multinomial(pWorkspace, k, n, w.begin(), z);

   Rcpp::IntegerVector X(k);
   for (unsigned int i=0; i<k; i++) X[i] = z[i];

   gsl_rng_free(pWorkspace);
   delete[] z;

   return X;
')

fr <- cxxfunction(signature(ns="integer", ws="numeric"), plugin="Rcpp", body='

   RNGScope tmp;
   unsigned int n = Rcpp::as<int>(ns);
   Rcpp::NumericVector w(ws);
   unsigned int k = w.size();
   Rcpp::IntegerVector x(k);

   rmultinom(n, w.begin(), k, x.begin());

   return x;
')



fgsl(30000, c(0.25, 0.25, 0.25, 0.25))
fr(30000, c(0.25, 0.25, 0.25, 0.25))
