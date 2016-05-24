
library(inline)

f1 <- cxxfunction(signature(xs="numeric"), plugin="Rcpp", body='
   Rcpp::NumericVector x(xs);
   return Rcpp::wrap(Rcpp::dnorm(x));
')

f2 <- cxxfunction(signature(xs="numeric"), plugin="RcppGSL",
                  include="#include <gsl/gsl_randist.h>",  body='
   double x = Rcpp::as<double>(xs);
   double y = gsl_ran_ugaussian_pdf(x);
   return Rcpp::wrap(y);
')

for (i in seq(-3,3,by=0.25)) print(c(f1(i), f2(i)))
