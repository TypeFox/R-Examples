### R code from vignette source 'RcppGSL-intro.Rnw'

###################################################
### code chunk number 1: setup
###################################################
require(inline)
library(RcppGSL)
options("width"=65)
rcppgsl.version <- packageDescription( "RcppGSL" )$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 7: inlineex1
###################################################
fx <- Rcpp::cppFunction("int sum_gsl_vector_int(RcppGSL::vector<int> vec) {
    int res = std::accumulate( vec.begin(), vec.end(), 0);
    return res;
}", depends="RcppGSL")


###################################################
### code chunk number 8: callinlineex1
###################################################
sum_gsl_vector_int(1:10)


###################################################
### code chunk number 10: inlinexex2
###################################################
Rcpp::cppFunction("double gsl_vector_sum_2(Rcpp::List data) {
    RcppGSL::vector<double> x = data[\"x\"];
    RcppGSL::vector<int> y = data[\"y\"];
    double res = 0.0;
    for (size_t i=0; i< x->size; i++) {
        res += x[i] * y[i];
    }
    return res;
}", depends= "RcppGSL")


###################################################
### code chunk number 11: callinlineex2
###################################################
data <- list( x = seq(0,1,length=10), y = 1:10 )
gsl_vector_sum_2(data)


###################################################
### code chunk number 23: RcppGSL-intro.Rnw:821-847
###################################################
require(inline)

inctxt='
   #include <gsl/gsl_matrix.h>
   #include <gsl/gsl_blas.h>
'

bodytxt='
  RcppGSL::matrix<double> M = sM;     // create data structures from SEXP
  int k = M.ncol();
  Rcpp::NumericVector n(k);           // to store results

  for (int j = 0; j < k; j++) {
    RcppGSL::vector_view<double> colview = gsl_matrix_column (M, j);
    n[j] = gsl_blas_dnrm2(colview);
  }
  return n;                           // return vector
'

foo <- cxxfunction(
    signature(sM="numeric"),
    body=bodytxt, inc=inctxt, plugin="RcppGSL")

## see Section 8.4.13 of the GSL manual: create M as a sum of two outer products
M <- outer(sin(0:9), rep(1,10), "*") + outer(rep(1, 10), cos(0:9), "*")
foo(M)


###################################################
### code chunk number 24: RcppGSL-intro.Rnw:853-854 (eval = FALSE)
###################################################
## package.skeleton("mypackage", foo)


###################################################
### code chunk number 26: RcppGSL-intro.Rnw:909-910 (eval = FALSE)
###################################################
## sourceCpp("gslNorm.cpp")


