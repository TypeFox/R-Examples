## TODO -- test that there is some variance across a set of bootstrap
##         samples

## TODO -- for ratios, test also that there is some variance in
##         numerator and denominator

## TODO -- test paired ego / alter datasets

## TODO -- test that calling bootstrap.estimates works when
##         total.popn.size is an argument and not an attribute of
##         the data frame (had to use parent.frame(2)) to fix
##         a bug about this

## TODO -- test cases where estimates should never be negative

## TODO -- look at
## http://stackoverflow.com/questions/8898469/is-it-possible-to-use-r-package-data-in-testthat-tests-or-run-examples
## to try and figure out the real way to include package data in
## unit tests...

set.seed(12345)

#########################################
## setup
## NB: see this for possible alternatives
## http://stackoverflow.com/questions/8898469/is-it-possible-to-use-r-package-data-in-testthat-tests-or-run-examples
load("mu284.RData")

#########################################
## rescaled (Rao / Wu) bootstrap

## TODO -- think about how many bootstrap reps
## we should use in the unit tests...
##M <- 1000
M <- 2000

context(paste0("variance estimators - rescaled bootstrap - correctness (M=", M, ")"))

rbsfn <- functional::Curry(bootstrap.estimates,
                           survey.design= ~ CL,
                           num.reps=M,
                           estimator.fn=MU284.estimator.fn,
                           weights="sample_weight",
                           bootstrap.fn="rescaled.bootstrap.sample")

test.boot <- llply(MU284.surveys,
                   function(svy) { do.call("rbind", rbsfn(survey.data=svy)) })

test.boot.summ <- ldply(test.boot,
                        summarize,
                        mean.TS82.hat=mean(TS82.hat),
                        mean.R.RMT85.P85.hat=mean(R.RMT85.P85.hat),
                        sd.TS82.hat=sd(TS82.hat),
                        sd.R.RMT85.P85.hat=sd(R.RMT85.P85.hat))

## TODO -- how to figure out what tolerance to use?
## for now, using .1 for everything, but this is pretty big; also,
## we'd expect different tolerances for different qois, i think.
qoi <- colnames(test.boot.summ)

l_ply(qoi,
      function(this.qoi) {
          l_ply(1:nrow(test.boot.summ),
                function(idx) {
                    expect_that(test.boot.summ[idx,this.qoi],
                                equals(MU284.boot.res.summ[idx,this.qoi],
                                       tolerance=.05),
                                label=paste0("qty: ", this.qoi,
                                             "; MU284 survey #", idx))
                })
      })

