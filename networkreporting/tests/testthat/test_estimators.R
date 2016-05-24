##########################################################
## test-estimators.R
##
## unit tests for functions that compute estimators
## from a dataset
##
## TODO -- eventually, develop a catalog of simple networks
##         that we can hand-compute estimator values for,
##         and that can be part of these tests

## NB: for each estimator, we should test the following:
##     * for one or more specific datasets, test that the
##       estimator return the (known) right answer
##
##     * [NOT IMPLEMENTED YET] test the accuracy of the
##       leave-one-out validation procedure
##
##     * [NOT IMPLEMENTED YET] for the same datasets,
##       test that the weighted versions of the estimator
##       returns the (known) right answer
##
##     * [NOT IMPLEMENTED YET] test that NAs are handled
##       in weighted and unweighted estimators
##
##     * test error handling in passing in and using
##       column names
##
##     * [NOT IMPLEMENTED YET] test the other ways
##       of passing in and using column names
##
##     * [NOT IMPLEMENTED YET] test the variance
##       estimators (this will likely be several tests)

## these tests use the toy networks that come
## packaged with the networksampling package
## TODO -- I don't understand why the package
## data aren't available without having to
## specify package=...
## (this could be a devtools thing?)
load("toynetworks.RData")
load("toynrnetworks.RData")

####################################
## multiplicity estimator
context("estimators - multiplicity")

ests <- ldply(toy.networks,
              function(tn) {
                return(c("est"=multiplicity.estimator(tn,
                                                      mult.col="mult.response")))
              })
truth <- ldply(toy.networks,
               function(tn) { c("truth"=attr(tn, "multiplicity.estimate")) })

tocheck <- merge(ests,
                 truth,
                 by=".id",
                 all=TRUE)

d_ply(tocheck,
      .(.id),
      function(x) {
        expect_that(x$est, equals(x$truth),
                    info=paste("estimate is", x$est,
                               "but it should be", x$truth),
                    label=paste("multiplicity estimate on toy network ",
                                x$".id"))
      })


####################################
## nsum estimator
context("estimators - nsum")

## TODO - add test for Killworth estimate of the
##        standard error (both for proportions and for totals)

ests <- ldply(toy.networks,
              function(tn) {
                return(data.frame(nsum.estimator(tn,
                                                 d.hat.col="d",
                                                 y.col="y",
                                                 total.popn.size=NA)))
              })
truth <- ldply(toy.networks,
               function(tn) { c("truth"=attr(tn, "nsum.estimate")) })

tocheck <- merge(ests,
                 truth,
                 by=".id",
                 all=TRUE)

d_ply(tocheck,
      .(.id),
      function(x) {
        expect_that(x$estimate, equals(x$truth),
                    info=paste("estimate is", x$est,
                               "but it should be", x$truth),
                    label=paste("nsum estimate on toy network ",
                                x$".id"))
      })

## TODO -- we should also check the numerators and the denominators
## TODO -- add weights and check that weighted estimates work properly
##         as well


####################################
## gnsum estimator
context("estimators - gnsum")

ests <- ldply(toy.networks,
              function(tn) {
                d.T.bar <- mean(subset(tn,hidden)$d)
                d.bar <- mean(tn$d)
                delta <- d.T.bar/d.bar

                return(data.frame(nsum.estimator(tn,
                                                 d.hat.col="d",
                                                 y.col="y",
                                                 deg.ratio=delta,
                                                 total.popn.size=NA)))
              })
truth <- ldply(toy.networks,
               function(tn) { c("truth"=attr(tn, "gnsum.estimate")) })

tocheck <- merge(ests,
                 truth,
                 by=".id",
                 all=TRUE)

d_ply(tocheck,
      .(.id),
      function(x) {
        expect_that(x$est, equals(x$truth),
                    info=paste("estimate is", x$est,
                               "but it should be", x$truth),
                    label=paste("gnsum estimate on toy network ",
                                x$".id"))
      })

####################################
## network reporting estimator
context("estimators - network survival")

## TODO (NB: see README from non-versioned networksampling folder)

## NB: for now, we're just testing the first network, since
##     the others have missing values, which we're not handling yet
#tocheck <- llply(1:length(toy.nr.networks),
tocheck <- plyr::llply(c(1),
              function(this.nrnet.idx) {

                this.nrnet <- toy.nr.networks[[this.nrnet.idx]]

                ## no weights in toy network data, so give everyone a weight
                ## of 1 for now
                this.nrnet$weight <- 1
                this.nrnet$nrid <- this.nrnet.idx

                this.attrib <- toy.nr.long.networks[[this.nrnet.idx]]
                this.attrib$ego.weight <- 1

                netsurv.est <- network.survival.estimator_(resp.data=this.nrnet,
                                                           attribute.data=this.attrib,
                                                           attribute.names=c("age", "sex"),
                                                           known.populations="d",
                                                           weights="weight",
                                                           attribute.weights="ego.weight",
                                                           total.kp.size=nrow(this.nrnet))

                truth <- attr(this.nrnet, "ns.estimate")
                truth <- plyr::rename(truth, c('est'='truth'))

                netsurv.est <- dplyr::left_join(netsurv.est, truth,
                                                by=c('age', 'sex'))

                return(netsurv.est)
              })

l_ply(tocheck,
      function(x) {
        expect_that(x$asdr.hat, equals(x$truth),
                    info=paste("estimate is", paste(x$asdr.hat, collapse=""),
                               "but it should be", paste(x$truth, collapse="")),
                    label=paste("network survival estimate on toy network ",
                                x$nrid[1]))
      })



