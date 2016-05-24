context("implemented reweighters")

reweighters <- list(adaboost=adaboostReweighter, arcfs=arcfsReweighter,
                 arcx4=arcx4Reweighter, vanillaBagger=vanillaBagger)

# see testReweighter.R for testReweighter, preds, response, wts and metadata

reweighterArgs <- list(prediction=preds, response=response,
                    m=0, weights=rep.int(1, length(preds)))

test_that("reweighters produce weights", {
  lapply(reweighters, function(reweighter) {
    out <- do.call(reweighter, reweighterArgs)
    weights <- out$weights
    expect_that(length(weights), equals(length(preds)))
  })
})

# test_that("vanilla bagger producers an estimator with 'oobPerf' attribute", {
#   out <- do.call(vanillaBagger, reweighterArgs)
#   estimator <- out$estimator
#   expect_that(is.null(estimator), is_false())
#   expect_that(is.null(attr(estimator, "oobPerf")), is_false())
# })