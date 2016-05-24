context("hhh4() model with shared overdispersion parameters")

## use a small subset of districts from the fluBYBW data
data("fluBYBW")
fluBWsub <- fluBYBW[, substr(colnames(fluBYBW), 1, 2) %in% "81"]
## stsplot_space(fluBWsub, labels = TRUE)

## set "neighbourhood" to order of adjacency + 1
neighbourhood(fluBWsub) <-  # nbOrder(neighbourhood(fluBWsub), maxlag = 5) + 1
    structure(
        c(1, 4, 3, 2, 2, 4, 2, 4, 3, 3, 4, 4, 5, 4, 1, 2, 3, 
          4, 5, 4, 2, 3, 4, 3, 4, 4, 3, 2, 1, 2, 3, 4, 3, 2, 2, 3, 3, 4, 
          4, 2, 3, 2, 1, 2, 4, 3, 3, 2, 3, 3, 4, 4, 2, 4, 3, 2, 1, 4, 2, 
          3, 2, 3, 3, 4, 4, 4, 5, 4, 4, 4, 1, 3, 4, 3, 2, 3, 3, 4, 2, 4, 
          3, 3, 2, 3, 1, 3, 2, 2, 3, 3, 4, 4, 2, 2, 3, 3, 4, 3, 1, 2, 3, 
          2, 3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 1, 2, 2, 3, 3, 3, 4, 3, 3, 3, 
          2, 2, 3, 2, 1, 2, 2, 3, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 2, 2, 
          4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 1, 2, 5, 4, 4, 4, 4, 4, 4, 3, 
          3, 3, 2, 2, 1),
        .Dim = c(13L, 13L),
        .Dimnames = list(
            c("8115", "8135", "8117", "8116", "8111", "8121", "8118", "8136",
              "8119", "8125", "8127", "8126", "8128"),
            c("8115", "8135", "8117", "8116", "8111", "8121", "8118", "8136",
              "8119", "8125", "8127", "8126", "8128")))

## a crazy model base
fluModel <- list(
    end = list(f = addSeason2formula(~0 + ri(type="iid"))),
    ne = list(f = ~0 + fe(1, unitSpecific = TRUE),
        weights = W_powerlaw(maxlag = 3)),
    start = list(random = rep.int(0, ncol(fluBWsub)))
)


if (FALSE) { # check derivatives
    fluDeriv <- hhh4(stsObj = fluBWsub,
                     control = c(fluModel, list(family = "NegBinM")),
                     check.analyticals = TRUE)
    ana <- fluDeriv$pen$fisher$analytic
    num <- fluDeriv$pen$fisher$numeric
    equal <- mapply(function (...) isTRUE(all.equal.numeric(...)),
                    ana, num, tolerance = 1e-4)
    dim(equal) <- dim(ana)
    Matrix::image(Matrix::Matrix(equal))
}


test_that("\"NegBinM\" fit is invariant to the ordering of the overdispersion parameters", {
    fluFitM <- hhh4(stsObj = fluBWsub, control = c(fluModel, list(
        family = "NegBinM"))) # identical to factor(colnames(fluBWsub), levels=colnames(fluBWsub))
    fluFitM_reordered <- hhh4(stsObj = fluBWsub, control = c(fluModel, list(
        family = factor(colnames(fluBWsub), levels=rev(colnames(fluBWsub))))))
    expect_equal(fluFitM_reordered$loglikelihood,
                 fluFitM$loglikelihood)
    expect_equal(fluFitM_reordered$margll,
                 fluFitM$margll)
    expect_equal(fluFitM_reordered$coefficients[names(fluFitM$coefficients)],
                 fluFitM$coefficients)
})


## fit a model with shared overdispersion parameters
fluFitShared <- hhh4(stsObj = fluBWsub, control = c(fluModel, list(
    family = factor(substr(colnames(fluBWsub), 3, 3) == "1",
                    levels = c(TRUE, FALSE), labels = c("region1", "elsewhere")))))

test_that("estimates with shared overdispersion are reproducible", {
    ## dput(coef(fluFitShared, se = TRUE))
    orig <- structure(
        c(0.0172448275799737, -2.29936227176632, -0.311391919170833, 
          0.0173369590386396, 0.242634649538434, -0.73402605050834, -0.0411427686831543, 
          -0.917845995715638, -0.324146451650439, -0.252506337389155, 0.153202205413176, 
          -0.857813219848051, -1.00758863915022, 2.01735387997105, 2.38047570484809, 
          -4.38317074697181, 2.46949727973784, 0.549903756338196, 1.12432744953686, 
          0.647372578569298, 0.21388842588635, -0.437822769909503, 0.255185408180267, 
          0.92949604237045, -1.09633602928844, 0.298117843865811, -0.68452091605681, 
          0.23456335139387, 0.162259631408099, 0.209619606465627, -0.10216429396362, 
          -0.629658878921399, 0.114133112372732, 0.823887580788133, 0.12141926111051, 
          0.113879127629599, 0.109816278251024, 0.221038616887962, 0.115707006557826, 
          0.187260599970159, 0.121830940397345, 0.172070355414403, 0.157444513096506, 
          0.254811666726125, 0.268571254537371, 0.215202234247305, 0.212970632033808, 
          0.262762514629277, 0.205440489731246, 0.0567461846032841, 0.154168532075271, 
          0.320248263514015, 0.309517737483193, 0.366585194306804, 0.370748971125027, 
          0.304859567470968, 0.397763842736319, 0.357894067104384, 0.380956131344983, 
          0.344676554711052, 0.37300484854814, 0.378382126329053, 0.342270280546076, 
          0.359489843015429),
        .Dim = c(32L, 2L),
        .Dimnames = list(
            c("ne.1.8115", "ne.1.8135", "ne.1.8117", "ne.1.8116", "ne.1.8111", "ne.1.8121", 
              "ne.1.8118", "ne.1.8136", "ne.1.8119", "ne.1.8125", "ne.1.8127", 
              "ne.1.8126", "ne.1.8128", "end.sin(2 * pi * t/52)", "end.cos(2 * pi * t/52)", 
              "end.ri(iid)", "neweights.d", "overdisp.region1", "overdisp.elsewhere", 
              "end.ri(iid).8115", "end.ri(iid).8135", "end.ri(iid).8117", "end.ri(iid).8116", 
              "end.ri(iid).8111", "end.ri(iid).8121", "end.ri(iid).8118", "end.ri(iid).8136", 
              "end.ri(iid).8119", "end.ri(iid).8125", "end.ri(iid).8127", "end.ri(iid).8126", 
              "end.ri(iid).8128"),
            c("Estimate", "Std. Error"))
        )
    
    expect_equal(coef(fluFitShared, se = TRUE), orig)
})

test_that("calibrationTest.oneStepAhead() works and \"final\" is equivalent to fit", {
    mysubset <- tail(fluFitShared$control$subset, 16)
    osa_final <- oneStepAhead(fluFitShared, tp = mysubset[1L]-1L,
                              type = "final", verbose = FALSE)
    idx <- 3:5  # ignore "method" and "data.name" in calibrationTest() output
    expect_equal(calibrationTest(osa_final, which = "dss")[idx], 
                 calibrationTest(fluFitShared, which = "dss", subset = mysubset)[idx])
})
