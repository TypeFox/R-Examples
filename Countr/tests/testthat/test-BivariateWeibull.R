test_that("dWeibullInter ... reproduces football results results", {
    print("~~~~~~ Testing football results ~~~~~~~~")
    fn <- system.file("extdata", "estimationParams.RDS", package = "Countr")
    parObj <- readRDS(fn) # readRDS("estimationParams.RDS")

    ## extract relevant information
    data <- parObj$dataUsed
    par <- parObj$par
    forecastDate <- parObj$forecastDate
    ExpPar <- 0.0065

    ## Collect score and form the covariate data.frame
    formHome <- FTHG  ~ HomeTeam - 1
    formAway <- FTAG  ~ AwayTeam - 1
    HG <- as.numeric(data$FTHG)
    AG <- as.numeric(data$FTAG)

    ## Computing the likelihood function
    HomeAdv <- par["HomeAdv"]
    betaAtt <- par[grepl("Attack_", names(par))]
    betaDef <- par[grepl("Defense_", names(par))]

    XHome <- model.matrix(formHome, model.frame(formHome, data))
    XAway <- model.matrix(formAway, model.frame(formAway, data))

    lambdaHome <- exp(XHome %*% betaAtt + XAway %*% betaDef + HomeAdv)
    lambdaAway <- exp(XAway %*% betaAtt + XHome %*% betaDef)

    ## apply the weights
    weights <-exp(-ExpPar *
                  as.numeric(forecastDate -
                             as.Date(as.character(data$Date),
                                     format ="%d/%m/%y")
                             ) / 3.5
                  )

    gam_weiH <- par["gam_weiH"]
    gam_weiA <- par["gam_weiA"]
    theta <- par["theta"]

    ## --------------------------- Weibull count -------------------------------
    loglikweiFrank0 <-
        sum(
           dBivariateWeibullCountFrankCopula(HG, AG,
                                             gam_weiH, lambdaHome,
                                             gam_weiA, lambdaAway,
                                             theta, "series_acc",
                                             1, TRUE) * weights)
    loglikweiFrank1 <-
        sum(
            dBivariateWeibullCountFrankCopula(HG, AG,
                                              gam_weiH, lambdaHome,
                                              gam_weiA, lambdaAway,
                                              theta, "conv_dePril",
                                              1, TRUE, conv_extrap = T) * weights)
    loglikweiFrank2 <-
       dBivariateWeibullCountFrankCopula_loglik(HG, AG,
                                                gam_weiH, lambdaHome,
                                                gam_weiA, lambdaAway,
                                                theta, "conv_dePril",
                                                1, TRUE, conv_extrap = T,
                                                weights = weights)
    loglikweiFrank3 <-
        dBivariateWeibullCountFrankCopula_loglik(HG, AG,
                                                 gam_weiH, lambdaHome,
                                                 gam_weiA, lambdaAway,
                                                 theta, "series_acc",
                                                 1, TRUE, weights = weights)

    expect_equal(object = loglikweiFrank0, expected = parObj$logLik,
                 tolerance = 0.5)
    expect_equal(object = loglikweiFrank1, expected = parObj$logLik,
                 tolerance = 0.5)
    expect_equal(object = loglikweiFrank2, expected = parObj$logLik,
                 tolerance = 0.5)
    expect_equal(object = loglikweiFrank3, expected = parObj$logLik,
                 tolerance = 0.5)

})
