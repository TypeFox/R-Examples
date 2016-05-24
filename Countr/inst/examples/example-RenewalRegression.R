test_that("test renewal -- McShane results", {
    print("~~~~~~ renewal regression-- McShane results ~~~~~~~~")
                # res <- readRDS("McShane_paperResults.RDS")
    fn <- system.file("extdata", "McShane_paperResults.RDS", package = "Countr")
    res <- readRDS(fn)

    y <- res$y
    data <- res$data
    form <-
        Y ~ GERMAN + EDU + VOC + UNI + CATH + PROT + MUSL + RURAL + YEAR_OF_B + AGEMARR
    ## =========================== gamma =====================================
    print("............ gamma ............")
    res <- renewal(formula = form, data = data, dist = "gamma",
                   computeHessian = FALSE,
                   control = renewal.control(trace = 0)
                   )
    ll <-  as.numeric(logLik(res))
    expect_less_than(abs(ll - (-2078)), 0.5)
})

test_that("test renewal -- McShane data --- prediction", {
    print("~~~~~~ renewal prediction-- McShane results ~~~~~~~~")
                # res <- readRDS("McShane_paperResults.RDS")
    fn <- system.file("extdata", "McShane_paperResults.RDS", package = "Countr")
    res <- readRDS(fn)

    y <- res$y
    data <- res$data
    form <-
        Y ~ GERMAN + EDU + VOC + UNI + CATH + PROT + MUSL + RURAL + YEAR_OF_B + AGEMARR
    ## =========================== weibull =====================================
    object <- renewal(formula = form, data = data, dist = "weibull",
                      computeHessian = TRUE, weiMethod = "series_acc",
                      control = renewal.control(trace = 0)
                      )

    predOld.response <- predict(object, type = "response", se.fit = TRUE)
    predOld.prob <- predict(object, type = "prob", se.fit = TRUE)

    ## newData (extracted from old Data)
    newData <- head(data)
    predNew.response <- predict(object, newdata = newData,
                                type = "response", se.fit = TRUE)
    predNew.prob <- predict(object, newdata = newData,
                            type = "prob", se.fit = TRUE)

    expect_equal(head(predOld.response$values),
                 predNew.response$values,
                 tolerance = 1e-3)

    expect_equal(head(predOld.response$se$scale),
                 predNew.response$se$scale,
                 tolerance = 1e-3)

    expect_equal(head(predOld.response$se$shape),
                 predNew.response$se$shape,
                 tolerance = 1e-3)

    expect_equal(head(predOld.prob$values),
                 predNew.prob$values,
                 tolerance = 1e-3)

    expect_equal(head(predOld.prob$se$scale),
                 predNew.prob$se$scale,
                 tolerance = 1e-3)

    expect_equal(head(predOld.prob$se$shape),
                 predNew.prob$se$shape,
                 tolerance = 1e-3)
})


## -------------------------- very slow due to conversion from R to C++
## test_that("test renewal -- McShane results -- user passed", {
##     print("~~~~~~ renewal regression-- McShane results user passed ~~~~~~~~")
##     res <- readRDS("McShane_paperResults.RDS")
##     y <- res$y
##     data <- res$data
##     form <-
##         Y ~ GERMAN + EDU + VOC + UNI + CATH + PROT + MUSL + RURAL + YEAR_OF_B + AGEMARR

##     ## =========================== weibull =====================================
##     print("............ weibull ............")
##     parNames <- c("scale", "shape")
##     sWei <- function(tt, distP) {
##         exp( -distP[["scale"]] * tt ^ distP[["shape"]])
##     }

##     .getExtrapol <- function(distP) {
##         c(2, distP[["shape"]])
##     }

##     customPars <- list(parNames = parNames,
##                        survivalFct = sWei,
##                        extrapolFct = .getExtrapol)

##     link <- list(scale = "log", shape = "log")
##     ## starting values
##     par0 <- coef(glm(form, family = "poisson", data = data))
##     names(par0) <- paste0("scale_", names(par0))
##     start <- c(par0, shape_ = 1)

##     ## control Parameters
##     control <- renewal.control(start = start, trace = 0)

##     res <- renewal(formula = form, data = data, dist = "custom", link = link,
##                    control = control, customPars = customPars,
##                    computeHessian = FALSE)

##     ll <-  as.numeric(logLik(res))
##     expect_less_than(abs(ll - (-2077)), 0.1)
## })
