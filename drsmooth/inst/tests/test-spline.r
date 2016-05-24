#' @name test_spline
#' @title Test spline function

test_that("Spline Correct", {
    
    # test scripts do not seem to recognize loaded libraries:
#     library(car)
#     library(clinfun)
#     library(multcomp)
#     library(pgirmess)
#     library(DTK)
#     library(mgcv)
#     library(segmented)

    dose <- DRdata[,"dose"]
    targetvariable <- DRdata[,"MF_Log"]
    spline <- mgcv::gam(targetvariable~s(dose, k=4), data=DRdata)
    
    expect_that (round(spline$coefficients[1], digits = 6), is_equivalent_to(2.142013))
    
    }
)