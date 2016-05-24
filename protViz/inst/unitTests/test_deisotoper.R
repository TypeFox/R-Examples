#R


test_deisotoper <-
function(){
    x <- list(mZ = c(726.068, 726.337, 726.589, 726.842, 727.343, 727.846, 728.346, 
        728.846, 729.348, 730.248, 730.336, 730.581, 730.836),
        intensity = c(6.77850e+03, 2.81688e+04, 6.66884e+04, 1.22032e+07, 
        9.90405e+06, 4.61409e+06, 1.50973e+06, 3.33996e+05, 5.09421e+04, 
        1.15869e+03, 2.14788e+03, 5.37853e+03, 5.79094e+02))


    out1 <- .Call("deisotoper_main", x$mZ, x$intensity, Z=1:4, averagine, 
        massError=0.01, PACKAGE="protViz")

    checkEqualsNumeric(out1$result[[1]][[1]], c(1, 4, 6, 8), tolerance=1.0e-6)
    checkEqualsNumeric(out1$result[[2]][[1]], c(1, 3, 4, 5, 6, 7, 8), tolerance=1.0e-6)
}

test_deisotoper()
