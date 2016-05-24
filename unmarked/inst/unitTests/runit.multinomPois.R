

test.removal <- function() {

    y <- matrix(c(
        5, 3, 2,
        3, 3, 1,
        2, 0, 0,
        0, 0, 0,
        0, 0, 0), nrow=5, ncol=3, byrow=TRUE)
    
    sc <- data.frame(x1 = c(NA, 2, 3, 4, 3))
    oc <- list(x2 = matrix(c(
        1, 1, 1,
        3, NA, 1,
        0, 0, 1,
        NA, NA, NA,
        NA, 1, 0), nrow=5, ncol=3, byrow=TRUE)) 
    
    umf1 <- unmarkedFrameMPois(y = y, siteCovs = sc, obsCovs = oc, 
        type="removal") 

    o2y <- diag(ncol(y))
    o2y[upper.tri(o2y)] <- 1
    checkEquals(obsToY(umf1), o2y)

    m1 <- multinomPois(~1 ~1, umf1)
    checkEqualsNumeric(coef(m1), c(1.5257743, -0.2328092), tol=1e-5)


    m2 <- multinomPois(~x2 ~1, umf1)
    checkEqualsNumeric(coef(m2), c(1.9159845, 0.2248897, -0.1808144), tol=1e-5)
    checkEquals(m2@sitesRemoved, 4:5)
    
    m3 <- multinomPois(~x2 ~x1, umf1)
    checkEqualsNumeric(m3@sitesRemoved, c(1, 4:5))
    checkEqualsNumeric(coef(m3), 
        c(1.9118525, -0.4071202, 8.3569943, 0.3232485), tol=1e-5)    
    
    }
    
    
    
test.double <- function() {
    y <- matrix(c(
        1, 0, 0,
        2, 1, 0,
        1, 0, 1,
        2, 1, 2,
        1, 0, 3,
        1, 1, 1), nrow=6, ncol=3, byrow=TRUE)
    oc <- matrix(c(
        1, 0,
        2, 1,
        1, 1,
        NA, 0,
        1, NA,
        NA, NA), nrow=6, ncol=2, byrow=TRUE)
        
    umf <- unmarkedFrameMPois(y = y, obsCovs = list(x=oc), type="double")    
    
    m1 <- multinomPois(~1 ~1, umf)
    checkEqualsNumeric(coef(m1), c(1.3137876, 0.2411609), tol=1e-5)
    
    m2 <- multinomPois(~x ~1, umf, starts=c(1.3, 0, 0.2))
    checkEquals(m2@sitesRemoved, 4:6)    
    }

      