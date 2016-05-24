boundary.fn <-
function(distribution="bernoulli", type1=0.05, type2=0.20, h0, h1) {
    
    # Wald boundaries
    boundary.A <- waldBoundary(type1 = type1, type2 = type2, boundary = "A", log = TRUE)
    boundary.B <- waldBoundary(type1 = type1, type2 = type2, boundary = "B", log = TRUE)
    
    # Intercept
    h1.intercept <- boundary.A / (C.fn(distribution, h1) - C.fn(distribution, h0))
    h0.intercept <- boundary.B / (C.fn(distribution, h1) - C.fn(distribution, h0))
    
    # Slope
    slope <- (D.fn(distribution, h1) - D.fn(distribution, h0)) / (C.fn(distribution, h1) - C.fn(distribution, h0))
    
    # Output matrix
    matrix(data = c(h0.intercept, slope, h1.intercept, slope),
           nrow = 2,
           byrow = TRUE,
           dimnames = list(c("h0", "h1"),
                           c("intercept", "slope")))
}
