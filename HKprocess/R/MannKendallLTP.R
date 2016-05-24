MannKendallLTP <- function(data) {
    
    nx <- as.integer(length(data)) # length of data
    nx_double_sum <- as.integer((nx*(nx-1))/2)
    
    # The equation numbers correspond to equation numbers in Hamed (2008) unless
    # mentioned otherwise

    # Call the score.c from the C library
    out1<-.C("score",nx,data,tr = as.integer(array(0,dim = c(1,3))),PACKAGE =
    "HKprocess")
    
    S <- (out1$tr)[1] # Eq(1)
    V0Sminus <- (out1$tr)[2]
    denom_ties <- (out1$tr)[3]

    rm(out1)
    
    V0S <- (nx * (nx - 1) * (2*nx + 5))/18 - V0Sminus/18 # Eq(4)
    u <- ifelse(S > 0,(S - 1)/sqrt(V0S),ifelse(S==0,0,(S + 1)/sqrt(V0S))) #Eq(6)
    pvalue_mann_kendall <- 2 * pnorm(-abs(u))
    
    rm(u)
    
    denominator <- sqrt(nx_double_sum - 0.5 * denom_ties)* sqrt(nx_double_sum)
    # The denominator is calculated according to eq(23.3.4) in Hipel and
    # McLeod (1994).
    
    tau <- S/denominator
    
    # Call the score0.c from the C library
    out2<-.C("score0",nx,data,tr = array(0,dim = c(1,nx_double_sum)),PACKAGE =
    "HKprocess")
    
    s0 <- median(out2$tr) # Eq(17)

    rm(out2)
    
    y <- data - s0 * (1:nx)
    ranky <- rank(y,ties.method = c("average"))
    z <- qnorm(ranky/(nx + 1)) # Eq(18)
    Hest <- mleHK(z)[3] # Eq(21) but it uses the ML estimator in Tyralis and
    # Koutsoyiannis (2011) instead.
    
    mHsign <- 0.5 - 2.87 * (nx^(-0.9067)) # Eq(22)
    sHsign <- 0.77654 * (nx^(-0.5)) - 0.0062 # Eq(23)
    pvalue_H <- 2 * pnorm(-abs((Hest - mHsign)/sHsign))
    
    a0 <- (1.0024 * nx - 2.5681)/(nx + 18.6693)
    a1 <- (-2.2510 * nx + 157.2075)/(nx + 9.2245)
    a2 <- (15.3402 * nx - 188.6140)/(nx + 5.8917)
    a3 <- (-31.4258 * nx + 549.8599)/(nx - 1.1040)
    a4 <- (20.7988 * nx - 419.0402)/(nx - 1.9248)
    B <- a0 + a1 * Hest + a2 * Hest^2 + a3 * Hest^3 + a4 * Hest^4 # Eq(24)
    
    # Call the VstarSfunction.c from the C library
    out3<-.C("VstarSfunction",nx,Hest,tr = array(0,dim = c(1,1)),PACKAGE =
                "HKprocess")
    V_star_S <- (c((out3$tr)[1]))
    
    rm(out3)
    
    VS <- B * V_star_S
    u <- ifelse(S > 0,(S - 1)/sqrt(VS),ifelse(S==0,0,(S + 1)/sqrt(VS))) # Eq(6)
    pvalue_mann_kendall_ltp <- 2 * pnorm(-abs(u))
    
    rm(u)
    
    Mann_Kendall <- setNames(c(tau,S,V0S,s0,denominator,pvalue_mann_kendall),
    c("Kendall_s_tau_statistic","Score","V0Score","Sen_slope","denominator",
    "2_sided_pvalue"))
    
    H_significance <- setNames(c(Hest,pvalue_H),c("Hest","2_sided_pvalue"))
    Mann_Kendall_LTP <- setNames(c(VS,pvalue_mann_kendall_ltp),c("VScore",
    "2_sided_pvalue"))

    results_names <- c("Mann_Kendall","Significance_of_H","Mann_Kendall_LTP")
    results <- list(Mann_Kendall,H_significance,Mann_Kendall_LTP)
    names(results) <- results_names
    
    return(results)
}