Chi2Dist <-
function(X){#compute the chi2 distance 
	# between the rows of a matrix
    # send back the distance and the vector of mass	
    
    # 1. transform X into Row profiles
    xip = rowSums(X)
    R <- X / xip
    # 2. Masses, weights
    xpp = sum(xip)             # grand total
    m <- xip / xpp             # masses 
    c <- colSums(X) / xpp      # row barycenter
    w = 1/c                    # columns weights
    # Preprocess R
    Rc = t(t(R) - c)           # deviations to barycenter
    Rtilde = t(t(Rc)*sqrt(w))  # weighted R
    S = Rtilde%*%t(Rtilde)     # covariance
    s =diag(S) # diag of
    D = (s - S) + t(s-S)       # Chi2 distance matrix
    return(list(Distance=D, masses=m))
}
