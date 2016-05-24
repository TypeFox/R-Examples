Hals.ses <-
function(Z, Hs, Ht, Hst.ls, rho, reg, b.lag, test.rng) {
    tau <- nrow(Z)
    xALS <- H.als.b(Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, rho=rho, reg=reg, b.lag=b.lag, Hs0=NULL, Ht0=NULL, Hst0.ls=NULL)
    rmse <- sqrt( mean( ( Z[ test.rng, ] - xALS$Z.hat[ test.rng, ] )^2 ) ) ; rmse
    als.se <- rmse * sqrt( xALS$ALS.g ) * sqrt( diag( xALS$inv.LHH ) )
    return( list( "estimates"=cbind( xALS$B[ tau, ], als.se ), "inv.LHH"=xALS$inv.LHH, "ALS.g"=xALS$ALS.g ) )
}
