getCutpointList <-
function(dataset, base.form, covs){
    cutpoints <- lapply(covs, function(cov)
        getCutpoint(dataset, base.form, cov)
                        )
    names(cutpoints) <- covs
    return(cutpoints)
}
