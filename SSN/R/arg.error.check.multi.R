arg.error.check.multi <- function(CorModels = c("Exponential.taildown"),
                                  use.nugget = TRUE,
                                  use.anisotropy = FALSE,
                                  addfunccol = NULL,
                                  family = "gaussian",
                                  EstMeth = "REML",
                                  ssn = ssn)
{
    if(!is.logical(use.nugget))
		return(list(Err = 1, message =
			"use.nugget argument must be TRUE or FALSE"))
    if(!is.logical(use.anisotropy))
        return(list(Err = 1, message =
			"use.anisotropy argument must be TRUE or FALSE"))
    if(EstMeth != "ML" & EstMeth != "REML")
        return(list(Err = 1, message =
                    "EstMeth must be ML (Max. Likelihood) or REML (Restricted Max. Likelihood)"))
    if(!(family == "gaussian" | family == "binomial" |
         family == "poisson"))
        return(list(Err = 1, message =
			"family must be gaussian or binomial or poisson"))

    Models <- c("Exponential.taildown","Spherical.taildown",
		"LinearSill.taildown","Mariah.taildown","Epanech.taildown",
		"Exponential.tailup","Spherical.tailup",
		"LinearSill.tailup","Mariah.tailup", "Epanech.tailup",
		"Exponential.Euclid","Spherical.Euclid","Gaussian.Euclid",
		"Cauchy.Euclid")

    ## Add variables that are factors in case a random effects model is specified by user
    cl <- sapply(ssn@obspoints@SSNPoints[[1]]@point.data,class)
    ind <- which(cl=="factor")
    if(length(ind)) Models <- c(Models,names(cl[ind]))

    ## check to make sure we have necessary arguments for stream and
    ## Euclidean distance models
    if(!(is.null(CorModels) || is.na(CorModels))) {
        if(length(setdiff(CorModels,Models))) {
            return(list(Err = 1, message = "Not a valid autocorrelation model model for CorModels"))
        }
        if(length(grep("tailup",CorModels)) & is.null(addfunccol))
            return(list(Err = 1, message = "Need an additive function column for tailup models"))
    }

    Err <- list(Err = 0)
    Err
}

