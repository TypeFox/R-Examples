


.fonctionMCMC <- function(data, x) {

return(info.nests(parameters=x, temperatures=data$data, 
                    derivate=data$derivate, weight=data$weight,
                    test=data$test, M0=data$M0, out="Likelihood", 
                    fixed.parameters=data$fixed.parameters, progress=FALSE, warnings=FALSE))

}
