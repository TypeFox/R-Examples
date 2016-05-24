##
## Function for creating 'CTRL' objects 
ctrl <- function(maxiters = 100L, abstol = 1e-7, reltol = 1e-6, feastol = 1e-7,
                 stepadj = 0.95, beta = 0.5, trace = TRUE){

    if(!is.integer(maxiters)){
        stop("\nThe count of maximal iterations must be an integer.\n")
    }
    if(maxiters < 1){
        stop("\nThe count of maximal iterations must be positive and greater or equal to one.\n")
    }
    if(!is.null(dim(abstol)) | length(abstol) > 1){
        stop("\nThe absolute tolerance for convergence must be a real scalar.\n")
    }
    if(!is.null(dim(reltol)) | length(reltol) > 1){
        stop("\nThe relative tolerance for convergence must be a real scalar.\n")
    }
    if(!is.null(dim(feastol)) | length(feastol) > 1){
        stop("\nThe feasabile tolerance for convergence must be a real scalar.\n")
    }
    if(abstol < 0 & reltol < 0){
        stop("\nAt least one of 'reltol' and 'abstol' must be positive.\n")
    }
    if(feastol <= 0){
        stop("\nThe convergence criteria for feasability must be positive.\n")
    }
    if(stepadj <= 0 || stepadj > 1.0){
        stop("\nStep-size adjustment must be in the interval: (0, 1].\n")
    }
    if(beta <= 0 || beta >= 1.0){
        stop("\nBacktracking parameter for domain of non-linear constraints\nmust be in the interval: (0, 1).\n")
    }
    
    new(CTRL, list(
        maxiters = maxiters,
        abstol = abstol,
        reltol = reltol,
        feastol = feastol,
        stepadj = stepadj,
        beta = beta,
        trace = as.logical(trace)[1])
        )
}
