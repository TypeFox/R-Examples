## frair_fit
# Wrapper function to fit functional response curves
# The one of the core functions of the frair package
frair_fit <- function(formula, data, response, start=list(), fixed=NULL){
    # Parse call, can check formula...
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf_list <- as.list(mf)
    expandmod <- terms(formula(mf_list$formula), data=data)
    expandform <- formula(expandmod)
    leftside <- all.vars(expandform[[2]])
    rightside <- all.vars(expandform[[3]])
    if(length(leftside)!=1 || length(rightside)!=1) {
        stop('Only simple formulae (e.g. y ~ x) are supported.')
    }
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    # moddata is the data we need for the fitting...
    moddata <- eval.parent(mf)
    names(moddata) <- c('Y', 'X')
    
    # Plausibly, 'response' might be the function itself (if the user hasen't provided a quoted string as they should have). Rather than bitch about it, we deal.
    if(!is.character(response)){
        response <- as.character(mf_list$response)
    }
   
    # Check we can deal with the requested response
    resp_known <- names(frair_responses(show=FALSE))
    resp_check <- match(response, resp_known, 0L)
    if(resp_check==0){
        stop(paste0(deparse(substitute(response)), ' is not a recognised response.\n   Use frair_responses(show=T) to see what has been implemented.'))
    }
    
    # Check start
    fr_checkstart(start, deparse(substitute(start)))
    
    # Check fixed
    if(!is.null(fixed)){
        if(!is.list(fixed) | is.null(names(fixed))){
            stop(paste0(deparse(substitute(fixed)), " must be a list containing single, named numeric values."))
        }
        if(length(fixed)>0 && any(lapply(fixed, length)!=1)){
            stop(paste0("The items in ", deparse(substitute(fixed)), " must be single, named numeric values."))
        }
        if(length(fixed)>0 && !(all(is.numeric(unlist(fixed))))){
            stop(paste0("The items in ", deparse(substitute(fixed)), " must be single, named numeric values."))
        }
    }
    
    # Check we have everything we need
    req_input <- names(formals(response))
    req_input  <- req_input[req_input!='X']
    input_matches <- match(req_input, c(names(start), names(fixed)), NA)
    if(any(is.na(input_matches))){
        missing_input <- req_input[is.na(input_matches)]
        if(length(missing_input)>1){
            stop(paste("Your requested response function requires input: ", paste(req_input, collapse=', '), ".\n",
                       "   The following items are missing: ", paste(missing_input, collapse=', '), ".\n",
                       "   Please provide them via 'start' or 'fixed', as appropriate.", sep=''))
        } else {
            stop(paste("Your requested response function requires input: ", paste(req_input, collapse=', '), ".\n",
                       "   The following item is missing: ", missing_input, ".\n",
                       "   Please provide them via 'start' or 'fixed', as appropriate.", sep=''))
        }
    }
    # Check we don't have superfluious input (bbmle::mle2 will bitch it's arse off otherwise [and rightly too!])
    input_matches <- match(c(names(start), names(fixed)), req_input, NA) # Inverse of above matching
    if(any(is.na(input_matches))){
        missing_input <- c(names(start), names(fixed))[is.na(input_matches)]
        if(length(missing_input)>1){
            stop(paste("Your requested response function requires input: ", paste(req_input, collapse=', '), ".\n",
                       "   The following items are not needed: ", paste(missing_input, collapse=', '), ".\n",
                       "   Please remove them from 'start' or 'fixed'.", sep=''))
        } else {
            stop(paste("Your requested response function requires input: ", paste(req_input, collapse=', '), ".\n",
                       "   The following item is not needed: ", missing_input, ".\n",
                       "   Please remove it from 'start' or 'fixed'.", sep=''))
        }
    }
    
    ## Go time!
    # Setup output
    out <- list('call' = call, 'x' = moddata$X, 'y'=moddata$Y, 'response'=response, 'xvar' = rightside, 'yvar' = leftside, optimvars=names(start), fixedvars=names(fixed))
    class(out) <- c('frfit', class(out))
    # In this instance, the sample is just the data itself...
    samp=c(1:nrow(moddata))
    ## Case specific fitting...
    frfunc <- get(unlist(frair_responses(show=FALSE)[[response]])[1], pos = "package:frair")
    frout <- frfunc(data=moddata, samp=c(1:nrow(moddata)), start=start, fixed=fixed)
    ## End case specific fitting...
    
    # Add the coefficients to the output
    out[['coefficients']] <- frout$out[req_input]
    # Add the sample, in this case, just a simple vector
    out[['sample']] <- samp
    # Finally, add the fit to the output (in case people need it)
    out[['fit']] <- frout$fit
    # Finally return our object!
    return(out)
}