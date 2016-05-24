run.cec.tests <- function()
{
    errors <- 0
    
    # test files
    tests <- list.files(system.file("cec_tests", package="CEC"))
    for (test in tests)
    {
        if (grepl(".R", test, perl=T))
        {
            testenv <- new.env()
            local(
                {
                    # just to trick R CMD check...
                    testname <- NULL
                    setup <- NULL
                },
                testenv
            )      
            source(system.file("cec_tests", test, package="CEC"), local=testenv)
            errors <- errors + local(
                {
                    local.errors <- 0
                    cat(paste("Test:",testname, "\n"))        
                    fs <- lsf.str()        
                    
                    # execute setup function if exists
                    if ("setup" %in% fs) eval(expr=body(setup), envir=testenv)        
                    for (fn in fs) 
                    {
                        # test cases 
                        if (grepl("test.", fn))
                        {      
                            cat(paste("---- ", fn))              
                            fbody = body(eval(parse(text=fn)))
                            
                            # evaluate test case function and catch (and count) errors
                            local.errors <- local.errors + tryCatch(
                                {                  
                                    eval(expr=fbody, envir=testenv)
                                    cat(": OK\n")
                                    0
                                },
                                error = function(er) {                  
                                    cat(": FAILED\n")
                                    warning(er$message, immediate.=T, call.=F)
                                    1
                                }
                            )             
                        }
                    }
                    local.errors
                },
                envir=testenv)
        }
    }
    if (errors > 0)
    {
        stop("One or more tests failed.")
    }
}
printmsg <- function(msg)
{    
    if (!is.null(msg))
        paste(msg, ":")
    else ""      
}

checkNumericVectorEquals <- function(ex, ac, msg=NULL, tolerance = .Machine$double.eps ^ 0.5)  
{  
    if (length(ex) != length(ac)) stop (paste(printmsg(msg),"Vectors have different length."))  
    for (i in seq(1, length(ex)))      
        if (!isTRUE(all.equal.numeric(ex[i], ac[i], tolerance=tolerance)))
            stop (paste(printmsg(msg), "Vectors differ at index:", i, ", expected:", ex[i], ", actuall:",ac[i]))
    
}

checkNumericEquals <- function(ex, ac, msg=NULL, tolerance = .Machine$double.eps ^ 0.5)
{
    if(!is.numeric(ex)) stop(paste(printmsg(msg), "Expression:",ex,"is not numeric type."))
    if(!is.numeric(ac)) stop(paste(printmsg(msg), "Expression:",ac,"is not numeric type."))
    if (!isTRUE(all.equal.numeric(ex, ac, tolerance=tolerance)))
        stop (paste(printmsg(msg), "Numeric values are different: expected:", ex, ", actuall:",ac, ", difference:", abs(ex - ac)))
}

checkEquals <- function(ex, ac, msg=NULL)
{
    if (!isTRUE(identical(ex, ac)))
        stop (paste(printmsg(msg), "Values are not identical: expected:", ex, ", actuall:",ac))
}

checkTrue <- function(exp, msg=NULL)
{
    if (!is.logical(exp))
    {
        stop(paste(printmsg(msg), "Expression has not logical type."))
    }
    if (!isTRUE(exp))
    {
        stop(paste(printmsg(msg), "Expression is not TRUE."))
    }
}

checkNumericMatrixEquals <- function(ex, ac, msg=NULL, tolerance = .Machine$double.eps ^ 0.5)  
{
    
    if (nrow(ex) != nrow(ac)) stop (paste(printmsg(msg),"Matrices have different dimensions."))
    if (ncol(ex) != ncol(ac)) stop (paste(printmsg(msg),"Matrices have different dimensions."))
    
    for (i in seq(1, nrow(ex)))  
        for (j in seq(1, ncol(ex)))    
            if (!isTRUE(all.equal.numeric(ex[i, j], ac[i, j], tolerance=tolerance)))
                stop (paste(printmsg(msg), "Matrices differ at row:", i, " col:", j, ": expected:", ex[i, j], ", actuall:",ac[i, j]))
    
}

# Maximum likelihood estimate of covariance matrix.
cov.mle <- function(M)
{
    mean <- colMeans(M)
    mat <- matrix(0,ncol(M),ncol(M))
    for (i in seq(1, nrow(M)))
    {
        v <- M[i,]   
        mat <- mat + (t(t(v - mean)) %*% t(v - mean))
    }
    mat <- mat / nrow(M)
    mat
}
