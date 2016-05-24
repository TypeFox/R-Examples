################################################################################
# Helper function: trim output values to the requested verbosity
trim.output <- function(res, verbosity = "terse"){

    # Check for a valid verbosity argument; if not valid, 
    # default to terse and warn
    if(!verbosity %in% c("onechar", "terse", "verbose")){
        warning("'", verbosity, "' is not a valid choice for result verbosity; defaulting to 'terse'.")
        verbosity <- "terse"
    }

    
    # If verbosity is "onechar", output just the first character of result
    if(verbosity %in% "onechar"){substr(res, 1, 1)} else

        # If verbosity is "terse", just the first word
        if(verbosity %in% "terse"){gsub(x = res, 
                                        pattern = " .*$", 
                                        replacement = "")} else 

            # If verbose, indeterminates indicate which criteria
            if(verbosity %in% "verbose"){res}
}

