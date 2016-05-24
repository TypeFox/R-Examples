
# Get and set global 'gridSVG' options

# Initial settings
assign("gridSVGoptions",
       list(
           id.sep=".",
           gPath.sep="::",
           vpPath.sep="::"
           ),
       .gridSVGEnv)

checkOptions <- function(options) {
    optionNames <- names(options)
    validOption <- sapply(options, is.character)
    if (any(! validOption))
        stop(paste("Invalid option for: ",
                   paste(dQuote(optionNames[! validOption]), collapse = ", "),
                   sep = ""))
}

# Get/set options
getSVGoption <- function(name) {
    oldOptions <- get("gridSVGoptions", .gridSVGEnv)
    optionNames <- names(oldOptions)
    if (name %in% optionNames) {
        oldOptions[[name]]
    }
}

getSVGoptions <- function() {
    get("gridSVGoptions", .gridSVGEnv)
}

setSVGoptions <- function(...) {
    oldOptions <- get("gridSVGoptions", .gridSVGEnv)
    options <- list(...)
    if (length(options)) {
        names <- names(options)
        optionNames <- names(oldOptions)
        names <- names[nchar(names) > 0 &
                       names %in% optionNames]
        if (length(options[names])) {
            newOptions <- oldOptions
            newOptions[names] <- options[names]
            checkOptions(newOptions)
            assign("gridSVGoptions", newOptions, .gridSVGEnv)
            invisible(oldOptions[names])
        } 
    } 
}
