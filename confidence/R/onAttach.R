.onAttach <-
function(libname, pkgname) {
    packageStartupMessage(
        "\nFor more information on the usage of the Confidence tool, type:\n", 
        'vignette("confidence")\n',
        "Examples of input files can be found in:\n", 
        sQuote(system.file("extdata", package = "confidence"))
    )
}  
