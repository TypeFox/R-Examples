.onAttach <-
function(libname, pkgname) {
    javaVersion <- utils::packageDescription(pkg = pkgname, fields = "SystemRequirements") 
    packageStartupMessage("Note: ", sQuote(pkgname), " requires ",
        javaVersion, ", which is available at www.java.com")
}  
