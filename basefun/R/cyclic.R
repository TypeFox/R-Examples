
### <FIXME> deriv? </FIXME>
cyclic_basis <- function(order, frequency, varname) {
    S <- 1:order * 2
    fm <- as.formula(paste("~ I(", rep(c("sin", "cos"), length(S)), "(",
                     rep(S, rep(2, length(S))), "* pi *", varname, "/", 
                     frequency, "))", collapse = "+"))
    return(as.basis(fm, remove_intercept = TRUE))
}
