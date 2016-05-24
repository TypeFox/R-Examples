bool <-
function(x=T, na=NA) {
    b <- (x!=0)
    ifelse(is.na(b), na, b)
}
