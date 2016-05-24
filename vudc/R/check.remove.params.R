check.remove.params <-
function (remove.absolute, remove.ratio) 
{
    if (!is.na(remove.absolute) && !is.na(remove.ratio)) {
        stop("At least one of remove.absolute and remove.ratio should be NA.")
    }
    if (!is.na(remove.ratio) && (remove.ratio < 0 || remove.ratio >= 
        0.5)) {
        stop("The remove.ratio should take value in range [0, 0.5).")
    }
    if (!is.na(remove.absolute) && (remove.absolute <= 0)) {
        stop("The remove.absolute should take value greater than 0.0.")
    }
}
