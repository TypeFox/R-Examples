### Test for whole number, with tolerance for representation
### From post by Tony Plate <tplate_at_acm.org>
is.wholenumber <- function(x, tolerance = .Machine$double.eps^0.5){
    if (!is.numeric(x)){
        return(FALSE)
    } else {
        return(isTRUE(all(abs(x - round(x)) < tolerance)))
    }
}
