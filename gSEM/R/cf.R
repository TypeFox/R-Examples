# Format a number for better printing
# Print "- 0.3" instead of "+ -0.3"

cf <- function(x, noPlus = FALSE){
    if(noPlus)
        ifelse(x > 0, paste0(x), paste0("-", -x))
    else
        ifelse(x > 0, paste0(" + ", x), paste0(" - ", -x))
}
