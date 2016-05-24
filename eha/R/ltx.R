ltx <- function(x,
                caption = NULL,
                label = NULL,
                dr = NULL,
                digits = max(options()$digits - 4, 3), ...){
    UseMethod("ltx")
}
