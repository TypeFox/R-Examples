"print.power.grouped" <-
function (x, ...) {
    cat("\n Two-sample Wald's test power calculation \n \t for grouped data \n")
    cat(paste(format(names(x), width = 15, justify = "right"), format(x), sep = " = "), sep = "\n")
    cat("\n")
}

