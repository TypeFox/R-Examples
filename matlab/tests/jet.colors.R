###
### $Id: jet.colors.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.jet.colors <- function(input, expected) {
    output <- do.call(getFromNamespace("jet.colors", "matlab"), input)
    identical(output, expected)
}

jet.expected.m0 <- character(0)
jet.expected.m1 <- "#00FFFF"
jet.expected.m8 <- c("#0000FF",
                     "#0080FF",
                     "#00FFFF",
                     "#80FF80",
                     "#FFFF00",
                     "#FF8000",
                     "#FF0000",
                     "#800000")

test.jet.colors(list(n = 0), jet.expected.m0)
test.jet.colors(list(n = 1), jet.expected.m1)
test.jet.colors(list(n = 8), jet.expected.m8)

