TIK <- function(Y, U = U, k = 1)
    {
    RES <- .Call( "TIK", Y, U, k, PACKAGE = "tsBSS")
    RES$Tik
    }


