CCK <- function(Y, k = 0:12)
    {
    RES <- .Call( "CCK", Y, k, PACKAGE = "tsBSS")
    RES$CCK
    }


