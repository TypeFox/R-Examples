dilation <- function(y,span)
  {
    y <- na.omit(y)
    nmes <- names(y)
    res <- .C("dilation",
              y = as.double(y),
              nym = as.integer(length(y)),
              sp = as.integer(span),
              ey = double(length(y)),
              PACKAGE = "VPdtw")$ey
    res <- rev(.C("dilation",
                  y = as.double(rev(res)),
                  nym = as.integer(length(res)),
                  sp = as.integer(span),
                  ey = double(length(res)),
                  PACKAGE = "VPdtw")$ey)
    names(res) <- nmes
    res
  }

