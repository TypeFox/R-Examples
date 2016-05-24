eHOF.modelnames <- factor(c("I", "II", "III", "IV", "V", "VI", "VII"), ordered = TRUE)

HOF.fun <- function (x, model, p, M = 1, ...) {
    model <- match.arg(model, eHOF.modelnames)
    a <- p[1]; b <- p[2]; c <- p[3]; d <- p[4]; f <- p[5]
    #v <- list(M = M, a = p[1], b = p[2], c = p[3], d = p[4], e = p[5])
    x <- scale01(x, ...)
    if (length(M) == 1)  M <- rep(M, length(x))
    fv <- switch(as.character(model), 
           I = M/(1 + exp(a)),
          II = M/(1 + exp(a + b * x)),
         III = M/(1 + exp(a + b * x)) * 1/(1 + exp(c)),
          IV = M/(1 + exp(a + b * x)) * 1/(1 + exp(c - b * x)),
           V = M/(1 + exp(a + b * x)) * 1/(1 + exp(c - d * x)),
          VI = M/(1 + exp(a + b * x)) * 1/(1 + exp(c - b * x)) + M/(1 + exp(a + b * (x - d))) * 1/(1 + exp(c - b * (x - d))),
         VII = M/(1 + exp(a + b * x)) * 1/(1 + exp(c - b * x)) + M/(1 + exp(a + b * (x - d))) * 1/(1 + exp(c - f * (x - d)))
     )
    return(fv)
}

