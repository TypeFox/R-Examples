prepanel.default.bwplot2 <-
function (x, y, horizontal = TRUE, nlevels, 
                                      origin = NULL, stack = TRUE,  ...) 
{
   if (any(!is.na(x) & !is.na(y))) {
        if (horizontal) {
            if (!is.factor(y)) {
                if (missing(nlevels)) 
                  nlevels <- length(unique(y))
                y <- factor(y, levels = 1:nlevels)
            }
            list(xlim = if (stack) {
                foo1 <- if (any(x[!is.na(x)] > 0)) range(tapply(x[x > 0], 
                  y[x > 0, drop = TRUE], sum, na.rm = TRUE), 
                  finite = TRUE) else 0
                foo2 <- if (any(x[!is.na(x)] < 0)) range(tapply(x[x < 0], 
                  y[x < 0, drop = TRUE], sum, na.rm = TRUE), 
                  finite = TRUE) else 0
                range(foo1, foo2)
            } else lattice:::scale.limits(c(x, origin)), ylim = levels(y), 
                yat = sort(unique(as.numeric(y))), dx = 1, dy = 1)
        }
        else {
            if (!is.factor(x)) {
                if (missing(nlevels)) 
                  nlevels <- length(unique(x))
                x <- factor(x, levels = 1:nlevels)
            }
            list(xlim = levels(x), xat = sort(unique(as.numeric(x))), 
                ylim = if (stack) {
                  foo1 <- if (any(y > 0)) range(tapply(y[y > 
                    0], x[y > 0], sum, na.rm = TRUE), finite = TRUE) else 0
                  foo2 <- if (any(y < 0)) range(tapply(y[y < 
                    0], x[y < 0], sum, na.rm = TRUE), finite = TRUE) else 0
                  range(foo1, foo2)
                } else lattice:::scale.limits(c(y, origin)), dx = 1, dy = 1)
        }
    }
    else lattice:::prepanel.null()
}
