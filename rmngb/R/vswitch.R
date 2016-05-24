vswitch <- function(x, ...)
    UseMethod("vswitch")

vswitch.factor <- function(x, ...) {
    x <- levels(x)[x]
    
    vswitch(x, ...)
}

vswitch.character <- function(x, ...) {
    listRes <- list(...)
    nRes <- names(listRes)
    if (is.null(nRes))
        stop("'...' should be named when 'x' is a character vector.")
    if (any(!is.na(x) & x %out% nRes))
        stop("Some values of 'x' do not correspond to any name in '...'.")
    x <- match(x, nRes)
    names(listRes) <- NULL
    
    do.call(vswitch, c(list(x = x), listRes))
}

vswitch.default <- function(x, ...) {
    listRes <- lapply(list(...),
                      function(x) rep(x, length.out = length(x)))
    if (! is.null(names(listRes)))
        warning("Named '...' with non character/factor 'x'. Names will be ignored.")
    if (any(!is.na(x) & (x > length(listRes) | x < 1)))
        stop("Some values of 'x' are out of range.")
    
    tabRes <- as.data.frame(listRes)
    iRows <- seq(length.out = nrow(tabRes))
    
    as.matrix(tabRes)[cbind(iRows, x)]
}
