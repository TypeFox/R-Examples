ScaleAdv=function (x, center = mean, scale = sd)
{
    if (class(x) != "matrix") {
        if (class(x) == "data.frame")
            x = as.matrix(x)
        else x = matrix(x, ncol = 1)
    }
    n = nrow(x)
    p = ncol(x)
    m = array(0, p)
    if (is.character(center))
        center = eval(parse(text = center))
    if (is.function(center)) {
        m = center(x)
        if (length(m) != p)
            m = apply(x, 2, center)
        x = x - matrix(1, nrow = n) %*% m
    }
    else if (length(center) == p & is.numeric(center)) {
        m = center
        x = x - matrix(1, nrow = n) %*% m
    }
    else if (!is.null(center))
        warning("Unknown center specification! Centering will be omitted.\n")
    s = array(1, p)
    if (is.character(scale))
        scale = eval(parse(text = scale))
    if (is.function(scale)) {
        s = scale(x)
        if (length(s) != p)
            s = apply(x, 2, scale)
        x = x/matrix(1, nrow = n) %*% s
    }
    else if (length(scale) == p & is.numeric(scale)) {
        s = scale
        x = x/matrix(1, nrow = n) %*% s
    }
    else if (!is.null(scale))
        warning("Unknown scale specification! Scaling will be omitted.\n")
    return(list(x = x, center = m, scale = s))
}
