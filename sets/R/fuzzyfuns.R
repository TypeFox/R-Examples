### fuzzy characteristic functions

### * general constructor

charfun_generator <-
function(FUN, height = 1)
{
    ret <- if (is.null(height)) {
        function(...)
            function(x) pmax(0, pmin(FUN(x, ...), 1))
    } else {
        if (height < 0 || height > 1)
            stop("Height must be in the unit interval.")
        function(...)
            function(x) {
                ret <- pmax(0, pmin(FUN(x, ...), 1))
                height * ret / max(ret)
            }
    }
    class(ret) <- "charfun_generator"
    ret
}

is.charfun_generator <-
function(x)
    inherits(x, "charfun_generator")

### * special functions

fuzzy_normal <-
function(mean = NULL, sd = 1, log = FALSE, height = 1, chop = 0)
{
    if (!is.null(height) && (height < 0 || height > 1))
        stop("Height must be in the unit interval.")
    function(x) {
        if (is.null(mean))
            mean <- mean(range(x))
        ret <- dnorm(x, mean = mean, sd = sd, log = log)
        if (!is.null(height))
            ret <- height * ret / max(ret)
        ret[ret <= chop] <- 0
        ret
    }
}
class(fuzzy_normal) <- "charfun_generator"

fuzzy_two_normals <-
function(mean = NULL, sd = c(1,1),
         log = c(FALSE, FALSE), height = 1, chop = 0)
{
    if (!is.null(mean)) {
        if(length(mean) != 2L) stop("Need two mean values.")
        if (mean[2] < mean[1]) {
            sd <- rev(sd)
            log <- rev(log)
        }
    }
    if (!is.null(height) && (height < 0 || height > 1))
        stop("Height must be in the unit interval.")

    sd <- rep(sd, length.out = 2)
    log <- rep(log, length.out = 2)

    function(x) {
        if (is.null(mean))
            mean <- x[trunc(seq(1, length(x), length.out = 4))[2:3]]
        ret <- rep(height, length.out = length(x))

        tmp <- dnorm(x[x <= mean[1]], mean = mean[1], sd = sd[1], log = log[1])
        if (!is.null(height))
            tmp <- height * tmp / max(tmp)
        ret[x <= mean[1]] <- tmp

        tmp <- dnorm(x[x >= mean[2]], mean = mean[2], sd = sd[2], log = log[2])
        if (!is.null(height))
            tmp <- height * tmp / max(tmp)
        ret[x >= mean[2]] <- tmp

        ret[ret <= chop] <- 0
        ret
    }
}
class(fuzzy_two_normals) <- "charfun_generator"

fuzzy_bell <-
function(center = NULL, cross = NULL, slope = 4, height = 1, chop = 0)
{
    if (!is.null(height) && (height < 0 || height > 1))
        stop("Height must be in the unit interval.")
    function(x) {
        if (is.null(center))
            center <- mean(range(x))
        if (is.null(cross))
            cross <- trunc(diff(range(x)) / 5)
        ret <- 1 / (1 + abs((x - center) / cross) ^ (2 * slope))
        if (!is.null(height))
            ret <- height * ret / max(ret)
        ret[ret <= chop] <- 0
        ret
    }
}
class(fuzzy_bell) <- "charfun_generator"

fuzzy_sigmoid <-
function(cross = NULL, slope = 0.5, height = 1, chop = 0)
{
    if (!is.null(height) && (height < 0 || height > 1))
        stop("Height must be in the unit interval.")
    function(x) {
        if (is.null(cross))
            cross <- mean(range(x))
        ret <- 1 / (1 + exp(-slope * (x - cross)))
        if (!is.null(height))
            ret <- height * ret / max(ret)
        ret[ret <= chop] <- 0
        ret
    }
}
class(fuzzy_sigmoid) <- "charfun_generator"

fuzzy_trapezoid <-
function(corners = NULL, height = c(1,1), return_base_corners = TRUE)
{
    height <- rep(height, length.out = 2)
    if (any(height > 1) || any(height < 0))
        stop("Height(s) must be in the unit interval.")
    if (!is.null(corners) && length(corners) != 4L)
        stop("Need four corner values.")
    function(x) {
        if (is.null(corners))
            corners <- x[trunc(seq(1, length(x), length.out = 6))[2:5]]
        ret <- approxfun(corners, c(0, height, 0), rule = 2)(x)
        if (return_base_corners)
            ret[match(corners[c(1,4)], x)] <- .Machine$double.eps
        ret
    }
}
class(fuzzy_trapezoid) <- "charfun_generator"

fuzzy_triangular <-
function(corners = NULL, height = 1, return_base_corners = TRUE)
{
    if (height > 1 || height < 0)
        stop("Height must be in the unit interval.")
    if (!is.null(corners) && length(corners) != 3L)
        stop("Need three corner values.")
    function(x) {
        if (is.null(corners))
            corners <- x[trunc(seq(1, length(x), length.out = 5))[2:4]]
        ret <- approxfun(corners, c(0, height, 0), rule = 2)(x)
        if (return_base_corners)
            ret[match(corners[-2], x)] <- .Machine$double.eps
        ret
    }
}
class(fuzzy_triangular) <- "charfun_generator"

fuzzy_cone <-
function(center = NULL, radius = 2, height = 1, return_base_corners = TRUE)
{
    if (height > 1 || height < 0)
        stop("Height must be in the unit interval.")
    function(x) {
        if (is.null(center))
            center <- x[trunc((1 + length(x)) / 2)]
        ret <- approxfun(center + c(-radius, 0, radius),
                         c(0, height, 0), rule = 2)(x)
        if (return_base_corners)
            ret[match(center + c(-radius, radius), x)] <- .Machine$double.eps
        ret
    }
}
class(fuzzy_cone) <- "charfun_generator"

fuzzy_pi3 <-
function(mid = NULL, min = NULL, max = NULL, height = 1,
         return_base_corners = TRUE)
{
    if (height > 1 || height < 0)
        stop("Height must be in the unit interval.")
    function(x) {
        if (is.null(mid))
            mid <- x[trunc((1 + length(x)) / 2)]
        if (is.null(min)) min <- mid - 2
        if (is.null(max)) max <- mid + 2
        ret <- ifelse(x < min | x > max,
               0,
               ifelse(x < mid,
                      height * (1 - (x - mid) ^ 2 / (min - mid) ^ 2),
                      height * (1 - (x - mid) ^ 2 / (max - mid) ^ 2)
                      )
               )
        if (return_base_corners)
            ret[match(c(min, max), x)] <- .Machine$double.eps
        ret
    }
}
class(fuzzy_pi3) <- "charfun_generator"

fuzzy_pi4 <- function (knots, height = 1, return_base_corners = TRUE)
{
    if (height > 1 || height < 0)
        stop("Height must be in the unit interval.")
    if (length(knots) != 4L)
        stop("Need four knots.")
    function(x) {
        ret <- ifelse(x <= knots[1] | x >= knots[4],
                      0,
                      ifelse(x > knots[1] & x <= ((knots[1] + knots[2]) / 2),
                             2  * ((x - knots[1]) / (knots[2] - knots[1]))^2,
                             ifelse(x > ((knots[1] + knots[2])/2) & x < knots[2],
                                    1 - 2*((x - knots[2]) / (knots[2] - knots[1]))^2,
                                    ifelse(x >= knots[2] & x <= knots[3],
                                           height,
                                           ifelse(x > knots[3] & x <= ((knots[3] + knots[4]) / 2),
                                                  1 - 2 * ((x - knots[3]) / (knots[4] - knots[3]))^2,
                                                  2 * ((x - knots[4]) / (knots[4] - knots[3]))^2)))))
        if (return_base_corners)
            ret[match(c(knots[1], knots[4]), x)] <- .Machine$double.eps
        ret
    }
}
class(fuzzy_pi4)<- "charfun_generator"

## * fuzzy set generators for convenience

.expand <-
function(universe = NULL)
{
    if (is.null(universe))
        universe <- sets_options("universe")
    if (is.null(universe))
        universe <- seq(0,20,0.1)
    as.set(eval(universe))
}

fuzzy_normal_gset <-
function(mean = NULL, sd = 1, log = FALSE, height = 1, chop = 0,
         universe = NULL)
    gset(charfun = fuzzy_normal(mean = mean, sd = sd, log = log,
                                height = height, chop = chop),
         universe = .expand(universe))

fuzzy_two_normals_gset <-
function(mean = NULL, sd = c(1,1), log = c(FALSE, FALSE),
         height = 1, chop = 0, universe = NULL)
    gset(charfun = fuzzy_two_normals(mean = mean, sd = sd, log = log,
                                     height = height, chop = chop),
         universe = .expand(universe))

fuzzy_bell_gset <-
function(center = NULL, cross = NULL, slope = 4, height = 1, chop = 0,
         universe = NULL)
    gset(charfun = fuzzy_bell(center = center, cross = cross, slope = slope,
                              height = height, chop = chop),
         universe = .expand(universe))

fuzzy_sigmoid_gset <-
function(cross = NULL, slope = 0.5, height = 1, chop = 0,
         universe = NULL)
    gset(charfun = fuzzy_sigmoid(cross = cross, slope = slope,
                                 height = height, chop = chop),
         universe = .expand(universe))

fuzzy_trapezoid_gset <-
function(corners = NULL, height = c(1, 1), universe = NULL,
         return_base_corners = TRUE)
    gset(charfun = fuzzy_trapezoid(corners = corners, height = height,
                                   return_base_corners = return_base_corners),
         universe = .expand(universe))

fuzzy_triangular_gset <-
function(corners = NULL, height = 1, universe = NULL,
         return_base_corners = TRUE)
    gset(charfun = fuzzy_triangular(corners = corners, height = height,
                                    return_base_corners = return_base_corners),
         universe = .expand(universe))

fuzzy_cone_gset <-
function(center = NULL, radius = 2, height = 1, universe = NULL,
         return_base_corners = TRUE)
    gset(charfun = fuzzy_cone(center = center, radius = radius,
                              height = height,
                              return_base_corners = return_base_corners),
         universe = .expand(universe))

fuzzy_pi3_gset <-
function(mid = NULL, min = NULL, max = NULL, height = 1, universe = NULL,
         return_base_corners = TRUE)
    gset(charfun = fuzzy_pi3(mid = mid, min = min, max = max, height = height,
                             return_base_corners = return_base_corners),
         universe = .expand(universe))

fuzzy_pi4_gset <-
function(knots, height = 1, universe = NULL, return_base_corners = TRUE)
    gset(charfun = fuzzy_pi4(knots = knots, height = height, return_base_corners = return_base_corners),
         universe = .expand(universe))

### * tuple generator

fuzzy_tuple <-
function(FUN = fuzzy_normal, n = 5, ..., universe = NULL, names = NULL)
{
    universe <- .expand(universe)
    F <- if (is.charfun_generator(FUN))
        function(i) gset(charfun = FUN(i, ...), universe = universe)
    else
        function(i) FUN(i, universe = universe, ...)

    if (length(n) == 1L)
        n <- .get_support(universe)[seq(from = 1, to = length(universe),
                                        length.out = n)]
    .structure(as.tuple(lapply(n, F)), names = names)
}

