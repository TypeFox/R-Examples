## Support for fuzzy sets needs the specification of negation,
## conjunction and disjunction via the functions N, T (t-norm) and S
## (t-conorm), respectively.
##
## We try to make our code "work" with arbitrary such triples by using
## functions .N., .T. and .S., respectively, and allow setting the fuzzy
## logic system via a dynamic variable.  If this turns out to be too
## cumbersome, we could also move back to hard-wiring the "usual" Zadeh
## connectives 1 minus, min and max, and provide ways to unlock the
## locked bindings as needed.

## .N. <- function(x) 1 - x
## .T. <- function(x, y) pmin(x, y)
## .S. <- function(x, y) pmax(x, y)
## .I. <- function(x, y) ifelse(x <= y, 1, y)

.N. <- function(x) fuzzy_logic()$N(x)

## ensure that .T.(NA, 0) == 0
.T. <- function(x, y)
    `[<-`(fuzzy_logic()$T(x, y),
          xor(is.na(y), is.na(x)) & pmin.int(x, y, na.rm = TRUE) == 0,
          0)

## ensure that .S.(NA, 1) == 1
.S. <- function(x, y)
    `[<-`(fuzzy_logic()$S(x, y),
          xor(is.na(y), is.na(x)) & pmax.int(x, y, na.rm = TRUE) == 1,
          1)

## ensure that .I.(0, NA) == 1
.I. <- function(x, y)
    `[<-`(fuzzy_logic()$I(x, y),
          is.na(y) & !is.na(x) & x == 0,
          1)

## Use dynamic variables for the fuzzy connectives.
## Note that there are also parametric fuzzy logic families, which we
## can set via
##   fuzzy_logic(NAME, PARAMS)

## One might think that using an environment for storing the fuzzy logic
## connectives, along the lines of
##   fuzzy_logic_db <- new.env()
## and adding
##   for(nm in names(family))
##     assign(nm, family[[nm]], envir = fuzzy_logic_db)
## and using
##   .N. <- function(...) fuzzy_logic_db$N(...)
## would be somewhat faster than
##   .N. <- function(...) fuzzy_logic()$N(...)
## but we could not find conclusive evidence for performance gains
## either way.

fuzzy_logic <-
local({
    family <- NULL
    function(new, ...) {
        if(!missing(new))
            family <<- fuzzy_logic_families[[new]](...)
        else
            family
    }
})

fuzzy_logic_family <-
function(name, T, S, N, I = NULL, params = NULL, meta = NULL)
{
    if(is.null(I))
        I <- function(x, y)
            stop("Implication not available.", call. = FALSE)
    .structure(list(name = name, T = T, S = S, N = N, I = I,
                    params = params, meta = meta),
               class = "fuzzy_logic_family")
}

fuzzy_logic_predicates <-
function()
    fuzzy_logic()$meta

fuzzy_logic_family_Zadeh <-
function()
    fuzzy_logic_family(name = "Zadeh",
                       N = function(x) 1 - x,
                       T = function(x, y) pmin(x, y),
                       S = function(x, y) pmax(x, y),
                       I = function(x, y) ifelse(x <= y, 1, y),
                       meta =
                       list(is_de_Morgan_triple = TRUE,
                            N_is_standard = TRUE,
                            T_is_continuous = TRUE,
                            T_is_Archimedean = FALSE)
                       )

fuzzy_logic_family_drastic <-
function()
    fuzzy_logic_family(name = "drastic",
                       N = function(x) 1 - x,
                       T = function(x, y)
                       ifelse(pmax(x, y) == 1, pmin(x, y), 0),
                       S = function(x, y)
                       ifelse(pmin(x, y) == 0, pmax(x, y), 1),
                       meta =
                       list(is_de_Morgan_triple = TRUE,
                            N_is_standard = TRUE,
                            T_is_continuous = FALSE,
                            T_is_Archimedean = TRUE,
                            T_generator = function(x) {
                                ifelse(x < 1, 2 - x, 0)
                            }
                            )
                       )

fuzzy_logic_family_product <-
function()
    fuzzy_logic_family(name = "product",
                       N = function(x) 1 - x,
                       T = function(x, y) x * y,
                       S = function(x, y) x + y - x * y,
                       I = function(x, y) pmin(y / x, 1),
                       meta =
                       list(is_de_Morgan_triple = TRUE,
                            N_is_standard = TRUE,
                            T_is_continuous = TRUE,
                            T_is_Archimedean = TRUE,
                            T_generator = function(x) { - log(x) }
                            )
                       )

fuzzy_logic_family_Lukasiewicz <-
function()
    fuzzy_logic_family(name = "Lukasiewicz",
                       N = function(x) 1 - x,
                       T = function(x, y) pmax(x + y - 1, 0),
                       S = function(x, y) pmin(x + y, 1),
                       I = function(x, y) pmin(1 - x + y, 1),
                       meta =
                       list(is_de_Morgan_triple = TRUE,
                            N_is_standard = TRUE,
                            T_is_continuous = TRUE,
                            T_is_Archimedean = TRUE,
                            T_generator = function(x) { 1 - x }
                            )
                       )

## Nilpotent minimum (\min_0).
fuzzy_logic_family_Fodor <-
function()
    fuzzy_logic_family(name = "Fodor",
                       N = function(x) 1 - x,
                       T = function(x, y)
                       ifelse(x + y > 1, pmin(x, y), 0),
                       S = function(x, y)
                       ifelse(x + y < 1, pmax(x, y), 1),
                       I = function(x, y)
                       ifelse(x <= y, 1, pmax(1 - x, y)),
                       meta =
                       list(is_de_Morgan_triple = TRUE,
                            N_is_standard = TRUE,
                            T_is_continuous = FALSE,
                            T_is_Archimedean = FALSE
                            )
                       )

## Frank family (Fodor & Roubens, page 20).
## One parameter family with 0 <= s <= \infty, where 0, 1 and \infty
## give Zadeh, product and Lukasiewicz, respectively.
fuzzy_logic_family_Frank <-
function(s)
{
    if(s < 0)
        stop("Invalid parameter.")
    if(s == 0) fuzzy_logic_family_Zadeh()
    else if(s == 1) fuzzy_logic_family_product()
    else if(s == Inf) fuzzy_logic_family_Lukasiewicz()
    else {
        T <- function(x, y)
            log(1 + (s^x - 1) * (s^y - 1) / (s - 1)) / log(s)
        fuzzy_logic_family(name = "Frank",
                           N = function(x) 1 - x,
                           T = T,
                           S = function(x, y) 1 - T(1 - x, 1 - y),
                           I = function(x, y)
                           ifelse(x <= y, 1,
                                  log(1 + (s - 1) * (s^y - 1) /
                                                    (s^x - 1)) / log(s)),
                           params = list(s = s),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE,
                                T_is_continuous = TRUE,
                                T_is_Archimedean = TRUE,
                                T_generator = function(x) {
                                    log((s - 1) / (s^x - 1))
                                }
                                ## Note:
                                ##   f^{-1}(x) = \log_s(1 + (s-1)e^{-x})
                                )
                           )
    }
}

## Hamacher family (Fodor & Roubens, page 21).
## This is a three parameter family of connectives T_\alpha, S_\beta,
## N_\gamma with \alpha >= 0, \beta, \gamma >= -1, which is a de Morgan
## triple iff \gamma > -1 and \alpha = (1 + \beta) / (1 + \gamma).
## Be nice and provide some default values ...
fuzzy_logic_family_Hamacher <-
function(alpha = NULL, beta = 0, gamma = 0)
{
    if(is.null(alpha))
        alpha <- (1 + beta) / (1 + gamma)
    else if((alpha < 0) || (beta < -1) || (gamma < -1))
        stop("Invalid parameter.")
    fuzzy_logic_family(name = "Hamacher",
                       N = function(x) (1 - x) / (1 + gamma * x),
                       T = function(x, y)
                       ifelse(x * y == 0,
                              0,
                              x * y /
                              (alpha + (1 - alpha) * (x + y - x * y))),
                       S = function(x, y)
                       (x + y + (beta - 1) * x * y) /
                        (1 + beta * x * y),
                       I = function(x, y)
                       ifelse(x <= y, 1,
                              ((y * (alpha + (1 - alpha) * x)) /
                               (x + (1 - alpha) * y * (x - 1)))),
                       params =
                       list(alpha = alpha, beta = beta, gamma = gamma),
                       meta =
                       list(is_de_Morgan_triple =
                            alpha == (1 + beta) / (1 + gamma),
                            N_is_standard = (gamma == 0),
                            T_is_continuous = TRUE,
                            T_is_Archimedean = TRUE,
                            T_generator = function(x) {
                                log((alpha + (1 - alpha) * x) / x)
                            }
                            ## Note:
                            ##   f^{-1}(x) = p / (e^x + p - 1)
                            )
                       )
}

## The following "families" are really families of t-norms, which we
## leverage into fuzzy logic families using standard negation and de
## Morgan triplet associated t-conorm.

## Schweizer-Sklar.
fuzzy_logic_family_Schweizer_Sklar <-
function(p)
{
    if(p == -Inf) fuzzy_logic_family_Zadeh()
    else if(p == 0) fuzzy_logic_family_product()
    else if(p == Inf) fuzzy_logic_family_drastic()
    else {
        T <- if(p < 0)
            function(x, y) (x^p + y^p - 1) ^ (1/p)
        else
            function(x, y) pmax(0, (x^p + y^p - 1) ^ (1/p))
        fuzzy_logic_family(name = "Schweizer-Sklar",
                           N = function(x) 1 - x,
                           T = T,
                           S = function(x, y) 1 - T(1 - x, 1 - y),
                           I = function(x, y)
                           ifelse(x <= y, 1, (1 - x^p + y^p) ^ (1/p)),
                           params = list(p = p),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE,
                                T_is_continuous = TRUE,
                                T_is_Archimedean = TRUE,
                                T_generator = function(x) {
                                    (1 - x^p) / p
                                }
                                ## Note:
                                ##   f^{-1}(x) = (1 - px)^{1/p}
                                )
                           )
    }
}

## Yager t-norm.
fuzzy_logic_family_Yager <-
function(p)
{
    if(p < 0)
        stop("Invalid parameter.")
    if(p == 0) fuzzy_logic_family_drastic()
    else if(p == Inf) fuzzy_logic_family_Zadeh()
    else {
        fuzzy_logic_family(name = "Yager",
                           N = function(x) 1 - x,
                           T = function(x, y)
                           pmax(0, 1 - ((1 - x)^p + (1 - y)^p)^(1/p)),
                           S = function(x, y)
                           pmin(1, (x^p + y^p) ^ (1/p)),
                           I = function(x, y)
                           ifelse(x <= y, 1,
                                  1 - ((1 - y)^p - (1 - x)^p)^(1/p)),
                           params = list(p = p),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE,
                                T_is_continuous = TRUE,
                                T_is_Archimedean = TRUE,
                                T_generator = function(x) {
                                    1 - x^p
                                }
                                ## Note:
                                ##   f^{-1}(x) = 1 - x^{1/p}
                                )
                           )
    }
}

## Dombi t-norm.
fuzzy_logic_family_Dombi <-
function(p)
{
    if(p < 0)
        stop("Invalid parameter.")
    if(p == 0) fuzzy_logic_family_drastic()
    else if(p == Inf) fuzzy_logic_family_Zadeh()
    else {
        T <- function(x, y)
            ifelse(x * y == 0,
                   0,
                   1 / (1 + ((1 / x - 1)^p + (1 / y - 1)^p)^(1 / p)))
        fuzzy_logic_family(name = "Dombi",
                           N = function(x) 1 - x,
                           T = T,
                           S = function(x, y) 1 - T(1 - x, 1 - y),
                           I = function(x, y)
                           ifelse(x <= y, 1,
                                  1 / (1 + ((1 / y - 1)^p -
                                            (1 / x - 1)^p)^(1/p))),
                           params = list(p = p),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE,
                                T_is_continuous = TRUE,
                                T_is_Archimedean = TRUE,
                                T_generator = function(x) {
                                    (1 / x - 1)^p
                                }
                                ## Note:
                                ##   f^{-1}(x) = 1/(1 + x^{1/p})
                                )
                           )
    }
}

## Aczel-Alsina t-norm.
fuzzy_logic_family_Aczel_Alsina <-
function(p)
{
    if(p < 0)
        stop("Invalid parameter.")
    if(p == 0) fuzzy_logic_family_drastic()
    else if(p == Inf) fuzzy_logic_family_Zadeh()
    else {
        T <- function(x, y) exp(- (abs(log(x))^p + abs(log(y))^p))
        fuzzy_logic_family(name = "Aczel-Alsina",
                           N = function(x) 1 - x,
                           T = T,
                           S = function(x, y) 1 - T(1 - x, 1 - y),
                           I = function(x, y)
                           ifelse(x <= y, 1,
                                  exp(-((abs(log(y))^p -
                                         abs(log(x))^p))^(1/p))),
                           params = list(p = p),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE,
                                T_is_continuous = TRUE,
                                T_is_Archimedean = TRUE,
                                T_generator = function(x) {
                                    (- log(x))^p
                                }
                                ## Note:
                                ##   f^{-1}(x) = \exp(-x^{1/p})
                                )
                           )
    }
}

## Sugeno-Weber t-norm.
fuzzy_logic_family_Sugeno_Weber <-
function(p)
{
    if(p < -1)
        stop("Invalid parameter.")
    if(p == -1) fuzzy_logic_family_drastic()
    else if(p == Inf) fuzzy_logic_family_product()
    else
        fuzzy_logic_family(name = "Sugeno-Weber",
                           N = function(x) 1 - x,
                           T = function(x, y)
                           pmax(0, (x + y - 1 + p * x * y) / (1 + p)),
                           S = function(x, y)
                           pmin(1, x + y - p * x * y / (1 + p)),
                           I = function(x, y)
                           ifelse(x <= y, 1,
                                  (1 + (1 + p) * y - x) / (1 + p * x)),
                           params = list(p = p),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE,
                                T_is_continuous = TRUE,
                                T_is_Archimedean = TRUE,
                                T_generator = function(x) {
                                    1 - log(1 + p * x) / log(1 + p)
                                }
                                ## Note
                                ##   f^{-1}(x) = ((1 + p)^{1-x} - 1) / p
                                )
                           )
}

## Dubois-Prade t-norm and t-conorm.

fuzzy_logic_family_Dubois_Prade <-
function(p)
{
    if((p < 0) || (p > 1))
        stop("Invalid parameter.")
    if(p == 0) fuzzy_logic_family_Zadeh()
    else if(p == 1) fuzzy_logic_family_product()
    else {
        T <- function(x, y) x * y / pmax(x, y, p)
        fuzzy_logic_family(name = "Dubois-Prade",
                           N = function(x) 1 - x,
                           T = T,
                           S = function(x, y) 1 - T(1 - x, 1 - y),
                           ## Solve T(x, z) = y for x >= y ...
                           I = function(x, y)
                           ifelse(x <= y, 1, pmax(p / x, 1) * y),
                           params = list(p = p),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE
                                )
                           )
    }
}

## Yu t-norm and t-conorm.

fuzzy_logic_family_Yu <-
function(p)
{
    if(p < -1)
        stop("Invalid parameter.")
    if(p == -1) fuzzy_logic_family_product()
    else if(p == Inf) fuzzy_logic_family_drastic()
    else {
        fuzzy_logic_family(name = "Yu",
                           N = function(x) 1 - x,
                           T = function(x, y)
                           pmax(0, (1 + p) * (x + y - 1) - p * x * y),
                           S = function(x, y)
                           pmin(1, x + y + p * x * y),
                           params = list(p = p),
                           meta =
                           list(is_de_Morgan_triple = TRUE,
                                N_is_standard = TRUE
                                )
                           )
    }
}


fuzzy_logic_families <-
    list("Zadeh" = fuzzy_logic_family_Zadeh,
         "drastic" = fuzzy_logic_family_drastic,
         "product" = fuzzy_logic_family_product,
         "Lukasiewicz" = fuzzy_logic_family_Lukasiewicz,
         "Fodor" = fuzzy_logic_family_Fodor,
         "Frank" = fuzzy_logic_family_Frank,
         "Hamacher" = fuzzy_logic_family_Hamacher,
         "Schweizer-Sklar" = fuzzy_logic_family_Schweizer_Sklar,
         "Yager" = fuzzy_logic_family_Yager,
         "Dombi" = fuzzy_logic_family_Dombi,
         "Aczel-Alsina" = fuzzy_logic_family_Aczel_Alsina,
         "Sugeno-Weber" = fuzzy_logic_family_Sugeno_Weber,
         "Dubois-Prade" = fuzzy_logic_family_Dubois_Prade,
         "Yu" = fuzzy_logic_family_Yu
         )

## See also e.g. http://en.wikipedia.org/wiki/Construction_of_t-norms
## for more information.

## This has the generators for most Archimedean t-norms, but typically
## not the residual implications: so we use the result that
##   I(x, y) = f^{(-1)}(\max(f(y) - f(x), 0))
## where f is the (additive) generator and f^{(-1)} its pseudoinverse,
## defined as
##   f^{(-1)}(x) = f^{-1}(x) if x \le f(0) and 0 otherwise.
## For the implication, note that x \le y implies f(x) \ge f(y) and
## hence I(x, y) = f^{(-1)}(0) = 1; otherwise,
##   I(x, y) = f^{(-1)}(f(y) - f(x)), x \ge y.
