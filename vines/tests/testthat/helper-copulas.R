# Arguments and copula parameters used by test-h.R and test-hinverse.R.

n <- 20  # Number of values of each variable.
np <- 10  # Number of values of each parameter.
tol <- 0.01  # Ignore differences smaller than tol.

X <- seq(from = 0, to = 1, length = n)
V <- seq(from = 0, to = 1, length = n)
XV <- merge(X, V)
colnames(XV) <- c("X", "V")

copulas <- list(
    indep = list(indepCopula()),
    normal = lapply(seq(from = -1, to = 1, length = np),
                    function (rho) normalCopula(rho)),
    t = apply(merge(seq(from = -1, to = 1, length = np),
                    seq(from = 1, to = 30, length = np)), 1,
              function (p) tCopula(p[1], df = as.integer(p[2]), df.fixed = TRUE)),
    clayton = lapply(seq(from = .Machine$double.eps^0.25, to = 75, length = np),
                     function (theta) claytonCopula(theta, use.indepC = "FALSE")),
    gumbel = lapply(seq(from = 1, to = 100, length = np),
                    function (theta) gumbelCopula(theta, use.indepC = "FALSE")),
    fgm = lapply(seq(from = -1, to = 1, length = np),
                 function (theta) fgmCopula(theta)),
    galambos = lapply(seq(from = 1, to = 25, length = np),
                      function (theta) galambosCopula(theta)),
    frank = lapply(seq(from = -45, to = 45, length = np),
                   function (theta) frankCopula(theta, use.indepC = "FALSE"))
)

test_that_for_each_copula <- function(desc, code) {
    for (copula_name in names(copulas)) {
        new_desc <- paste(desc, "for the", copula_name, "copula")
        for (copula in copulas[[copula_name]]) {
            env <- parent.frame()
            assign("copula", copula, env)
            testthat:::test_code(new_desc, substitute(code), env = env)
        }
    }
}
