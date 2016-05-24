### <======================================================================>
setClass("ghyp",
         representation(call = "call",
                        lambda = "numeric",
                        alpha.bar = "numeric",
                        chi = "numeric",
                        psi = "numeric",
                        mu = "numeric",
                        sigma = "matrix",
                        gamma = "numeric",
                        model = "character",
                        dimension = "numeric",
                        expected.value = "numeric",
                        variance = "matrix",
                        parametrization = "character",
                        data = "matrix"),
          prototype(call = call("ghyp"),
                    lambda = 1,
                    alpha.bar = 1,
                    chi = 1,
                    psi = 1,
                    mu = 0,
                    sigma = matrix(0),
                    gamma = 0,
                    model = "Symmetric Generalized Hyperbolic",
                    dimension = 1,
                    expected.value = numeric(0),
                    variance = matrix(0),
                    parametrization = "alpha.bar",
                    data = matrix(0))
)
### <---------------------------------------------------------------------->

### <======================================================================>
setClass("mle.ghyp",
         representation(n.iter = "numeric",
                        llh = "numeric",
                        converged = "logical",
                        error.code = "numeric",
                        error.message = "character",
                        fitted.params = "logical",
                        aic = "numeric",
                        parameter.variance = "matrix",
                        trace.pars = "list"),

         prototype(n.iter = numeric(0),
                   llh = numeric(0),
                   converged = FALSE,
                   error.code = 0,
                   error.message = character(0),
                   fitted.params = logical(0),
                   aic = numeric(0),
                   parameter.variance = matrix(0),
                   trace.pars = list()),

         contains = "ghyp"
)
### <---------------------------------------------------------------------->
