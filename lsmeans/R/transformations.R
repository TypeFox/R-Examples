# Code to implement transformations my way

### My modification/expansion of stats:make.link()
### Also, if not found, returns make.link("identity") modified with
##     unknown = TRUE, name = link
## In addition, I make all links truly monotone on (-Inf, Inf) in
##     lieu of valideta
##
## Extensions to make.link results:
##     unknown: set to TRUE if link is unknown
##     mult: scalar multiple of transformation
##
.make.link = function(link) {
    if (link %in% c("logit", "probit", "cauchit", "cloglog", "identity", "log"))
        result = stats::make.link(link)
    else result = switch(link,
         sqrt = { tmp = make.link("sqrt") 
             tmp$linkinv = function(eta) pmax(0, eta)^2
             tmp$mu.eta = function(eta) 2*pmax(0, eta)
             tmp },
         `1/mu^2` = { tmp = make.link("1/mu^2") 
             tmp$linkinv = function(eta) 1/sqrt(pmax(0, eta))
             tmp$mu.eta = function(eta) -1/(2*pmax(0, eta)^1.5)
             tmp },
         inverse = { tmp = make.link("inverse") 
             tmp$linkinv = function(eta) 1/pmax(0, eta)
             tmp$mu.eta = function(eta) -1/pmax(0, eta)^2
             tmp },
         `/` = .make.link("inverse"),
         reciprocal = .make.link("inverse"),
         log10 = list(
             linkinv = function(eta) 10^eta,
             mu.eta = function(eta) 10^eta * log(10),
             name = "log10"
         ),
         log2 = list(
             linkinv = function(eta) 2^eta,
             mu.eta = function(eta) 2^eta * log(2),
             name = "log2"
         ),
         asin.sqrt = make.tran("asin.sqrt"),
         `asin.sqrt./` = make.tran("asin.sqrt", 100),
         asinh.sqrt = list(
             linkinv = function(eta) sinh(eta)^2,
             mu.eta = function(eta) sinh(2 * eta),
             name = "asinh(sqrt(mu))"
         ),
         `+.sqrt` = {
             tmp = .make.link("sqrt")
             tmp$mult = 2
             tmp
         },
         
         { # default if not included, flags it as unknown
             tmp = stats::make.link("identity")
             tmp$unknown = TRUE
             tmp$name = link
             tmp
         }
    )
    result
}

# Implementation of additional transformations, typically ones with parameters
# Returns a list like stats::make.link, but often with an additional "param" member
# types:
#       glog: log(mu + param)
make.tran = function(type = c("genlog", "power", "boxcox", "sympower", "asin.sqrt"), param = 1) {
    type = match.arg(type)
    origin = 0
    mu.lbl = "mu"
    if (length(param) > 1) {
        origin = param[2]
        param = param[1]
        mu.lbl = paste0("(mu - ", round(origin, 3), ")")
    }
    switch(type,
        genlog = {
            if((origin < 0) || (origin == 1))
                stop('"genlog" transformation must have a positive base != 1')
            logbase = ifelse(origin == 0, 1, log(origin))
            xlab = ifelse(origin == 0, "", paste0(" (base ", round(origin, 3), ")"))
            list(linkfun = function(mu) log(pmax(mu + param, 0)) / logbase,
                 linkinv = function(eta) pmax(exp(logbase * eta), .Machine$double.eps) - param,
                 mu.eta = function(eta) logbase * pmax(exp(logbase * eta), .Machine$double.eps),
                 valideta = function(eta) TRUE,
                 param = c(param, origin),
                 name = paste0("log(mu + ", round(param,3), ")", xlab)
            )
        },
        power = {
            if (param == 0) {
                if(origin == 0) make.link("log")
                else make.tran("genlog", -origin)
            }
            else list(
                linkfun = function(mu) pmax(mu - origin, 0)^param,
                linkinv = function(eta) origin + pmax(eta, 0)^(1/param),
                mu.eta = function(eta) pmax(eta, 0)^(1/param - 1) / param,
                valideta = function(eta) all(eta > 0),
                param = c(param, origin),
                name = ifelse(param > 0, 
                              paste0(mu.lbl, "^", round(param,3)),
                              paste0(mu.lbl, "^(", round(param,3), ")"))
            )
        },
        boxcox = {
            if (param == 0) {
                result = if(origin == 0) make.link("log")
                         else make.tran("genlog", -origin)
                return (result)
            }
            min.eta = ifelse(param > 0, -1 / param, -Inf)
            xlab = ifelse(origin == 0, "", paste0(" with origin at ", round(origin, 3)))
            list(
                linkfun = function(mu) ((mu - origin)^param - 1) / param,
                linkinv = function(eta) origin + (1 + param * pmax(eta, min.eta))^(1/param),
                mu.eta = function(eta) (1 + param * pmax(eta, min.eta))^(1/param - 1),
                valideta = function(eta) all(eta > min.eta),
                param = c(param, origin),
                name = paste0("Box-Cox (lambda = ", round(param, 3), ")", xlab)
            )
        },
        sympower = {
            if (param <= 0) 
                stop('"sympower" transformation requires positive param')
            if (origin == 0) 
                mu.lbl = paste0("(", mu.lbl, ")")
            absmu.lbl = gsub("\\(|\\)", "|", mu.lbl)
            list(linkfun = function(mu) sign(mu - origin) * abs(mu - origin)^param,
                 linkinv = function(eta) origin + sign(eta) * abs(eta)^(1/param),
                 mu.eta = function(eta) (abs(eta))^(1/param - 1),
                 valideta = function(eta) all(eta > min.eta),
                 param = c(param, origin),
                 name = paste0(absmu.lbl, "^", round(param,3), " * sign", mu.lbl)
            )
        },
        asin.sqrt = {
            mu.lbl = ifelse(param == 1, "mu", paste0("mu/", round(param,3)))
            list(linkfun = function(mu) asin(sqrt(mu/param)),
                 linkinv = function(eta) param * sin(pmax(pmin(eta, pi/2), 0))^2,
                 mu.eta = function(eta) param * sin(2*pmax(pmin(eta, pi/2), 0)),
                 valideta = function(eta) all(eta <= pi/2) && all(eta >= 0),
                 name = paste0("asin(sqrt(", mu.lbl, "))")
            )
        }
    )
}
