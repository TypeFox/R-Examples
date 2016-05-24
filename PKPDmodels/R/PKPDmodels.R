##' Substitute the expression \code{sub} for the name \code{nm} in
##' \code{expr} by walking the tree.
##'
##' @title "Sub[stitute] expr[ession]"
##'
##' @param expr an expression
##' @param nm a name for which to substitute
##' @param sub the expression to substitute for name nm
##'
##' @return the expression with all occurrences of nm replaced by sub
##' @note this function is called recursively
subexpr <- function(expr, nm, sub) 
{
    if (length(expr) == 1) {
        if (is.name(expr) && expr == nm) return(sub[[1]])
        return(expr)
    }
    for (j in 2:length(expr)) expr[[j]] <- subexpr(expr[[j]], nm, sub)
    expr
}

##' Return a formula for the PK model with linear
##' elimination according to the number of compartments, the administration form and the dosage pattern.
##' 
##' @title Expressions for PK models with linear elimination
##'
##' @param admin form of administration of the drug, one of
##'    \code{"bolus"}, \code{"infusion"} or \code{"oral"}.  Defaults to
##'    \code{"bolus"}. 
##' @param dosage form of dosage, one of \code{"sd"} (single dose),
##'    \code{"md"} (multiple, equally-spaced doses) and \code{"ss"}
##'    (steady-state).  Defaults to \code{"sd"}.
##' @param subst a list of formulas of substitutions to perform
##' @param cpt scalar integer - the number of model compartments.
##' @return a formula
##' @examples
##' ## single-dose oral administration
##' PKexpr("oral", "sd")
PKexpr <- function(admin=c("bolus", "infusion", "oral"),
                   dosage=c("sd", "md", "ss"),
                   subst=list(),
                   cpt=1L) {
    stopifnot((cpt <- as.integer(cpt)[1]) > 0, cpt < 4)
    frm <- list(                        # one compartment
                list(bolus =
                     list(sd = ~ dose * exp(-k * t) / V,

                          md = ~(dose/V) * ((1-exp(-N*k*tau))/(1-exp(-k*tau))) * exp(-k*(t-(N-1)*tau)),

                          ss = ~(dose/V)/(1-exp(-k*tau))*(exp(-k*(t-(TimeSS))))
                          ),
                     infusion =
                     list(sd = ~(dose/TInf) / (k*V)(1-exp(-k*TInf))*(exp(-k*(t-TInf))),
                          
                          md = ~(dose/TInf) / (k*V) * (1-exp(-k*TInf))*
                          (1-exp(-N * k * tau)) / (1-exp(-k * tau)) *
                          (exp(-k * (t - (N-1)*tau - TInf))),
                          
                          ss = ~(dose/Tinf) /(k*V)*(1-exp(-k * TInf)) *
                          exp(-k * (t-TimeSS-TInf))/(1-exp(-k*tau))
                          ),
                     oral =
                     list(sd = ~(dose/V) * (ka/(ka-k)) * (exp(-k*t)-exp(-ka*t)),

                          md = ~(dose/V) * (ka/(ka-k)) * (exp(-k*(t-(N-1)*tau))*(1-exp(-N*k*tau)) /
                                                          (1-exp(-k*tau))-exp(-ka*(t-(N-1)*tau)) *
                                                          (1-exp(-N*ka*tau))/(1-exp(-ka*tau))),
                          
                          ss = ~(dose/V) * (ka/(ka-k)) * (exp(-k*(t-TimeSS))/(1-exp(-k*tau)) -
                                                          exp(-ka*(t-TimeSS))/(1-exp(-ka*tau)))
                          )
                     ),
                list(                   # 2 compartment
                     bolus =
                     list(sd = ~(dose/V) * (A*exp(-alpha1 * t) + B*exp(-alpha2 * t)),
                          
                          md = ~(dose/V) * (((1-A*exp(-N*alpha1*tau))/(1-A*exp(-alpha1*tau))) *
                                            A*exp(-alpha1*(t-(N-1)*tau))
                                            +((1-B*exp(-N*alpha2*tau))/(1-B*exp(-alpha2*tau))) *
                                            B*exp(-alpha2*(t-(N-1)*tau))),
                          
                          ss = ~(dose/V)*((A*exp(-alpha1*(t-(TimeSS))))/(1-A*exp(-alpha1*tau))+
                                          (B*exp(-alpha2*(t-(TimeSS))))/(1-B*exp(-alpha2*tau)))
                          ),
                     infusion =list(),
                     oral =list()
                     ),
                list()                  # 3 compartment (not yet available)
                )[[cpt]][[match.arg(admin)]][[match.arg(dosage)]]
                    
    for (i in seq_along(subst)) {
        stopifnot(class(subfrm <- eval(subst[[i]])) == "formula",
                  is.name(subfrm[[2]]))
        frm <- subexpr(frm, subfrm[[2]], as.expression(subfrm[[3]]))
    }
    frm
}    

##' Create a model function with gradient evaluation (and, optionally,
##' Hessian evaluation) for a model according to the number of compartments,
##' the form of administration and dosage of
##' the drug after performing any substitutions given. 
##' 
##'
##' The substitutions are given as a list of formulas, such as 
##' \code{list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV))}.  They are applied left
##' to right.
##' 
##' @title PK models with linear elimination
##'
##' @param admin form of administration of the drug, one of
##'    \code{"bolus"}, \code{"infusion"} or \code{"oral"}.  Defaults to
##'     \code{"bolus"}. 
##' @param dosage type of dosage of the drug, one of \code{"sd"}
##'    (single dose), \code{"md"} (multiple dose) or \code{"ss"}
##'    (steady-state).  Defaults to \code{"sd"}.
##' @param subst a list of formulas of substitutions to perform
##' @param cpt scalar integer - the number of model compartments.
##' @param hessian a logical value indicating whether the second
##'    derivatives should be calculated and incorporated in the return
##'    value. 
##' @return a byte-compiled model function with gradient evaluation 
##'
##' @examples
##' ## return a function with substitutions
##' PKmod("bolus", "sd", list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV)))
##'
PKmod <- function(admin  = c("bolus", "infusion", "oral"),
                  dosage = c("sd", "md", "ss"),
                  subst  = list(),
                  cpt    = 1L,
                  hessian= FALSE) {
    stopifnot((cpt <- as.integer(cpt)[1]) > 0L, cpt < 4L)
    admin  <- match.arg(admin)
    dosage <- match.arg(dosage)
    frm    <- PKexpr(admin, dosage, subst, cpt)
    covariates <- list(list(bolus=
                            list(sd=c("dose", "t"),
                                 md=c("dose","t","TInf"),
                                 ss=c("dose", "t")),
                            infusion=
                            list(sd=character(0),
                                 md=c("N", "tau"),
                                 ss=c("TimeSS", "N", "tau")),
                            oral=
                            list(sd=c("dose", "t"),
                                 md=c("dose","t","TInf"),
                                 ss=c("dose", "t"))))[[cpt]][[admin]][[dosage]]
    pnms <- setdiff(all.vars(frm), covariates)
    cmpfun(deriv(frm, pnms, c(covariates, pnms), hessian=hessian))
}
