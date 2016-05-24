NRFUopt <- function(Ctot=NULL, c1, c2, theta, CV0=NULL, CVpop=NULL, N=Inf, type.sw){
    if (!(type.sw %in% c("cost", "cv")))
        stop("type.sw must be cost or cv.\n")
    if (is.null(c1) | is.null(c2) | is.null(theta))
        stop("c1, c2, and theta cannot be NULL.\n")
    if (c1 < 0 | c2 < 0)
        stop("Unit costs, c1 and c2, must be positive.\n")
    if (c1 < 0 | c2 < 0)
        stop("Unit costs c1 and c2 must be positive.\n")
    if (theta <= 0 | theta > 1)
        stop("Response probability must be in (0,1].\n")

    if (type.sw == "cost"){
        if (is.null(Ctot)) stop("Ctot must be specified for fixed cost allocation.\n")
        allocation <- "fixed cost"
    }
    if (type.sw == "cv"){
        if (is.null(CV0)) stop("CV0 must be specified for fixed CV allocation.\n")
        if (CV0 <0 )
            stop("CV0 must be positive.\n")
        allocation <- "fixed CV"
    }

    v.opt <- sqrt(c1/c2/theta)

    if (type.sw=="cv"){
        n1.opt <- (1/v.opt) * (1-theta*(1-v.opt)) / ((CV0/CVpop)^2 + 1/N)
    }
    if (type.sw=="cost"){
        n1.opt <- Ctot / (c1 + c2*v.opt*(1-theta))
    }

    n2 <- v.opt * (1-theta) * n1.opt
    Ctot.chk <- c1*n1.opt + c2*v.opt*(1-theta)*n1.opt

    if (!is.null(CVpop)){
        CV0.chk <- sqrt(CVpop^2/n1.opt * (1-n1.opt/N + (1-v.opt)/v.opt * (1-theta)))
        CV0.chk <- round(CV0.chk,4)
    }
            else {CV0.chk <- NULL}

    n.srs <- n1.opt* 1/(theta + (1-theta)/v.opt - n1.opt/N)  / theta
    c.ratio <- n1.opt/n.srs * (1 + c2/c1 *v.opt * (1-theta))

    if (v.opt > 1) warning("v.opt > 1: Solution is not feasible.\n")

    list("allocation" = allocation,
         "Total variable cost" = Ctot.chk,
         "Response rate" = theta,
         "CV" = CV0.chk,
         "v.opt" = round(v.opt,4),
         "n1.opt" = n1.opt,
         "Expected n2" = n2,
         "Expected total cases (2-phase)" = n1.opt + n2,
         "srs sample for same cv" = n.srs,
         "Cost Ratio: Two phase to srs" = round(c.ratio,3))
}
