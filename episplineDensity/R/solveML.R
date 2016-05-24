solveML <- function (data, softinfo, epiparameters, caseinfo, constraints, 
outer.trace = FALSE, postproc.controls = postproc.control(),
return.constraints = TRUE, opts, local.opts, ...)
{
#
# Do ML estimating.
#
less.const.function <- function (r, softinfo, epiparameters, caseinfo, c.out,                                   constraints) {
#
# If there's a moment constraint, call the nonlinear.constraints function.
# 
#
    nonlin <- FALSE
    linear.part <- as.vector (constraints$A %*% r) - constraints$b
    linear.jac <- constraints$A
    if (any (names (softinfo) == "upperbound1moment"))
        nonlin <- TRUE
    if (any (names (softinfo) == "upperbound2moment")) # ensure positivity?
        nonlin <- TRUE
    if (nonlin == FALSE)
        return (list (constraints = linear.part, jacobian = linear.jac))
# Nonlinear needed.
    nonlin.out <- nonlinear.constraints (epiparameters, softinfo, r) 
    return (list (constraints = c(linear.part, nonlin.out$g),
                  jacobian = rbind (linear.jac, t(nonlin.out$grad_g))))
}

equal.const.function <- function (r, softinfo, epiparameters, caseinfo, c.out, 
                                       constraints) {
#
# Call the integrate.to.one function if needed.
#
linear.part <- as.vector (constraints$Aeq %*% r) - constraints$beq
linear.jac <-  constraints$Aeq

if (any (names (softinfo) == "integrateToOne") && !is.null (softinfo$integrateToOne) && softinfo$integrateToOne == TRUE)
{
    nonlinear.part <- integrate.to.one (epiparameters, softinfo, r)
    return (list (constraints = c(linear.part, nonlinear.part$nonlinconst.eq.val),
                  jacobian = rbind(linear.jac, t(nonlinear.part$nonlinconst.eq.grad))))
}
    else return (list (constraints = linear.part, jacobian = linear.jac))
}

# This is for debugging for now
if (any (names (constraints) == "A")) {
    A <- constraints$A
    b <- constraints$b
}
else {
    A <- NULL; b <- NULL
}
if (any (names (constraints) == "Aeq")) {
    Aeq <- constraints$Aeq
    beq <- constraints$beq
}
else {
    Aeq <- NULL; beq <- NULL
}

# "c.out" plays the role of "c" in Matlab. "c" is not a good R name.
# This item holds the sum of the vectors of coefficients for each 
# data point. We should probably compute it as a matrix and carry
# that matrix around, rather than recomputing it point by point.
#
c.out <- numeric ( (epiparameters$order + 2) * epiparameters$N + 1)
for (i in 1:length(data))
    c.out <- c.out + coeff(data[i], epiparameters)

#
# You know what else should work?
# c.out <- rowSums (sapply (data, function (i) coeff (i, epiparameters)))))
#

# Vector of starting values. Start at 1's for right now.
r_in <- rep (1, (epiparameters$order + 2)*epiparameters$Ndiscr + 1)


lower_bounds <- constraints$r_lower
upper_bounds <- constraints$r_upper

#
# Surgery on Aeq, beq to produce a subset of full rank.
#
if (any (Aeq.dupes <- duplicated (Aeq)))
{
cat ("Warning: deleting", sum (Aeq.dupes), "redundant constraints.\n")
    beq <- beq[!Aeq.dupes]
    Aeq <- Aeq[!Aeq.dupes,]
    constraints$beq <- beq
    constraints$Aeq <- Aeq
}
#
# Set up the list of arguments, then use "do.call" to run auglag. The
# advantage of building the list is that we need not supply eval_g_eq or 
# eval_g_ineq if there are no corresponding constraints.
#
list.of.args <- list (x0 = r_in,  eval_f = obj.nloptr, 
        lb = lower_bounds, ub = upper_bounds, 
        softinfo = softinfo, epiparameters = epiparameters, caseinfo = caseinfo, c.out = c.out, 
        constraints = constraints, ...)

if (!is.null (A)){
    list.of.args$eval_g_ineq   <- less.const.function
}
if (!is.null (Aeq)) {
    list.of.args$eval_g_eq     <- equal.const.function
}
#
# Options for the optimizer are set up in setup.optargs.
#
list.of.args$opts <- setup.optargs (length (r_in), opts, local.opts)
#
# The big call.
#
result <- do.call ("nloptr", list.of.args)

#
# Add stuff for return value
#
result$softinfo <- softinfo
result$epiparameters <- epiparameters
result$caseinfo <- caseinfo
result$c.out <- c.out
result$x <- data
result$opts <- opts

if (return.constraints) result$constraints <- constraints
#
# Call postproc() to compute the density estimates (and do other 
# stuff, like drawing pictures). Add class info and return.
# 
result <- postproc(result, postproc.controls) # add some stuff
class (result) <- c("episplineDensity", "nloptr")
return (result)
}
