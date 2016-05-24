##############################################################################
# A set of functions to help specify the model                               #
##############################################################################
# 
# Initial values shall be specified as named lists, with names:
# beta, tau, alpha, theta
#
# Constraints list may have only 2 names: lower and upper. Sublists of 
# these shall be named list with names: beta, alpha 
#
# Constraints & Initial values must be matrix objects & have:
#   beta  :: rownames = grp.labels, colnames = random effects vars 
#   tau   :: names = grp.labels
#   alpha :: names = fixed effects vars 
#   theta :: rownames = upper level covariates, colnames = rand eff covars
#
# Missing (unspecified) constraints will be taken +/- Inf (upper/lower)
# Missing (unspecified) initial values will be automatically assigned
###############################################################################

###############################################################################
# Specify number of sampling & burn-in iterations; and suggest sampler to use
###############################################################################
hbglm.sampler.control <- function(num.samples = -1, samp.factor = 50,
    sampler = c("slice"))
{
    if (samp.factor < 1)
        stop("Invalid 'samp.factor' arg: must be an integer >= 1")
    if (length(sampler) == 2 && sampler == c("slice", "sns"))
        sampler = "slice"
    if (! sampler %in% c("slice", "sns")) {
        warning("Unsupported MCMC sampler requested, using slice sampler")
        sampler = "slice"
    }
    list(num.samples = num.samples, samp.factor = samp.factor,
         sampler = sampler)
}

###########################################
# Specify constraints and initial values
###########################################
hbglm.model.control <- function(constraints = NULL,
    initializer = c("user", "regression"), user.init.val = NULL)
{
    # Handle constraints
    if (!is.null(constraints)) {
      if (!is.list(constraints))
        stop("Arg constraints to hbglm.model.control() must be a list.")
      if (all(is.null(constraints$upper), is.null(constraints$lower)))
        stop("Arg constraints to hbglm.model.control() missing lower & upper")
    }
  
    # Handle initializer
    if (!is.null(user.init.val))
         initializer = "user"  
    if (length(initializer) == 2 && initializer == c("user", "regression"))
          initializer <- "regression"  # default choice
    if (initializer == "user") { 
        if (is.null(user.init.val))
            stop("user.init.val can't be NULL when initializer = user.")
        if (!is.list(user.init.val))
            stop("user.init.val must be a list when initializer = user.")
    }
    else if (initializer != "regression") 
        stop("hbglm.model.control(): Arg initializer not recognized")

    scale.data = FALSE  # to be implemented
    list(scale.data    = scale.data,
         constraints   = constraints,
         initializer   = initializer,
         user.init.val = user.init.val)
}

###########################################################################
# Main function for setting initial values                              
###########################################################################
# 
# Shall only be called by hbglm()
# Args: 
#   family      - see details of family arg to hbglm()
#   model       - internal data structure of hbglm()
#   bounds      - the object returned by get.box.bounds()
#   init.vals   - data.frame in the format given in hbglm.model.control()
#                 (contains suggested initial values)
#   default     - a default value used to initialize all params not specified 
# Returns:
#   If suggested initial values satisfy constraints then they are accepted
#   else a suitable constraint satisfying initial value is chosen
#   
# Note:
#   1. Ensures that all constraints (in bounds) are satisfied
#   2. POLICY for choosing actual initial values: 
#          see the source code of function - choose.init()
#   3. NA is accepted as a missing entry 
###########################################################################
get.init.vals <- function(family, model, bounds, model.control, default=0.01)
{
    # User speficified or auto-generated initial values
    init.vals <- default.init.val(model.control, family, model, 
                                  def.val = default)
    tau <- init.vals$tau
    Sigma <- init.vals$Sigma
    
    # Make initial values data structure for beta, alpha, theta
    beta <- matrix(rep(NA, model$J * model$K), nrow = model$J)
    colnames(beta) <- model$rand.cov; rownames(beta) <- model$grp.labels
    alpha <- if (model$has.fixed) rep(NA, model$M) else NULL
    if (!is.null(tau)) names(alpha) <- model$fixed.cov
    theta <- NULL
    if (model$has.upper.level) {
        theta <- matrix(rep(NA, model$L * model$K), nrow = model$L)
        colnames(theta) <- model$rand.cov
        if (model$L > 1) rownames(theta) <- model$upper.cov
    }
   
    # Sort rows of data.frame according to vector
    # Returns NULL if some entry of vec isn't in rownames(df)
    sort.df.rows <- function(df, vec)
    {
        if (!all(vec %in% rownames(df))) return(NULL)
        od <- sapply(vec, function(x) which(x == rownames(df)))
        df <- df[od, , drop = FALSE]
        return(df)
    }

    # Replace defaults with suggested initial values 
    if (!is.null(init.vals)) {
        # Handle beta
        if (!is.null(init.vals$beta)) {
            init <- sort.df.rows(init.vals$beta, model$grp.labels)
            if (is.null(init)) stop("Incorrect initial val spec: missing rows")
            init.vars <- colnames(init)
            for (i in 1:length(init.vars)) {
                var <- init.vars[i]
                j <- which(model$rand.cov == var)
                beta[ , j] <- init[[var]]
            }     
        }
        # Handle alpha
        if (model$has.fixed && !is.null(init.vals$alpha)) {
            init <- init.vals$alpha
            if (!is.vector(init)) stop("Incorrect alpha init val: not a vector")
            init.vars <- names(init)
            if (!all(init.vars %in% model$fixed.cov)) 
                stop("Incorrect alpha init val: fixed effect not recognized")
            for (i in 1:length(init.vars)) {
                var <- init.vars[i]
                j <- which(model$fixed.cov == var)
                alpha[j] <- init[i]
            }
        }
        # Handle tau
        if (family$has.tau && !is.null(init.vals$tau)) {
            init <- init.vals$tau
            if (!is.vector(init)) stop("Incorrect tau init val: not a vector")
            init.vars <- names(init)
            if (!all(init.vars %in% model$grp.labels)) 
                stop("Incorrect tau init val spec: grp label not recognized")
            for (i in 1:length(init.vars)) {
                var <- init.vars[i]
                j <- which(model$grp.labels == var)
                tau[j] <- init[i]
            }
        }
        # Handle theta
        if (!is.null(init.vals$theta)) {
            if (model$upper.reg) {
              init <- sort.df.rows(init.vals$theta, model$upper.cov)
              if (is.null(init)) stop("Incorrect theta init val: missing rows")
            } else init <- init.vals$theta
            init.vars <- colnames(init)
            for (i in 1:length(init.vars)) {
                var <- init.vars[i]
                j <- which(model$rand.cov == var) 
                theta[ , j] <- init[[var]]
            }
        }
    } # Finished inserting suggested intial values
    
    # Function to select a single value within 2 bounds
    # This implements the 'POLICY' give in 'Note' in docs of get.init.val()
    choose.init <- function(suggest = NA, high = Inf, low = -Inf, default=0.01) 
    {
        eps <- min(default, (high - low)/2) # small +ve quantity  

        # No suggestion given - choose a value close to 0.01
        if (is.na(suggest) || is.null(suggest))
            return(min(high - eps, max(default, low + eps)))
        if (abs(suggest) == Inf)
            stop("Bad initial value suggested: found infinity")

        # If 'suggest' in bounds
        if (low <= suggest && suggest <= high) return(suggest)

        # If not in bounds 
        if (suggest < low) return(low + eps)
        else { # suggest > high
            val <- if (low <= 0) min(default, high - eps)
                   else (high + low)/2
            return(val)
        }
    }   
 
    # Reconcile initial values with constraints (given in 'bounds')
    for (j in 1:ncol(beta)) { # beta
      beta[ , j] <- sapply(1:nrow(beta), function(i) {
                                           choose.init(suggest=beta[i,j], 
                                      high = bounds$beta.hi[i,j], 
                                      low = bounds$beta.lo[i,j])})
    }
    if (family$has.tau) {
      for (i in 1:model$J) {# tau, constraints maybe provided by family
        tau[i] <- choose.init(suggest=tau[i], high=bounds$tau.hi[i], 
                                              low=bounds$tau.lo[i])
      }
    } 
    if (model$has.fixed) { # alha
      for (m in 1:model$M) { 
        alpha[m] <- choose.init(suggest=alpha[m], high = bounds$alpha.hi[m],
                                                  low = bounds$alpha.lo[m])
      }
    }
                  
    start <- list(
        beta  = beta,
        alpha = alpha,
        tau   = tau,
        theta = theta,
        Sigma = Sigma
    )
    return(start) 
}

##############################################################################
# Main function for handling constraints.
# Shall only be called by hbglm().
##############################################################################
# Users may only specify constraints on: beta, alpha
# family may specify constraints only on: tau
# Args: 
#   family      - see details of family arg to hbglm()
#   model       - internal data structure of hbglm()
#   user.constraints - data.frame in the format given in hbglm.model.control()
# Returns:
#   A list with hi/lo slots for each of: beta, alpha, tau
#   Makes sure constraints are stricter of user and family specifications
##############################################################################
get.box.bounds <- function(family, model,  user.constraints = NULL)
{
    nparam <-(model$J + model$L + 1) * model$K + model$J
    # Default constraints 
    beta  = list(lo = -Inf, hi = +Inf)
    alpha = list(lo = -Inf, hi = +Inf)
    tau   = list(lo = -Inf, hi = +Inf)

    # Apply family constraints
    if (!is.null(family$constraint)) {
        fc <- family$constraint
        #if (!is.null(fc$beta)) {
        #  if (!is.null(fc$beta$upper)) beta$hi = min(beta$hi, fc$beta$upper)
        #  if (!is.null(fc$beta$lower)) beta$lo = max(beta$lo, fc$beta$lower)
        #}     
        if (!is.null(fc$tau)) {
          if (!is.null(fc$tau$upper)) tau$hi = min(tau$hi, fc$tau$upper)
          if (!is.null(fc$tau$lower)) tau$lo = max(tau$lo, fc$tau$lower)
        }     
    }
    # Constraint data structure
    beta.lo <- matrix(rep(beta$lo, model$J * model$K), nrow = model$J) 
    beta.hi <- matrix(rep(beta$hi, model$J * model$K), nrow = model$J) 
    rownames(beta.lo) <- model$grp.labels; colnames(beta.lo) <- model$rand.cov
    rownames(beta.hi) <- model$grp.labels; colnames(beta.hi) <- model$rand.cov

    alpha.lo <- NULL; alpha.hi <- NULL
    if (model$has.fixed) {
        alpha.lo <- rep(alpha$lo, model$M); names(alpha.lo) <- model$fixed.cov 
        alpha.hi <- rep(alpha$hi, model$M); names(alpha.hi) <- model$fixed.cov
    }

    tau.lo <- NULL; tau.hi <- NULL
    if (family$has.tau) {
        tau.lo <- rep(tau$lo, model$J); names(tau.lo) <- model$grp.labels 
        tau.hi <- rep(tau$hi, model$J); names(tau.hi) <- model$grp.labels
    }

    # Sort rows of data.frame according to vector
    # Returns NULL if some entry of vec isn't in rownames(df)
    sort.df.rows <- function(df, vec) 
    {
        if (!all(vec %in% rownames(df))) return(NULL)
        od <- sapply(vec, function(x) which(x == rownames(df)))
        df <- df[od, , drop = FALSE]
        return(df)
    }
 
    # Helper function to handle user specified constraints
    insert.constraints <- function(cons.list, hi = TRUE) 
    {
        # Handle beta
        if (!is.null(cons.list$beta)) {
            cons <- sort.df.rows(cons.list$beta, model$grp.labels)
            if (is.null(cons)) stop("Incorrect beta constraint: missing rows")
            cons.vars <- colnames(cons)
            for (i in 1:length(cons.vars)) {
                var <- cons.vars[i]
                j <- which(model$rand.cov == var)
                if (hi) beta.hi[ , j] <- min(cons[[var]], beta.hi[ ,j])
                else    beta.lo[ , j] <- max(cons[[var]], beta.lo[ ,j])
            }
        }
        # Handle alpha 
        if (model$has.fixed && !is.null(cons.list$alpha)) {
            cons <- cons.list$alpha
            if (!is.vector(cons)) stop("Incorrect alpha constraint: not a vector")
            cons.vars <- names(cons)
            if (!all(cons.vars %in% model$fixed.cov))
                stop("Incorrect alpha constraint: fixed cov not recognized")
            for (i in 1:length(cons.vars)) {
                var <- cons.vars[i]
                j <- which(model$fixed.cov == var)
                if (hi) alpha.hi[j] <- min(cons[i], alpha.hi[j])
                else    alpha.lo[j] <- max(cons[i], alpha.lo[j])
            }
        }
    } # end helper function function 
    
    # Apply user constraints 
    if (!is.null(user.constraints)) {
        uc <- user.constraints
        if (!is.null(uc$lower)) 
          insert.constraints(uc$lower, hi = FALSE)
        if (!is.null(uc$upper)) 
          insert.constraints(uc$upper, hi = TRUE)
    }
    # Form the bounds into a vector in the order: beta, theta, tau, Sigma
    bounds <- list(
        user.const = user.constraints,
#        upper = c(as.vector(beta.hi), if(family$has.tau) tau.hi else NULL, if(model$has.fixed) alpha.hi else NULL),
#        lower = c(as.vector(beta.lo), if(family$has.tau) tau.lo else NULL, if(model$has.fixed) alpha.lo else NULL),
        beta.hi  = beta.hi,  beta.lo  = beta.lo,
        tau.hi   = tau.hi,   tau.lo   = tau.lo,
        alpha.hi   = alpha.hi,   alpha.lo   = alpha.lo
    )
    if (any(bounds$upper < bounds$lower))
        stop("Incorrect constraint spec: some upper bound < upper bound")
    return(bounds)
}

##############################################################################
# Generate default initial values (these maybe user specified)
##############################################################################
# Shall only be called by hbglm()
#
# Args: 
#   model         - internal data structure of hbglm()
#   model.control - list, see docs of hbglm() for more information
# Returns:
#   A list with keys 'beta', 'alpha', 'tau', 'theta', 'Sigma' & values as
#   data.frames specifying initial values. May have 'NA' in the data.frames  
##############################################################################
default.init.val <- function(model.control, family, model, def.val=0.01)
{
   mc <- model.control
   if (mc$initializer == "user") {
       if (is.null(mc$user.init.val))
           stop("user.init.val can't be NULL when initializer = user.")
       if (!is.list(mc$user.init.val))
           stop("user.init.val must be a list when initializer = user.")
       return(mc$user.init.val)
   } else if (mc$initializer == "regression") {
       return(compute.init.guess(model, family, def.val = def.val))
   } else
       stop("hbglm.model.control(): Arg initializer not recognized")
}

# Get initial values by doing a regression
compute.init.guess <- function(model, family, def.val=0.01)
{
  default <- def.val
  # Initialize 
  beta  = matrix(rep(def.val, model$J * model$K), nrow = model$J)
  tau   = if(family$has.tau) rep(def.val, model$J) else NULL
  alpha = if(model$has.fixed) rep(def.val, model$M) else NULL
  theta = if(model$has.upper.level) 
            matrix(rep(def.val, model$L * model$K), nrow = model$L) else NULL
  Sigma = if(model$has.upper.level) def.val * diag(model$K) else NULL
   
  # Perform regressions to compute parameter guesses
  if (F) {
    # beta
    glm.fam <- NULL
    if (is.character(family$fam.name)) 
      glm.fam <- get(family$fam.name, mode = "function")
    else 
      glm.fam <- family
    
    for (j in 1:model$J) {
      beta.j <- glm.fit(model$mat.lower.split[[j]]$X, 
                        model$mat.lower.split[[j]]$y,
                        family = glm.fam())$coefficients
      beta.j[is.na(beta.j)] <- default
      beta[j, ] <- beta.j
    }

    # theta
    for (k in 1:model$K) {
      theta.k <- lm.fit(model$mat.upper, beta[ , k])$coefficients
      theta.k[is.na(theta.k)] <- default
      theta[k, ] <- theta.k
    }
  }
  
  # Return appropriate initializer data structure 
  beta <- data.frame(beta)
  rownames(beta) <- model$grp.labels; colnames(beta) <- model$rand.cov
  if (!is.null(theta)) {
      theta <- data.frame(theta)
      if (model$L > 1) rownames(theta) <- model$upper.cov 
      colnames(theta) <- model$rand.cov
  }
  if (!is.null(tau)) names(tau) <- model$grp.labels 
  if (!is.null(alpha)) names(alpha) <- model$fixed.cov 
  if (!is.null(Sigma)) {
      rownames(Sigma) <- model$rand.cov 
      colnames(Sigma) <- model$rand.cov
  } 
  return(list(beta=beta, tau=tau, alpha=alpha, theta=theta, Sigma=Sigma))
}
##############################################################################
