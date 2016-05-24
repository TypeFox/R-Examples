JAGScontrol <- function(variables, n.iter = 1000, thin = 1, burn.in = 0, seed,
                        rng = c("base::Wichmann-Hill", "base::Marsaglia-Multicarry",
                          "base::Super-Duper", "base::Mersenne-Twister")) {
  rng <- match.arg(rng)
  c(list(variables = variables, n.iter = n.iter, thin = thin, burn.in = burn.in, 
         RNG = c(list(".RNG.name" = rng),
           if (!missing(seed)) list(".RNG.seed" = as.integer(seed)))))
}

JAGScall <- function(model, y, prefix, control, ...) UseMethod("JAGScall")

JAGScall.BMMsetup <- function(model, y, prefix, control, ...) {
  dummy <- model
  model <- list(k = 2, priors = BMMpriors(y = y), inits = "initsFS",
                aprioriWeights = 1, restrict = "none", no.empty.classes = FALSE)
  n <- names(dummy)
  s <- names(model)
  p <- pmatch(n, s)
  if(any(is.na(p)))
    stop(paste("\nInvalid name(s) in model :", paste(n[is.na(p)], collapse=" ")))
  names(dummy) <- s[p]
  for (i in names(dummy)) {
    model[[i]] <- dummy[[i]]
  }
  model <- BMMmodel(y, model$k, model$priors, model$inits,
                    model$aprioriWeights, model$no.empty.classes, model$restrict, ...)
  if (!inherits(model, "BMMmodel")) stop("Model not specified correctly")
  JAGScall(model, y, prefix, control)
}

JAGScall.default <- function(model, y, prefix, control, ...) {
  if (!inherits(model, "JAGSmodel")) stop("Only for use with 'JAGSmodel' objects!")
  if (!is.null(control$RNG)) model$inits <- c(model$inits, control$RNG)
  if (length(model$bugs) > 1) model$bugs <- paste(model$bugs, collapse = prefix)
  FILE <- paste(prefix, "bug", sep = ".")
  write(model$bugs, file = FILE)
  JAGSmodel <- jags.model(FILE, inits = model$inits, data = model$data)
  if (control$burn.in > 0) jags.samples(JAGSmodel, control$variables, control$burn.in, control$burn.in)
  results <- as.mcmc(coda.samples(JAGSmodel, control$variables, control$n.iter, control$thin))
  index <- grep("tau", colnames(results))
  variables <- unique(sapply(colnames(results), function(x)
                             strsplit(x, "\\[")[[1]][1]))
  results[, index] <- 1/results[, index]
  colnames(results) <- sub("tau", "sigma2", colnames(results))
  variables <- sub("tau", "sigma2", variables)
  list(results = results, model = model, variables = variables)
}




