# Generate and fit all 1-way, 2-way, ... k-way terms in a glm

Kway <- function(formula, family=poisson, data, ..., order=nt, prefix="kway") {

   if (is.character(family)) 
       family <- get(family, mode = "function", envir = parent.frame())
   if (is.function(family)) 
       family <- family()
   if (is.null(family$family)) {
       print(family)
       stop("'family' not recognized")
   }
   if (missing(data)) 
        data <- environment(formula)

   models <- list()
   mod <- glm(formula, family=family, data, ...)
   mod$call$formula <- formula
   terms <- terms(formula)
   tl <- attr(terms, "term.labels")
   nt <- length(tl)
   models[[1]] <- mod
   for(i in 2:order) {
       models[[i]] <- update(mod, substitute(.~.^p, list(p = i)))
   }      # null model
   mod0 <- update(mod, .~1)
   models <- c(list(mod0), models)
   names(models) <- paste(prefix, 0:order, sep = ".")
   class(models) <- "glmlist"
   models
}

