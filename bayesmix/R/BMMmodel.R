BMMmodel <-
function(y, k, priors, inits = "initsFS", aprioriWeights = 1, no.empty.classes = FALSE, restrict = "none", ...) {
  if (missing(y)) {
    call <- match.call(expand.dots = TRUE)
    model <- as.list(call[-1])
    model <- lapply(model, eval)
    class(model) <- "BMMsetup"
  }
  else {
    if (!is.null(dim(y))) {
      if (dim(y)[1] == 1) y <- y[1, ]
      else if (dim(y)[2] == 1) y <- y[ ,1]
      else stop("Only univariate data allowed")
    }
    y <- as.numeric(y)
    N <- as.numeric(length(y))
    
    if (!inherits(priors, "BMMpriors") & is.list(priors)) priors <- BMMpriors(priors, y = y)
    if (!inherits(priors, "BMMpriors")) stop("Priors not specified correctly")
    if (!all(sapply(priors$var, length) %in% c(0, 1, k))) stop("Priors not specified correctly - dimension differ")
    if (is.character(inits)) {
      inits <- get(inits)(y, k, restrict, ...)
    }
    for (i in names(inits)[names(inits) %in% names(priors$var)]) {
      if (length(inits[[i]]) != length(priors$var[[i]])) stop("Priors and inits dimension differ!")
    }
    if (length(aprioriWeights) != k) e <- rep(aprioriWeights[1], k) else e <- aprioriWeights
    model <- list()
    model$inits <- priors$inits <- inits
    if (priors$name[1] == "condconjugate") priors$var$B <- rep(NA, length(inits$tau*priors$var$B0inv))
    priors$inits <- priors$inits[names(priors$inits) %in% c("mu", "eta", "tau")]
    priors$var <- priors$var[!names(priors$var) %in% names(priors$inits)]
  
    index <- sapply(priors$var, function(x) (length(x) > 0) && !is.na(x))
    const <-  priors$var[index]
    const <- c(const, k = k, N = N)
    var <- c(priors$inits, priors$var[!index], list(e = e), list(y = y), list(S = y))  
    if (no.empty.classes) {
      var <- c(list(ind = matrix(0, nrow = N, ncol = k)),
               list(tot = vector(length = k)),
               list(seg = diag(k)),
               var)
      const <- c(const, list(Itot = rep(1, k)))
    }
    varlist <- varSpec(c(const, var))
    
    bugs <- paste("var \n", varlist, sep = "")
    restrict <- match.arg(restrict, c("none", "mu", "tau"))
    if (restrict == "mu") {
      bugs <- paste(bugs,
                    "model\t{\n\tfor (i in 1:N) {\n\t\ty[i] ~ dnorm(mu,tau[S[i]]);\n\t\tS[i] ~ dcat(eta[]);\n\t}\n",
                    sep = "")
    }
    else if (restrict == "tau") {
      bugs <- paste(bugs,
                    "model\t{\n\tfor (i in 1:N) {\n\t\ty[i] ~ dnorm(mu[S[i]],tau);\n\t\tS[i] ~ dcat(eta[]);\n\t}\n",
                    sep = "")
    }
    else {
      bugs <- paste(bugs,
                    "model\t{\n\tfor (i in 1:N) {\n\t\ty[i] ~ dnorm(mu[S[i]],tau[S[i]]);\n",
                    "\t\tS[i] ~ dcat(eta[]);\n\t}\n", sep = "")
    }
    bugs <- paste(bugs, modelPriors(priors, restrict), sep = "")
    model$data <- c(const, list(e = e), list(y = y))

    if (no.empty.classes) {
      bugs <- paste(bugs, "\tfor (i in 1:N) {\n\t\tind[i,] <- seg[S[i],];\n\t}\n",
                    "\tfor (j in 1:k){\n\t\ttot[j] <- sum(ind[,j]);\n",
                    "\t\tItot[j] ~ dinterval(tot[j], 0);\n\t}\n", sep = "")
      model$data$seg <- diag(k)
      if (!"S" %in% names(model$inits)) {
        posterior <- matrix(model$inits$eta * stats::dnorm(rep(y, each = k), model$inits$mu, sqrt(1/model$inits$tau)),
                            ncol = k, byrow = TRUE)
        S <- max.col(posterior)
        model$inits$S <- S
      }
      if (any(tabulate(model$inits$S, k) == 0)) stop("Please provide a valid initialization of S.")
    }
    bugs <- paste(bugs,paste("\teta[] ~ ddirch(e[]);\n}\n"), sep = "")
    if (length(priors$name) > 1) {
      if (priors$name[2] == "tau") {
        if(!all(c(model$data$g0G0Half,model$data$g0Half) > 0) & !("S0" %in% names(model$inits)))
          stop("Priors not specified correctly: Need an initial value for S0 with improper hierarchical prior.")
      }
      else stop("Should not be possible to have an hierarchical prior other than tau")
    }
    
    model$bugs <- bugs
    class(model) <- c("BMMmodel", "JAGSmodel")
  }
  model
}

modelParameters <-
function(priors) {
  parlist <- NULL
  for (i in names(priors$var)) {
    if (!any(is.na(priors$var[[i]]))) {
      cc <- priors$var[[i]]
      if (length(cc) > 1) {
        for (j in seq_along(cc)) {
          parlist <- paste(parlist, "\t", i,"[",j,"] <- ",cc[j],";\n ", sep = "")
        }
      }
    }
  }
  parlist
}

modelPriors <-
function(priors, restrict) {
  mu <- b0 <- B0inv <- B <- tau <- nu0Half <- nu0S0Half <- S0 <- g0Half <- g0G0Half <- numeric(0)
  variants <- c("independence", "condconjugate")
  variant <- match.arg(tolower(priors$name[1]), variants)
  var <- c(priors$var, priors$inits)
  if (restrict == "tau") var$tau <- rep(NA, 1)
  if (restrict == "mu") var$mu <- rep(NA, 1)
  for (i in seq_along(var)) {
    if (length(var[[i]]) > 1) assign(names(var[i]), paste(names(var[i]),"[j]", sep = ""))
    else assign(names(var[i]), names(var[i]))
  }
  pr <- NULL
  if (variant == "independence") {
    pr <- c(pr, paste(mu, " ~ dnorm(", b0,",", B0inv, ");\n", sep = ""))
  }
  else if (variant == "condconjugate") {
    pr <- c(pr, paste(mu, " ~ dnorm(", b0,",", B, ");\n", sep = ""))
    pr <- c(pr, paste(B, " <- ", B0inv , "*", tau, ";\n", sep = ""))
  }
  pr <- c(pr, paste(tau ," ~ dgamma(", nu0Half, ",", nu0S0Half,");\n", sep = ""))

  if(length(priors$name) > 1 && priors$name[2] == "tau") {
    pr <- c(pr, paste(S0," ~ dgamma(", g0Half, ",",g0G0Half,");\n", sep = ""))
    pr <- c(pr, paste(nu0S0Half , " <- ", nu0Half ," * ", S0, ";\n", sep = ""))
  }

  priorsSpec <- paste("\tfor (j in 1:k) {\n\t\t",
                      paste(pr[grep("j", pr)], collapse = "\t\t"), "\t}\n\t",
                      paste(pr[-grep("j", pr)], collapse = "\t"), "\n",
                      sep = "")
  priorsSpec
}

varSpec <-
function(var)  {
  var <- var[c(which(names(var) != "S"), which(names(var) == "S"))]
  varlist = NULL
  for (i in names(var)) {
    if (!is.null(dim(var[[i]]))) {
      cc <- dim(var[[i]])
      varlist = paste(varlist, "\t", i, "[", cc[1],",", cc[2], "],\n", sep = "")
    }
    else {
      cc <- length(var[[i]])
      if (i == "S") {
        varlist = paste(varlist, "\t", i, "[", cc,"]; \n\n", sep = "")
      }
      else if (cc <= 1) {
        varlist = paste(varlist, "\t", i,",\n ", sep = "")
      }
      else{
        varlist = paste(varlist, "\t", i, "[",cc,"], \n", sep = "")
      }
    }
  }
  varlist
}

print.JAGSmodel <- function(x, ...) {
  cat("Data for nodes: ", paste(names(x$data), collapse = ", "), "\n", sep ="")
  cat("Initial values for nodes: ", paste(names(x$inits), collapse = ", "), "\n\n", sep ="")
  cat("Model specification in BUGS language:\n\n")
  cat(x$bugs)
}

