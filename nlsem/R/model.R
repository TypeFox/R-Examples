# model.R
#
# created Sep/23/2014, NU
# last mod Sep/04/2015, NU

#--------------- main functions ---------------

# Define model specification for different SEMs with nonlinear effects;
# possible objects classes are 'singleClass', 'semm', 'nsemm'; exported function
specify_sem <- function(num.x, num.y, num.xi, num.eta, xi, eta,
                        constraints=c("indirect", "direct1", "direct2"), num.classes=1,
                        rel.lat="default", interaction="none"){

  # check arguments
  if (!is.numeric(num.x) || !is.numeric(num.y) || !is.numeric(num.xi)
      || !is.numeric(num.eta) || !is.numeric(num.classes)) {
    stop("Number of variables or classes must be numeric.")
  } else if (num.x < num.xi || num.y < num.eta) {
    stop("The model contains not enough observed variables.")
  }
  stopifnot(num.x > 0, num.y >= 0, num.xi > 0, num.eta >= 0, num.classes > 0)

  # check if only defined xi's are in the interaction and eta is in there
  # for num.eta > 1
  if (interaction != "none") {
    interact.matrix <- calc_interaction_matrix(unlist(strsplit(interaction, ",")))
    if (max(interact.matrix) > num.xi) {
      stop("Interaction effects contain more xi's than defined.")
    }
    if (num.eta > 1) {
      if (!grepl("eta", interaction)) {
        stop("For more than one eta, specification for interaction must be something like eta1~xi1:xi2.")
      }
    }
  }

  # class of model
  model.class <- get_model_class(num.classes, interaction)

  # check latent variables and get indices for Lambda.x and Lambda.y
  xi.s <- unlist(strsplit(xi, ","))
  if (length(xi.s) != num.xi) {
    stop("Number of xi's and assignation of x's to xi's does not match. See ?specify_sem.")
  }
  xi.ind <- list()
  for (i in seq_len(num.xi)) {
    xi.ind[[i]] <- grep_ind(xi.s[i])
    if (max(xi.ind[[i]]) > num.x) {
      stop("Number of x's assinged to xi exceeds x's specified. See ?specify_sem.")
    }
  }

  eta.s <- unlist(strsplit(eta, ","))
  if (length(eta.s) != num.eta) {
    stop("Number of eta's and assignation of y's to eta's does not match. See ?specify_sem.")
  }
  eta.ind <- list()
  for (i in seq_len(num.eta)) {
    eta.ind[[i]] <- grep_ind(eta.s[i])
    if (max(eta.ind[[i]]) > num.y) {
      stop("Number of y's assinged to eta exceeds y's specified. See ?specify_sem.")
    }
  }

  # create matrices with default constraints
  # Lambda.x
  Lambda.x <- matrix(0, nrow=num.x, ncol=num.xi)
  for (i in seq_len(num.xi)){
    Lambda.x[xi.ind[[i]], i] <- c(1, rep(NA, length(xi.ind[[i]]) - 1))
  }
  # Lambda.y
  Lambda.y <- matrix(0, nrow=num.y, ncol=num.eta)
  for (i in seq_len(num.eta)){
    Lambda.y[eta.ind[[i]], i] <- c(1, rep(NA, length(eta.ind[[i]]) - 1))
  }
  # Gamma and Beta
  if (rel.lat == "default"){
    Gamma <- matrix(nrow=num.eta, ncol=num.xi)
    Beta  <- diag(num.eta)
  }
  else {
    GB    <- rel_lat(rel.lat, num.eta=num.eta, num.xi=num.xi)
    Gamma <- tryCatch({ GB[[grep("G", names(GB))]] },
                        error=function(e) matrix(nrow=num.eta, ncol=num.xi) )
    Beta  <- tryCatch({ GB[[grep("B", names(GB))]] },
                        error=function(e) diag(num.eta))
  }
  # Theta.d
  Theta.d <- diag(NA, nrow=num.x)
  # Theta.e
  Theta.e <- diag(NA, nrow=num.y)
  # Psi
  Psi <- matrix(NA, nrow=num.eta, ncol=num.eta)
  # Psi must be symmetrical, upper.tri = lower.tri -> in fill_model
  Psi[upper.tri(Psi)] <- 0
  # Phi
  Phi <- matrix(NA, nrow=num.xi, num.xi)
  # Phi must be symmetrical, upper.tri = lower.tri -> in fill_model
  Phi[upper.tri(Phi)] <- 0
  # nu's
  nu.x <- matrix(NA, nrow=num.x, ncol=1)
  for (i in seq_len(num.xi)) nu.x[xi.ind[[i]][1]] <- 0
  nu.y <- matrix(NA, nrow=num.y, ncol=1)
  for (i in seq_len(num.eta)) nu.y[eta.ind[[i]][1]] <- 0
  # alpha
  alpha <- matrix(NA, nrow=num.eta, ncol=1)
  # tau
  tau <- matrix(NA, nrow=num.xi, ncol=1)
  # Omega
  Omega <- matrix(0, nrow=num.xi, ncol=num.xi)
  if (interaction != "none"){

    interaction.s <- unlist(strsplit(interaction, ","))
    eta.logical <- matrix(nrow=num.eta, ncol=length(interaction.s))

    for (i in seq_len(num.eta)){
      eta.logical[i,] <- grepl(paste0("eta[",i,"]"), interaction.s)
    }

    which.eta <- apply(eta.logical, 1, which)

    if (nrow(eta.logical) == 1) {
      ind <- calc_interaction_matrix(interaction.s)
      Omega[ind] <- NA
      # check if Omega has row echelon form
      test_omega(Omega)

    } else {
      Omega <- array(0, dim=c(num.xi, num.xi, num.eta))
      for (i in seq_len(num.eta)) {
        eta.row <- which(eta.logical[i,])
        ind <- calc_interaction_matrix(interaction.s[eta.row])
        Omega[,,i][ind] <- NA
        test_omega(Omega[,,i])
      }
    }
  }

  # make a list of the matrices for each class
  matrices <- list()
  for (c in seq_len(num.classes)) {
    if (model.class == "singleClass") {
      matrices[[c]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                            Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                            Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                            nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                            Omega=Omega)
    } else if (model.class == "semm") {
      matrices[[c]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                            Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                            Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                            nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau)
    } else {
      matrices[[c]] <- list(Lambda.x=Lambda.x, Lambda.y=Lambda.y,
                            Gamma=Gamma, Beta=Beta, Theta.d=Theta.d,
                            Theta.e=Theta.e, Psi=Psi, Phi=Phi,
                            nu.x=nu.x, nu.y=nu.y, alpha=alpha, tau=tau,
                            Omega=Omega)
    }
  }
  names(matrices) <- paste0("class", seq_len(num.classes))

  # class weights w
  w <- matrix(1/num.classes, nrow=num.classes, ncol=1)

  constraints <- match.arg(constraints)

  # create model with matrices and info
  model <- list(matrices=matrices, info=list(num.xi=num.xi, num.eta=num.eta,
                                             num.x=num.x, num.y=num.y,
                                             constraints=constraints,
                                             num.classes=num.classes,
                                             par.names=list(), w=w))

  class(model) <- model.class
  # add parameter names to model
  model$info$par.names <- get_parnames(model=model, constraints=constraints)
  # bounds for parameters (variances > 0)
  model$info$bounds <- bounds(model)

  model
}

# Create model matrices from a dataframe with columns label (for parameter
# labels) and class 1 to class n; only needed when user wants to have full
# control over constraints, etc.; exported function
create_sem <- function(dat){

  stopifnot(is.data.frame(dat))

  Lambda.x <- as.character(dat$label[grep("Lambda.x", dat$label)])
  Lambda.y <- as.character(dat$label[grep("Lambda.y", dat$label)])
  Gamma    <- as.character(dat$label[grep("Gamma", dat$label)])
  Beta     <- as.character(dat$label[grep("Beta", dat$label)])
  Theta.d  <- as.character(dat$label[grep("Theta.d", dat$label)])
  Theta.e  <- as.character(dat$label[grep("Theta.e", dat$label)])
  Psi      <- as.character(dat$label[grep("Psi", dat$label)])
  Phi      <- as.character(dat$label[grep("Phi", dat$label)])
  nu.x     <- as.character(dat$label[grep("nu.x", dat$label)])
  nu.y     <- as.character(dat$label[grep("nu.y", dat$label)])
  alpha    <- as.character(dat$label[grep("alpha", dat$label)])
  tau      <- as.character(dat$label[grep("tau", dat$label)])
  Omega    <- as.character(dat$label[grep("Omega", dat$label)])

  # number of latent and indicator variables and classes
  num.x       <- length(nu.x)
  num.y       <- length(nu.y)
  num.xi      <- length(tau)
  num.eta     <- length(alpha)
  num.classes <- ncol(dat) - 1

  # create matrices
  matrices <- list()
  for (c in seq_len(num.classes)) {
    Lambda.x.matrix <- matrix(dat[dat$label %in% Lambda.x, paste0("class",c)],
                              nrow=num.x, ncol=num.xi)
    Lambda.y.matrix <- matrix(dat[dat$label %in% Lambda.y, paste0("class",c)],
                              nrow=num.y, ncol=num.eta)
    Gamma.matrix    <- matrix(dat[dat$label %in% Gamma, paste0("class",c)],
                              nrow=num.eta, ncol=num.xi)
    Beta.matrix     <- matrix(dat[dat$label %in% Beta, paste0("class",c)],
                              nrow=num.eta, ncol=num.eta)
    Theta.d.matrix  <- matrix(dat[dat$label %in% Theta.d, paste0("class",c)],
                              nrow=num.x, ncol=num.x)
    Theta.e.matrix  <- matrix(dat[dat$label %in% Theta.e, paste0("class",c)],
                              nrow=num.y, ncol=num.y)
    Psi.matrix      <- matrix(dat[dat$label %in% Psi, paste0("class",c)],
                              nrow=num.eta, ncol=num.eta)
    Phi.matrix      <- matrix(dat[dat$label %in% Phi, paste0("class",c)],
                              nrow=num.xi, ncol=num.xi)
    nu.x.matrix     <- matrix(dat[dat$label %in% nu.x, paste0("class",c)],
                              nrow=num.x, ncol=1)
    nu.y.matrix     <- matrix(dat[dat$label %in% nu.y, paste0("class",c)],
                              nrow=num.y, ncol=1)
    alpha.matrix    <- matrix(dat[dat$label %in% alpha, paste0("class",c)],
                              nrow=num.eta, ncol=1)
    tau.matrix      <- matrix(dat[dat$label %in% tau, paste0("class",c)],
                              nrow=num.xi, ncol=1)
    Omega.matrix    <- matrix(dat[dat$label %in% Omega, paste0("class",c)],
                              nrow=num.xi, ncol=num.xi)

    if (all(is.na(Omega.matrix))) {
      matrices[[c]] <- list(Lambda.x=Lambda.x.matrix,
                            Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                            Beta=Beta.matrix, Theta.d=Theta.d.matrix,
                            Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                            Phi=Phi.matrix, nu.x=nu.x.matrix,
                            nu.y=nu.y.matrix, alpha=alpha.matrix,
                            tau=tau.matrix)
    } else {
      matrices[[c]] <- list(Lambda.x=Lambda.x.matrix,
                            Lambda.y=Lambda.y.matrix, Gamma=Gamma.matrix,
                            Beta=Beta.matrix, Theta.d=Theta.d.matrix,
                            Theta.e=Theta.e.matrix, Psi=Psi.matrix,
                            Phi=Phi.matrix, nu.x=nu.x.matrix,
                            nu.y=nu.y.matrix, alpha=alpha.matrix,
                            tau=tau.matrix, Omega=Omega.matrix)
    }
  }
  names(matrices) <- paste0("class", 1:num.classes)

  w <- matrix(1/num.classes, nrow=num.classes, ncol=1)

  constraints <- attr(dat, "constraints")

  info <- list(num.xi=num.xi, num.eta=num.eta, num.x=num.x, num.y=num.y,
               constraints=constraints, num.classes=num.classes,
               w=w)

  model <- list(matrices=matrices, info=info)

  if (all(is.na(Omega.matrix))) {
    interaction <- "none"
  } else interaction <- "not_empty"

  class(model) <- get_model_class(num.classes, interaction)

  model$info$par.names <- get_parnames(model, constraints=constraints)
  model$info$bounds <- bounds(model)

  model
}

# Count free parameters of a model created with specify_sem (i.e. NAs in
# the model are counted); exported function
count_free_parameters <- function(model) {

  if (class(model) == "singleClass") {
    constraints <- "direct1"
  } else {
    constraints <- model$info$constraints
  }
  switch(EXPR = constraints,

    indirect = {

      res <- sum(unlist(lapply(model$matrices$class1, is.na)))
      parnames <- c("Phi", "tau")
      num.c <- 0
      for (class in names(model$matrices[-1])) {
        for (parname in parnames) {
          num.c <- num.c +
            sum(unlist(lapply(model$matrices[[class]][[parname]], is.na)))
        }
      }
      res <- res + num.c
    },

    direct1 = {

      res <- 0
      for (class in names(model$matrices))
        res <- res + sum(unlist(lapply(model$matrices[[class]], is.na)))

    },

    direct2 = {

      res <- sum(unlist(lapply(model$matrices$class1, is.na)))

      if (class(model) == "semm") {
        parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau")
      } else {
        parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau", "Omega")
      }

      num.c <- 0
      for (class in names(model$matrices[-1])) {
        for (parname in parnames) {
          num.c <- num.c +
            sum(unlist(lapply(model$matrices[[class]][[parname]], is.na)))
        }
      }
      res <- res + num.c
    }
  )     # end of switch

  res
}

# Fill a model created with specify_sem with parameters given as a vector;
# mostly needed to simulate data from a prespecified model; exported function
fill_model <- function(model, parameters) {

  stopifnot(class(model) == "singleClass" || class(model) == "semm"
            || class(model) == "nsemm")

  stopifnot(count_free_parameters(model) == length(parameters))

  matrices <- model$matrices

  if (class(model) == "singleClass") {
    constraints <- "direct1"
  } else {
    constraints <- model$info$constraints
  }
  switch(EXPR = constraints,

    indirect = {

      parnames <- c("Phi", "tau")
      ind.list <- list()
      num.list <- list()
      for (class in names(matrices)) {
        for (parname in parnames) {
          ind.list[[class]][[parname]] <- is.na(matrices[[class]][[parname]])
          num.list[[class]][[parname]] <-
            length(matrices[[class]][[parname]][is.na(matrices[[class]][[parname]])])
          }
      }

      for (i in seq_along(matrices$class1)) {
        matrix.i <- matrices$class1[[i]]
        # number of NA's in matrix
        num.na <- length(matrix.i[is.na(matrix.i)])
          if (num.na > 0) {
            matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
            parameters <- parameters[-(1:num.na)]
            for (class in names(matrices)) {
              matrices[[class]][[i]] <- matrix.i
            }
          }
      }
      for (class in names(matrices)[-1]) {
        for (parname in parnames) {
          ind <- ind.list[[class]][[parname]]
          num <- num.list[[class]][[parname]]
          if (num > 0) {
            matrices[[class]][[parname]][ind] <- parameters[1:num]
            parameters <- parameters[-(1:num)]
          }
        }
      }
    },

    direct1 = {

      for (class in names(matrices)) {
        for (i in seq_along(matrices[[class]])) {
          matrix.i <- matrices[[class]][[i]]
          # number of NA's in matrix
          num.na <- length(matrix.i[is.na(matrix.i)])
          if (num.na > 0) {
            matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
            parameters <- parameters[-(1:num.na)]
            matrices[[class]][[i]] <- matrix.i
          }
        }
      }

    },

    direct2 = {

      if (class(model) == "semm") {
        parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau")
      } else {
        parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau", "Omega")
      }
      ind.list <- list()
      num.list <- list()
      for (class in names(matrices)) {
        for (parname in parnames) {
          ind.list[[class]][[parname]] <- is.na(matrices[[class]][[parname]])
          num.list[[class]][[parname]] <-
            length(matrices[[class]][[parname]][is.na(matrices[[class]][[parname]])])
        }
      }

      for (i in seq_along(matrices$class1)) {
        matrix.i <- matrices$class1[[i]]
        # number of NA's in matrix
        num.na <- length(matrix.i[is.na(matrix.i)])
          if (num.na > 0) {
            matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
            parameters <- parameters[-(1:num.na)]
            for (class in names(matrices)) {
              matrices[[class]][[i]] <- matrix.i
            }
          }
      }
      for (class in names(matrices)[-1]){
        for (parname in parnames) {
          ind <- ind.list[[class]][[parname]]
          num <- num.list[[class]][[parname]]
          if (num > 0) {
            matrices[[class]][[parname]][ind] <- parameters[1:num]
            parameters <- parameters[-(1:num)]
          }
        }
      }
    }
  )     # end of switch

  # fill symmetric matrices
  for (class in names(matrices)) {
    tryCatch({matrices[[class]]$Phi <-
      fill_symmetric(matrices[[class]]$Phi)}, error=function(e) e,
      warning=function(w) w)
    tryCatch({matrices[[class]]$Psi <-
      fill_symmetric(matrices[[class]]$Psi)}, error=function(e) e,
      warning=function(w) w)
  }

  out <- list(matrices=matrices, info=model$info)
  class(out) <- class(model)
  out
}

#--------------- helper functions ---------------

# all NOT exported

# Check if model or matrices are filled
check_filled <- function(x) {
    if (anyNA(unlist(x))) stop("Model is not filled.")
}
# TODO What's going on with this function? Why does it not work with
# indirect=TRUE??

# Grep indices for Lambda matrices from input that defines which indicators
# are asociated with which latent variable
grep_ind <- function(x){
  tryCatch({
    if (length(unlist(strsplit(x, "-"))) > 1){
        as.numeric(gsub("^.*[x,y]([0-9]+).*[x,y]([0-9]+)$", "\\1",
        x)):as.numeric(gsub("^.*[x,y]([0-9]+).*[x,y]([0-9]+)$", "\\2", x))
    } else {
        as.numeric(gsub("^.([0-9]+).*$", "\\1", x))
    }
  }, warning = function(war) {
    stop("Wrong input for specifying exogenous or endogonous latent variables (xi or etas). See ?specify_sem.")
    # this might never be evaluated, error catched earlier
  })
}

# Returns matrix which specifies which latent variables interact with each
# other
calc_interaction_matrix <- function(x){
  tryCatch({
    rows <- as.numeric(gsub("^.*xi([0-9]+):xi[0-9]+$", "\\1", x))
    cols <- as.numeric(gsub("^.*xi.*:xi([0-9]+)$", "\\1", x))
    mat  <- sort_interaction_effects(rows, cols)
    mat
  }, warning = function(war) {
    stop("Wrong input for interaction. See ?specify_sem.")
  }, error = function(err) { # perhaps error catching is unnecessary
    stop("Wrong input for interaction. See ?specify_sem.")
  })
}

# Ensures that interaction effects are in the correct order when passed to
# Omega
sort_interaction_effects <- function(rows, cols){

  for (i in seq_along(rows)){
    if (rows[i] > cols[i]){
      rows_i <- rows[i]
      cols_i <- cols[i]
      rows[i] <- cols_i
      cols[i] <- rows_i
    }
  }
  cbind(rows, cols)
}

# Tests if input for Omega is in the correct format; Omega needs to be in
# row echelon form; returns nothing if Omega has the correct form
test_omega <- function(Omega){

  if (anyNA(Omega)){

    interactions <- Omega
    #diag(interactions) <- 0
    ind <- which(is.na(interactions), arr.ind=TRUE)
    dim <- nrow(interactions)
    msg <- "Interactions are not well-defined. Please change order of xi's. See ?specify_sem for details."
    # TODO This error message does not make sense for # eta > 1. Maybe get
    # rid of test_omega all together? It's more a debugging tool than
    # anything else anyway.

    # test if any rows are 0 in between interaction effects
    na.r <- rowSums(interactions)
    for (i in seq_along(na.r)){
      if (!is.na(na.r[i]))
        if (is.na(sum(na.r[-c(1:i)])))
          stop(msg)
    }

    if (max(ind[,1]) > 1){
      for (i in 2:max(ind[,1])){
        if (min(which(is.na(interactions[i-1,]))) >=
        min(which(is.na(interactions[i,])))) {
          stop(msg)
        }
      }
    }
    invisible(NULL)
  } else {
      invisible(NULL)
  }
}

# Creates Beta and Gamma matrices according to the input obtained by
# rel.lat; matrices define relationships for latent variables except for
# interaction effects
rel_lat <- function(x, num.eta, num.xi){

  error.msg <- "Latent variables misspecified. Must be of the form 'eta1~xi1' or 'eta2~eta1'. See ?specify_sem for details."

  x.s       <- unlist(strsplit(x, ","))
  which.xi  <- which(grepl("xi", x.s))
  which.eta <- which(!grepl("xi", x.s))

  if (length(which.xi) == 0){
      G <- matrix(NA, nrow=num.eta, ncol=num.xi)
  } else {
    G <- matrix(0, nrow=num.eta, ncol=num.xi)

    for (i in which.xi){
      xi.s <- unlist(strsplit(x.s[i], "~"))
      if (length(xi.s) < 2) stop(error.msg)

      etas <- unlist(strsplit(xi.s[1], "[+]"))
      xis  <- unlist(strsplit(xi.s[2], "[+]"))
      tryCatch({
        ind.xi <- as.numeric(gsub("^.*xi([0-9]+).*$", "\\1", xis))
        ind.eta <- as.numeric(gsub("^.*eta([0-9]+).*$", "\\1", etas))
        G[ind.eta, ind.xi] <- NA
      }, error = function(e) stop(error.msg)
      , warning = function(w) stop(error.msg)
      )
    }
  }

  B <- diag(1, num.eta)
  for (i in which.eta){
    eta.s <- unlist(strsplit(x.s[i], "~"))
    if (length(eta.s) < 2) stop(error.msg)

    eta.rows <- unlist(strsplit(eta.s[1], "[+]"))
    eta.cols <- unlist(strsplit(eta.s[2], "[+]"))
    tryCatch({
      ind.rows <- as.numeric(gsub("^.*eta([0-9]+).*$", "\\1", eta.rows))
      ind.cols <- as.numeric(gsub("^.*eta([0-9]+).*$", "\\1", eta.cols))
      if (any(sapply(ind.rows, function(x) {is.element(x, ind.cols)})))
        stop(error.msg)
      B[ind.rows, ind.cols] <- NA
    }, error = function(e) stop(error.msg), warning = function(w) stop(error.msg)
    )
  }
  out <- list(Gamma=G, Beta=B)
  out
}

# Defines model class of a given specification; possible output:
# 'singleClass', 'semm', 'nsemm'
get_model_class <- function(num.classes, interaction) {
  if (num.classes == 1) {
    model.class <- "singleClass"
  } else {
    if (interaction != "none") {
      model.class <- "nsemm"
    } else model.class <- "semm"
  }
  model.class
}

# Obtains parameter names from a given model; used in specify_sem
get_parnames <- function(model, constraints=c("indirect", "direct1",
                         "direct2")) {

  constraints <- match.arg(constraints)
  switch(EXPR = constraints,

    indirect = {
      parnames <- c("Phi", "tau")
    },

    direct1 = {
      lst <- unlist(lapply(model$matrices$class1, is.na))
      parnames <- names(lst[lst])
    },

    direct2 = {
      if (class(model) == "semm") {
        parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau")
      } else {
        parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau", "Omega")
      }
    }
  )     # end of switch

  par.names <- list()
  lst <- unlist(lapply(model$matrices$class1, is.na))
  par.names[["class1"]] <- names(lst[lst])
  if (class(model) == "nsemm" || class(model) == "semm") {
    for (class in names(model$matrices)[-1]) {
      if (constraints != "direct1") {
        ind.parnames <- unlist(lapply(parnames, FUN=function(x) grep(x,
          par.names$class1)))
        par.names[[class]] <- names(lst[lst])[ind.parnames]
      } else {
        lst <- unlist(lapply(model$matrices[[class]], is.na))
        par.names[[class]] <- names(lst[lst])
      }
    }
  } else {
    par.names <- par.names$class1
  }
  par.names
}

# Fills upper.tri of a (filled) matrix which should be symmetric
fill_symmetric <- function(mat) {
  for (i in seq_len(nrow(mat))) {
    for (j in i:ncol(mat)) {
        mat[i,j] <- mat[j,i]
    }
  }
  mat
}

# Vector for diagonal indices
diag_ind <- function(num) diag(matrix(seq_len(num^2), num))

# Set bounds for parameters to (0, Inf)
bounds <- function(model, constraints=c("indirect", "direct1", "direct2")) {

  if (model$info$num.x > 1){
      t.d <- paste0("Theta.d", diag_ind(model$info$num.x))
  } else t.d <- "Theta.d"
  if (model$info$num.y > 1){
      t.e <- paste0("Theta.e", diag_ind(model$info$num.y))
  } else t.e <- "Theta.e"
  if (model$info$num.eta > 1){
      psi <- paste0("Psi", diag_ind(model$info$num.eta))
  } else psi <- "Psi"
  if (model$info$num.xi > 1){
      phi <- paste0("Phi", diag_ind(model$info$num.eta))
  } else phi <- "Phi"

  if (class(model) == "singleClass") {

    lower <- rep(-Inf, count_free_parameters(model))
    upper <- rep(Inf, count_free_parameters(model))

    lower[model$info$par.names %in% c(t.d, t.e, psi, phi)] <- 0
    out <- list(upper=upper, lower=lower)

  } else if (class(model) == "semm" || class(model) == "nsemm") {

    lower.class <- list()
    upper.class <- list()

    for (class in names(model$matrices)) {

      lower.class[[class]] <- rep(-Inf, length(model$info$par.names[[class]]))
      upper.class[[class]] <- rep(Inf, length(model$info$par.names[[class]]))

      lower.class[[class]][model$info$par.names[[class]] %in% c(t.d, t.e, psi, phi)] <- 0
    }
    out <- list(upper=upper.class, lower=lower.class)
  }
  out
}

