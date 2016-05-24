# em.R
#
# last mod: Aug/27/2015, NU

# Performs EM-algorithm for different models of class 'singleClass',
# 'semm', and 'nsemm'
em <- function(model, data, start, qml=FALSE, verbose=FALSE,
               convergence=1e-02, max.iter=100, m=16, optimizer=c("nlminb",
               "optim"), max.mstep=1, max.singleClass=1, neg.hessian=TRUE,
               ...) {

  stopifnot(class(model) == "singleClass" || class(model) == "semm" ||
            class(model) == "nsemm")

  if (anyNA(data)) stop("Data contains NAs. Please remove.")

  if (class(model) == "nsemm" & neg.hessian == TRUE) {
    neg.hessian = FALSE
    warning("Negative Hessian cannot be computed for model of class 'nsemm'. neg.hessian will be set to FALSE.\n")
  }

  if (is.matrix(data)) {
    data <- data
  } else if (is.data.frame(data)) {
    data <- as.matrix(data)
  } else {
    stop("data need to be a matrix or a data frame.")
  }

  if (!count_free_parameters(model) == length(start)){
    stop("Number of starting parameters is not equal to number of free parameters in model.")
  }

  if (ncol(data) != (model$info$num.x + model$info$num.y)) {
    stop("Number of columns in data does not match number of x's and y's.")
  }

  # if (class(model) == "singleClass" || class(model) == "nsemm"){
  #     n.na <- length(which(is.na(model$matrices$class1$Omega)))
  #     if (any(start[-c(1:(length(start) - n.na))] == 0)){
  #         stop("Starting parameters for Omega should not be 0.")
  #     }
  # }

  if (anyNA(model$matrices$class1$Omega) && model$info$num.eta > 1){
    stop("Model with interaction effects and num.eta > 1 cannot be fitted (yet).")
  }

  if(verbose == TRUE) {
    cat("-----------------------------------\n")
    cat("Starting EM-algorithm for", class(model), "\n")
    cat(paste("Convergence: ", convergence, "\n"))
    cat("-----------------------------------\n")
    cat("-----------------------------------\n")
  }

  ll.ret   <- NULL
  num.iter <- 0     # number of iterations
  if (class(model) == "semm" || class(model) == "nsemm") {
    par.new <- start
  } else {
    par.new <- convert_parameters_singleClass(model, start)
  }
  ll.new <- 0

  run <- TRUE
  while(run) { # as long as no convergence is reached

    if (num.iter > 3){
      if (ll.new - ll.old > 0) {
        warning("Loglikelihood should be increasing.")
      }
    }

    if(verbose == TRUE) {
      cat(paste("Iteration", num.iter+1, "\n"))
      cat("Doing expectation-step \n")
    }

    # Update loglikelihood
    ll.old <- ll.new
    par.old <- par.new

    # E-step
    switch(class(model),
       "singleClass" = {
            names(model$matrices$class1)[grep("Phi", names(model$matrices$class1))] <- "A"
            # rename Phi to A, since LMS algorithm estimates A
            P <- estep_lms(model=model, parameters=par.old, dat=data,
                           m=m, ...)
        },
       "semm" = {
            P <- estep_semm(model=model, parameters=par.old, data=data)
            model$info$w <- colSums(P) / nrow(data)
            if (verbose == TRUE) {
                cat("Class weights: ", round(model$info$w, digits=4), "\n")
            }
        },
        "nsemm" = {
            res <- estep_nsemm(model=model, parameters=par.old,
                               data=data, max.singleClass=max.singleClass,
                               qml=qml, convergence=convergence, ...)
            P            <- res$P
            model$info$w <- res$w.c
            par.old      <- res$par.old
            if (verbose == TRUE) {
                cat("Class weights: ", round(model$info$w, digits=4), "\n")
            }
        }
    )

    if(verbose == TRUE){
        cat("Doing maximization-step \n")
    }
    
    # M-step
    switch(class(model),
      "singleClass" = {
        m.step <- mstep_lms(model=model, P=P, dat=data,
                            parameters=par.old, m=m,
                            optimizer=optimizer,
                            max.mstep=max.mstep, ...)
      },
      "semm" = {
        m.step <- mstep_semm(model=model, parameters=par.old, P=P,
                             data=data, optimizer=optimizer,
                             max.mstep=max.mstep, ...) },
      "nsemm" = {
        m.step <- mstep_nsemm(model=model, parameters=par.old, P=P,
                              data=data, optimizer=optimizer,
                              max.mstep=max.mstep, ...) }
    )

    if(verbose == TRUE) {
      cat("Results of maximization \n")
      cat(paste0("Loglikelihood: ", round(-m.step$objective, 3), "\n"))
      cat(paste0("Convergence message: ", m.step$convergence[1], "\n"))
      cat(paste0("Number of iterations: ", m.step$iterations, "\n"))
      cat("----------------------------------- \n")
    }
  
    ll.new   <- m.step$objective
    ll.ret   <- c(ll.ret, ll.new)
    par.new  <- unlist(m.step$par)
    num.iter <- num.iter + 1

    if(num.iter == max.iter){
      warning("Maximum number of iterations was reached. EM algorithm might not have converged.")
      break
    }
    if (abs(ll.old - ll.new) < convergence) run <- FALSE
  }

  
  if(verbose == TRUE) {
    cat("-----------------------------------\n")
    cat("EM completed \n")
    #cat(paste0("Previous loglikelihood: ", round(-ll.old, 3), "\n"))
    #cat(paste0("Final loglikelihood: ", round(-ll.new, 3),"\n"))
    cat("-----------------------------------\n")

    cat("-----------------------------------\n")
    if (neg.hessian == TRUE) {
        cat("Computing negative Hessian \n")
    } else {
        cat("Computing final model \n")
    }
    cat("-----------------------------------\n")
  }

  switch(class(model),
    "singleClass" = {
      final <- mstep_lms(model=model, P=P, dat=data,
                         parameters=par.new, neg.hessian=neg.hessian,
                         m=m, optimizer=optimizer,
                         max.mstep=max.mstep, ...)
      coefficients <- final$par
      names(coefficients) <- model$info$par.names
      # Transform parameters back to Phi
      A <- matrix(0, nrow=model$info$num.xi, ncol=model$info$num.xi)
      A[lower.tri(A, diag=TRUE)] <- coefficients[grep("Phi",
        names(coefficients))]
      Phi <- A %*% t(A)
      coefficients[grep("Phi", names(coefficients))] <- Phi[lower.tri(Phi, diag=TRUE)]
    },
    "semm" = {
      final <- mstep_semm(model=model, parameters=par.old, P=P,
                           data=data, neg.hessian=neg.hessian,
                           optimizer=optimizer, max.mstep=max.mstep, ...)
      if (is.numeric(final$par)) {
        coefficients <- get_class_parameters(model, final$par)
        for (class in names(model$matrices)) {
          names(coefficients[[class]]) <-
            get_class_parnames(model)[[class]]
        }
      } else {
        coefficients <- final$par
        for (class in names(model$matrices)) {
          names(coefficients[[class]]) <- model$info$par.names[[class]]
        }
      }
    },
    "nsemm" = {
      final <- mstep_nsemm(model=model, parameters=par.old, P=P,
                           data=data, neg.hessian=neg.hessian,
                           optimizer=optimizer, max.mstep=max.mstep, ...)
      if (is.numeric(final$par)) {
        coefficients <- get_class_parameters(model, final$par)
        for (class in names(model$matrices)) {
          names(coefficients[[class]]) <-
            get_class_parnames(model)[[class]]
        }
      } else {
        coefficients <- final$par
        for (class in names(model$matrices)) {
          names(coefficients[[class]]) <- model$info$par.names[[class]]
        }
      }
    }
  )

  # convergence of em
  if (num.iter == max.iter) {
    em_convergence <- "no"
  } else {em_convergence <- "yes"}

  info <- model$info[c("num.xi", "num.eta", "num.x", "num.y",
    "constraints", "num.classes")]
  info$n <- nrow(data)

  out <- list(model.class=class(model), coefficients=coefficients,
              objective=-final$objective, em.convergence=em_convergence,
              neg.hessian=final$hessian, loglikelihoods=-ll.ret, info=info)

  # attach w for semm and nsemm
  if (class(model) == "semm" || class(model) == "nsemm") out$info$w <-
    model$info$w

  class(out) <- "emEst"
  out
}


