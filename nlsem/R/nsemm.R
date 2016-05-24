# nsemm.R
#
# created: Nov/14/2014, KN
# last mod: Aug/27/2015, NU

#--------------- main functions ---------------

estep_nsemm <- function(model, parameters, data, max.singleClass, qml,
                        convergence, verbose=FALSE, ...) {

  num.classes <- model$info$num.classes

  class.parameters <- get_class_parameters(model, parameters)
  par.new <- NULL

  # lms or qml for each class
  # Note that B is not estimated
  for (c in seq_len(num.classes)) {
    lms.model <- lms_ify(model, c)

    if (qml == FALSE) {
      # em for lms
      suppressWarnings(
      est <- em(model=lms.model, data=data, start=class.parameters[[c]],
                verbose=verbose, neg.hessian=FALSE,
                max.iter=max.singleClass, convergence=convergence, ...)
      )
      # suppress warnings since they are non-informative in this
      # intermediate step
  
      pars <- est$coefficients[model$info$par.names[[c]]]

      par.new <- c(par.new, pars)
    } else {
      est <- mstep_qml(model=lms.model, data=data, parameters=class.parameters[[c]],
                       neg.hessian=FALSE, max.iter=max.singleClass, ...)

      pars <- est$par[model$info$par.names[[c]]]
      par.new <- c(par.new, pars)
    }
  }

  # e-step for semm
  # Note that Omega and A are not estimated
  P <- estep_semm(model=model, parameters=par.new, data=data)
  w.c <- colSums(P) / nrow(data)

  res <- list(P=P, w.c=w.c, par.old=par.new)
  res
}

mstep_nsemm <- function(model, parameters, P, data, optimizer, max.mstep,
                      control=list(), ...) {

  est <- mstep_semm(model=model, parameters=parameters, P=P,
                              data=data, optimizer=optimizer,
                              max.mstep=max.mstep, control=control, ...)

  est
}

#--------------- helper functions ---------------

# create singleClass model for a specific class of an nsemm model
lms_ify <- function(model, c) {
    lms.model <- list(matrices=list(class1=model$matrices[[c]]),
                      info=model$info)
    if (model$info$constraints == "direct1") {
      lms.model$info$par.names <- model$info$par.names[[c]]
    } else {
      lms.model$info$par.names <- model$info$par.names$class1
    }
    lms.model$info$bounds$upper <- model$info$bounds$upper[[c]]
    lms.model$info$bounds$lower <- model$info$bounds$lower[[c]]
    lms.model$info$num.classes <- 1
    class(lms.model) = "singleClass"
    lms.model
}
