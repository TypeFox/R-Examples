#' Gaussian Class with args (mu, sigma) or (pi, tau)
#' @description object to represent normal distributions
#' normal with mean mu, and standard deviation sigma
#' precision is pi, 1 / sigma squared
#' tau is mu * pi, the precision adjusted mean
#' @param mu mean skill
#' @param sigma std dev of skill
#' @param pi precision ( 1 / sigma^2)
#' @param tau precision adjusted mean (mu / sigma^2) 
Gaussian <- setRefClass('Gaussian',
  fields = list(pi = "numeric", tau = "numeric"),
  methods = list(
    initialize = function(mu = NULL, sigma = NULL, pi = NULL, tau = NULL) {
      if(!is.null(pi) & !is.null(tau)) {
        .self$pi <- pi
        .self$tau <- tau
      }
      else if(!is.null(mu) & !is.null(sigma)) {
        .self$pi <- sigma ^ -2
        .self$tau <- .self$pi * mu
      }
      else {
        .self$pi <- 0
        .self$tau <- 0
      }
      .self 
    },
    MuSigma = function() {
      if(pi == 0) { 
        mu = 0
        sigma = Inf 
      }
      else { 
        mu = tau / pi
        sigma = sqrt(1 / pi) 
      }
      return(c(mu, sigma))
    },
    show = function() {
      print(sprintf("Guassian [(mu, sigma), (pi, tau)]: [(%s, %s), (%s, %s)]", 
        round(MuSigma()[1], 3), round(MuSigma()[2], 3), round(pi, 3), round(tau, 3)))
    }, 
    mu = function() return(MuSigma()[1]),
    sigma = function() return(MuSigma()[2])
  )                                                                 
)

#' multiply Gaussians
#' @param x a Gaussian
#' @param y a Gaussian
Multiply <- function(x, y) {
    z <- Gaussian$new()
    z$pi  <- x$pi + y$pi
    z$tau <- x$tau + y$tau
    z
}

#' divide Gaussians
#' @param x a Gaussian
#' @param y a Gaussian
Divide <- function(x, y) {
    z <- Gaussian$new()
    z$pi  <- x$pi - y$pi
    z$tau <- x$tau - y$tau
    z
}

#' multiply Gaussians
#' @param x a Gaussian
#' @param y a Gaussian
setMethod("*", signature(e1 = "Gaussian", e2 = "Gaussian"), function(e1, e2) Multiply(e1, e2))
#' divide Gaussians
#' @param x a Gaussian
#' @param y a Gaussian
setMethod("/", signature(e1 = "Gaussian", e2 = "Gaussian"), function(e1, e2) Divide(e1, e2))


GetName <- function(x) { return(x$name) }
GetNames <- function(list) { return(Map(GetName, list)) }

PrettyPrint <- function(x) { 
  string <- c("[\"")
  string <- cbind(string, list(x[1]))
  if (length(x) > 1) {
    for(i in 2:length(x)) {
      string <- cbind(string, c("\", \""))
      string <- cbind(string, list(x[i]))
    }
  }
  string <- cbind(string, c("\"]"))
  return(string)
}


# ... and callSuper is for potential classes that inhererit from Gaussian
# http://stat.ethz.ch/R-manual/R-devel/library/methods/html/refClass.html

# "messages from factors", called factors in Doug Zongker's version 
# it stores messages in Gaussian form
# and are referenced by factor name

#' Variable Class has inputs value (skill), name and messages, a list of messages which are referenced by factor name 
#' @param value value of Gaussian,
#' @param name name of Variable in factor graph, e.g. skill, performance, team and team difference
#' @param messages list of messages which are referenced by factor name
Variable <- setRefClass("Variable", 
  fields = list(value = "Gaussian", name = "character", messages = "list"),
  methods = list(
    initialize = function(value = Gaussian$new(), name = "", messages = list(), ...) {
      .self$name <- name
      .self$value <- value
      .self$messages <- messages
      callSuper(...)
      .self
    },
    # update message to be new message
    # save old message
    # update value to be (value / old message) * new message
    UpdateMessage = function(factor = Factor$new(), message = Gaussian$new()) {
      old_message <- messages[[factor$name]]
      value <<- value / old_message * message
      messages[factor$name] <<- message
      .self
    },  
    # .self$variable cannot have assignment <<- but must have <-
    UpdateValue = function(factor = Factor$new(), new_value = Gaussian$new()) {
      old_message <- messages[[factor$name]]
      old_value <- value
      messages[factor$name] <<- (new_value * old_message) / old_value
      .self$value <- new_value
      .self
    },  
    GetMessage = function(factor) {
      messages[[factor$name]]
    },
    show = function() {
      mu <- value$tau / value$pi
      sigma <- sqrt(1 / value$pi)
      print(sprintf("[Variable] %s with value [(mu, sigma), (pi, tau)]: [(%s, %s), (%s, %s)]",
        name, round(mu, 3), round(sigma, 3), round(value$pi, 3), round(value$tau, 3)))
      if (length(messages) == 0) {
         print("Messages from [Factors] have not been assigned")
      }
      else {
        # factors is not a list of factors, but a list of Gaussians with factor names as keys 
        factor_names <- names(messages)
        string <- PrettyPrint(factor_names)
        string <- cbind("[Factors]: ", string)
        do.call("cat", string)
      }
    }
  ) 
)

#' Base Factor class from which other Factor Classes inherit
#' @param variables list of variables that factor is connected to
#' @param name to keep track of factor objects and display purposes
#' @param ... Reference Class inheritance
Factor <- setRefClass("Factor",
  fields = list(variables = "list", name = "character"),
  methods = list(
    initialize = function(variables = list(), name = "", ...) {
      .self$name <- name
      .self$variables <- variables
      if(length(variables) > 0) {
        for (i in 1:length(variables)) {
           AttachFactor(variables[[i]])
        }                  
      }
      .self
    },
    AttachFactor = function(variable) {
      variable$messages[name] <- Gaussian$new()
    },
    show = function() {
      variable_names <- unlist(GetNames(variables))
      string <- cbind("[Factor]:", name, "with [Variable(s)]:",  PrettyPrint(variable_names))
      do.call("cat", string)
    }
  )
)


#' PriorFactor
#' @description
#' Factor classes that implement the 5 update equationsin Herbrich
#' PriorFactor - That is, the starting skill of a player, push the parameters to the variable.
#' @param ... Reference Class inheritance
#' @param Variable 
#' @param param the starting skill value (Gaussian object) to update the Variable
PriorFactor <- setRefClass("PriorFactor",
  contains = "Factor",
  fields = list(variable = "Variable", param = "Gaussian"),
  methods = list(
    initialize = function(variable = Variable$new(), param = Gaussian$new(), ...) {
      variables <- list(variable)
      callSuper(variables, ...)
      .self$param <- param
    },
    Start = function() {
      variables[[1]]$UpdateValue(.self, param)
    },
    show = function() {
      variable_names <- unlist(GetNames(variables))
      string <- cbind("[PriorFactor]:", name, "with [Variable]:",  PrettyPrint(variable_names))
      do.call("cat", string)
    }
  )
)

#' Likelihood Factor Class
#' @description Factor classes that implement the 5 update equationsin Herbrich                       
#' Connects two variables, the value of one being the mean of the message sent to the other.
#' @param ... Reference Class inheritance
# default variance is BETA^2
LikelihoodFactor <- setRefClass("LikelihoodFactor",
  contains = "Factor",
  fields = list(mean = "Variable", value = "Variable", variance = "numeric"),
  methods = list(
    initialize = function(mean_variable = Variable$new(), value_variable = Variable$new(), variance = BETA^2, ...) {
      callSuper(variables = list(mean_variable, value_variable), ...)
      .self$mean = mean_variable
      .self$value = value_variable
      .self$variance = variance
    },
    UpdateValue = function() {
    "Update the value after a change in the mean (going down) in the TrueSkill factor graph."
      y = (.self$mean)$value
      fy = (.self$mean)$GetMessage(.self)
      a = 1.0 / (1.0 + .self$variance * (y$pi - fy$pi))
      message <- Gaussian$new(pi = a * (y$pi - fy$pi), tau = a * (y$tau - fy$tau))
      # print(sprintf("message from [Likelihood Factor]: %s", .self$name))
      # print(message)
      value$UpdateMessage(.self, message)
    },
    UpdateMean = function() {
    "Update the mean after a change in the value (going up in the TrueSkill factor graph. "
      # Note this is the same as UpdateValue, with self.mean and self.value interchanged.
      x = (.self$value)$value
      fx = (.self$value)$GetMessage(.self)
      a = 1.0 / (1.0 + .self$variance * (x$pi - fx$pi))
      mean$UpdateMessage(.self, Gaussian$new(pi = a * (x$pi - fx$pi), tau = a * (x$tau - fx$tau)))
    }
  )
)

#' SumFactor
#' @description
#' Factor classes that implement the 5 update equationsin Herbrich
#' @details
#' A factor that connects a sum variable with 1 or more terms, which are summed after being multiplied by fixed (real) coefficients.
#' SumFactor$UpdateTerm()
#' Swap the coefficients around to make the term we want to update
#' be the 'sum' of the other terms and the factor's sum, eg.,
#' change:
#'
#'    x = y_1 + y_2 + y_3
#'
#' to
#'
#'    y_2 = x - y_1 - y_3
#'
#' then use the same update equation as for UpdateSum.
SumFactor <- setRefClass("SumFactor",
  contains = "Factor",
  fields = list(sum = "Variable", terms = "list", coeffs = "list"),
  methods = list(
    initialize = function(sum_variable = Variable$new(), term_variables = list(Variable$new()), coeffs = list(1), ...) {
      stopifnot(length(coeffs) == length(term_variables))
      .self$sum = sum_variable
      .self$terms = term_variables
      .self$coeffs = coeffs
      callSuper(variables = append(list(sum_variable), term_variables), ...)
    },
    # var is a Variable, y is a list of Values, fy is a list of messages, a is a list of numerics
    InternalUpdate = function(var = Variable$new(), y = list(), fy = list(), a = list()) {
      fn = function(a, y, fy) return(a^2 / (y$pi - fy$pi))
      fn2 = function(a, y, fy) return(a * (y$tau - fy$tau) / (y$pi - fy$pi))
      
      new_pi = 1.0 / sum(unlist(mapply(fn, a, y ,fy, SIMPLIFY = F)))
      new_tau <- new_pi * sum(unlist(mapply(fn2, a, y, fy, SIMPLIFY = F)))
      
      new_msg <- Gaussian(pi = new_pi, tau = new_tau)
      #print(new_msg)
      var$UpdateMessage(.self, new_msg)
    },
    UpdateSum = function() {
      "Update the sum value (moving down in the factor graph)."
      y <- Map(function(x) return(x$value), terms)  	    
      fy <- Map(function(x) return(x$GetMessage(.self)), terms)
      a <- coeffs
      InternalUpdate(sum, y, fy, a)
    },
    UpdateTerm = function(index) {
      " Update one of the term values (moving up in the factor graph). "
      b <- coeffs
      a <- list()
      for(i in 1:length(b)) {
      	a[[i]] <- -b[[i]] / b[[index]]
      }
      
      a[index] <- 1 / b[[index]]
      
      v <- .self$terms
      v[index] <- .self$sum
      
      y <- Map(function(x) return(x$value), v)
      fy <- Map(function(x) return(x$GetMessage(.self)), v)
      
      InternalUpdate(terms[[index]], y , fy, a)
    }                                                    
  )
)

#' @title Vwin
#' @description update rules for approximate marginals for the win and draw cases
#' required by truncate factor
#' @param e argument from TruncateFactor: e <- .self$epsilon * sqrt(div$pi)  
#' @param t argument from TruncateFactor: t <- div$tau / sqrt(div$pi)
Vwin <- function(t, e) {
  # dnorm is PDF
  # pnorm is CDF
  # qnorm is inverse CDF
  x <- t - e
  return(dnorm(x) / pnorm(x))
}

#' @title Wwin
#' @description update rules for approximate marginals for the win and draw cases
#' required by truncate factor
#' @param e argument from TruncateFactor: e <- .self$epsilon * sqrt(div$pi)  
#' @param t argument from TruncateFactor: t <- div$tau / sqrt(div$pi)
Wwin <- function(t, e) {
  return(Vwin(t, e) * (Vwin(t, e) + t - e))
}

#' @title Vdraw
#' @description update rules for approximate marginals for the win and draw cases
#' required by truncate factor
#' @param e argument from TruncateFactor: e <- .self$epsilon * sqrt(div$pi)  
#' @param t argument from TruncateFactor: t <- div$tau / sqrt(div$pi)
Vdraw <- function(t, e) {
  return((dnorm(-e - t) - dnorm(e - t)) / (pnorm(e - t) - pnorm(-e - t)))
}

#' @title Wdraw
#' @description update rules for approximate marginals for the win and draw cases
#' required by truncate factor
#' @param e argument from TruncateFactor: e <- .self$epsilon * sqrt(div$pi)  
#' @param t argument from TruncateFactor: t <- div$tau / sqrt(div$pi)
Wdraw <- function(t, e) {
  return((Vdraw(t, e) ^ 2) + ((e - t) * dnorm(e - t) + (e + t) * dnorm(e + t)) / (pnorm(e - t) - pnorm( -e - t)))		
}

#' TruncateFactor
#' @description Factor classes that implement the 5 update equationsin Herbrich
#' @details A factor for (approximately) truncating the team difference
#' distribution based on a win or a draw (the choice of which is
#' determined by the functions you pass as V and W). 
TruncateFactor <- setRefClass("TruncateFactor",
  contains = "Factor",
  fields = list(variable = "Variable", V = "function", W = "function", epsilon = "numeric"),
  methods = list(
    initialize = function(variable = Variable$new(), V = Vwin, W = Wwin, epsilon = 0.10, ...) {
      callSuper(variables = list(variable), ...)
      .self$variable = variable
      .self$V = V
      .self$W = W
      .self$epsilon = epsilon
    },
    Update = function() {
      x = (.self$variable)$value
      fx = (.self$variable)$GetMessage(.self)
                 
      div <- x / fx
      
      t <- div$tau / sqrt(div$pi)
      e <- .self$epsilon * sqrt(div$pi) 
      
      v <- V(t, e)
      w <- W(t, e)
 
      new_pi <- (div$pi / (1.0 - w))
      new_tau <- (div$tau + sqrt(div$pi) * v) / (1.0 - w)
      
      new_val <- Gaussian(pi = new_pi, tau = new_tau)
      variable$UpdateValue(.self, new_val)
    },
    show = function() {
      variable_names <- unlist(GetNames(variables))
      string <- cbind("[Truncate Factor]:", name, "with [Variable(s)]:", PrettyPrint(variable_names))
      do.call("cat", string)
    }
  )
)
