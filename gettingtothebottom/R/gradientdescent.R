#' Gradient Descent Algorithm
#' 
#' \code{gdescent} Performs gradient descent algorithm given 
#' an objective function and a gradient for the objective function
#' 
#' @param f objective function as a function of X, y, and b
#' @param grad_f gradient of f as a function of X,y, and b
#' @param X matrix of independent variables
#' @param y vector containing dependent variable
#' @param alpha (optional) step size for the algorithm
#' @param iter (optional) the number of iterations to include in the algorithm
#' @param liveupdates (optional) if TRUE, the function will print live updates showing the norm of the gradient vector in each iteration
#' @param tol (optional) tolerance for determining convergence
#' @param intercept (optional) if TRUE, the model includes an estimate for the intercept
#' @param autoscaling (optional) if TRUE, the function will automatically rescale the columns of X (divides each element in X by the maximal element in that column)
#' 
#' @import grid
#' 
#' @export
#' 
#' @examples
#' #--------------------------------------------------------
#' # EXAMPLE 1 - A Simple Example
#' #--------------------------------------------------------
#' 
#' # Generate some data for a simple bivariate example
#' set.seed(12345)
#' x <- sample(seq(from = -1, to = 1, by = 0.1), size = 50, replace = TRUE)
#' y <- 2*x + rnorm(50)
#' plot(x,y)
#' 
#' # Setting up for gradient descent
#' X <- as.matrix(x)
#' y <- as.vector(y)
#' f <- function(X,y,b) {
#'    (1/2)*norm(y-X%*%b,"F")^{2}
#' }
#' grad_f <- function(X,y,b) {
#'    t(X)%*%(X%*%b - y)
#' }
#' 
#' # Run a simple gradient descent example
#' simple_ex <- gdescent(f,grad_f,X,y,0.01)
#' 
#' # We can compare our gradient descent results with what we get if we use the lm function
#' lm(y~X)
#' 
#' # Notice that the algorithm may diverge if the step size (alpha) is not small enough
#' # THE FOLLOWING NOT RUN
#' # simple_ex2 <- gdescent(f,grad_f,X,y,alpha=0.05,liveupdates=TRUE)
#' # The live updates show the norm of the gradient in each iteration.  
#' # We notice that the norm of the gradient diverges when alpha is not small enough.
#' 
#' #--------------------------------------------------------
#' # EXAMPLE 2 - Linear Regression & Feature Scaling
#' #--------------------------------------------------------
#' 
#' f <- function(X,y,b) {
#'   (1/2)*norm(y-X%*%b,"F")^{2}
#' }
#' grad_f <- function(X,y,b) {
#'   t(X)%*%(X%*%b - y)
#' }
#' 
#' data(moviebudgets)
#' X <- as.matrix(moviebudgets$budget)
#' y <- as.vector(moviebudgets$rating)
#' # THE FOLLOWING NOT RUN
#' # movies1 <- gdescent(f,grad_f,X,y,1e-4,5000)
#' 
#' # We can compare our gradient descent results with what we get if we use the lm function
#' # THE FOLLOWING NOT RUN
#' # lm(y~X)
#' 
#' # Compare the above result with what we get without feature scaling 
#' # Not run:
#' # movies2 <- gdescent(f,grad_f,X,y,alpha=1e-19,iter=10000,liveupdates=TRUE,autoscaling=FALSE)
#' ## Note that running the gradient descent algorithm on unscaled column vectors 
#' ## requires a much smaller step size and many more iterations.
#' 
#' #--------------------------------------------------------
#' # EXAMPLE 3 - Multivariate Linear Regression
#' #--------------------------------------------------------
#' 
#' f <- function(X,y,b) {
#'   (1/2)*norm(y-X%*%b,"F")^{2}
#' }
#' grad_f <- function(X,y,b) {
#'   t(X)%*%(X%*%b - y)
#' }
#' 
#' data(baltimoreyouth)
#' B <- baltimoreyouth
#' X <- matrix(c(B$farms11,B$susp11,B$sclemp11,B$abshs11), nrow=nrow(B), byrow=FALSE)
#' y <- as.vector(B$compl11)
#' # THE FOLLOWING NOT RUN
#' # meals_graduations <- gdescent(f,grad_f,X,y,0.01,12000)
#' 
#' # We can compare our gradient descent results with what we get if we use the lm function
#' # THE FOLLOWING NOT RUN
#' # lm(y~X)
#' 
#' #--------------------------------------------------------
#' # EXAMPLE 4 - Logistic Regression
#' #--------------------------------------------------------
#' 
#' set.seed(12345)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n*p),n,p)
#' b <- matrix(rnorm(p),p,1)
#' e <- 0.5*matrix(rnorm(n),n,1)
#' z <- X%*%b + e
#' y <- as.vector((plogis(z) <= runif(n)) + 0)
#' 
#' l <- function(X,y,b) {
#'   -t(y)%*%(X%*%b) + sum(log(1+exp(X%*%b)))
#' }
#' grad_l <- function(X,y,b) {
#'   -t(X)%*%(y-plogis(X%*%b))
#' }
#' alpha = 1/(0.25*svd(cbind(1,X))$d[1]**2)
#' 
#' # Use gradient descent algorithm to solve logistic regression problem 
#' # THE FOLLOWING NOT RUN
#' # logistic_ex <- gdescent(l,grad_l,X,y,alpha=alpha,iter=15000)
#' 
#' # Use glm function to solve logistic regression problem
#' # THE FOLLOWING NOT RUN
#' # glm(y~X, family=binomial)
#' 
#' @author Jocelyn T. Chi
#' 
gdescent <- function(f,grad_f,X,y,alpha=1e-6,iter=3000,liveupdates=FALSE,tol=1e-6,intercept=TRUE,autoscaling=TRUE) {
  if (alpha <= 0) {
    return(cat('Please select a positive value for alpha.'))
  }
  if (intercept == TRUE) {
    cols <- ncol(X)
    len <- nrow(X)
    X <- matrix(c(rep(1,len),X),ncol=cols+1,byrow=FALSE)
  }
  if (autoscaling == TRUE) {
    maxcol_values <- numeric(length = ncol(X))
    for (i in 1:ncol(X)) {
      maxcol_values[i] <- max(X[,i])
      X[,i] <- X[,i]/max(X[,i])
    } 
  }
  b = as.vector(rep(0, ncol(X)))
  values.f <- numeric(length = iter)
  values.b <- matrix(NA,ncol(X),iter)
  values.bnorm <- matrix(NA,1,iter)
  values.g <- matrix(NA,ncol(X),iter)
  values.gnorm <- matrix(NA,1,iter)
  for (i in 1:iter) {
    values.f[i] <- f(X,y,b)
    values.g[,i] <- grad_f(X,y,b)
    values.gnorm[i] <- norm(as.matrix(values.g[,i]),"F")
    if (values.f[i] > 1e300 | values.gnorm[i] > 1e300) {
      values.f <- values.f[1:(i-1)]
      values.b <- values.b[,1:(i-1)]
      values.bnorm <- values.bnorm[1:(i-1)]
      values.g <- values.g[,1:(i-1)]
      values.gnorm <- values.gnorm[1:(i-1)]
      return(cat('The step size is too large.  Please choose a smaller step size.'))
    }
    if (liveupdates == TRUE) {
      cat(c(i,values.gnorm[i])) 
      cat("\n")
    }
    if (values.gnorm[i] < tol){
      digits <- max(3L, getOption("digits") - 3L)
      if (autoscaling == TRUE) {
        b <- b/maxcol_values 
      }
      cat(c('Minimum function value:\n', format(values.f[i],digits)))
      cat("\n\n")
      if (intercept==TRUE) {
        cat(c('Intercept:\n',format(b[1],digits))) 
        cat("\n\n")
      }
      if(intercept==TRUE) {
        cat(c('Coefficient(s):\n',format(b[2:length(b)],digits))) 
      } else {
        cat(c('Coefficient(s):\n',format(b,digits))) 
      }
      cat("\n\n")
      values.f <- values.f[1:i]
      values.b <- values.b[,1:i-1]
      values.bnorm <- values.bnorm[1:i-1]
      values.g <- values.g[,1:i]
      values.gnorm <- values.gnorm[1:i]
      break
    }
    b <- b - alpha*grad_f(X,y,b)
    if (autoscaling == TRUE) {
      values.b[,i] <- b/maxcol_values
    } else {
      values.b[,i] <- b
    }
    values.bnorm[i] <- norm(as.matrix(values.b[,i]),"F")
  }
  if (values.gnorm[length(values.gnorm)] > tol){
    cat('Minimum function value not attained. A better result might be obtained by decreasing the step size or increasing the number of iterations.')
    values.f <- values.f[1:i-1]
    values.b <- values.b[,1:i-1]
    values.bnorm <- values.bnorm[1:i-1]
    values.g <- values.g[,1:i-1]
    values.gnorm <- values.gnorm[1:i-1]
  }
  return(list(f=values.f,b=values.b,iterates=values.bnorm,gradient=values.g,gradient_norm=values.gnorm,iter=i))
}

#' Gradient Descent Algorithm - Plotting the Loss Function
#' 
#' \code{plot_loss} Plots the loss function of an object containing the results of a gradient descent object implementation
#' 
#' @param obj Object containing the results of a gradient descent implementation
#' 
#' @export
#' 
#' @examples
#' # Generate some data for a simple bivariate example
#' set.seed(12345)
#' x <- sample(seq(from = -1, to = 1, by = 0.1), size = 50, replace = TRUE)
#' y <- 2*x + rnorm(50)
#' 
#' # Components required for gradient descent
#' X <- as.matrix(x)
#' y <- as.vector(y)
#' f <- function(X,y,b) {
#'    (1/2)*norm(y-X%*%b,"F")^{2}
#' }
#' grad_f <- function(X,y,b) {
#'    t(X)%*%(X%*%b - y)
#' }
#' 
#' # Run a simple gradient descent example
#' simple_ex <- gdescent(f,grad_f,X,y,alpha=0.01)
#' 
#' # Plot the loss function
#' plot_loss(simple_ex)
#' 
#'   
#' @author Jocelyn T. Chi
#' 
plot_loss <- function(obj) {
  iter <- loss <- NULL
  df1 <- data.frame(iter = c(1:length(obj$f)), loss = obj$f)
  plot <- ggplot(df1, aes(x=iter, y=loss))
  plot + geom_point(colour="black") + theme_bw(base_size=14) + xlab("Iterations") + ylab("Value of the objective function") 
}


#' Gradient Descent Algorithm - Plotting the Iterates
#' 
#' \code{plot_iterates} Plots the iterates of an object containing the results of a gradient descent object implementation
#' 
#' @param obj Object containing the results of a gradient descent implementation
#' 
#' @export
#' @examples
#' # Generate some data for a simple bivariate example
#' set.seed(12345)
#' x <- sample(seq(from = -1, to = 1, by = 0.1), size = 50, replace = TRUE)
#' y <- 2*x + rnorm(50)
#' 
#' # Components required for gradient descent
#' X <- as.matrix(x)
#' y <- as.vector(y)
#' f <- function(X,y,b) {
#'    (1/2)*norm(y-X%*%b,"F")^{2}
#' }
#' grad_f <- function(X,y,b) {
#'    t(X)%*%(X%*%b - y)
#' }
#' 
#' # Run a simple gradient descent example
#' simple_ex <- gdescent(f,grad_f,X,y,alpha=0.01)
#' 
#' # Plot the iterates
#' plot_iterates(simple_ex)
#'  
#' @author Jocelyn T. Chi
#'
plot_iterates <- function(obj) {
  iter <- change_in_iterates <- NULL
  df1 <- data.frame(iter = c(1:length(obj$iterates)), change_in_iterates = apply(obj$b-obj$b[,ncol(obj$b)],2,FUN=function(x){norm(as.matrix(x),'f')}))
  plot <- ggplot(df1, aes(x=iter, y=change_in_iterates))
  plot + geom_point(colour="black") + theme_bw(base_size=14) + xlab("Iterations") + ylab("Distance of each iterate to the solution")
}

#' Gradient Descent Algorithm - Plotting the Gradient Function
#' 
#' \code{plot_gradient} Plots the norm of the gradient function of an object containing the results of a gradient descent object implementation
#' 
#' @param obj Object containing the results of a gradient descent implementation
#' 
#' @export
#' @examples
#' # Generate some data for a simple bivariate example
#' set.seed(12345)
#' x <- sample(seq(from = -1, to = 1, by = 0.1), size = 50, replace = TRUE)
#' y <- 2*x + rnorm(50)
#' 
#' # Components required for gradient descent
#' X <- as.matrix(x)
#' y <- as.vector(y)
#' f <- function(X,y,b) {
#'    (1/2)*norm(y-X%*%b,"F")^{2}
#' }
#' grad_f <- function(X,y,b) {
#'    t(X)%*%(X%*%b - y)
#' }
#' 
#' # Run a simple gradient descent example
#' simple_ex <- gdescent(f,grad_f,X,y,alpha=0.01)
#' 
#' # Plot the norm of the gradient function
#' plot_gradient(simple_ex)
#'  
#' @author Jocelyn T. Chi
#'
plot_gradient <- function(obj) {
  iter <- gradient <- NULL
  df1 <- data.frame(iter = c(1:length(obj$gradient_norm)), gradient = abs(obj$gradient_norm))
  plot <- ggplot(df1, aes(x=iter, y=gradient))
  plot + geom_point(colour="black") + theme_bw(base_size=14) + xlab("Iterations") + ylab("Norm of the norm of the gradient function") 
}

#' Gradient Descent Algorithm - Plots Depicting Gradient Descent Results in Example 1 Using Different Choices for the Step Size
#' 
#' \code{example.alpha} Plots the function values for gradient descent in Example 1 (without the intercept) given a particular value for alpha.
#' 
#' @param alpha Step-size for the gradient descent algorithm.
#' 
#' @export
#' @examples
#' example.alpha(0.01)
#' example.alpha(0.12)
#'  
#' @author Jocelyn T. Chi
#'
example.alpha <- function(alpha) {
  if (alpha > 0.13) {
    return(cat('Please select a smaller value for alpha.'))
  }
  # Run Example 1 with the given alpha
  set.seed(12345)
  x <- as.vector(sample(seq(from = -1, to = 1, by = 0.1), size = 50, replace = TRUE))
  y <- 2*x + rnorm(50)
  X <- as.matrix(x)
  y <- as.vector(y)
  f <- function(X,y,b) {
    (1/2)*norm(y-X%*%b,"F")^{2}
  }
  grad_f <- function(X,y,b) {
    return(t(X)%*%(X%*%b - y))
  }
  alpha_ex <- gdescent(f,grad_f,X,y,alpha,500,intercept=FALSE)
  # Plot the function values
  obj <- alpha_ex
  ols <- function(b) {
    (1/2)*norm(y-X%*%b,"F")^{2}
  }
  f_val <- numeric(length=length(alpha_ex$b))
  for (i in 1:length(alpha_ex$b)) {
    f_val[i] <- f(X,y,alpha_ex$b[i])
  }
  loss <- iterates <- NULL
  df1 <- data.frame(loss = f_val, iterates = alpha_ex$b)
  plot <- ggplot(df1, aes(x=iterates, y=loss))
  b <- as.matrix(seq(-abs(max(alpha_ex$b))+2,abs(max(alpha_ex$b))+2,length.out=1e3))
  df <- data.frame(b=b,f=apply(b,1,FUN=ols))
  plot + geom_point(colour="red") + geom_path(arrow = arrow(length = unit(0.5, "cm")), colour="red") + geom_line(data=df,aes(x=b,y=f)) + theme_bw(base_size=14) + xlab("b") + ylab("Value of the Loss Function") 
}


#' Gradient Descent Algorithm - Plots Depicting How Different Choices of Alpha Result in Differing Quadratic Approximations
#' 
#' \code{example.quadratic.approx} Shows how the quadratic approximations for the function in Example 1 change with choice of alpha.
#' 
#' @param alpha1 A smaller step-size for the gradient descent algorithm.
#' @param alpha2 A larger step-size for the gradient descent algorithm.  (alpha2 > alpha1)
#' 
#' @export
#' @examples
#' example.quadratic.approx(alpha1=0.01, alpha2=0.12)
#'  
#' @author Jocelyn T. Chi
#'
example.quadratic.approx <- function(alpha1=0.01, alpha2=0.12) {
  anchor <- 5
  if (alpha1 > 0.13 | alpha2 > 0.13) {
    return(cat('Please select a smaller value for alpha.'))
  }
  if (alpha2 <= alpha1) {
    return(cat('Please select a value for alpha2 that is larger than alpha1.'))
  }
  if (alpha2 < 0.06) {
    return(cat('Please select a larger value for alpha2.'))
  }
  set.seed(12345)
  x <- as.vector(sample(seq(from = -1, to = 1, by = 0.1), size = 50, replace = TRUE))
  y <- 2*x + rnorm(50)
  X <- as.matrix(x)
  y <- as.vector(y)
  # Plot the quadratic approximations for each
  ols <- function(b,X,y) {
    return((1/2)*norm(y-X%*%b,"F")^2)
  }
  grad_ols <- function(b,X,y) {
    return(t(X)%*%(X%*%b - y))
  }
  f <- NULL
  b <- as.matrix(seq(-5+2.12,5+2.12,length.out=1e3))
  df_orig <- data.frame(b=b,f=apply(b,1,FUN=function(b) {ols(b,X=X,y=y)}))
  quadratic_approx <- function(b,b0=anchor,alpha) {
    approx <- ols(b0,X,y) + (grad_ols(b0,X,y)*(b-b0)) + (1/(2*alpha))*(b-b0)^2
    return(approx)
  }
  df_init <- data.frame(df_orig[which(round(df_orig$b,2)==anchor),])
  df_alpha_small <- data.frame(b=b,f=apply(b,1,FUN=function(b) {quadratic_approx(b,alpha=alpha1)}))
  df_alpha_large <- data.frame(b=b,f=apply(b,1,FUN=function(b) {quadratic_approx(b,alpha=alpha2)}))
  df_points <- data.frame(rbind(df_orig[which.min(df_orig$f),],df_alpha_small[which.min(df_alpha_small$f),],df_alpha_large[which.min(df_alpha_large$f),]))
  df_minima <- data.frame(rbind(df_orig[which.min(df_alpha_large$f),],df_orig[which.min(df_alpha_small$f),]))
  plot <- ggplot(data=df_orig,aes(x=b,y=f)) + geom_line(colour="black") + geom_line(data=df_alpha_small,aes(x=b,y=f),colour='red') + geom_line(data=df_alpha_large,aes(x=b,y=f),colour='blue') + ylim(c(min(df_alpha_large$f)-10,max(df_orig$f))) + xlab("b") + ylab("Value of the Loss Function") + theme_bw(base_size=14)
  plot_vlines <- plot + geom_vline(xintercept=df_points$b,colour="gray",linetype = "longdash")
  plot_vlines_points <- plot_vlines + geom_point(data=df_points, aes(x=b, y=f,colour = factor(round(f,2)))) + scale_colour_manual(values=c('blue','black','red'),guide='none')
  plot_vlines_points_minima <- plot_vlines_points + geom_point(data=df_minima, aes(x=b, y=f), color="black") + geom_point(data=df_init, aes(x=b, y=f), color="green4")
  suppressWarnings(print(plot_vlines_points_minima))
}

