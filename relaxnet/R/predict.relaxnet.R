
## predict function for relaxnet objects -- currently uses colnames
## of beta matrix in order to subset newx -- may need to change this

predict.relaxnet <- function(object,
                             newx,
                             which.model,
                             s = NULL,
                             type = c("link", "response", "coefficients",
                                     "nonzero", "class"),
                             exact = FALSE,
                             ##offset,
                             ...) {

  type = match.arg(type)

  if(!identical(exact, FALSE))
    stop("For the current version, the only supported value of\n",
         "the exact argument is FALSE")
  
  if(which.model == "main") {
    return(predict(object$main.glmnet.fit,
                   newx = newx,
                   s = s,
                   type = type,
                   exact = exact,
                   ## offset = offset,
                   ...))
  }

  ## otherwise we want predictions from one of the relaxed models
  

  ## intercept only model
  
  if(inherits(object$relax.glmnet.fits[[which.model]],
              "relaxnet.intercept.only")) {

    int.model <- object$relax.glmnet.fits[[which.model]]

    if(type == "coefficients") {

      if(length(s) != 1) stop("for type = coefficients, s must have length 1")

      return(predict(object$main.glmnet.fit,
                     s = object$main.glmnet.fit$lambda[1],
                     type = "coefficients",
                     ...))
    }

    return(switch(class(object$main.glmnet.fit)[1],
                  elnet =
                    if(type %in% c("link", "response")) {
                      return(rep(int.model,
                                 nrow(newx)))
                    } else {
                           stop("type = ", type,
                                "has not been implemented yet\n",
                                "for intercept only elnet models")
                    },

                  lognet =
                    switch(type,
                           link = rep(int.model,
                             nrow(newx)),
                           
                           response = 1 / (1 +
                             exp(-rep(int.model,
                                      nrow(newx)))),
                           ## check to make sure I'm returning the
                           ## class the right way (i.e. as 0's/1's)
                           class = rep(ifelse(int.model>=0, 1, 0),
                                       nrow(newx)),
                           stop("type = ", type,
                                "has not been implemented yet\n",
                                "for intercept only lognet models"))))

  }

  ## if type is coefs or nonzero, we need to make sure that the result is in
  ## agreement with the full predictor matrix from the original relaxnet call
  ## as opposed to just the predictors that were in the relaxed model
  
  if(type == "coefficients") {

    if(length(s) != 1) stop('for type = "coefficients", s must have length 1')

    ## just do this to get a single column matrix of the right size

    result <- predict(object$main.glmnet.fit, type = "coefficients",
                      s = object$main.glmnet.fit$lambda[1])

    ## now get the actual coefs

    ## subset the columns of newx to conform with the relaxed model

    actual.coefs <- predict(object$relax.glmnet.fits[[which.model]],
                    newx = newx[,
                      rownames(object$relax.glmnet.fits[[which.model]]$beta)],
                    s = s,
                    type = type,
                    exact = exact,
                    ##offset = offset,
                    ...)

    ## set all to zero, keeping dimnames
    
    result[] <- 0

    ## replace relevant ones with actual
    
    result[rownames(actual.coefs), ] <- actual.coefs

    return(result)
  }

  ## type == nonzero returns column indices for x corresponding to predictors
  ## who's coefficient was not set to zero. Need to convert the indices from
  ## indices of the x used for the relaxed model back to the full x

  ## note: seems like when s has length one, predict.glmnet(type = nonzero)
  ## returns a single-column data frame. Keeping this behavior for now
  
  if(type == "nonzero") {

    if(length(s) != 1) stop('for type = "nonzero", s must have length 1')

    relax.x.indices <- predict(object$relax.glmnet.fits[[which.model]],
                                     s = s,
                                     type = type,
                                     exact = exact)[[1]]


    nonzero.col.names <-
      rownames(object$relax.glmnet.fits[[which.model]]$beta)[relax.x.indices]

    full.x.indices <-
      data.frame(X1 = which(rownames(object$main.glmnet.fit$beta)
                       %in% nonzero.col.names))

    return(full.x.indices)
  }
  
  ## otherwise type is not coefs or nonzero and length will be OK
  ## (i.e. not different from main)
  
  predict(object$relax.glmnet.fits[[which.model]],
          newx = newx[,
            rownames(object$relax.glmnet.fits[[which.model]]$beta)],
          s = s,
          type = type,
          exact = exact,
          ...)  
}
