setClass("FLXMRcondlogit",
         representation(strata="ANY"),
         contains = "FLXMRglm")

FLXMRcondlogit <- function(formula=.~., strata) {
  z <- new("FLXMRcondlogit", weighted=TRUE, formula=formula, strata=strata,
           family="multinomial", name=paste("FLXMRcondlogit"))

  
  z@defineComponent <- expression({
    predict <- function(x, ...) 
      tcrossprod(x, t(coef))
    
    logLik <- function(x, y, strata) {
      llh_all <- vector("numeric", length = length(y))
      eta <- predict(x)
      llh_all[as.logical(y)] <- eta[as.logical(y)]
      ((tapply(llh_all, strata, sum) - tapply(exp(eta), strata, function(z) log(sum(z))))/tabulate(strata))[strata]
    }
    
    new("FLXcomponent",
        parameters=list(coef=coef),
        logLik=logLik, predict=predict,
        df=df)
  })
  
  z@fit <- function(x, y, w, strata){
    index <- w > 0
    fit <- survival::coxph.fit(x[index,,drop=FALSE], survival::Surv(1-y, y)[index], strata[index], weights=w[index], control = survival::coxph.control(),
                               method = "exact", rownames = seq_len(nrow(y))[index])
    coef <- coef(fit)
    df <- length(coef)
    with(list(coef = coef, df = df), eval(z@defineComponent))
  }
  z
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRcondlogit"),
          function(model, data, formula, lhs=TRUE, ...)
{
  formula <- RemoveGrouping(formula)
  if(is.null(model@formula))
    model@formula = formula
  model@fullformula = update(terms(formula, data=data), model@formula)
  ## Ensure that an intercept is included
  model@fullformula <- update(model@fullformula, ~ . + 1)
  if (lhs) {
    mf <- model.frame(model@fullformula, data=data, na.action = NULL)
    model@x <- model.matrix(attr(mf, "terms"), data=mf)
    response <- as.matrix(model.response(mf))
    model@y <- model@preproc.y(response)
  }
  else {
    mt1 <- terms(model@fullformula, data=data)
    mf <- model.frame(delete.response(mt1), data=data, na.action = NULL)
    mt <- attr(mf, "terms")
    model@x <- model.matrix(mt, data=mf)
  }
  strata <- update(model@strata, ~ . + 0)
  mf <- model.frame(strata, data=data, na.action=NULL)
  model@strata <- as.integer(model.matrix(attr(mf, "terms"), data=mf))
  ## Omit the intercept for identifiability
  model@x <- model@x[,attr(model@x, "assign") != 0, drop=FALSE]
  model@x <- model@preproc.x(model@x)
  model
})

setMethod("FLXmstep", signature(model = "FLXMRcondlogit"), function(model, weights, ...) {
  apply(weights, 2, function(w) model@fit(model@x, model@y,  w, model@strata))
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRcondlogit"), function(model, components, ...) {
  sapply(components, function(x) x@logLik(model@x, model@y, model@strata))
})

setMethod("existGradient", signature(object = "FLXMRcondlogit"),
          function(object) FALSE)

