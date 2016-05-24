FlimFormula <-
function(f) {
  if(is(f, "formula")) {
    responses  <- all.vars(f[[2]])
    predictors <- all.vars(f[[3]])
    covariates <- predictors[!(predictors %in% responses)]
    no.f <- length(responses)
    reg.fmlas <- list()
    pred.terms <- attr(terms(f), "term.labels")
    for(k in 1:no.f) {
      resp <- paste0(responses[k], ".inc")
      cvrt <- pred.terms
      reg.fmlas[[k]] <- reformulate(cvrt, resp, attr(terms(f), "intercept"))
    } 
    ret <- list(input = f,
                responses = responses,
                covariates = covariates,
                predictors = predictors,
                predictor.terms = pred.terms,
                reg.fmlas = reg.fmlas)
  }
  if(is.list(f)) {
    ret <- list(input = f,
                responses = c(),
                covariates = c(),
                predictors = c(),
                predictor.terms = list(),
                reg.fmlas = list())
    no.f <- length(f)
    for(k in 1:no.f) {
      fk <- f[[k]]
      ret$responses <- c(ret$responses, all.vars(fk[[2]]))
      ret$predictors <- c(ret$predictors, all.vars(fk[[3]]))
      pred.terms <- attr(terms(fk), "term.labels")
      ret$predictor.terms[[k]] <- pred.terms
      resp <- paste0(all.vars(fk[[2]]), ".inc")
      cvrt <- pred.terms
      ret$reg.fmlas[[k]] <- reformulate(cvrt, resp, attr(terms(fk), "intercept"))
    }
    ret$predictors <- unique(ret$predictors)
    ret$covariates <- ret$predictors[!(ret$predictors %in% ret$responses)]
  }
  ret
}
