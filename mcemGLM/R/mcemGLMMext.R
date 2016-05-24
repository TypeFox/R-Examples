mcemGLMMext <- function(object, it = 20, controlEM) {
  if (class(object) != "mcemGLMM") {
    stop("Wrong object type.")
  }
  
  if (missing(controlEM)) {
    ctrl <- list(     EMit = it,
                      MCit = nrow(object$randeff),
                       MCf = ifelse(      is.null(object$call$controlEM$MCf),  1.25,       object$call$controlEM$MCf),
                      verb = ifelse(     is.null(object$call$controlEM$verb), FALSE,      object$call$controlEM$verb),
                      MCsd = object$MCsd, 
                   EMdelta = ifelse(  is.null(object$call$controlEM$EMdelta),  0.025,   object$call$controlEM$EMdelta),
                 EMepsilon = ifelse(is.null(object$call$controlEM$EMepsilon), 0.001, object$call$controlEM$EMepsilon))
    
    fit0 <- mcemGLMM(       fixed = eval(object$call$fixed),
                           random = eval(object$call$random), 
                             data = eval(object$call$data), 
                           family = eval(object$call$family), 
                           vcDist = eval(object$call$vcDist), 
                               df = eval(object$call$df), 
                        controlEM = ctrl, 
                     controlTrust = eval(object$call$controlTrust), 
                          initial = tail(eval(object$mcemEST), 1))
  } else {
    fit0 <- mcemGLMM(       fixed = eval(object$call$fixed),
                           random = eval(object$call$random),
                             data = eval(object$call$data),
                           family = eval(object$call$family),
                           vcDist = eval(object$call$vcDist),
                               df = eval(object$call$df),
                        controlEM = controlEM,
                     controlTrust = eval(object$call$controlTrust),
                          initial = tail(eval(object$mcemEST), 1))
  }
  fit0$call <- object$call
  return(fit0)
}