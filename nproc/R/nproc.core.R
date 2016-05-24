nproc.core <- function(x = NULL, y, method = c("logistic", "penlog", "svm", "randomforest",
                                          "lda", "nb", "ada", "custom"), score = NULL, conf = FALSE,
                  alphalist = seq(from = 0, to = 1,by = 0.01), delta = 0.05,
                  split = TRUE, cv = FALSE, fold = 5, loc.prob = NULL, n.cores = 1) {
  method = match.arg(method)

    if (!cv) {
      fit = npc(x, y, method, score = score, alpha = alphalist, delta = delta,
                split = split, n.cores = n.cores)
      # if(is.null(loc.prob)) loc.prob = fit$loc.prob obj = npc.core(fit$y, fit$prob,
      # alpha = alphalist, delta = delta, loc.prob = loc.prob)
      pred = predict(fit, x, score)
      v = getroc(pred$pred.label, y)
    } else {
      n = length(y)
      n.fold = ceiling(n/fold)
      ind = sample(1:n)
     # rocmat = array(0, c(fold, length(alphalist), 2))

      rocmat = mcmapply(1:fold,FUN = function(i){
        cind = (i - 1) * n.fold + 1:n.fold
        cind = intersect(1:n, cind)
        te.ind = ind[cind]
        tr.ind = setdiff(1:n, te.ind)
        fit = npc(x[tr.ind, ], y[tr.ind], method, score = score[tr.ind], alpha = alphalist,
                  delta = delta, split = split, loc.prob = loc.prob, n.cores = 1)
        pred = predict(fit, x[te.ind, ], score[te.ind])
        getroc(pred$pred.label, y[te.ind])
      }, mc.cores = n.cores)
      rocmat = array(rocmat, c(length(alphalist),2,5))

#       for(i in 1:fold) {
#         cind = (i - 1) * n.fold + 1:n.fold
#         cind = intersect(1:n, cind)
#         te.ind = ind[cind]
#         tr.ind = setdiff(1:n, te.ind)
#         fit = npc(x[tr.ind, ], y[tr.ind], method, score = score[tr.ind], alpha = alphalist,
#                   delta = delta, split = split, loc.prob = loc.prob, n.cores = n.cores)
#         pred = predict(fit, x[te.ind, ], score[te.ind])
#         rocmat[i, , ] = getroc(pred$pred.label, y[te.ind])
#       }
      v = apply(rocmat, c(1, 2), mean)

    }

  if (is.null(loc.prob)&(!cv))
    loc.prob = fit$loc.prob

  return(list(v = v, loc.prob = loc.prob))
}


pred.npc.core<-function(object, newx = NULL){
  if (object$method == "logistic") {
    pred.score = predict(object$fit, data.frame(newx), type = "response")
  } else if (object$method == "penlog") {
    pred.score = predict(object$fit$glmnet.fit, newx = newx, type = "response",
                         s = object$fit$lambda.min)
  } else if (object$method == "svm") {
    pred.score = attr(predict(object$fit, newx, decision.values = TRUE), "decision.values")[,
                                                                                            1]
  } else if (object$method == "randomforest") {
    pred.score = predict(object$fit, newx, type = "prob")[, 2]
  } else if (object$method == "lda") {
    pred.score = predict(object$fit, data.frame(newx))$posterior[, 2]
  } else if (object$method == "nb") {
    pred.score = predict(object$fit, data.frame(newx), type = "raw")[, 2]
  } else if (object$method == "ada") {
    pred.score = predict(object$fit, data.frame(newx), type = "probs")[, 2]
  }
  pred.score = as.vector(pred.score)
if (object$sign == TRUE) {
  pred.label = outer(pred.score, object$cutoff, ">")
} else {
  pred.label = outer(pred.score, object$cutoff, "<=")
}
list(pred.label = pred.label, pred.score = pred.score)

}
