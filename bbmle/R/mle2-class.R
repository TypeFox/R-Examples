## must go before setAs to avoid warnings
setClass("mle2", representation(call = "language",
                                call.orig = "language",
                                coef = "numeric",
                                fullcoef = "numeric",
                                vcov = "matrix",
                                min = "numeric",
                                details = "list",
                                minuslogl = "function",
                                method = "character",
                                data="list",
                                formula="character",
                                optimizer="character"))

setAs("mle","mle2", function(from,to) {
  new("mle2",
      call=from@call,
      call.orig=from@call,
      coef=from@coef,
      fullcoef=from@fullcoef,
      vcov=from@vcov,
      min=from@min,
      details=from@details,
      minuslogl=from@minuslogl,
      method=from@method,
      data=list(),
      formula="",
      optimizer="optim")
})
                

setClass("summary.mle2", representation(call = "language",
                               coef = "matrix",
                               m2logL = "numeric"))

setClass("profile.mle2", representation(profile="list",
                                       summary="summary.mle2"))


setClass("slice.mle2", representation(profile="list",
                                       summary="summary.mle2"))

setIs("profile.mle2", "slice.mle2")

