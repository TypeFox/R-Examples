compareModels <-
function (groups, estprobs = list(lda = NULL, rf = NULL),
            gpnames = NULL, robust = TRUE, print = TRUE)
{
  tab <- table(groups)
  checknam <- sapply(estprobs, function(x) all(names(x) ==
                                               names(tab)))
  if (!all(checknam))
    stop(c(paste("Levels of 'groups' are:", names(tab)),
           paste("List elements", paste(names(tab)[!checknam],
                                        collapse = " & "), "do not match these levels")))
  models <- factor(names(estprobs), levels = names(estprobs))
  if (is.null(models))
    stop("Elements of the list 'estprobs' must be named")
  g <- length(levels(groups))
  n <- length(groups)
  m <- length(estprobs)
  selmat <- cbind(1:n, unclass(groups))
  probs <- as.vector(sapply(estprobs, function(x) x[selmat]))
  df <- data.frame(gp = rep(groups, m), prob = probs, model = rep(models,
                                                        rep(n, m)), obs = factor(rep(1:n, m)))
  if (robust)
    mod <- rlm(prob ~ model + obs, data = df)
  else mod <- lm(prob ~ model + obs, data = df)
  pred <- predict(mod, type = "terms", terms = c("model", "obs"))
  bmod <- pred[match(models, df$model), "model"] + attr(pred,
                                                        "constant")
  gpmod <- lm(pred[, "obs"] ~ -1 + gp, data = df)
  gptab <- summary(gpmod)$coef
  bgp <- gptab[, 1] + attr(pred, "constant")
  if(!is.null(gpnames))names(bgp) <- gpnames
  avsegp <- sqrt(mean(gptab[, 2]^2))
  names(bmod) <- levels(models)
  coeff <- summary(mod)$coef[, 1:2]
  cnam <- rownames(coeff)
  modlab <- paste("model", levels(models), sep = "")
  modrows <- match(modlab[-1], cnam)
  semod <- coeff[modrows, 2]
  avsemod <- sqrt(mean(semod^2))
  resrows <- resid(gpmod)
  if (print) {
    print("Average accuracies for groups:")
    print(round(bgp, 4))
    print(c(`Approx sed` = round(avsegp, 4)))
    print("Average accuracies for methods:")
    print(round(bmod, 4))
    print(c(`Approx sed` = round(avsemod, 4)))
  }
  invisible(list(modelAVS = bmod, modelSE = avsemod, gpAVS = bgp,
                 gpSE = avsegp, obsEff = resrows))
}

