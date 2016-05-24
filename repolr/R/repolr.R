repolr <-
function(formula, subjects, data, times, categories, corr.mod = "independence",
           alpha = 0.5, po.test = FALSE, fixed = FALSE, poly = NULL, space = NULL,
           diffmeth = "analytic", fit.opt = rep(NA, 5)){


 # start
 call <- match.call()

 # set-up
 corr.mods <- c("ar1", "uniform", "independence")
 icorr.mod <- as.integer(match(corr.mod, corr.mods, -1))
 if (icorr.mod < 1){stop("corr.mod: set to be independence, ar1 or uniform")}
 rcorr.mod <- corr.mod
 diffmeths <- c("analytic", "numeric")
 idiffmeth <- as.integer(match(diffmeth, diffmeths, -1))
 if (idiffmeth < 1){stop("diffmeth: set to be numeric or analytic")}
 if(times[1] != 1){stop("times: times should be vector with first value set to 1")}
 alpha <- as.double(alpha)
 po.test <- as.logical(po.test)
 fixed <- as.logical(fixed)
 set.fit.opt <- c(cmaxit = 10, omaxit = 5, ctol = 0.001, otol = 0.00001, h = 0.01)
 set.fit.opt[which(is.na(fit.opt) == FALSE)] <- fit.opt[which(is.na(fit.opt) == FALSE)]
 if(set.fit.opt[1] < 4){set.fit.opt[1] <- 4}
 if(alpha < 0.05 | alpha > 0.95){stop("alpha: invalid correlation parameter")}
 categories <- as.integer(categories)
 categories1 <- categories - 1
 subjects <- as.character(subjects)
 isubject <- as.integer(match(subjects, names(as.data.frame(data)), -1))
 if (isubject < 1){stop("subjects: unknown subject name")}
 orig.formula <- as.formula(formula)
 if(corr.mod == "independence" | length(times) == 1){alpha <- 0;   fixed <- TRUE;   corr.mod <- "uniform"}
 if(length(times) == 1){diffmeth <- "analytic"}
 if(is.null(poly) == FALSE){po.test <- FALSE; poly <- as.integer(poly)}
 if(sum(diff(times) > 0) != (length(times) - 1)){stop("times: invalid vector of times")}
 if(is.null(space) != TRUE){if(sum(diff(space) > 0) != categories1){stop("space: invalid vector of spacings")}}
 if(is.null(space) != TRUE){if(space[1] != 1){stop("space: space should be vector with first value set to 1")}}

 # model matrix
 exdata <- ord.expand(space = space, formula = formula, times = times, poly = poly,
                                   data = data, subjects = subjects, categories = categories)
 formula <- exdata$formula
 Xmat <- list(design = Matrix::Matrix(model.matrix(formula, data = exdata$data), sparse = TRUE))

 # initialize
 var.names <- Xmat$design@Dimnames[[2]]
 mod.glm <- glm(formula, family = binomial(), data = exdata$data)
 if(is.null(poly) == FALSE){
  poly1 <- poly + 1
  var.names[1:poly1] <- c("Intercept", paste("poly(cuts, ",categories - 1, ")", 1:poly, sep = ""))
  polydesign <- as(Xmat$design[1:(categories1), 1:(poly1)], "CsparseMatrix")
  polycuts <- as.numeric(polydesign %*% mod.glm$coefficients[1:(poly1)])
  xsmat <- smat(coeff=polycuts)
 } else {
  xsmat <- smat(coeff = mod.glm$coefficients[1:categories1])
 }
 xcmat <- cmat(ctimes = times, alpha = alpha, corrmod = corr.mod, 
                               diffmeth = diffmeth, h = set.fit.opt[5])
 xhgmat <- hgmat(mod = mod.glm, cmat = xcmat, smat = xsmat, X = Xmat, 
                                         modtype = "glm", diffmeth = diffmeth)
 xicormat <- list(irmat = xhgmat$icormat)
 mod.ordgee <- ordgee(mod = mod.glm, icormat = xicormat, X = Xmat, modtype = "glm",
                    ctimes = times, categories = categories,
                    omaxit = as.integer(set.fit.opt[2]), otol = as.double(set.fit.opt[4]))
 coeffs <- mod.ordgee$coefficients

 # fit model
 ogstop <- 1; iter <- 0; convergence <- FALSE
 while (ogstop > as.double(set.fit.opt[3]) & iter <= as.integer(set.fit.opt[1])){

  iter <- iter+1
  if(is.null(poly) == FALSE){
   polycuts <- as.numeric(polydesign %*% coeffs[1:(poly1)])
   xsmat <- smat(coeff = polycuts)
  } else {
   polycuts <- NA
   xsmat <- smat(coeff = coeffs[1:categories1])
  }
  xcmat <- cmat(ctimes = times, alpha = alpha, corrmod = corr.mod, 
                                  diffmeth = diffmeth, h = set.fit.opt[5])
  if(fixed == FALSE){
   xhgmat <- hgmat(mod = mod.ordgee, cmat = xcmat, smat = xsmat, 
                                  X = Xmat, modtype = "gee", diffmeth = diffmeth)
   xupalpha <- upalpha(hgmat = xhgmat, alpha = alpha, diffmeth = diffmeth, h = set.fit.opt[5])
   xicormat <- list(irmat = xhgmat$icormat)
  } else {
   xicormat <- icormat(mod = mod.ordgee, smat = xsmat, cmat = xcmat, modtype = "gee")
   xupalpha <- list(gvb = NA, ggvb = NA, alpha = alpha)
  }
  mod.ordgee <- ordgee(mod = mod.ordgee, icormat = xicormat, X = Xmat, modtype = "gee",
                 ctimes = times, categories = categories,
                 omaxit = as.integer(set.fit.opt[2]), otol = as.double(set.fit.opt[4]))
  crit <- abs(1 - sqrt(sum((coeffs / mod.ordgee$coefficients)^2) / length(coeffs)))
  if (iter <= 4){ogstop <- 1} else {ogstop <- crit}
  coeffs <- mod.ordgee$coefficients
  if(fixed == FALSE){alpha <- xupalpha$alpha}

 }
 if(ogstop <= as.double(set.fit.opt[3])){convergence <- TRUE}
 warn.message <- paste("Model did not converge: iter ", as.character(iter), 
                     " and crit ", as.character(crit), sep="")
 if(convergence == FALSE){warning(warn.message)}

 # variance matrices
 if(fixed != FALSE){
  xhgmat <- hgmat(mod = mod.ordgee, cmat = xcmat, smat = xsmat, 
                              X = Xmat, modtype = "gee", diffmeth = diffmeth)
 }
 robust.var <- Matrix::solve(xhgmat$hmat) %*% xhgmat$gmat %*% Matrix::solve(xhgmat$hmat)
 naive.var <- Matrix::solve(xhgmat$hmat)
 robust.var@Dimnames[[1]] <- robust.var@Dimnames[[2]] <- var.names
 naive.var@Dimnames[[1]] <- naive.var@Dimnames[[2]] <- var.names
 coeffs <- as.numeric(coeffs); names(coeffs) <- var.names
 if(is.null(poly) == FALSE){
  polycuts.robust <- sqrt(as.numeric(Matrix::diag(polydesign %*% robust.var[1:(poly1), 1:(poly1)] 
                                                  %*% Matrix::t(polydesign))))
  polycuts.naive <- sqrt(as.numeric(Matrix::diag(polydesign %*% naive.var[1:(poly1), 1:(poly1)] 
                                                 %*% Matrix::t(polydesign))))
 } else {
  polycuts.robust <- polycuts.naive <- NA
 }

 # po test
 if(po.test == TRUE){
  form.vars <- attr(terms.formula(as.formula(formula)), "variables")
  resp.var <- attr(terms.formula(as.formula(formula)), "response")
  term.labels <- attr(terms.formula(as.formula(formula)), "term.labels")
  resp.label <- as.character(form.vars[[resp.var + 1]])
  exformula <- as.formula(paste(resp.label, "~",
            "cuts - 1 + cuts / (", paste(term.labels[2:length(term.labels)], collapse="+"),
                  paste(")", sep = ""), sep = ""))
  exX_mat <- model.matrix(exformula, data = exdata$data)
  exXmat <- list(design = as(exX_mat, "CsparseMatrix"))
  xpotest <- potest(mod = mod.ordgee, hgmat = xhgmat, X = exXmat, 
                                      categories = categories, ctimes = times)
  testchi <- pchisq(xpotest$teststat, xpotest$testdf, lower.tail = FALSE)
  potest.out <- list(po.stat = xpotest$teststat, po.df = xpotest$testdf, po.chi = testchi)
 } else {
  potest.out <- list(po.stat = NA, po.df = NA, po.chi = NA)
 }

 # output
 ordgee.mod <- list(title = "repolr: 2016-02-26 version 3.4",
                          call = call,
                          data = call[["data"]],
                          subjects = subjects,
                          formula = formula, 
                          orig.formula = orig.formula,
                          corr.mod = rcorr.mod,
                          times = times,
                          categories = categories,
                          poly.mod = list(poly = poly, polycuts = list(coeff = polycuts, 
                               robust.se = polycuts.robust, naive.se = polycuts.naive), space = space),
                          max.id = as.numeric(mod.ordgee$max.id),
                          id = as.numeric(mod.ordgee$id),
                          y = as.numeric(mod.ordgee$y),
                          linear.predictors = as.numeric(mod.ordgee$linear.predictors),
                          fitted.values = as.numeric(mod.ordgee$fitted.values),
                          coefficients = coeffs,
                          robust.var = as.matrix(robust.var),
                          naive.var = as.matrix(naive.var),
                          fixed = fixed,
                          alpha = as.numeric(alpha),
                          convergence = convergence,
                          iter = iter,
                          fit.opt = set.fit.opt,
                          diffmeth = diffmeth,
                          grad1 = as.numeric(xupalpha$gvb),
                          grad2 = as.numeric(xupalpha$ggvb),
                          crit = crit,
                          po.test = potest.out)
  
 # set-up class for summary function
 class(ordgee.mod) <- "repolr"
 
 # end
 return(ordgee.mod)

}
