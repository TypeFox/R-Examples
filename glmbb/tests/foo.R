
 library(glmbb)
 data(crabs)

 gout <- glm(satell ~ color * spine * width * weight, family = poisson,
     data = crabs)
 sum(! is.na(gout$coefficients))
 gout <- glm(satell ~ 1, family = poisson, data = crabs)
 sum(! is.na(gout$coefficients))

 # default criterion AIC

 gout <- glmbb(satell ~ color * spine * width * weight, data = crabs)

 fits <- ls(envir = gout$envir, pattern = "^sha1")
 length(fits)
 criteria <- Map(function(x) get(x, envir = gout$envir)$criterion, fits)
 formulae <- Map(function(x) get(x, envir = gout$envir)$formula, fits)
 names(criteria) <- NULL
 names(formulae) <- NULL
 criteria <- unlist(criteria)
 formulae <- sapply(formulae, tidy.formula.hierarchical)
 fred <- data.frame(criteria, formulae, stringsAsFactors = FALSE)
 fred <- fred[order(criteria), ]
 fred <- fred[fred$criteria <= min(fred$criteria) + gout$cutoff, ]
 w <- fred$criteria
 w <- w - w[1]
 w <- exp(- w / 2)
 w <- w / sum(w)
 fred <- data.frame(criterion = fred$criteria, weight = w,
     formula = fred$formulae, stringsAsFactors = FALSE)
 opt <- options(width = 132)
 print(fred, right = FALSE, row.names = FALSE, print.gap = 2)
 options(opt)

 # check criteria
 criteria.too <- Map(function(x) get(x, envir = gout$envir)$aic, fits)
 names(criteria.too) <- NULL
 criteria.too <- unlist(criteria.too)
 identical(criteria, criteria.too)

 # check we do indeed have all less than cutoff

 gout.full <- glmbb(satell ~ color * spine * width * weight,
     data = crabs, cutoff = Inf)

 fits <- ls(envir = gout$envir, pattern = "^sha1")
 criteria <- Map(function(x) get(x, envir = gout$envir)$criterion, fits)
 criteria <- unlist(criteria)
 fits.full <- ls(envir = gout.full$envir, pattern = "^sha1")
 criteria.full <- Map(function(x)
     get(x, envir = gout.full$envir)$criterion, fits)
 criteria.full <- unlist(criteria.full)
 length(fits)
 length(fits.full)
 min(criteria) == min(criteria.full)
 inies <- which(criteria.full <= min(criteria.full) + gout$cutoff)
 idx <- match(names(criteria.full)[inies], names(criteria))
 all(! is.na(idx))
 all(criteria.full[inies] == criteria[idx])

 # now BIC

 gout <- glmbb(satell ~ color * spine * width * weight,
     family = poisson, data = crabs, criterion = "BIC")

 fits <- ls(envir = gout$envir, pattern = "^sha1")
 length(fits)
 criteria <- Map(function(x) get(x, envir = gout$envir)$criterion, fits)
 formulae <- Map(function(x) get(x, envir = gout$envir)$formula, fits)
 names(criteria) <- NULL
 names(formulae) <- NULL
 criteria <- unlist(criteria)
 formulae <- sapply(formulae, tidy.formula.hierarchical)
 fred <- data.frame(criteria, formulae, stringsAsFactors = FALSE)
 fred <- fred[order(criteria), ]
 fred <- fred[fred$criteria <= min(fred$criteria) + gout$cutoff, ]
 w <- fred$criteria
 w <- w - w[1]
 w <- exp(- w / 2)
 w <- w / sum(w)
 fred <- data.frame(criterion = fred$criteria, weight = w,
     formula = fred$formulae, stringsAsFactors = FALSE)
 print(fred, right = FALSE, row.names = FALSE, print.gap = 2)

 # check criteria
 criteria.too <- Map(function(x) BIC(get(x, envir = gout$envir)), fits)
 names(criteria.too) <- NULL
 criteria.too <- unlist(criteria.too)
 identical(criteria, criteria.too)

 # now AICc

 gout <- glmbb(satell ~ color * spine * width * weight,
     family = poisson, data = crabs, criterion = "AICc", cutoff = 5)
 fits <- ls(envir = gout$envir, pattern = "^sha1")
 criteria <- Map(function(x) get(x, envir = gout$envir)$criterion, fits)
 criteria.too <- Map(function(x) get(x, envir = gout$envir)$aic, fits)
 p.too <- Map(function(x)
     sum(! is.na(get(x, envir = gout$envir)$coefficients)), fits)
 n <- nrow(crabs)
 criteria.too <- Map(function(x, p) x + 2 * p * (p + 1) / (n - p - 1),
     criteria.too, p.too)
 all.equal(criteria, criteria.too)

 opt <- options(width = 132)
 summary(gout)
 summary(gout, cutoff = 2)
 summary(gout, cutoff = 8)

