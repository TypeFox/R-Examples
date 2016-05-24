
 library(glmbb)
 data(crabs)
 options(width = 132)

 gout <- glmbb(satell ~ color * spine * width * weight,
     data = crabs, criterion = "BIC", graphical = TRUE)
 sout <- summary(gout)
 sout
 
 # is this correct?

 gout.full <- glmbb(satell ~ color * spine * width * weight,
     data = crabs, criterion = "BIC", cutoff = Inf)

 e <- gout.full$envir
 fits <- ls(envir = e, pattern = "^sha1")
 criteria <- Map(function(x) get(x, envir = e)$criterion, fits)
 formulae <- Map(function(x) get(x, envir = e)$formula, fits)
 names(criteria) <- NULL
 names(formulae) <- NULL
 criteria <- unlist(criteria)
 isgraphical <- sapply(formulae, isGraphical)
 formulae <- sapply(formulae, tidy.formula.hierarchical)

 criteria <- criteria[isgraphical]
 formulae <- formulae[isgraphical]
 inies <- criteria <= min(criteria) + gout$cutoff
 criteria <- criteria[inies]
 formulae <- formulae[inies]
 formulae <- formulae[order(criteria)]
 criteria <- criteria[order(criteria)]
 w <- criteria - min(criteria)
 w <- exp(- w / 2)
 w <- w / sum(w)
 all.equal(criteria, sout$results$criterion)
 all.equal(formulae, sout$results$formula)
 all.equal(w, sout$results$weight)

