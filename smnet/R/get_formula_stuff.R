	get_formula_stuff<-function(formula, data){
	  gp         <- interpret.formula(formula.object = formula)
	  mf         <- match.call(expand.dots = FALSE)
	  m          <- match(c("formula", "data"), names(mf), 0L)
	  mf         <- mf[c(1L, m)]
	  mf$formula <- gp$fake.formula
	  mf$drop.unused.levels <- TRUE
	  mf[[1]]    <- as.name("model.frame")
	  mf         <- eval(mf)
	  list(gp = gp, formula = formula, mf = mf)
	}