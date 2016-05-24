# uses lme and groupedData from nlme
# uses as.sets from sets for power set of X

splitvars <- function(fixed) {
  if (length(fixed) == 0 || fixed == "") {
    return(fixed)
  }
  fixed.terms = fixed
  seps = c(':','\\*','\\|','/')
  for (sep in seps) {
    tokens = strsplit(fixed.terms, sep)
    fixed.terms = c()
    for (tok in tokens) {
      fixed.terms = c(fixed.terms, sub("-","",unlist(tok)))
    }
    fixed.terms = unlist(fixed.terms)
  }
  return(unique(fixed.terms))
}


getIC <- function(fit, ictype) {
  if (substr(ictype,1,3) == "AIC") {
    icval = AIC(fit)
    if (ictype == "AICc") {
      # add the correction
      npar = length(fit$coefficients) + 1
      n = nobs(fit)
      icval = icval + 2 * npar * (npar + 1) / (n - npar - 1)
    }
  } else if (ictype == "BIC") {
    icval = BIC(fit)
  } else {
    warning(paste("Unrecognized IC option", ictype, "using AIC instead"))
    icval = AIC(fit)
  }
  icval
}

fmi <- function(coly, ...) UseMethod("fmi")

fmi.default <- function(coly, candidates=c(""), fixed=c(""), data=list(), modeltype="lm", random=~1, ic="AIC", ...) {
  return(FindMinIC.default(coly, candidates, fixed, data, modeltype, random, ic, ...))
}

fmi.formula <- function(formula, data=list(), na.action=na.omit, fixed=c(""), random=~1, ...) {
  return(FindMinIC.formula(formula, data, na.action, fixed, random, ...))
}

FindMinIC <- function(coly, ...) UseMethod("FindMinIC")

FindMinIC.default <- function(coly, candidates=c(""), fixed=c(""), data=list(), modeltype="lm", random=~1, ic = "AIC", ...){
  fixed.vars = splitvars(fixed)
  cand.vars = splitvars(candidates)
  random.terms = paste(random)[-1]
  group.vars = splitvars(random.terms)
  do.lme = FALSE
  if (modeltype == "lme") {
    do.lme = TRUE
    comb.vars = unique(c(fixed.vars, cand.vars, coly, group.vars))
  } else {
    comb.vars = unique(c(fixed.vars, cand.vars, coly))
  }
  # strip out whitespace
  comb.vars = gsub(" ","", comb.vars, fixed=TRUE)

  comb.vars = comb.vars[which(comb.vars!="1")]
  comb.vars = comb.vars[which(comb.vars!=".")]

  # This subsetting should work now, but it's not clear that it's necessary or useful
  # so for now passing in entire data frame
  # tmp.df = subset(data, select=comb.vars)
  tmp.df = data

  if (do.lme) {
    # allowing full formula specification for groupedData
    # note that using groupedData instead of passing random into lme directly
    # keeps the user from being able to use certain syntaxes for the specification
    # of the random effects.  Specifically, crossed random effects specified by something like 
    #    random=pdBlocked(list(pdIdent(~1), pdIdent(~sample-1), pdIdent(~dilut-1)))
    # are not allowed.  If these become necessary, groupedData use here will need to be rethought
    gdcall = parse(text = paste("groupedData(",
                     coly, " ~ ", random.terms,
                     ", data=tmp.df)"))
    tmp.gds = eval(gdcall)
    
  } else {
    tmp.gds = tmp.df
  }
  fixed.model = ""
  firstelem = TRUE
  # for the rest of the fixed elements
  for (f in fixed) {
    if (f != "") {
      if (firstelem) {
        fixed.model = paste(fixed.model, f, sep="")
        firstelem = FALSE
      } else {
        fixed.model = paste(fixed.model, "+", f, sep="")
      }
    }
  }

  results = list()

#
# define model space for group analysis
#
  xpnames.set = as.set(candidates)
  Xp = parse(text = 2^xpnames.set) # power set
  M = length(Xp)

  for(m in 1:M){ # all possible models given covariates in candidates
    if(m == 1) model = "1"
    if(m != 1){
      model = NULL
      op = "+"
      xtext = as.character(Xp[[m] ] )[-1] # remove string "list()"
      P = length(xtext)
#
# the following 5 lines could be more well thought out, I think...
#
      if(P == 1) op = ""
      for(p in 1:(P - 1) ){ # weird for P == 1, but it works...
        model = paste(model, xtext[p], op, sep = "")
      }
      if(P > 1) model = paste(model, xtext[P], sep = "")
    }
    if (do.lme) {
      call = parse(text = paste("withRestarts(", modeltype, "(",
                     coly,
                     " ~ ",
                     fixed.model,
                     " + ",
                     model,
                     ", method = 'ML'",
                     ", random = random",
                     ", data = tmp.gds",
                     ",...)",
                     ")",
                     sep = "") )
      fit = try(eval(call), silent = TRUE)
      if(class(fit) == "try-error") {
        warning(paste("Received error <", attr(fit, "condition"), "> while calling ", call))
        next
      }

      res = list(call = fit$call, IC = getIC(fit, ic), ictype = ic, formula = formula(fit))
      class(res) = "cm" # for candidate model
      results[[length(results)+1]] = res

    } else {
      # TODO ok to convert modeltype directly into R command?
      call = parse(text = paste(modeltype,
                     "(",
                     coly,
                     " ~ ",
                     fixed.model,
                     " + ",
                     model,
                     ", data = tmp.gds",
                     ",...)",
                     sep = "") )
      fit = try(eval(call), silent = TRUE);
      if(class(fit) == "try-error") {
        warning(paste("Received error <", attr(fit, "condition"), "> while calling ", call))
        next
      }

      res = list(call = fit$call, IC = getIC(fit, ic), ictype = ic, formula = formula(fit))
      class(res) = "cm" # for candidate model
      results[[length(results)+1]] = res
    }
  }

  if (length(results) > 1) {
    results = results[order(as.numeric(lapply(results, IC.cm)))]
  }
  first = results[[1]]
  if (do.lme) {
    # for lme, need to do REML instead of ML here
    call = parse(text = paste(modeltype,"(",
                   deparse(first$call$fixed, width.cutoff=500),
                   ", data = tmp.gds",
                   ", random = random",
                   ", method = 'REML'",
                   ",...)",
                   sep = ""
                   ) )
  } else {
    call = first$call
  }
  fit = eval(call)
  first$fit = fit
  answer = list(results=results, data=tmp.gds, first=first,
                modeltype=modeltype, random=random)
  class(answer) = "cmList"

  return(answer)
} # end FindMinIC.default

FindMinIC.formula <- function(formula, data=list(), na.action=na.omit, fixed=c(""), random=~1, ...)
{
    mf <- model.frame(formula=formula, data=data, na.action=na.action)
    nms <- names(mf)
    trms <- unlist(strsplit(paste(formula)[-1], '\\+'))
    trms <- gsub(" ","", trms, fixed=TRUE)
    coly <- trms[[1]]
    candidates <- unique(c(trms[-1], nms[-1]))
    est <- FindMinIC.default(coly, candidates, data=data, fixed=fixed, random=random, ...)
    est$call <- match.call()
    est$formula <- formula
    est$terms <- attr(mf, "terms")

    est
}
