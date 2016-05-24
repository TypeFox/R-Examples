"do.formula.trellis.xysplom" <-
function(formula, data, na.action=na.pass) {
  ## based on S-Plus do.formula.trellis
  deparen <- function(expr) {
    ## removes all parentheses from an expression, may be overkill here, 
    ## deparse-parse cycle (as in vi()) changes x~y|z to x~(y|z) and the
    ## extra parens in the parse tree surprise do.formula.trellis.
    ## 'Proper' fix may be to change precendence of tilde operator.
    if(mode(expr) == "(") expr <- expr[[2]]
    if(is.recursive(expr))
      for(i in seq(along = expr))
        if(mode(expr[[i]]) != "missing") expr[[i]] <- Recall(expr[[i]])
    expr
  }
  formula <- deparen(formula)
  if (length(formula) == 2) {
    formula <- formula[c(1,2,2)]
    if (length(formula[[3]]) == 3)
      formula[[2]] <- formula[[3]][[2]]
    y.formula <- NULL
  }
  else
    y.formula <- formula[[2]]

  tmp.formula <- formula[1:2]
  
  bar.loc <- if.R(r= if(length(strsplit (deparse(formula), "\\|")[[1]])==2) 2,
##                   if(names(attr(terms(formula, "|"),"specials"))=="|") 2,
                  s=attr(terms(formula, "|"), "specials")$"|")
  if(!is.null(bar.loc)) {
    if(bar.loc == 2) {
      g.formula <- formula[[3]][[3]]
      tmp.formula[[2]] <- g.formula
      g.formula <- tmp.formula
      g <- as.data.frame(model.frame(g.formula, data, na.action=na.action))
      x.formula <- formula[[3]][[2]]
    }
    else if(bar.loc == 1)
      stop("bar.loc == 1")
    ##formula[[2]][[1]] <- as.name("~")
  }
  else {
    g.formula <- NULL
    g <- NULL
    x.formula <- formula[[3]]
  }
  tmp.formula[[2]] <- x.formula
  x.formula <- tmp.formula
  ## acxf2 is a hack to permit xysplom(~data) 
  acxf2 <- as.character(x.formula[[2]])
  if (length(acxf2) == 1)
    if.R(r={x <-
              if(is.null(dim(data))) get(acxf2)
              else data[,acxf2,drop=FALSE]},
         s= x <- get(acxf2, data))
    else
    x <- FALSE
  if (class(x) != "data.frame")
    x <- as.data.frame(model.frame(x.formula, data, na.action=na.action))
  
  if (is.null(y.formula)) {
    y.formula <- x.formula
    y <- x
  }
  else {
    tmp.formula[[2]] <- y.formula
    y.formula <- tmp.formula
    y <- as.data.frame(model.frame(y.formula, data, na.action=na.action))
  }
  
  list(x=x, y=y, group=g,
       x.formula=x.formula,
       g.formula=g.formula,
       y.formula=y.formula)
}

