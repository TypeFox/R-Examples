# vegan : adonis

adonis.II <- function(formula,data=NULL,...) {
  f <- formula
  left <- formula[[2]]
  left <- eval(left,data,parent.frame())
  formula[[2]] <- NULL
  right.frame <- model.frame(formula,data,drop.unused.levels=TRUE)
  right <- attr(terms(formula),"term.labels")
  right.formula <- paste(right,collapse="+")
  form <- as.formula(paste0("left~",right.formula))
  aov.tab <- vegan::adonis(form,data=right.frame,...)$aov.tab
  n.terms <- length(right)
  SS <- MS <- F <- p <- rep(0,n.terms)
  for (i in 1:n.terms) {
    names.i <- c(right[-i],right[i])
    right.formula.i <- paste(names.i,collapse="+")
    form.i <- as.formula(paste0("left~",right.formula.i))
    mod.i <- vegan::adonis(form.i,data=right.frame,...)
    SS[i] <- mod.i$aov.tab[right[i],"SumsOfSqs"]
    MS[i] <- mod.i$aov.tab[right[i],"MeanSqs"]
    F[i] <- mod.i$aov.tab[right[i],"F.Model"]
    p[i] <- mod.i$aov.tab[right[i],"Pr(>F)"]
  }
  result <- data.frame(aov.tab$SumsOfSqs,aov.tab$MeanSqs,aov.tab$Df,aov.tab$F.Model,aov.tab$"Pr(>F)")
  rownames(result) <- rownames(aov.tab)
  colnames(result) <- c("Sum Sq","Mean Sq","Df","F","Pr(>F)")
  result[1:n.terms,"Sum Sq"] <- SS
  result[1:n.terms,"Mean Sq"] <- MS
  result[1:n.terms,"F"] <- F
  result[1:n.terms,"Pr(>F)"] <- p
  class(result) <- c("anova","data.frame")
  attr(result,"heading") <- attr(aov.tab,"heading")
  attr(result,"heading")[2] <- "Type II tests"
  attr(result,"heading")[3] <- paste("Response:",deparse(attr(terms(f),"variables")[[2]]))
  return(result)
}
