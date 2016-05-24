MVA.anova <- function(object,...) {
  right <- attr(terms(object),"term.labels")
  right <- right[right %in% attr(terms(object$terminfo),"term.labels")]
  aov.tab <- anova(object,by="terms",...)
  n.terms <- length(right)
  Df <- Var <- F <- p <- rep(0,n.terms)
  for (i in 1:n.terms) {
    names.i <- c(right[-i],right[i])
    right.formula.i <- paste(names.i,collapse="+")
    form.i <- as.formula(paste0(".~",right.formula.i))
    mod.i <- update(object,form.i)
    aov.tab.i <- anova(mod.i,by="terms",...)
    Df[i] <- aov.tab.i[right[i],1]
    Var[i] <- aov.tab.i[right[i],2]
    F[i] <- aov.tab.i[right[i],3]
    p[i] <- aov.tab.i[right[i],4]
  }
  result <- aov.tab
  result[1:n.terms,1] <- Df
  result[1:n.terms,2] <- Var
  result[1:n.terms,3] <- F
  result[1:n.terms,4] <- p
  h <- attr(aov.tab,"heading")[1]
  hsplit <- strsplit(h,split="\n")[[1]]
  hsplit[2] <- "Type II tests"
  attr(result,"heading")[1] <- paste0(paste(hsplit,collapse="\n"),"\n")
  return(result)
}
