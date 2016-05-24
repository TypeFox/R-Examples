pcor <-
function(x,y,z,semi=FALSE,use="complete.obs",method=c("pearson","kendall","spearman")) {
  if (is.list(z)) {
    z <- as.data.frame(do.call("cbind",z))
  }
  tab <- data.frame(x=x,y=y,z)
  form.x <- as.formula(paste0("x~",paste(colnames(tab)[-c(1,2)],collapse="+")))
  form.y <- as.formula(paste0("y~",paste(colnames(tab)[-c(1,2)],collapse="+")))
  m.x <- lm(form.x,data=tab,na.action="na.exclude")
  m.y <- lm(form.y,data=tab,na.action="na.exclude")
  value <- if (!semi) {
    cor(resid(m.x),resid(m.y),use=use,method=method)
  } else {
    cor(tab$x,resid(m.y),use=use,method=method)
  }
  result <- value
  return(result)
}
