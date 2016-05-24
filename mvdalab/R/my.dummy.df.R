my.dummy.df <- function(data, contr = "contr.niets") {
  dat <- data
  dat$x <- 1:nrow(dat)
  options(contrasts = c(contr, "contr.poly"))
  options(na.action='na.pass')
  new.data <- data.frame(model.matrix(x ~., data = dat)[, -1])
  options(na.action = 'na.omit')
  options(contrasts = c("contr.treatment", "contr.poly"))
  return(new.data)
}