fit.saery <-
function(X, ydi, D, md, sigma2edi, model=c("INDEP","AR1","MA1"), conf.level = 0.95){
  model <- toupper(model)
  model <- match.arg(model)
  switch(model,
         INDEP = fit.saery.indep(X, ydi, D, md, sigma2edi, conf.level),
         AR1 = fit.saery.AR1(X, ydi, D, md, sigma2edi, conf.level),
         MA1 = fit.saery.MA1(X, ydi, D, md, sigma2edi, conf.level)
  )
}
