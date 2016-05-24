eblup.saery <-
function(X, ydi, D, md, sigma2edi, model = c("INDEP","AR1","MA1"), plot = FALSE, type = "I", B=100){
  model <- toupper(model)
  model <- match.arg(model)
  switch(model,
         INDEP = eblup.saery.indep(X, ydi, D, md, sigma2edi, plot),
         AR1 = eblup.saery.AR1(X, ydi, D, md, sigma2edi, plot, type, B),
         MA1 = eblup.saery.MA1(X, ydi, D, md, sigma2edi, plot, type, B)
  )
}
