#######  examples from User's Guide section 6

  require("dse")
  
  data("eg1.DSE.data", package = "dse") 

  cat("Estimation...\n")

  model.eg1.ls <- estVARXls(trimNA(eg1.DSE.data), warn=F)
  subsample.data <- tfwindow(eg1.DSE.data,start=c(1972,1),end=c(1992,12),
                             warn=FALSE)

  summary(model.eg1.ls)
  model.eg1.ls # or print(model.eg1.ls)
  
  tfplot(model.eg1.ls)
  tfplot(model.eg1.ls, start=c(1990,1))

  checkResiduals(model.eg1.ls, plot.=F, pac=T)


  model.eg1.ss <- estSSfromVARX(trimNA(eg1.DSE.data)) 
  summary(model.eg1.ss)
  model.eg1.ss # or print(model.eg1.ss)
 
  informationTests(model.eg1.ls, model.eg1.ss)
 
