

vred = function(a1, a2=NULL, vredTyp = "hyp", variogramModel, 
         pdf1=NULL, pdf2=NULL, aover=NULL, dist=NULL, inner = 0, resol=5) {
#  print(asdf)
#  return(0)
  model = variogramModel$model
  imod = imodel(model)
  param = variogramModel$params
  ci = 0
  if (vredTyp == "ind"){
    if (is(a1,"list")) {
      a2 = a1[[2]]
      a1 = a1[[1]]
    } else if (is.null(a2)) a2 = a1
    a1 = coordinates(a1)
    a2 = coordinates(a2)
    ip1 = dim(a1)[1]
    ip2 = dim(a2)[1]
    vreda = .Fortran("vredind", ci, ip1, ip2, a1, a2, length(param), param, imod)
#  } else if (vredTyp == "pdf") {
#    vreda = .Fortran("vredpdf",ci,c1,c2,ip1,ip2,ipb,pdf1,pdf2,pdfb,length(param),param,model)
  } else if (vredTyp == "hyp") {
    vreda = .Fortran("vredhyp", ci, a1, a2, dist, length(param), param,
          as.integer(resol), imod)
  }
###### Nugget needs to be implemented
  if (!is.null(aover)) {
    nug = 0
  } else {
    nug = 0
  }

  vreda[[1]]
}
