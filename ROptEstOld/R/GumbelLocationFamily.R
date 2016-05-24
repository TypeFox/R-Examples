##################################################################
## Gumbel location family
##################################################################
GumbelLocationFamily <- function(loc = 0, scale = 1, trafo){
    if(missing(trafo)) trafo <- matrix(1, dimnames = list("loc","loc"))
    modParam <- function(theta){}
    body(modParam) <- substitute({ Gumbel(loc = theta, scale = sd) },
                                 list(sd = scale))
    res <- L2LocationFamily(loc = loc,
                     name = "Gumbel location family",
                     locname = c("loc"="loc"),
                     centraldistribution = Gumbel(loc = 0, scale = scale),
                     modParam = modParam,
                     LogDeriv = function(x) (1 - exp(-x/scale))/scale,
                     L2derivDistr.0 = (1 - Exp(rate = 1))/scale,
                     FisherInfo.0 = matrix(1/scale^2,
                                    dimnames = list("loc","loc")),
                     distrSymm = NoSymmetry(),
                     L2derivSymm = FunSymmList(NonSymmetric()),
                     L2derivDistrSymm = DistrSymmList(NoSymmetry()),
                     trafo = trafo, .returnClsName = "GumbelLocationFamily")
    if(!is.function(trafo))
       f.call <- substitute(GumbelLocationFamily(loc = l, scale = s,
                          trafo = matrix(Tr, dimnames = list("loc","loc"))),
  	                     list(l = loc, s = scale, Tr = trafo))
  	else
       f.call <- substitute(GumbelLocationFamily(loc = l, scale = s, trafo = Tr),
  	                     list(l = loc, s = scale, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}
