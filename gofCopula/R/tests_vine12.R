gofKendallKS = function(copula, x, param = 0, param.est = T, margins = "ranks", M = 100, execute.times.comp = T){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofKendallKS", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    x = .margins(x, margins)
  }
  if (param.est == T){param = 0}
  
  if (copula == "independence"){fam = 0}
  else if (copula == "gaussian"){fam = 1}
  # else if (copula == "t"){fam = 2}
  else if (copula == "clayton"){fam = 3}
  else if (copula == "gumbel"){fam = 4}
  else if (copula == "frank"){fam = 5}
#   else if (copula == "joe"){fam = 6}
#   else if (copula == "bb1"){fam = 7}
#   else if (copula == "bb6"){fam = 8}
#   else if (copula == "bb7"){fam = 9}
#   else if (copula == "bb8"){fam = 10}
#   else if (copula == "survival clayton"){fam =13}
#   else if (copula == "survival gumbel"){fam =14}
#   else if (copula == "survival joe"){fam =16}
#   else if (copula == "survival bb1"){fam =17}
#   else if (copula == "survival bb6"){fam =18}
#   else if (copula == "survival bb7"){fam =19}
#   else if (copula == "survival bb8"){fam =20}
#   else if (copula == "90 clayton"){fam =23}
#   else if (copula == "90 gumbel"){fam =24}
#   else if (copula == "90 joe"){fam =26}
#   else if (copula == "90 bb1"){fam =27}
#   else if (copula == "90 bb6"){fam =28}
#   else if (copula == "90 bb7"){fam =29}
#   else if (copula == "90 bb8"){fam =30}
#   else if (copula == "270 clayton"){fam =33}
#   else if (copula == "270 gumbel"){fam =34}
#   else if (copula == "270 joe"){fam =36}
#   else if (copula == "270 bb1"){fam =37}
#   else if (copula == "270 bb6"){fam =38}
#   else if (copula == "270 bb7"){fam =39}
#   else if (copula == "270 bb8"){fam =40}
  else if (copula == "t"){stop("t-copula is not implemented for the GoF test based on Kendall's process")}
  else {stop("This copula is not implemented for gofKendallKS.")}
  res = BiCopGofTest(u1 = x[,1], u2 = x[,2], family = fam, par = param, par2 = 0, method = "kendall", B = M)
  
  structure(class = "gofCOP", 
            list(method = sprintf("Parametric bootstrap goodness-of-fit test (Kolmogorov-Smirnof) based on Kendall's process"), 
                 erg.tests = matrix(c(res$p.value.KS, res$statistic.KS), ncol = 2, 
                                    dimnames = list("KendallKS", c("p.value", "test statistic")))))
}
