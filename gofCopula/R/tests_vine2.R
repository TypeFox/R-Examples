########################################################################
gofWhite = function(copula, x, M = 1000, param = 0, param.est = T, df = 0, df.est = T, margins = "ranks", execute.times.comp = T){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofWhite", M = M, print.res = F)
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
  if (df.est == T){df = 0}
  
  if (copula == "independence"){fam = 0}
  else if (copula == "gaussian"){fam = 1}
  else if (copula == "t"){fam = 2}
  else if (copula == "clayton"){fam = 3}
  else if (copula == "gumbel"){fam = 4}
  else if (copula == "frank"){fam = 5}
#   else if (copula == "joe"){fam = 6}
#   else if (copula == "survival clayton"){fam =13}
#   else if (copula == "survival gumbel"){fam =14}
#   else if (copula == "survival joe"){fam =16}
#   else if (copula == "90 clayton"){fam =23}
#   else if (copula == "90 gumbel"){fam =24}
#   else if (copula == "90 joe"){fam =26}
#   else if (copula == "270 clayton"){fam =33}
#   else if (copula == "270 gumbel"){fam =34}
#   else if (copula == "270 joe"){fam =36}
#   else if (copula == "bb1" || copula == "survival bb1" || copula == "90 bb1" || copula == "270 bb1"){stop("BB1 copula is not implemented for the GoF test based on White's information matrix equality.")}
#   else if (copula == "bb6" || copula == "survival bb6" || copula == "90 bb6" || copula == "270 bb6"){stop("BB6 copula is not implemented for the GoF test based on White's information matrix equality.")}
#   else if (copula == "bb7" || copula == "survival bb7" || copula == "90 bb7" || copula == "270 bb7"){stop("BB7 copula is not implemented for the GoF test based on White's information matrix equality.")}
#   else if (copula == "bb8" || copula == "survival bb8" || copula == "90 bb8" || copula == "270 bb8"){stop("BB8 copula is not implemented for the GoF test based on White's information matrix equality.")}
  else {stop("This copula is not implemented for gofWhite.")}
  res = BiCopGofTest(u1 = x[,1], u2 = x[,2], family = fam, par = param, par2 = df, method = "white", max.df = dim(x) - 1, B = M)
  structure(class = "gofCOP", 
            list(method = sprintf("Goodness-of-fit test based on White's information equality matrix"), 
                 erg.tests = matrix(c(res$p.value, res$statistic), ncol = 2, 
                                    dimnames = list("White", c("p.value", "test statistic")))))
}
