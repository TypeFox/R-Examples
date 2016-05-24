opsff_compare_logic <- function(e1, e2, ops) {
  res <- ff(vmode="logical", length=max(c(length(e1), length(e2))))
  if(inherits(e1, "ff_vector") && inherits(e2, "ff_vector")){
    if(length(e1) != length(e2)) stop("operator requires equal length ff_vectors")
    for(i in chunk(e1, RECORDBYTES = .rambytes[vmode(e1)] + .rambytes[vmode(e2)])){
      res[i] <- switch(ops,                       
                       "==" = e1[i] == e2[i],
                       ">" = e1[i] > e2[i],
                       "<" = e1[i] < e2[i],
                       "!=" = e1[i] != e2[i],
                       "<=" = e1[i] <= e2[i],
                       ">=" = e1[i] >= e2[i],
                       "&" = e1[i] & e2[i],
                       "|" = e1[i] | e2[i]
                       )      
    }
  }else if(inherits(e1, "ff_vector")){
    if(length(e2) != 1 & ops != "!") stop("operator requires length 1 for e2, recycling not possible")
    for (i in chunk(e1)){
      res[i] <- switch(ops,                       
                       "==" = e1[i] == e2,
                       ">" = e1[i] > e2,
                       "<" = e1[i] < e2,
                       "!=" = e1[i] != e2,
                       "<=" = e1[i] <= e2,
                       ">=" = e1[i] >= e2,
                       "&" = e1[i] & e2,
                       "|" = e1[i] | e2,
                       "!" = !e1[i]
                       )      
    }    
  }else if(inherits(e2, "ff_vector")){
    if(length(e1) != 1) stop("operator requires length 1 for e1, recycling not possible")
    for (i in chunk(e2)){
      res[i] <- switch(ops,                       
                       "==" = e1 == e2[i],
                       ">" = e1 > e2[i],
                       "<" = e1 < e2[i],
                       "!=" = e1 != e2[i],
                       "<=" = e1 <= e2[i],
                       ">=" = e1 >= e2[i],
                       "&" = e1 & e2[i],
                       "|" = e1 | e2[i]
                       )      
    }
  }  
  res
}
mathff_math <- function(e1, math) {
  i1 <- i2 <- 0L
  opts <- options(warn=-1)
  on.exit(options(opts))
  expr <- switch(math,                       
                 "abs" = expression(abs(e1[i1:i2])),
                 "sign" = expression(sign(e1[i1:i2])),
                 "sqrt" = expression(sqrt(e1[i1:i2])),
                 "ceiling" = expression(ceiling(e1[i1:i2])),
                 "floor" = expression(floor(e1[i1:i2])),
                 "trunc" = expression(trunc(e1[i1:i2])),
                 "log10" = expression(log10(e1[i1:i2])),
                 "log2" = expression(log2(e1[i1:i2])),
                 "log1p" = expression(log1p(e1[i1:i2])),
                 "acos" = expression(acos(e1[i1:i2])),
                 "acosh" = expression(acosh(e1[i1:i2])),
                 "asin" = expression(asin(e1[i1:i2])),
                 "asinh" = expression(asinh(e1[i1:i2])),
                 "atan" = expression(atan(e1[i1:i2])),
                 "atanh" = expression(atanh(e1[i1:i2])),
                 "exp" = expression(exp(e1[i1:i2])),
                 "expm1" = expression(expm1(e1[i1:i2])),
                 "cos" = expression(cos(e1[i1:i2])),
                 "cosh" = expression(cosh(e1[i1:i2])),
                 "sin" = expression(sin(e1[i1:i2])),
                 "sinh" = expression(sinh(e1[i1:i2])),
                 "tan" = expression(tan(e1[i1:i2])),
                 "tanh" = expression(tanh(e1[i1:i2])),
                 "gamma" = expression(gamma(e1[i1:i2])),
                 "lgamma" = expression(lgamma(e1[i1:i2])),
                 "digamma" = expression(digamma(e1[i1:i2])),
                 "trigamma" = expression(trigamma(e1[i1:i2]))
                 ) 
  ffvecapply(eval(expr), X = e1, VMODE = NULL, VBYTES = NULL, 
             RETURN = TRUE, FF_RETURN = TRUE, FROM = "i1", TO = "i2")  
}
mathff_math_log <- function(e1, math, base) {
  i1 <- i2 <- 0L
  expr <- expression(log(e1[i1:i2], base=base))
  ffvecapply(eval(expr), X = e1, VMODE = NULL, VBYTES = NULL, 
             RETURN = TRUE, FF_RETURN = TRUE, FROM = "i1", TO = "i2")  
}
mathff_math2 <- function(e1, math, digits) {
  i1 <- i2 <- 0L
  if(math == "round"){
    ffvecapply(round(e1[i1:i2], digits=digits), X = e1, VMODE = NULL, VBYTES = NULL, 
               RETURN = TRUE, FF_RETURN = TRUE, FROM = "i1", TO = "i2")
  }else if(math == "signif"){
    ffvecapply(signif(e1[i1:i2], digits=digits), X = e1, VMODE = NULL, VBYTES = NULL, 
               RETURN = TRUE, FF_RETURN = TRUE, FROM = "i1", TO = "i2")
  }
}
opsff_arith <- function(e1, e2, ops) {
  res <- ff(vmode="double", length=max(c(length(e1), length(e2))))
  if(inherits(e1, "ff_vector") && inherits(e2, "ff_vector")){
    if(length(e1) != length(e2)) stop("requires equal length ff_vectors")
    for(i in chunk(e1, RECORDBYTES = .rambytes[vmode(e1)] + .rambytes[vmode(e2)])){
      res[i] <- switch(ops,                       
                       "+" = e1[i] + e2[i],
                       "-" = e1[i] - e2[i],
                       "*" = e1[i] * e2[i],
                       "^" = e1[i] ^ e2[i],
                       "%%" = e1[i] %% e2[i],
                       "%/%" = e1[i] %/% e2[i],
                       "/" = e1[i] / e2[i]
                       )      
    }
  }else if(inherits(e1, "ff_vector")){
    if(length(e2) != 1) stop("requires length 1 for e2, recycling not possible")
    for (i in chunk(e1)){
      res[i] <- switch(ops,                       
                       "+" = e1[i] + e2,
                       "-" = e1[i] - e2,
                       "*" = e1[i] * e2,
                       "^" = e1[i] ^ e2,
                       "%%" = e1[i] %% e2,
                       "%/%" = e1[i] %/% e2,
                       "/" = e1[i] / e2
                       )      
    }    
  }else if(inherits(e2, "ff_vector")){
    if(length(e1) != 1) stop("requires length 1 for e1, recycling not possible")
    for (i in chunk(e2)){
      res[i] <- switch(ops,                       
                       "+" = e1 + e2[i],
                       "-" = e1 - e2[i],
                       "*" = e1 * e2[i],
                       "^" = e1 ^ e2[i],
                       "%%" = e1 %% e2[i],
                       "%/%" = e1 %/% e2[i],
                       "/" = e1 / e2[i],
                       )      
    }
  }  
  res
}


#' Arithmetic Operators for ff vectors
#'
#' These binary operators perform arithmetic on numeric ff vectors. 
#' Arith family: \cr
#' \itemize{
#' \item{Arith: }{\code{"+", "-", "*", "/", "^", "\%\%", "\%/\%"}}
#' }
#' The operators require either \code{x} or \code{y} to be an \code{ff_vector} or both. In case either \code{x} or \code{y} is not an \code{ff_vector},
#' the other object needs to be of length 1. Recycling is not implemented.
#'
#' @rdname ff_arithmetic
#' @aliases Arithmetic +.ff_vector -.ff_vector *.ff_vector /.ff_vector ^.ff_vector %%.ff_vector %/%.ff_vector
#' @param x either a numeric \code{ff_vector} or a vector of length 1 in RAM in which case y should be an \code{ff_vector}
#' @param y either a numeric \code{ff_vector} or a vector of length 1 in RAM in which case x should be an \code{ff_vector}
#' @return an \code{ff_vector}. For the definition of the operators see the base package of R.

#' @rdname ff_arithmetic
#' @usage \method{+}{ff_vector} (x, y)
#' @export
"+.ff_vector" <- function(x, y) opsff_arith(x, y, "+")
#' @rdname ff_arithmetic
#' @usage \method{-}{ff_vector} (x, y)
#' @export
"-.ff_vector" <- function(x, y) opsff_arith(x, y, "-")
#' @rdname ff_arithmetic
#' @usage \method{*}{ff_vector} (x, y)
#' @export
"*.ff_vector" <- function(x, y) opsff_arith(x, y, "*")
#' @rdname ff_arithmetic
#' @usage \method{/}{ff_vector} (x, y)
#' @export
"/.ff_vector" <- function(x, y) opsff_arith(x, y, "/")
#' @rdname ff_arithmetic
#' @usage \method{^}{ff_vector} (x, y)
#' @export
"^.ff_vector" <- function(x, y) opsff_arith(x, y, "^")
#' @rdname ff_arithmetic
#' @usage \method{\%\%}{ff_vector} (x, y)
#' @export
"%%.ff_vector" <- function(x, y) opsff_arith(x, y, "%%")
#' @rdname ff_arithmetic
#' @usage \method{\%/\%}{ff_vector} (x, y)
#' @export
"%/%.ff_vector" <- function(x, y) opsff_arith(x, y, "%/%")

#' Ops for ff vectors
#'
#' These operators implement \code{ff_vector} specific operators and handle the following operators from the
#' Ops family: \cr
#' \itemize{
#' \item{Compare: }{\code{"==", "!=", "<", "<=", ">=", ">"}}
#' \item{Logic: }{\code{"&", "|", "!"}}
#' }
#' The operators require either \code{x} or \code{y} to be an \code{ff_vector} or both. In case either \code{x} or \code{y} is not an \code{ff_vector},
#' the other object needs to be of length 1. Recycling is not implemented.
#'
#' @rdname ff_ops
#' @aliases Operators >.ff_vector <.ff_vector ==.ff_vector !=.ff_vector <=.ff_vector >=.ff_vector &.ff_vector |.ff_vector !.ff_vector  
#' @param x either a numeric \code{ff_vector} or a vector of length 1 in RAM in which case y should be an \code{ff_vector}
#' @param y either a numeric \code{ff_vector} or a vector of length 1 in RAM in which case x should be an \code{ff_vector}
#' @return an \code{ff_vector}. For the definition of the operators see the base package of R.

#' @rdname ff_ops
#' @usage \method{>}{ff_vector} (x, y)
#' @export
">.ff_vector" <- function(x, y) opsff_compare_logic(x, y, ">")
#' @rdname ff_ops
#' @usage \method{<}{ff_vector} (x, y)
#' @export
"<.ff_vector" <- function(x, y) opsff_compare_logic(x, y, "<")
#' @rdname ff_ops
#' @usage \method{==}{ff_vector} (x, y)
#' @export
"==.ff_vector" <- function(x, y) opsff_compare_logic(x, y, "==")
#' @rdname ff_ops
#' @usage \method{!=}{ff_vector} (x, y)
#' @export
"!=.ff_vector" <- function(x, y) opsff_compare_logic(x, y, "!=")
#' @rdname ff_ops
#' @usage \method{<=}{ff_vector} (x, y)
#' @export
"<=.ff_vector" <- function(x, y) opsff_compare_logic(x, y, "<=")
#' @rdname ff_ops
#' @usage \method{>=}{ff_vector} (x, y)
#' @export
">=.ff_vector" <- function(x, y) opsff_compare_logic(x, y, ">=")
#' @rdname ff_ops
#' @usage \method{&}{ff_vector} (x, y)
#' @export
"&.ff_vector" <- function(x, y) opsff_compare_logic(x, y, "&")
#' @rdname ff_ops
#' @usage \method{|}{ff_vector} (x, y)
#' @export
"|.ff_vector" <- function(x, y) opsff_compare_logic(x, y, "|")
#' @rdname ff_ops
#' @usage \method{!}{ff_vector} (x)
#' @export
"!.ff_vector" <- function(x) opsff_compare_logic(x, NULL, "!")


#' Math for ff vectors
#'
#' These mathematical functions implement \code{ff_vector} specific math and handle the following functions from the
#' Math family: \cr
#' \itemize{
#' \item{Math: }{\code{"abs", "sign", "sqrt", "ceiling", "floor", "trunc", "log", "log10", "log2", "log1p", "acos", "acosh", "asin", "asinh", "atan", "atanh", "exp", "expm1", "cos", "cosh", "sin", "sinh", "tan", "tanh", "gamma", "lgamma", "digamma", "trigamma"}}
#' \item{Math2: }{\code{"round", "signif"}}
#' }
#' The operators require \code{x} to be an \code{ff_vector}.
#'
#' @rdname ff_math
#' @aliases Math abs.ff_vector sign.ff_vector sqrt.ff_vector ceiling.ff_vector floor.ff_vector trunc.ff_vector log.ff_vector log10.ff_vector log2.ff_vector log1p.ff_vector acos.ff_vector acosh.ff_vector asin.ff_vector asinh.ff_vector atan.ff_vector atanh.ff_vector exp.ff_vector expm1.ff_vector cos.ff_vector cosh.ff_vector sin.ff_vector sinh.ff_vector tan.ff_vector tanh.ff_vector gamma.ff_vector lgamma.ff_vector digamma.ff_vector trigamma.ff_vector round.ff_vector signif.ff_vector
#' @param x a numeric \code{ff_vector}
#' @param base base for \code{log}
#' @param digits digits for \code{round} and \code{signif}
#' @param ... for \code{trunc}, currently not used
#' @return an \code{ff_vector}. For the definition of the operators see the base package of R.

#' @rdname ff_math
#' @usage \method{abs}{ff_vector} (x)
#' @export
"abs.ff_vector" <- function(x) mathff_math(x, "abs")
#' @rdname ff_math
#' @usage \method{sign}{ff_vector} (x)
#' @export
"sign.ff_vector" <- function(x) mathff_math(x, "sign")
#' @rdname ff_math
#' @usage \method{sqrt}{ff_vector} (x)
#' @export
"sqrt.ff_vector" <- function(x) mathff_math(x, "sqrt")
#' @rdname ff_math
#' @usage \method{ceiling}{ff_vector} (x)
#' @export
"ceiling.ff_vector" <- function(x) mathff_math(x, "ceiling")
#' @rdname ff_math
#' @usage \method{floor}{ff_vector} (x)
#' @export
"floor.ff_vector" <- function(x) mathff_math(x, "floor")
#' @rdname ff_math
#' @usage \method{trunc}{ff_vector} (x, ...)
#' @export
"trunc.ff_vector" <- function(x, ...) mathff_math(x, "trunc")
#' @rdname ff_math
#' @usage \method{log10}{ff_vector} (x)
#' @export
"log10.ff_vector" <- function(x) mathff_math(x, "log10")
#' @rdname ff_math
#' @usage \method{log2}{ff_vector} (x)
#' @export
"log2.ff_vector" <- function(x) mathff_math(x, "log2")
#' @rdname ff_math
#' @usage \method{log1p}{ff_vector} (x)
#' @export
"log1p.ff_vector" <- function(x) mathff_math(x, "log1p")
#' @rdname ff_math
#' @usage \method{acos}{ff_vector} (x)
#' @export
"acos.ff_vector" <- function(x) mathff_math(x, "acos")
#' @rdname ff_math
#' @usage \method{acosh}{ff_vector} (x)
#' @export
"acosh.ff_vector" <- function(x) mathff_math(x, "acosh")
#' @rdname ff_math
#' @usage \method{asin}{ff_vector} (x)
#' @export
"asin.ff_vector" <- function(x) mathff_math(x, "asin")
#' @rdname ff_math
#' @usage \method{asinh}{ff_vector} (x)
#' @export
"asinh.ff_vector" <- function(x) mathff_math(x, "asinh")
#' @rdname ff_math
#' @usage \method{atan}{ff_vector} (x)
#' @export
"atan.ff_vector" <- function(x) mathff_math(x, "atan")
#' @rdname ff_math
#' @usage \method{atanh}{ff_vector} (x)
#' @export
"atanh.ff_vector" <- function(x) mathff_math(x, "atanh")
#' @rdname ff_math
#' @usage \method{exp}{ff_vector} (x)
#' @export
"exp.ff_vector" <- function(x) mathff_math(x, "exp")
#' @rdname ff_math
#' @usage \method{expm1}{ff_vector} (x)
#' @export
"expm1.ff_vector" <- function(x) mathff_math(x, "expm1")
#' @rdname ff_math
#' @usage \method{cos}{ff_vector} (x)
#' @export
"cos.ff_vector" <- function(x) mathff_math(x, "cos")
#' @rdname ff_math
#' @usage \method{cosh}{ff_vector} (x)
#' @export
"cosh.ff_vector" <- function(x) mathff_math(x, "cosh")
#' @rdname ff_math
#' @usage \method{sin}{ff_vector} (x)
#' @export
"sin.ff_vector" <- function(x) mathff_math(x, "sin")
#' @rdname ff_math
#' @usage \method{sinh}{ff_vector} (x)
#' @export
"sinh.ff_vector" <- function(x) mathff_math(x, "sinh")
#' @rdname ff_math
#' @usage \method{tan}{ff_vector} (x)
#' @export
"tan.ff_vector" <- function(x) mathff_math(x, "tan")
#' @rdname ff_math
#' @usage \method{tanh}{ff_vector} (x)
#' @export
"tanh.ff_vector" <- function(x) mathff_math(x, "tanh")
#' @rdname ff_math
#' @usage \method{gamma}{ff_vector} (x)
#' @export
"gamma.ff_vector" <- function(x) mathff_math(x, "gamma")
#' @rdname ff_math
#' @usage \method{lgamma}{ff_vector} (x)
#' @export
"lgamma.ff_vector" <- function(x) mathff_math(x, "lgamma")
#' @rdname ff_math
#' @usage \method{digamma}{ff_vector} (x)
#' @export
"digamma.ff_vector" <- function(x) mathff_math(x, "digamma")
#' @rdname ff_math
#' @usage \method{trigamma}{ff_vector} (x)
#' @export
"trigamma.ff_vector" <- function(x) mathff_math(x, "trigamma")

#' @rdname ff_math
#' @usage \method{log}{ff_vector} (x, base)
#' @export
"log.ff_vector" <- function(x, base = exp(1)) mathff_math_log(x, "log", base = base)

#' @rdname ff_math
#' @usage \method{round}{ff_vector} (x, digits)
#' @export
"round.ff_vector" <- function(x, digits = 0) mathff_math2(x, "round", digits = digits)
#' @rdname ff_math
#' @usage \method{signif}{ff_vector} (x, digits)
#' @export
"signif.ff_vector" <- function(x, digits = 6) mathff_math2(x, "signif", digits = digits)



