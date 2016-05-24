
setClass("S2",
slots = c(s2 = "numeric", y2 = "numeric", u2 = "numeric", p.value.y2="numeric",
   p.value.s2="numeric", p.value.u2="numeric", data.name="character"))


setGeneric("S2", function(object) standardGeneric("S2"))


setMethod("show",
signature = "S2",
definition = function(object) {
    cat("            Chi-squared type tests for Multivariate Normality", "\n\n", sep = " ")
  ##  cat("---------------------------------------------", "\n", sep = " ")
    cat("  data :", object@data.name, "\n\n", sep = " ")
    cat("  McCulloch (S2)            :", object@s2, "\n", sep = " ")
    cat("  p-value.S2                :", object@p.value.s2, "\n\n", sep = " ")
    cat("  Nikulin-Rao-Robson (Y2)   :", object@y2, "\n", sep = " ")
    cat("  p-value.Y2                :", object@p.value.y2, "\n\n", sep = " ")
    cat("  Dzhaparidze-Nikulin (U2)  :", object@u2, "\n", sep = " ")
    cat("  p-value.U2                :", object@p.value.u2, "\n\n", sep = " ")
    cat(if(object@p.value.s2 > 0.05){"  Result  : Data are multivariate normal (sig.level = 0.05)"}
        else {"  Result  : Data are not multivariate normal (sig.level = 0.05)"},"\n")
  ##  cat("---------------------------------------------", "\n\n", sep = " ")
    invisible(NULL)
})


## ------------------

setClass("ad",
slots = c(AD = "numeric", p.value="numeric", data.name="character"))

setGeneric("ad", function(object) standardGeneric("ad"))


setMethod("show",
signature = "ad",
definition = function(object) {
    cat("            Anderson-Darling test for Multivariate Normality", "\n\n", sep = " ")
  ##  cat("---------------------------------------------", "\n", sep = " ")
    cat("  data :", object@data.name, "\n\n", sep = " ")
    cat("  AD              :", object@AD, "\n", sep = " ")
    cat("  p-value         :", object@p.value, "\n\n", sep = " ")
    cat(if(object@p.value > 0.05){"  Result  : Data are multivariate normal (sig.level = 0.05)"}
        else {"  Result  : Data are not multivariate normal (sig.level = 0.05)"},"\n")
  ##  cat("---------------------------------------------", "\n\n", sep = " ")
    invisible(NULL)
})



## ------------------

setClass("cm",
slots = c(CM = "numeric", p.value="numeric", data.name="character"))

setGeneric("cm", function(object) standardGeneric("cm"))


setMethod("show",
signature = "cm",
definition = function(object) {
    cat("            Cramer-von Mises test for Multivariate Normality", "\n\n", sep = " ")
  ##  cat("---------------------------------------------", "\n", sep = " ")
    cat("  data :", object@data.name, "\n\n", sep = " ")
    cat("  CM              :", object@CM, "\n", sep = " ")
    cat("  p-value         :", object@p.value, "\n\n", sep = " ")
    cat(if(object@p.value > 0.05){"  Result  : Data are multivariate normal (sig.level = 0.05)"}
        else {"  Result  : Data are not multivariate normal (sig.level = 0.05)"},"\n")
  ##  cat("---------------------------------------------", "\n\n", sep = " ")
    invisible(NULL)
})

## ------------------

setClass("dh",
slots = c(DH = "numeric", p.value="numeric", data.name="character"))

setGeneric("dh", function(object) standardGeneric("dh"))


setMethod("show",
signature = "dh",
definition = function(object) {
    cat("            Doornik-Hansen test for Multivariate Normality", "\n\n", sep = " ")
  ##  cat("---------------------------------------------", "\n", sep = " ")
    cat("  data :", object@data.name, "\n\n", sep = " ")
    cat("  DH              :", object@DH, "\n", sep = " ")
    cat("  p-value         :", object@p.value, "\n\n", sep = " ")
    cat(if(object@p.value > 0.05){"  Result  : Data are multivariate normal (sig.level = 0.05)"}
        else {"  Result  : Data are not multivariate normal (sig.level = 0.05)"},"\n")
  ##  cat("---------------------------------------------", "\n\n", sep = " ")
    invisible(NULL)
})


## ------------------

setClass("r",
slots = c(R = "numeric", p.value="numeric", data.name="character"))

setGeneric("r", function(object) standardGeneric("r"))


setMethod("show",
signature = "r",
definition = function(object) {
    cat("            Royston test for Multivariate Normality", "\n\n", sep = " ")
  ##  cat("---------------------------------------------", "\n", sep = " ")
    cat("  data :", object@data.name, "\n\n", sep = " ")
    cat("  R               :", object@R, "\n", sep = " ")
    cat("  p-value         :", object@p.value, "\n\n", sep = " ")
    cat(if(object@p.value > 0.05){"  Result  : Data are multivariate normal (sig.level = 0.05)"}
        else {"  Result  : Data are not multivariate normal (sig.level = 0.05)"},"\n")
  ##  cat("---------------------------------------------", "\n\n", sep = " ")
    invisible(NULL)
})



## ------------------

setClass("hz",
slots = c(HZ = "numeric", p.value="numeric", data.name="character"))

setGeneric("hz", function(object) standardGeneric("hz"))


setMethod("show",
signature = "hz",
definition = function(object) {
    cat("            Henze-Zirkler test for Multivariate Normality", "\n\n", sep = " ")
  ##  cat("---------------------------------------------", "\n", sep = " ")
    cat("  data :", object@data.name, "\n\n", sep = " ")
    cat("  HZ              :", object@HZ, "\n", sep = " ")
    cat("  p-value         :", object@p.value, "\n\n", sep = " ")
    cat(if(object@p.value > 0.05){"  Result  : Data are multivariate normal (sig.level = 0.05)"}
        else {"  Result  : Data are not multivariate normal (sig.level = 0.05)"},"\n")
  ##  cat("---------------------------------------------", "\n\n", sep = " ")
    invisible(NULL)
})



