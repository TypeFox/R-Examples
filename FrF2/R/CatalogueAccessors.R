
"[.catlg" <- function(catlg,i){
       class(catlg) <- "list"
       aus <- catlg[i]
       if (!is.list(aus[1])) aus <- list(aus)
       class(aus) <- c("catlg","list")
       aus}
## same function for "[[.catlg" would destroy functionality of lapply

## define generics
res <- function(catlg) UseMethod("res")
nclear.2fis <- function(catlg) UseMethod("nclear.2fis")
clear.2fis <- function(catlg) UseMethod("clear.2fis")
#all.2fis.clear <- function(catlg) UseMethod("all.2fis.clear")
# not made generic because of conflict with generic all
nfac <- function(catlg) UseMethod("nfac")
nruns <- function(catlg) UseMethod("nruns")
WLP <- function(catlg) UseMethod("WLP")
dominating <- function(catlg) UseMethod("dominating")

## define defaults
#res.default <- function(catlg) stop(gettext("res works on class catlg objects only"))
#nclear.2fis.default <- function(catlg) stop(gettext("nclear.2fis works on class catlg objects only"))
#clear.2fis.default <- function(catlg) stop(gettext("clear.2fis works on class catlg objects only"))
#all.2fis.clear.default <- function(catlg) stop(gettext("all.2fis.clear works on class catlg objects only"))
#nfac.default <- function(catlg) stop(gettext("nfac works on class catlg objects only"))
#nruns.default <- function(catlg) stop(gettext("nruns works on class catlg objects only"))
#WLP.default <- function(catlg) stop(gettext("WLP works on class catlg objects only"))
#dominating.default <- function(catlg) stop(gettext("dominating works on class catlg objects only"))

res.catlg <- function(catlg) sapply(catlg, function(obj) obj$res)
nclear.2fis.catlg <- function(catlg) sapply(catlg, function(obj) obj$nclear.2fis)
clear.2fis.catlg <- function(catlg) lapply(catlg, function(obj) obj$clear.2fis)
all.2fis.clear.catlg <- function(catlg) lapply(catlg, function(obj) obj$all.2fis.clear)
nfac.catlg <- function(catlg) sapply(catlg, function(obj) obj$nfac)
nruns.catlg <- function(catlg) sapply(catlg, function(obj) obj$nruns)
WLP.catlg <- function(catlg) lapply(catlg, function(obj) obj$WLP)
dominating.catlg <- function(catlg) sapply(catlg, function(obj) {
      hilf <- obj$dominating
      if (is.null(hilf)) TRUE else hilf
      })

