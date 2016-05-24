## S4 object definitions and assigment/accessor functions for the slots.
##
## created  10.09.03 alexandros karatzoglou
## updated  23.08.05 


setClass("kernel",representation("function",kpar="list"))
setClass("kernelMatrix",representation("matrix"),prototype=structure(.Data=matrix()))

setClassUnion("listI", c("list","numeric","vector","integer","matrix"))
setClassUnion("output", c("matrix","factor","vector","logical","numeric","list","integer","NULL"))
setClassUnion("input", c("matrix","list"))
setClassUnion("kfunction", c("function","character"))
setClassUnion("mpinput", c("matrix","data.frame","missing"))
setClassUnion("lpinput", c("list","missing"))
setClassUnion("kpinput", c("kernelMatrix","missing"))



setClass("vm", representation(alpha = "listI", ## since setClassUnion is not  working
                              type = "character",
                              kernelf = "kfunction",
                              kpar = "list",
                              xmatrix = "input",
                              ymatrix = "output",
                              fitted = "output",  
                              lev = "vector",
                              nclass = "numeric",
                              error = "vector",
                              cross = "vector",
                              n.action= "ANY",
                              terms = "ANY",
                              kcall = "call"), contains= "VIRTUAL")
                                        #Generic Vector Machine object 


if(!isGeneric("type")){
  if (is.function("type"))
    fun <- type
  else fun <- function(object) standardGeneric("type")
  setGeneric("type", fun)
}
setMethod("type", "vm", function(object) object@type)
setGeneric("type<-", function(x, value) standardGeneric("type<-"))
setReplaceMethod("type", "vm", function(x, value) {
  x@type <- value
  x
})


if(!isGeneric("kernelf")){
  if (is.function("kernelf"))
    fun <- kernelf
  else fun <- function(object) standardGeneric("kernelf")
  setGeneric("kernelf", fun)
}
setMethod("kernelf", "vm", function(object) object@kernelf)
setGeneric("kernelf<-", function(x, value) standardGeneric("kernelf<-"))
setReplaceMethod("kernelf", "vm", function(x, value) {
  x@kernelf <- value
  x
})

if(!isGeneric("kpar")){
  if (is.function("kpar"))
    fun <- kpar
  else fun <- function(object) standardGeneric("kpar")
  setGeneric("kpar", fun)
}
setMethod("kpar", "vm", function(object) object@kpar)
setGeneric("kpar<-", function(x, value) standardGeneric("kpar<-"))
setReplaceMethod("kpar", "vm", function(x, value) {
  x@kpar <- value
  x
})

if(!isGeneric("kcall")){
  if (is.function("kcall"))
    fun <- kcall
  else fun <- function(object) standardGeneric("kcall")
  setGeneric("kcall", fun)
}
setMethod("kcall", "vm", function(object) object@kcall)
setGeneric("kcall<-", function(x, value) standardGeneric("kcall<-"))
setReplaceMethod("kcall", "vm", function(x, value) {
  x@kcall <- value
  x
})

setMethod("terms", "vm", function(x, ...) x@terms)
setGeneric("terms<-", function(x, value) standardGeneric("terms<-"))
setReplaceMethod("terms", "vm", function(x, value) {
  x@terms <- value
  x
})



if(!isGeneric("xmatrix")){
  if (is.function("xmatrix"))
    fun <- xmatrix
  else fun <- function(object) standardGeneric("xmatrix")
  setGeneric("xmatrix", fun)
}
setMethod("xmatrix", "vm", function(object) object@xmatrix)
setGeneric("xmatrix<-", function(x, value) standardGeneric("xmatrix<-"))
setReplaceMethod("xmatrix", "vm", function(x, value) {
  x@xmatrix <- value
  x
})

if(!isGeneric("ymatrix")){
  if (is.function("ymatrix"))
    fun <- ymatrix
  else fun <- function(object) standardGeneric("ymatrix")
  setGeneric("ymatrix", fun)
}
setMethod("ymatrix", "vm", function(object) object@ymatrix)
setGeneric("ymatrix<-", function(x, value) standardGeneric("ymatrix<-"))
setReplaceMethod("ymatrix", "vm", function(x, value) {
  x@ymatrix <- value
  x
})

setMethod("fitted", "vm", function(object, ...) object@fitted)
setGeneric("fitted<-", function(x, value) standardGeneric("fitted<-"))
setReplaceMethod("fitted", "vm", function(x, value) {
  x@fitted <- value
  x
})

if(!isGeneric("lev")){
  if (is.function("lev"))
    fun <- lev
  else fun <- function(object) standardGeneric("lev")
  setGeneric("lev", fun)
}
setMethod("lev", "vm", function(object) object@lev)
setGeneric("lev<-", function(x, value) standardGeneric("lev<-"))
setReplaceMethod("lev", "vm", function(x, value) {
  x@lev <- value
  x
})

if(!isGeneric("nclass")){
  if (is.function("nclass"))
    fun <- nclass
  else fun <- function(object) standardGeneric("nclass")
  setGeneric("nclass", fun)
}
setMethod("nclass", "vm", function(object) object@nclass)
setGeneric("nclass<-", function(x, value) standardGeneric("nclass<-"))
setReplaceMethod("nclass", "vm", function(x, value) {
  x@nclass <- value
  x
})

if(!isGeneric("alpha")){
  if (is.function("alpha"))
    fun <- alpha
  else fun <- function(object) standardGeneric("alpha")
  setGeneric("alpha", fun)
}
setMethod("alpha", "vm", function(object) object@alpha)
setGeneric("alpha<-", function(x, value) standardGeneric("alpha<-"))
setReplaceMethod("alpha", "vm", function(x, value) {
  x@alpha <- value
  x
})

if(!isGeneric("error")){
  if (is.function("error"))
    fun <- error
  else fun <- function(object) standardGeneric("error")
  setGeneric("error", fun)
}
setMethod("error", "vm", function(object) object@error)
setGeneric("error<-", function(x, value) standardGeneric("error<-"))
setReplaceMethod("error", "vm", function(x, value) {
  x@error <- value
  x
})

if(!isGeneric("cross")){
  if (is.function("cross"))
    fun <- cross
  else fun <- function(object) standardGeneric("cross")
  setGeneric("cross", fun)
}
setMethod("cross", "vm", function(object) object@cross)
setGeneric("cross<-", function(x, value) standardGeneric("cross<-"))
setReplaceMethod("cross", "vm", function(x, value) {
  x@cross <- value
  x
})

if(!isGeneric("n.action")){
  if (is.function("n.action"))
    fun <- n.action
  else fun <- function(object) standardGeneric("n.action")
  setGeneric("n.action", fun)
}
setMethod("n.action", "vm", function(object) object@n.action)
setGeneric("n.action<-", function(x, value) standardGeneric("n.action<-"))
setReplaceMethod("n.action", "vm", function(x, value) {
  x@n.action <- value
  x
})




setClass("ksvm", representation(param = "list",
                                scaling = "ANY",
                                coef = "ANY",
                                alphaindex = "ANY",
                                b = "numeric",
				obj = "vector",
                                SVindex = "vector",
                                nSV = "numeric",
                                prior = "list",
                                prob.model = "list"
                                ), contains="vm")

if(!isGeneric("param")){
  if (is.function("param"))
    fun <- param
  else fun <- function(object) standardGeneric("param")
  setGeneric("param", fun)
}
setMethod("param", "ksvm", function(object) object@param)
setGeneric("param<-", function(x, value) standardGeneric("param<-"))
setReplaceMethod("param", "ksvm", function(x, value) {
  x@param <- value
  x
})

if(!isGeneric("scaling")){
  if (is.function("scaling"))
    fun <- scaling
  else fun <- function(object) standardGeneric("scaling")
  setGeneric("scaling", fun)
}
setMethod("scaling", "ksvm", function(object) object@scaling)
setGeneric("scaling<-", function(x, value) standardGeneric("scaling<-"))
setReplaceMethod("scaling", "ksvm", function(x, value) {
  x@scaling<- value
  x
})

if(!isGeneric("obj")){
  if (is.function("obj"))
     fun <- obj
  else fun <- function(object) standardGeneric("obj")
  setGeneric("obj", fun)
}
setMethod("obj", "ksvm", function(object) object@obj)
setGeneric("obj<-", function(x, value) standardGeneric("obj<-"))
setReplaceMethod("obj", "ksvm", function(x, value) {
   x@obj<- value
   x
})


setMethod("coef", "ksvm", function(object, ...) object@coef)
setGeneric("coef<-", function(x, value) standardGeneric("coef<-"))
setReplaceMethod("coef", "ksvm", function(x, value) {
  x@coef <- value
  x
})

if(!isGeneric("alphaindex")){
  if (is.function("alphaindex"))
    fun <- alphaindex
  else fun <- function(object) standardGeneric("alphaindex")
  setGeneric("alphaindex", fun)
}
setMethod("alphaindex", "ksvm", function(object) object@alphaindex)
setGeneric("alphaindex<-", function(x, value) standardGeneric("alphaindex<-"))
setReplaceMethod("alphaindex", "ksvm", function(x, value) {
  x@alphaindex <- value
  x
})

if(!isGeneric("b")){
  if (is.function("b"))
    fun <- b
  else fun <- function(object) standardGeneric("b")
  setGeneric("b", fun)
}
setMethod("b", "ksvm", function(object) object@b)
setGeneric("b<-", function(x, value) standardGeneric("b<-"))
setReplaceMethod("b", "ksvm", function(x, value) {
  x@b <- value
  x
})

if(!isGeneric("SVindex")){
  if (is.function("SVindex"))
    fun <- SVindex
  else fun <- function(object) standardGeneric("SVindex")
  setGeneric("SVindex", fun)
}
setMethod("SVindex", "ksvm", function(object) object@SVindex)
setGeneric("SVindex<-", function(x, value) standardGeneric("SVindex<-"))
setReplaceMethod("SVindex", "ksvm", function(x, value) {
  x@SVindex <- value
  x
})

if(!isGeneric("nSV")){
  if (is.function("nSV"))
    fun <- nSV
  else fun <- function(object) standardGeneric("nSV")
  setGeneric("nSV", fun)
}
setMethod("nSV", "ksvm", function(object) object@nSV)
setGeneric("nSV<-", function(x, value) standardGeneric("nSV<-"))
setReplaceMethod("nSV", "ksvm", function(x, value) {
  x@nSV <- value
  x
})

if(!isGeneric("prior")){
  if (is.function("prior"))
    fun <- prior
  else fun <- function(object) standardGeneric("prior")
  setGeneric("prior", fun)
}
setMethod("prior", "ksvm", function(object) object@prior)
setGeneric("prior<-", function(x, value) standardGeneric("prior<-"))
setReplaceMethod("prior", "ksvm", function(x, value) {
  x@prior <- value
  x
})

if(!isGeneric("prob.model")){
  if (is.function("prob.model"))
    fun <- prob.model
  else fun <- function(object) standardGeneric("prob.model")
  setGeneric("prob.model", fun)
}
setMethod("prob.model", "ksvm", function(object) object@prob.model)
setGeneric("prob.model<-", function(x, value) standardGeneric("prob.model<-"))
setReplaceMethod("prob.model", "ksvm", function(x, value) {
  x@prob.model <- value
  x
})


setClass("lssvm", representation(param = "list",
                                scaling = "ANY",
                                coef = "ANY",
                                alphaindex = "ANY",
                             ##    prob.model = "list",
                                b = "numeric",
                                 nSV = "numeric"
                                 ), contains="vm")



##setMethod("prob.model", "lssvm", function(object) object@prob.model)
##setGeneric("prob.model<-", function(x, value) standardGeneric("prob.model<-"))
##setReplaceMethod("prob.model", "lssvm", function(x, value) {
##  x@prob.model <- value
##  x
##})

setMethod("param", "lssvm", function(object) object@param)
setReplaceMethod("param", "lssvm", function(x, value) {
  x@param <- value
  x
})

setMethod("scaling", "lssvm", function(object) object@scaling)
setReplaceMethod("scaling", "lssvm", function(x, value) {
  x@scaling<- value
  x
})

setMethod("coef", "lssvm", function(object, ...) object@coef)
setReplaceMethod("coef", "lssvm", function(x, value) {
  x@coef <- value
  x
})

setMethod("alphaindex", "lssvm", function(object) object@alphaindex)
setReplaceMethod("alphaindex", "lssvm", function(x, value) {
  x@alphaindex <- value
  x
})

setMethod("b", "lssvm", function(object) object@b)
setReplaceMethod("b", "lssvm", function(x, value) {
  x@b <- value
  x
})

setMethod("nSV", "lssvm", function(object) object@nSV)
setReplaceMethod("nSV", "lssvm", function(x, value) {
  x@nSV <- value
  x
})

setClass("kqr", representation(param = "list",
                               scaling = "ANY",
                               coef = "ANY",
                               b = "numeric"
                               ), contains="vm")


setMethod("b", "kqr", function(object) object@b)
setReplaceMethod("b", "kqr", function(x, value) {
  x@b <- value
  x
})

setMethod("scaling", "kqr", function(object) object@scaling)
setReplaceMethod("scaling", "kqr", function(x, value) {
  x@scaling <- value
  x
})

setMethod("coef", "kqr", function(object) object@coef)
setReplaceMethod("coef", "kqr", function(x, value) {
  x@coef <- value
  x
})

setMethod("param", "kqr", function(object) object@param)
setReplaceMethod("param", "kqr", function(x, value) {
  x@param <- value
  x
})

## failed attempt to get rid of all this above


## mkaccesfun <- function(cls)
#{
#  snames <- slotNames(cls)
## 
#
#  for(i in 1:length(snames))
#    { resF <- paste("\"",snames[i],"\"",sep="")
#      if(!isGeneric(snames[i]))
#        eval(parse(file="",text=paste("setGeneric(",resF,",function(object)","standardGeneric(",resF,")",")",sep=" ")))   
#    setGeneric(snames[i], function(object) standardGeneric(snames[i]))
# 
#  setMethod(snames[i], cls, function(object)  eval(parse(file="",text=paste("object@",snames[i],sep=""))))
#  resG  <-  paste("\"",snames[i],"<-","\"",sep="")
#eval(parse(file="",text=paste("setGeneric(",resG,",function(x, value)","standardGeneric(",resG,")",")",sep=" ")))
#  setReplaceMethod(snames[i], cls, function(x, value) {
#    eval(parse(file="",text=paste("x@",snames[i],"<-value",sep="")))
#    x
#  })                   
#  }
#}


setClass("prc", representation(pcv = "matrix",
                               eig = "vector",
                               kernelf = "kfunction",
                               kpar = "list",
                               xmatrix = "input",
                               kcall = "ANY",
                               terms = "ANY",
                               n.action = "ANY"),contains="VIRTUAL")
#accessor functions 
if(!isGeneric("pcv")){
  if (is.function("pcv"))
    fun <- pcv
  else fun <- function(object) standardGeneric("pcv")
  setGeneric("pcv", fun)
}
setMethod("pcv", "prc", function(object) object@pcv)
setGeneric("pcv<-", function(x, value) standardGeneric("pcv<-"))
setReplaceMethod("pcv", "prc", function(x, value) {
  x@pcv <- value
  x
})

if(!isGeneric("eig")){
  if (is.function("eig"))
    fun <- eig
  else fun <- function(object) standardGeneric("eig")
  setGeneric("eig", fun)
}
setMethod("eig", "prc", function(object) object@eig)
setGeneric("eig<-", function(x, value) standardGeneric("eig<-"))
setReplaceMethod("eig", "prc", function(x, value) {
  x@eig <- value
  x
})



setMethod("kernelf","prc", function(object) object@kernelf)
setReplaceMethod("kernelf","prc", function(x, value){
  x@kernelf <- value
  x
})

setMethod("xmatrix","prc", function(object) object@xmatrix)
setReplaceMethod("xmatrix","prc", function(x, value){
  x@xmatrix <- value
  x
})

setMethod("kcall","prc", function(object) object@kcall)
setReplaceMethod("kcall","prc", function(x, value){
  x@kcall <- value
  x
})

setMethod("terms","prc", function(x, ...) x@terms)
setReplaceMethod("terms","prc", function(x, value){
  x@terms <- value
  x
})

setMethod("n.action","prc", function(object) object@n.action)
setReplaceMethod("n.action","prc", function(x, value){
  x@n.action <- value
  x
})




##kernel principal components object
setClass("kpca", representation(rotated = "matrix"),contains="prc")
#accessor functions 

if(!isGeneric("rotated")){
  if (is.function("rotated"))
    fun <- rotated
  else fun <- function(object) standardGeneric("rotated")
  setGeneric("rotated", fun)
}
setMethod("rotated", "kpca", function(object) object@rotated)
setGeneric("rotated<-", function(x, value) standardGeneric("rotated<-"))
setReplaceMethod("rotated", "kpca", function(x, value) {
  x@rotated <- value
  x
})

## kernel maximum mean discrepancy

setClass("kmmd", representation(H0="logical", AsympH0 ="logical", kernelf = "kfunction", Asymbound="numeric", Radbound="numeric", xmatrix="input", mmdstats="vector"))

if(!isGeneric("mmdstats")){
  if (is.function("mmdstats"))
    fun <- mmdstats
  else fun <- function(object) standardGeneric("mmdstats")
  setGeneric("mmdstats", fun)
}
setMethod("mmdstats","kmmd", function(object) object@mmdstats)
setGeneric("mmdstats<-", function(x, value) standardGeneric("mmdstats<-"))
setReplaceMethod("mmdstats","kmmd", function(x, value){
  x@mmdstats <- value
  x
})


if(!isGeneric("Radbound")){
  if (is.function("Radbound"))
    fun <- Radbound
  else fun <- function(object) standardGeneric("Radbound")
  setGeneric("Radbound", fun)
}

setMethod("Radbound","kmmd", function(object) object@Radbound)
setGeneric("Radbound<-", function(x, value) standardGeneric("Radbound<-"))
setReplaceMethod("Radbound","kmmd", function(x, value){
  x@Radbound <- value
  x
})


if(!isGeneric("Asymbound")){
  if (is.function("Asymbound"))
    fun <- Asymbound
  else fun <- function(object) standardGeneric("Asymbound")
  setGeneric("Asymbound", fun)
}
setMethod("Asymbound","kmmd", function(object) object@Asymbound)
setGeneric("Asymbound<-", function(x, value) standardGeneric("Asymbound<-"))
setReplaceMethod("Asymbound","kmmd", function(x, value){
  x@Asymbound <- value
  x
})

if(!isGeneric("H0")){
  if (is.function("H0"))
    fun <- H0
  else fun <- function(object) standardGeneric("H0")
  setGeneric("H0", fun)
}
setMethod("H0","kmmd", function(object) object@H0)
setGeneric("H0<-", function(x, value) standardGeneric("H0<-"))
setReplaceMethod("H0","kmmd", function(x, value){
  x@H0 <- value
  x
})


if(!isGeneric("AsympH0")){
  if (is.function("AsympH0"))
    fun <- AsympH0
  else fun <- function(object) standardGeneric("AsympH0")
  setGeneric("AsympH0", fun)
}
setMethod("AsympH0","kmmd", function(object) object@AsympH0)
setGeneric("AsympH0<-", function(x, value) standardGeneric("AsympH0<-"))
setReplaceMethod("AsympH0","kmmd", function(x, value){
  x@AsympH0 <- value
  x
})

setMethod("kernelf","kmmd", function(object) object@kernelf)
setReplaceMethod("kernelf","kmmd", function(x, value){
  x@kernelf <- value
  x
})




setClass("ipop", representation(primal = "vector",
                                dual = "numeric",
                                how = "character"
                               ))

if(!isGeneric("primal")){
  if (is.function("primal"))
    fun <- primal
  else fun <- function(object) standardGeneric("primal")
  setGeneric("primal", fun)
}
setMethod("primal", "ipop", function(object) object@primal)
setGeneric("primal<-", function(x, value) standardGeneric("primal<-"))
setReplaceMethod("primal", "ipop", function(x, value) {
  x@primal <- value
  x
})

if(!isGeneric("dual")){
  if (is.function("dual"))
    fun <- dual
  else fun <- function(object) standardGeneric("dual")
  setGeneric("dual", fun)
}
setMethod("dual", "ipop", function(object) object@dual)
setGeneric("dual<-", function(x, value) standardGeneric("dual<-"))
setReplaceMethod("dual", "ipop", function(x, value) {
  x@dual <- value
  x
})

if(!isGeneric("how")){
  if (is.function("how"))
    fun <- how
  else fun <- function(object) standardGeneric("how")
  setGeneric("how", fun)
}
setMethod("how", "ipop", function(object) object@how)
setGeneric("how<-", function(x, value) standardGeneric("how<-"))
setReplaceMethod("how", "ipop", function(x, value) {
  x@how <- value
  x
})

# Kernel Canonical Correlation Analysis
setClass("kcca", representation(kcor = "vector",
                                xcoef = "matrix",
                                ycoef = "matrix"
                                ##xvar = "matrix",
                                ##yvar = "matrix"
                                ))


if(!isGeneric("kcor")){
  if (is.function("kcor"))
    fun <- kcor
  else fun <- function(object) standardGeneric("kcor")
  setGeneric("kcor", fun)
}
setMethod("kcor", "kcca", function(object) object@kcor)
setGeneric("kcor<-", function(x, value) standardGeneric("kcor<-"))
setReplaceMethod("kcor", "kcca", function(x, value) {
  x@kcor <- value
  x
})

if(!isGeneric("xcoef")){
  if (is.function("xcoef"))
    fun <- xcoef
  else fun <- function(object) standardGeneric("xcoef")
  setGeneric("xcoef", fun)
}
setMethod("xcoef", "kcca", function(object) object@xcoef)
setGeneric("xcoef<-", function(x, value) standardGeneric("xcoef<-"))
setReplaceMethod("xcoef", "kcca", function(x, value) {
  x@xcoef <- value
  x
})

if(!isGeneric("ycoef")){
  if (is.function("ycoef"))
    fun <- ycoef
  else fun <- function(object) standardGeneric("ycoef")
  setGeneric("ycoef", fun)
}
setMethod("ycoef", "kcca", function(object) object@ycoef)
setGeneric("ycoef<-", function(x, value) standardGeneric("ycoef<-"))
setReplaceMethod("ycoef", "kcca", function(x, value) {
  x@ycoef <- value
  x
})

##if(!isGeneric("xvar")){
##  if (is.function("xvar"))
##    fun <- xvar
##  else fun <- function(object) standardGeneric("xvar")
##  setGeneric("xvar", fun)
##}
##setMethod("xvar", "kcca", function(object) object@xvar)
##setGeneric("xvar<-", function(x, value) standardGeneric("xvar<-"))
##setReplaceMethod("xvar", "kcca", function(x, value) {
##  x@xvar <- value
##  x
##})

##if(!isGeneric("yvar")){
##  if (is.function("yvar"))
##    fun <- yvar
##  else fun <- function(object) standardGeneric("yvar")
##  setGeneric("yvar", fun)
##}
##setMethod("yvar", "kcca", function(object) object@yvar)
##setGeneric("yvar<-", function(x, value) standardGeneric("yvar<-"))
##setReplaceMethod("yvar", "kcca", function(x, value) {
##  x@yvar <- value
##  x
##})

## Gaussian Processes object
setClass("gausspr",representation(tol = "numeric",
                                  scaling = "ANY",
                                  sol = "matrix",
                                  alphaindex="list",
                                  nvar = "numeric"
                                  ),contains="vm")



setMethod("alphaindex","gausspr", function(object) object@alphaindex)
setReplaceMethod("alphaindex","gausspr", function(x, value){
  x@alphaindex <- value
  x
})

if(!isGeneric("sol")){
  if (is.function("sol"))
    fun <- sol
  else fun <- function(object) standardGeneric("sol")
  setGeneric("sol", fun)
}
setMethod("sol","gausspr", function(object) object@sol)
setGeneric("sol<-", function(x, value) standardGeneric("sol<-"))
setReplaceMethod("sol","gausspr", function(x, value){
  x@sol <- value
  x
})


setMethod("scaling","gausspr", function(object) object@scaling)
setReplaceMethod("scaling","gausspr", function(x, value){
  x@scaling <- value
  x
})


setMethod("coef", "gausspr", function(object, ...) object@alpha)


# Relevance Vector Machine object 
setClass("rvm", representation(tol = "numeric",
                               nvar = "numeric",
                               mlike = "numeric",
                               RVindex = "vector",
                               coef = "ANY",
                               nRV = "numeric"),contains ="vm")


if(!isGeneric("tol")){
  if (is.function("tol"))
    fun <- tol
  else fun <- function(object) standardGeneric("tol")
  setGeneric("tol", fun)
}
setMethod("tol", "rvm", function(object) object@tol)
setGeneric("tol<-", function(x, value) standardGeneric("tol<-"))
setReplaceMethod("tol", "rvm", function(x, value) {
  x@tol <- value
  x
})


setMethod("coef", "rvm", function(object, ...) object@coef)
setReplaceMethod("coef", "rvm", function(x, value) {
  x@coef <- value
  x
})

if(!isGeneric("RVindex")){
  if (is.function("RVindex"))
    fun <- RVindex
  else fun <- function(object) standardGeneric("RVindex")
  setGeneric("RVindex", fun)
}
setMethod("RVindex", "rvm", function(object) object@RVindex)
setGeneric("RVindex<-", function(x, value) standardGeneric("RVindex<-"))
setReplaceMethod("RVindex", "rvm", function(x, value) {
  x@RVindex <- value
  x
})

if(!isGeneric("nvar")){
  if (is.function("nvar"))
    fun <- nvar
  else fun <- function(object) standardGeneric("nvar")
  setGeneric("nvar", fun)
}
setMethod("nvar", "rvm", function(object) object@nvar)
setGeneric("nvar<-", function(x, value) standardGeneric("nvar<-"))
setReplaceMethod("nvar", "rvm", function(x, value) {
  x@nvar <- value
  x
})

if(!isGeneric("nRV")){
  if (is.function("nRV"))
    fun <- nRV
  else fun <- function(object) standardGeneric("nRV")
  setGeneric("nRV", fun)
}
setMethod("nRV", "rvm", function(object) object@nRV)
setGeneric("nRV<-", function(x, value) standardGeneric("nRV<-"))
setReplaceMethod("nRV", "rvm", function(x, value) {
  x@nRV <- value
  x
})

setMethod("coef", "rvm", function(object, ...) object@alpha)

if(!isGeneric("mlike")){
  if (is.function("mlike"))
    fun <- mlike
  else fun <- function(object) standardGeneric("mlike")
  setGeneric("mlike", fun)
}
setMethod("mlike", "rvm", function(object) object@mlike)
setGeneric("mlike<-", function(x, value) standardGeneric("mlike<-"))
setReplaceMethod("mlike", "rvm", function(x, value) {
  x@mlike <- value
  x
})


setClass("inchol",representation("matrix",
                                 pivots="vector",
                                 diagresidues="vector",
                                 maxresiduals="vector"),
         prototype=structure(.Data=matrix(),
           pivots=vector(),
           diagresidues=vector(), 
           maxresiduals=vector()))


if(!isGeneric("pivots")){
if (is.function("pivots"))
  fun <- pivots
else fun <- function(object) standardGeneric("pivots")
setGeneric("pivots", fun)
}
setMethod("pivots", "inchol", function(object) object@pivots)
setGeneric("pivots<-", function(x, value) standardGeneric("pivots<-"))
setReplaceMethod("pivots", "inchol", function(x, value) {
  x@pivots <- value
  x
})

if(!isGeneric("diagresidues")){
if (is.function("diagresidues"))
  fun <- diagresidues
else fun <- function(object) standardGeneric("diagresidues")
setGeneric("diagresidues", fun)
}
setMethod("diagresidues", "inchol", function(object) object@diagresidues)
setGeneric("diagresidues<-", function(x,value) standardGeneric("diagresidues<-"))
setReplaceMethod("diagresidues", "inchol", function(x, value) {
  x@diagresidues <- value
  x
})

if(!isGeneric("maxresiduals")){
if (is.function("maxresiduals"))
  fun <- maxresiduals
else fun <- function(object) standardGeneric("maxresiduals")
setGeneric("maxresiduals", fun)
}
setMethod("maxresiduals", "inchol", function(object) object@maxresiduals)
setGeneric("maxresiduals<-", function(x,value) standardGeneric("maxresiduals<-"))
setReplaceMethod("maxresiduals", "inchol", function(x, value) {
  x@maxresiduals <- value
  x
})


## csi object
setClass("csi",representation(Q = "matrix",
                              R = "matrix",
                              truegain = "vector",
                              predgain = "vector"),contains="inchol")

if(!isGeneric("Q")){
if (is.function("Q"))
  fun <- Q
else fun <- function(object) standardGeneric("Q")
setGeneric("Q", fun)
}
setMethod("Q", "csi", function(object) object@Q)
setGeneric("Q<-", function(x, value) standardGeneric("Q<-"))
setReplaceMethod("Q", "csi", function(x, value) {
  x@Q <- value
  x
})

if(!isGeneric("R")){
if (is.function("R"))
  fun <- R
else fun <- function(object) standardGeneric("R")
setGeneric("R", fun)
}
setMethod("R", "csi", function(object) object@R)
setGeneric("R<-", function(x, value) standardGeneric("R<-"))
setReplaceMethod("R", "csi", function(x, value) {
  x@R <- value
  x
})

if(!isGeneric("truegain")){
if (is.function("truegain"))
  fun <- truegain
else fun <- function(object) standardGeneric("truegain")
setGeneric("truegain", fun)
}
setMethod("truegain", "csi", function(object) object@truegain)
setGeneric("truegain<-", function(x, value) standardGeneric("truegain<-"))
setReplaceMethod("truegain", "csi", function(x, value) {
  x@truegain <- value
  x
})

if(!isGeneric("predgain")){
if (is.function("predgain"))
  fun <- predgain
else fun <- function(object) standardGeneric("predgain")
setGeneric("predgain", fun)
}
setMethod("predgain", "csi", function(object) object@predgain)
setGeneric("predgain<-", function(x, value) standardGeneric("predgain<-"))
setReplaceMethod("predgain", "csi", function(x, value) {
  x@predgain <- value
  x
})


setClass("specc",representation("vector",
                                centers="matrix",
                                size="vector",
                                kernelf="kfunction",
                                withinss = "vector"
                                ),prototype=structure(.Data=vector(),
                                    centers = matrix(),
                                    size=matrix(),
                                    kernelf = ls,
                                    withinss=vector()))


if(!isGeneric("centers")){
if (is.function("centers"))
  fun <- centers
else fun <- function(object) standardGeneric("centers")
setGeneric("centers", fun)
}
setMethod("centers", "specc", function(object) object@centers)
setGeneric("centers<-", function(x,value) standardGeneric("centers<-"))
setReplaceMethod("centers", "specc", function(x, value) {
  x@centers <- value
  x
})

if(!isGeneric("size")){
if (is.function("size"))
  fun <- size
else fun <- function(object) standardGeneric("size")
setGeneric("size", fun)
}
setMethod("size", "specc", function(object) object@size)
setGeneric("size<-", function(x,value) standardGeneric("size<-"))
setReplaceMethod("size", "specc", function(x, value) {
  x@size <- value
  x
})

if(!isGeneric("withinss")){
if (is.function("withinss"))
  fun <- withinss
else fun <- function(object) standardGeneric("withinss")
setGeneric("withinss", fun)
}
setMethod("withinss", "specc", function(object) object@withinss)
setGeneric("withinss<-", function(x,value) standardGeneric("withinss<-"))
setReplaceMethod("withinss", "specc", function(x, value) {
  x@withinss <- value
  x
})

setMethod("kernelf","specc", function(object) object@kernelf)
setReplaceMethod("kernelf","specc", function(x, value){
  x@kernelf <- value
  x
})



setClass("ranking",representation("matrix",
 				    convergence="matrix",
				    edgegraph="matrix"),
				    prototype=structure(.Data=matrix(),
                                    convergence=matrix(),
				    edgegraph=matrix()))

if(!isGeneric("convergence")){
if (is.function("convergence"))
  fun <- convergence
else fun <- function(object) standardGeneric("convergence")
setGeneric("convergence", fun)
}
setMethod("convergence", "ranking", function(object) object@convergence)
setGeneric("convergence<-", function(x,value) standardGeneric("convergence<-"))
setReplaceMethod("convergence", "ranking", function(x, value) {
  x@convergence <- value
  x
})

if(!isGeneric("edgegraph")){
if (is.function("edgegraph"))
  fun <- edgegraph
else fun <- function(object) standardGeneric("edgegraph")
setGeneric("edgegraph", fun)
}
setMethod("edgegraph", "ranking", function(object) object@edgegraph)
setGeneric("edgegraph<-", function(x,value) standardGeneric("edgegraph<-"))
setReplaceMethod("edgegraph", "ranking", function(x, value) {
  x@edgegraph <- value
  x
})

## online learning algorithms class

setClass("onlearn", representation(
                                kernelf = "kfunction",
                                buffer = "numeric",
                                kpar = "list",
                                xmatrix = "matrix",
                                fit = "numeric",
                                onstart = "numeric",
                                onstop = "numeric",
                                alpha = "ANY",
                                rho = "numeric",
                                b = "numeric",
                                pattern ="ANY",
                                type="character"
                               ))


if(!isGeneric("fit")){
  if (is.function("fit"))
    fun <- fit
  else fun <- function(object) standardGeneric("fit")
  setGeneric("fit", fun)
}
setMethod("fit","onlearn", function(object) object@fit)
setGeneric("fit<-", function(x, value) standardGeneric("fit<-"))
setReplaceMethod("fit","onlearn", function(x, value){
  x@fit <- value
  x
})

if(!isGeneric("onstart")){
  if (is.function("onstart"))
    fun <- onstart
  else fun <- function(object) standardGeneric("onstart")
  setGeneric("onstart", fun)
}
setMethod("onstart", "onlearn", function(object) object@onstart)
setGeneric("onstart<-", function(x, value) standardGeneric("onstart<-"))
setReplaceMethod("onstart", "onlearn", function(x, value) {
  x@onstart <- value
  x
})

if(!isGeneric("onstop")){
  if (is.function("onstop"))
    fun <- onstop
  else fun <- function(object) standardGeneric("onstop")
  setGeneric("onstop", fun)
}
setMethod("onstop", "onlearn", function(object) object@onstop)
setGeneric("onstop<-", function(x, value) standardGeneric("onstop<-"))
setReplaceMethod("onstop", "onlearn", function(x, value) {
  x@onstop <- value
  x
})

if(!isGeneric("buffer")){
  if (is.function("buffer"))
    fun <- buffer
  else fun <- function(object) standardGeneric("buffer")
  setGeneric("buffer", fun)
}
setMethod("buffer", "onlearn", function(object) object@buffer)
setGeneric("buffer<-", function(x, value) standardGeneric("buffer<-"))
setReplaceMethod("buffer", "onlearn", function(x, value) {
  x@buffer <- value
  x
})

setMethod("kernelf","onlearn", function(object) object@kernelf)
setReplaceMethod("kernelf","onlearn", function(x, value){
  x@kernelf <- value
  x
})

setMethod("kpar","onlearn", function(object) object@kpar)
setReplaceMethod("kpar","onlearn", function(x, value){
  x@kpar <- value
  x
})

setMethod("xmatrix","onlearn", function(object) object@xmatrix)
setReplaceMethod("xmatrix","onlearn", function(x, value){
  x@xmatrix <- value
  x
})


setMethod("alpha","onlearn", function(object) object@alpha)
setReplaceMethod("alpha","onlearn", function(x, value){
  x@alpha <- value
  x
})

setMethod("b","onlearn", function(object) object@b)
setReplaceMethod("b","onlearn", function(x, value){
  x@b <- value
  x
})

setMethod("type","onlearn", function(object) object@type)
setReplaceMethod("type","onlearn", function(x, value){
  x@type <- value
  x
})

if(!isGeneric("rho")){
  if (is.function("rho"))
    fun <- rho
  else fun <- function(object) standardGeneric("rho")
  setGeneric("rho", fun)
}
setMethod("rho", "onlearn", function(object) object@rho)
setGeneric("rho<-", function(x, value) standardGeneric("rho<-"))
setReplaceMethod("rho", "onlearn", function(x, value) {
  x@rho <- value
  x
})

if(!isGeneric("pattern")){
  if (is.function("pattern"))
    fun <- pattern
  else fun <- function(object) standardGeneric("pattern")
  setGeneric("pattern", fun)
}
setMethod("pattern", "onlearn", function(object) object@pattern)
setGeneric("pattern<-", function(x, value) standardGeneric("pattern<-"))
setReplaceMethod("pattern", "onlearn", function(x, value) {
  x@pattern <- value
  x
})





setClass("kfa",representation(alpha = "matrix",
                              alphaindex = "vector",
                              kernelf = "kfunction",
                              xmatrix = "matrix",
                              kcall = "call",
                              terms = "ANY" )) 


setMethod("coef", "kfa", function(object, ...) object@alpha)

setMethod("kernelf","kfa", function(object) object@kernelf)
setReplaceMethod("kernelf","kfa", function(x, value){
  x@kernelf <- value
  x
})

setMethod("alphaindex","kfa", function(object) object@alphaindex)
setReplaceMethod("alphaindex","kfa", function(x, value){
  x@alphaindex <- value
  x
})

setMethod("alpha","kfa", function(object) object@alpha)
setReplaceMethod("alpha","kfa", function(x, value){
  x@alpha <- value
  x
})

setMethod("xmatrix","kfa", function(object) object@xmatrix)
setReplaceMethod("xmatrix","kfa", function(x, value){
  x@xmatrix <- value
  x
})


setMethod("kcall","kfa", function(object) object@kcall)
setReplaceMethod("kcall","kfa", function(x, value){
  x@kcall <- value
  x
})


setMethod("terms","kfa", function(x, ...) x@terms)
setReplaceMethod("terms","kfa", function(x, value){
  x@terms <- value
  x
})


## kernel hebbian algorithm object
setClass("kha", representation(eskm ="vector"),contains="prc")

## accessor functions 

if(!isGeneric("eskm")){
  if (is.function("eskm"))
    fun <- eskm
  else fun <- function(object) standardGeneric("eskm")
  setGeneric("eskm", fun)
}
setMethod("eskm", "kha", function(object) object@eskm)
setGeneric("eskm<-", function(x, value) standardGeneric("eskm<-"))
setReplaceMethod("eskm", "kha", function(x, value) {
  x@eskm <- value
  x
})

