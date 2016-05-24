setClass("CGHoptions",
         representation=list(
           select    = 'character',
           calling   = 'logical',
           wavenorm  = 'character',
           GCnorm    = 'character',
           nblevels  = 'numeric',
           alpha     = 'numeric',
           beta      = 'numeric',
           nbprocs   = 'numeric',
           itermax   = 'numeric'),
         prototype(select    = "mBIC",
                   calling   = FALSE,
                   wavenorm  = "none",
                   GCnorm    = "none",
                   nblevels  = 3,
                   alpha     = 0.2,
                   beta      = 0.5,
                   nbprocs   = 1,
                   itermax   = Inf)
         )

setMethod(
          f = "[",
          signature = "CGHoptions",
          definition = function(x,i,j,drop){
            if (i=="select")     {return(x@select)}      else {}
            if (i=="calling")    {return(x@calling)}     else {}
            if (i=="wavenorm")   {return(x@wavenorm)}    else {}
            if (i=="GCnorm")     {return(x@GCnorm)}      else {}
            if (i=="nblevels")   {return(x@nblevels)}    else {}
            if (i=="alpha")      {return(x@alpha)}       else {}
            if (i=="beta")       {return(x@beta)}        else {}
            if (i=="nbprocs")    {return(x@nbprocs)}     else {}
            if (i=="itermax")    {return(x@itermax)}     else {}
          }          
          )

setMethod(
          f = "print",
          signature = "CGHoptions",
          definition = function(x){
            cat("****** CGHoption print ******\n")
            A = data.frame(options = rep(NA,9), value = rep(NA,9),row.names=NULL)
            A[1,1]  = "select"
            A[1,2]  = x@select
            A[2,1]  = "calling"
            A[2,2]  = x@calling
            A[3,1]  = "wavenorm"
            A[3,2]  = x@wavenorm
            A[4,1]  = "GCnorm"
            A[4,2]  = x@GCnorm
            A[5,1]  = "nblevels"
            A[5,2]  = x@nblevels
            A[6,1]  = "alpha"
            A[6,2]  = x@alpha
            A[7,1]  = "beta"
            A[7,2]  = x@beta
            A[8,1]  = "nbprocs"
            A[8,2]  = x@nbprocs            
            A[9,1]  = "itermax"
            A[9,2]  = x@itermax
            print(A)            
          }          
          )

setMethod(
          f = "show",
          signature = "CGHoptions",
          definition = function(object){
            cat("****** CGHoption show ******\n")
            A = data.frame(options = rep(NA,9), value = rep(NA,9))
            A[1,1]  = "select"
            A[1,2]  = object@select
            A[2,1]  = "calling"
            A[2,2]  = object@calling
            A[3,1]  = "wavenorm"
            A[3,2]  = object@wavenorm
            A[4,1]  = "GCnorm"
            A[4,2]  = object@GCnorm
            A[5,1]  = "nblevels"
            A[5,2]  = object@nblevels
            A[6,1]  = "alpha"
            A[6,2]  = object@alpha
            A[7,1]  = "beta"
            A[7,2]  = object@beta
            A[8,1]  = "nbprocs"
            A[8,2]  = object@nbprocs
            A[9,1]  = "itermax"
            A[9,2]  = object@itermax
            print(A)            
          }          
          )

setGeneric("select<-",function(object,value){standardGeneric("select<-")})
setReplaceMethod(
                 f="select",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@select =value
                   return (object)
                 }
                 )
setGeneric("calling<-",function(object,value){standardGeneric("calling<-")})
setReplaceMethod(
                 f="calling",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@calling =value
                   return (object)
                 }
                 )
setGeneric("wavenorm<-",function(object,value){standardGeneric("wavenorm<-")})
setReplaceMethod(
                 f="wavenorm",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@wavenorm =value
                   return (object)
                 }
                 )
setGeneric("GCnorm<-",function(object,value){standardGeneric("GCnorm<-")})
setReplaceMethod(
                 f="GCnorm",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@GCnorm =value
                   return (object)
                 }
                 )

setGeneric("nblevels<-",function(object,value){standardGeneric("nblevels<-")})
setReplaceMethod(
                 f="nblevels",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@nblevels =value
                   return (object)
                 }
                 )
setGeneric("alpha<-",function(object,value){standardGeneric("alpha<-")})
setReplaceMethod(
                 f="alpha",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@alpha =value
                   return (object)
                 }
                 )
setGeneric("beta<-",function(object,value){standardGeneric("beta<-")})
setReplaceMethod(
                 f="beta",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@beta =value
                   return (object)
                 }
                 )
setGeneric("nbprocs<-",function(object,value){standardGeneric("nbprocs<-")})
setReplaceMethod(
                 f="nbprocs",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@nbprocs =value
                   return (object)
                 }
                 )
setGeneric("itermax<-",function(object,value){standardGeneric("itermax<-")})
setReplaceMethod(
                 f="itermax",
                 signature="CGHoptions",
                 definition=function(object,value){
                   object@itermax =value
                   return (object)
                 }
                 )

setGeneric( name = "fun2run"   ,def = function(.Object){standardGeneric("fun2run")})

setMethod(
          f = "fun2run",
          signature = "CGHoptions",
          definition = function(.Object){
            if        (.Object@calling==TRUE & (.Object@wavenorm!="none" | .Object@GCnorm!="none") ){
              instr = parse(text = "Res = ILSclust(.Object,CGHo,uniKmax,multiKmax)")
              cat("[multiseg] ILSclust running \n")
            } else if (.Object@calling==TRUE & .Object@wavenorm=="none" & .Object@GCnorm=="none" ){
              instr = parse(text = "Res = multisegclust(.Object,CGHo,uniKmax,multiKmax)")
              cat("[multiseg] multisegclust running \n")    
            } else if (.Object@calling==FALSE & .Object@wavenorm=="none" & .Object@GCnorm=="none"){
              instr = parse(text = "Res = multisegmean(.Object,CGHo,uniKmax,multiKmax)")
              cat("[multiseg] multisegmean running \n")
            } else if (.Object@calling==FALSE & (.Object@wavenorm!="none" | .Object@GCnorm!="none") ){
              instr = parse(text = "Res = ILS(.Object,CGHo,uniKmax,multiKmax)")
              cat("[multiseg] ILS running \n")
            }
            return(instr)
          }          
          )



