setClassUnion("numOrNULL", c("numeric","NULL"))
setClassUnion("factorOrNULL", c("factor","NULL"))


setClass(Class = "CGHdata",representation(Y = "list",GCcontent="numOrNULL",genomic.position="numOrNULL",probeID="factorOrNULL"))

setMethod(
          f          = "initialize",
          signature  = "CGHdata",
          definition =          
          function(.Object,Y) {
            
            Y = data.frame(Y)            
            patient.status = ""

            if ("GCcontent" %in% colnames(Y)){
              j = which(colnames(Y)=="GCcontent")
              .Object@GCcontent = Y[,j]
              Y = Y[,-j]
            } else {
              .Object@GCcontent=NULL
            }
            if ("genomic.position" %in% colnames(Y)){
              if (is.unsorted(Y$genomic.position)){
                cat("[raw data: error] sort genomic.positions \n")
                stop()
              } else {              
                j = which(colnames(Y)=="genomic.position")
                .Object@genomic.position = Y[,j]
                Y = Y[,-j]
              }
            } else {
              .Object@genomic.position = NULL
            }
            
            if ("probeID" %in% colnames(Y)){
              if (sum(duplicated(Y$probeID)>0)){
                cat("[raw data: error] remove duplicated names in probeID \n")
                stop()
              } else {                
                j = which(colnames(Y)=="probeID")              
                .Object@probeID = Y[,j]
                Y = Y[,-j]
              }
            } else {
              .Object@probeID = NULL
            }
            if (dim(Y)[2]==1){
              patient.status = "U"
            } else {
              patient.status = "M"
            }
                        
            ########################################
            #                                      #
            # check data format on the dataframe   #
            #                                      #
            ########################################
            
            cat("[raw data: check] minimum number of positions per signal \n")
            n.com    = dim(Y)[1]
            min.ncom = 10
            
            if ( n.com<min.ncom){
              cat("[raw data: check] Problem with dataset \n")
              stop("[raw data: check] It must contain more than ", min.ncom , " points \n")
            }            
            
            if (patient.status=="M"){
              cat("[raw data: check] number of records per position \n")
              j        = which(apply((apply(Y,1,is.na)),2,sum)==dim(Y)[2])   
              if (length(j) > 0){
                cat("[raw data: check] Problem with non observed positions \n")
                cat("Position(s) number", j ," is/are missing in every patient(s) \n")
              }
            } else {
              cat("[raw data: check] number of records per position \n")
              if (sum(is.na(Y))>0){
                cat("[raw data: check] Problem with non observed positions \n")
                cat("Position(s) number", which(is.na(Y))," is/are missing \n")
              }
            }
            cat("[raw data: format] changing dataframe to compact format (list) \n")
            
            Res = lapply(names(Y), FUN = function(m){
              Y[,names(Y)==m]
            })
            names(Res) = names(Y)
            .Object@Y  = Res            
            return(.Object)
          }
          )
          
          
            


###############################################################
#                                                             #  
# GENERIC DEFINITION FOR CLASS CGHDATA                        #
#                                                             #
###############################################################


setGeneric( name = "getuniKmax"            ,def = function(.Object,CGHo,uniKmax=NULL){standardGeneric("getuniKmax")})
setGeneric( name = "getmultiKmax"          ,def = function(.Object,CGHo,uniKmax=NULL,multiKmax=NULL){standardGeneric("getmultiKmax")})
setGeneric( name = "multisegout"           ,def = function(.Object,seg.rep,Res,Kselect){standardGeneric("multisegout")})
setGeneric( name = "multisegmean"          ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("multisegmean")})
setGeneric( name = "multisegmixt"          ,def = function(.Object,CGHo,uniKmax,multiKmax,phi){standardGeneric("multisegmixt")})
setGeneric( name = "ILS"                   ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("ILS")})
setGeneric( name = "ILSclust.output"       ,def = function(.Object,mu,phi,tau){standardGeneric("ILSclust.output")})
setGeneric( name = "ILSclust"              ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("ILSclust")})
setGeneric( name = "multisegclust.output"  ,def = function(.Object,mu,phi,tau){standardGeneric("multisegclust.output")})
setGeneric( name = "multisegclust"         ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("multisegclust")})
setGeneric( name = "multiseg"              ,def = function(.Object,CGHo,uniKmax=NULL,multiKmax=NULL){standardGeneric("multiseg")})
setGeneric( name = "golden.search"         ,def = function(.Object,CGHo,uniKmax,multiKmax){standardGeneric("golden.search")})
setGeneric( name = "uniseg"                ,def = function(.Object,CGHo,uniKmax=NULL){standardGeneric("uniseg")})
setGeneric( name = "correctdata"           ,def = function(.Object,bias){standardGeneric("correctdata")})
setGeneric( name = "revertdata"            ,def = function(.Object,bias){standardGeneric("revertdata")})
setGeneric( name = "getbias"               ,def = function(.Object,CGHo,mu,bias,phi=c(),tau=c()){standardGeneric("getbias")})
setGeneric( name = "DP2EM"                 ,def = function(.Object,mu,theta=NULL){standardGeneric("DP2EM")})




setGeneric("removebias<-",function(object,value){standardGeneric("removebias<-")})
setReplaceMethod(
                 f="removebias",
                 signature="CGHdata",
                 definition=function(object,value){
                   object@Y = lapply(object@Y,FUN=function(y){y-value})
                   return (object)
                 }
                 )
setGeneric("revertbias<-",function(object,value){standardGeneric("revertbias<-")})
setReplaceMethod(
                 f="revertbias",
                 signature="CGHdata",
                 definition=function(object,value){
                   object@Y = lapply(object@Y,FUN=function(y){y+value})
                   return (object)
                 }
                 )
setMethod(
          f = "[",
          signature = "CGHdata",
          definition = function(x,i,j,drop){   
            if (i=="Y")                  {return(x@Y)}                else {}
            if (i=="GCcontent")          {return(x@GCcontent)}        else {}
            if (i=="genomic.position")   {return(x@genomic.position)} else {}
            if (i=="probeID")            {return(x@probeID)}          else {}            
          }          
          )
setMethod(
          f = "summary",
          signature = "CGHdata",
          definition = function(object){
            cat("****** Summary of CGHd object ******\n")
            cat("[CGHd summary] Patients ID \n");
            cat(names(object@Y),fill=TRUE,"\n")
            cat("\n")
            cat("[CGHd summary] number of points\n")
            print(length(object@Y[[1]]))
            cat("[CGHd summary] probeID records\n")
            print(!is.null(object@probeID))            
             cat("[CGHd summary] genomic position \n")
            print(!is.null(object@genomic.position))
            cat("[CGHd summary] GC content records\n")
            print(!is.null(object@GCcontent))
          }          
          )
setMethod(
          f = "show",
          signature = "CGHdata",
          definition = function(object){
            cat("****** CGHdata show ******\n")
            cat("[CGHd show] Data are in the list format [[patient]]\n")
            cat("[CGHd show] Data sample: \n")
            cat("Y[[",names(object@Y)[1],"]]\n",sep="")            
            print(object@Y[[1]][1:5])
            if (length(names(object@Y))>1){
              cat("Y[[",names(object@Y)[2],"]] \n",sep="")
              print(object@Y[[2]][1:5])              
            }
            cat("[CGHd show] probeID sample: \n")
            print(object@probeID[1:5])
            cat("[CGHd show] genomic positions sample: \n")
            print(object@genomic.position[1:5])
            cat("[CGHd show] GC content sample: \n")
            print(object@GCcontent[1:5])
          }          
          )

setMethod(
          f = "print",
          signature = "CGHdata",
          definition = function(x){
            cat("****** CGHdata print ******\n")
            cat("[CGHd print] Data are in the list format [[patient]]\n")
            cat("[CGHd print] Data sample: \n")
            cat("Y[[",names(x@Y)[1],"]]\n",sep="")            
            print(x@Y[[1]][1:5])
            if (length(names(x@Y))>1){
              cat("Y[[",names(x@Y)[2],"]] \n",sep="")            
              print(x@Y[[2]][1:5])              
            }
            cat("[CGHd print] probeID sample: \n")
            print(x@GCcontent[1:5])
            cat("[CGHd print] genomic positions sample: \n")
            print(x@genomic.position[1:5])
            cat("[CGHd print] GC content sample: \n")
            print(x@GCcontent[1:5])
          }          
          )

