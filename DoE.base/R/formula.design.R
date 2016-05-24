formula.design <- function(x, ..., response = NULL, degree = NULL, FUN=NULL, use.center=NULL, use.star=NULL, use.dummies=FALSE){
      ## open: aggregation of repeat.only designs
      ## open: aggregation of parameter designs
      ###### implement with function aggregate.design
      
      
      ## error checks
      xnam <- deparse(substitute(x))
      
      ## avoid NOTE about no global function definition
      ##if (!is.loaded("FrF2")) iscube <- function(design) NULL
      ## does not work

      
      if (!"design" %in% class(x)) stop("This function is applicable to class design objects only.")
      di <- design.info(x)
      ## capture designs from conf.design
          if (is.null(di)) stop("formula.design does not work for class design from package conf.design")
          if (is.null(use.center)) if (di$type=="ccd") use.center <- TRUE else use.center <- FALSE
          if (is.null(use.star)) if (di$type=="ccd") use.star <- TRUE else use.star <- FALSE

      if (is.null(di$response.names)) stop("formula.design needs at least one response in the design")
      if (!(is.null(degree) | is.numeric(degree))) stop("degree must be numeric or NULL")
      if (is.numeric(degree)){ 
           if(!degree > 0) stop("degree must be positive")
           if (!degree == round(degree)) stop("degree must be an integer number")
           }
      if (!is.null(response)) 
         if (length(response)>1) stop("formula.design can only handle one response at a time")

      ## prepare wide format design, if it has not yet been aggregated
      ## character responses only
      if(!is.null(di$responselist)){
      if (is.null(response)){
          if (di$response.names[1]==di$responselist[1,1]){
              if (is.null(FUN)) FUN <- "mean"
              x <- aggregate(x, FUN=FUN)
              di <- design.info(x)
              response <- di$response.names[1]
              assign(response, undesign(x)[,response])
          }}
      else{ 
      if (is.character(response) & !as.character(response) %in% di$response.names){
          if (!response %in% colnames(di$responselist)) stop("invalid response name")
              response <- which(colnames(di$responselist)==response)
              if (is.null(FUN)) FUN <- "mean"
              x <- aggregate(x, FUN=FUN)
              di <- design.info(x)
              response <- di$response.names[response]
              assign(response, undesign(x)[,response])
          }
          }
          } 

      ## check available responses
      respnam <-di$response.names
      respnamOK <- intersect(colnames(x),respnam)
      if (is.null(respnamOK) | length(respnamOK)==0)
          stop("For formula.design, the design requires at least one response to be available.")
      respnamOK <- respnamOK[which(sapply(x[,respnamOK], function(obj) all(!is.na(obj))))]
      if (length(respnamOK)==0) stop("the design does not contain any response variable with complete observations")
      respposOK <- which(di$response.names %in% respnamOK)
      
      ## check response given by user
      if (!is.null(response)){
         if (!(is.character(response) | is.numeric(response))) 
              stop("response must be a character string of the response name or a position number")
         if (is.numeric(response)) {
               if (response < 1 | response > length(respnam) | !response==round(response)) 
                     stop("if numeric, response must be an integer from 1 to ", length(respnam))
               response <- respnam[response]
               ## now response is character
            }
         if (is.character(response)){
               if (!response %in% colnames(x)) 
                     stop("response is not a column of x")
               if (!response %in% respnam) 
                     stop("response has not been declared a response variable")
               if (!response %in% respnamOK) 
                     stop("response has missing values, which precludes default analysis of the design")
            }
         }
      else response <- respnamOK[1]
      ## else: no response given by user
      
    ##  if (length(grep("center",di$type))>0 & !is.loaded("FrF2")) 
    ##      stop("For working with center point designs, package FrF2 must be loaded")
    ## is it possible to protect users from this error ?
      if ((length(grep("center",di$type))>0 | length(grep("ccd",di$type))>0) & !use.center){
          if (di$type=="ccd") x <- pickcube(x)
          x <- x[iscube(x),]
          fn <- names(di$factor.names)
          assign(response, x[,response])
          if (!is.null(di$block.name)){ if (nlevels(di$block.name)>1) assign(di$block.name, x[,di$block.name])
                                         else di$block.name <- NULL}
          for (i in 1:di$nfactors) 
            assign(fn[i], x[,fn[i]])
          class(x) <- c("design","data.frame")
          di$nruns <- di$ncube
          di$ncenter <- NULL
          di$replications <- di$replications[1]
          design.info(x) <- di
          #message("analysis without center points")
          # message removed, because it is more annoying than helpful 
          # particularly with repeated usage
      } 
          
      if (length(grep("ccd",di$type))>0 & use.center & !use.star){
          x <- pickcube(x)
#          fn <- names(di$factor.names)
#          assign(response, x[,response])
#          if (!is.null(di$block.name)) assign(di$block.name, x[,di$block.name])
#          for (i in 1:di$nfactors) 
#            assign(fn[i], x[,fn[i]])
#          class(x) <- c("design","data.frame")
#          di$nruns <- di$ncube
#          design.info(x) <- di
          #message("analysis without center points")
          # message removed, because it is more annoying than helpful 
          # particularly with repeated usage
      } 

     ## long format repeated measurement or parameter designs 
      if (di$repeat.only){ 
          if (is.null(FUN)) FUN <- "mean"
          x <- aggregate(reptowide(x),response=response, FUN=FUN)
          ## new response names!!
          response <- response.names(x)[1]
          di <- design.info(x)
          fn <- names(di$factor.names)
          assign(response, undesign(x)[,response])
          if (!is.null(di$block.name)) assign(block.name, undesign(x)[,block.name])
          for (i in 1:di$nfactors) 
            assign(fn[i], undesign(x)[,fn[i]])
          message("analysing repeated measurement ", FUN)
      }
      if (length(grep("param",di$type))>0 & length(grep("wide",di$type))==0){
          if (is.null(FUN)) FUN <- "mean"
          x <- aggregate(paramtowide(x),response=response, FUN=FUN)
          response <- response.names(x)[1]
          di <- design.info(x)
          fn <- names(di$factor.names)
          assign(response, undesign(x)[,response])
          for (i in 1:di$nfactors) 
            assign(fn[i], undesign(x)[,fn[i]])
          message("analysing outer array ", FUN)
      }
      
      ## identify and check response candidates
      type <- di$type
      
      ## default degrees: 1 for pb and oa, 2 for everything else
      if (is.null(degree) && substr(type,1,2) %in% c("pb","oa")) degree <- 1
      if (is.null(degree)) degree <- 2

      factor.names <- di$factor.names
      ## remove error effects from the formula for pb designs 
      ## for designs from package version 1.3 onwards
      if (length(grep("pb",type))>0){
          if (!(is.null(di$ndummies) | use.dummies)) factor.names <- factor.names[1:(di$nfactors-di$ndummies)]
      }
      
      block.name <- di$block.name
      if (!(is.null(block.name))){
      ## now degree is given
          if (degree==1){
              aus <- formula(paste(response, paste(c(block.name,names(factor.names)),collapse="+"),sep="~"))
              if (length(grep("center",type))>0){ 
                      if (use.center){
                      center <- !iscube(x)
                      aus <- formula(paste(response, paste(c(block.name, names(factor.names), "center"), collapse="+"),sep="~"))
                      }
                      else{
                      aus <- formula(paste(response, paste(c(block.name, names(factor.names)), collapse="+"),sep="~"))
                      }
                      }
              }
          if (degree > 1){ 
              if (substr(type,1,2) %in% c("pb","oa")) warning("degree > 1 is often inadequate with design types pb and oa")
              if (!type %in% c("bbd.blocked","ccd")){
                  if (length(grep("center",type))>0){
                    if (use.center){
                      center <- !iscube(x)
                      aus <- formula(paste(response, paste(c(block.name, 
                          paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""), "center"), collapse="+"),sep="~"))
                    }       
                    else 
                       aus <- formula(paste(response, 
                         paste(block.name, paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""),sep="+"),
                         sep="~"))
                  }
                  else{
                    aus <- formula(paste(response, 
                      paste(block.name,paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""),sep="+"),
                      sep="~"))
                  }
                  }
                  else{ 
                    ##bbd.blocked oder ccd
                    aus <- formula(paste(response, 
                      paste(block.name,paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""),
                      paste(paste("I(",names(factor.names),"^",degree,")",sep=""),collapse="+"),sep="+"),
                      sep="~"))
                    #else aus <- formula(paste(response, 
                    #  paste(block.name, paste("FO(",paste(names(factor.names),collapse=","),")",sep=""),
                    #      paste("TWI(",paste(names(factor.names),collapse=","),")",sep=""),
                    #      paste("PQ(",paste(names(factor.names),collapse=","),")",sep=""),sep="+"),
                    #  sep="~"))
                      }
          }
          }
      else {
          ## no blocks
          if (degree==1){
                 aus <- formula(paste(response, paste(names(factor.names),collapse="+"),sep="~"))
              if (length(grep("center",type))>0){ 
                   if (use.center){
                      center <- !iscube(x)
                      aus <- formula(paste(response, paste(c(names(factor.names), "center"), collapse="+"),sep="~"))
                      }
                   else
                      aus <- formula(paste(response, paste(c(names(factor.names)), collapse="+"),sep="~"))
              }
              }
          if (degree > 1){ 
              if (substr(type,1,2) %in% c("pb","oa")) warning("degree > 1 is often inadequate with design types pb and oa")
              aus <- formula(paste(response, paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""),sep="~"))
              if (type %in% c("bbd","ccd","lhs")){ 
                    aus <- formula(paste(response, 
                      paste(paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""),
                      paste(paste("I(",names(factor.names),"^",degree,")",sep=""),collapse="+"),sep="+"),
                      sep="~"))
                    #else aus <- formula(paste(response, 
                    #      paste("SO(",paste(names(factor.names),collapse=","),")",sep=""),
                    #      sep="~"))
                      }
              if (length(grep("Dopt",di$type)) > 0 ){
                  aus <- formula(paste(response, paste(as.character(di$formula), collapse=" ")))
              }
              if (length(grep("center",type))>0){ 
                   if (use.center){
                      center <- !iscube(x)
                      aus <- formula(paste(response, 
                      paste(c(paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""), "center"), collapse="+"),
                      sep="~"))
                      }
                   else aus <- formula(paste(response, 
                      paste(c(paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep="")), collapse="+"),
                      sep="~"))
                      }
              else {
              if (di$type=="ccd"){ 
                 aus <- formula(paste(response, paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""),sep="~"))
              if (use.center & !use.star) {
                 center <- !iscube(x)
                 aus <- formula(paste(response, 
                      paste(c(paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""), "center"), collapse="+"),
                      sep="~"))
                 }
              }
              else
              if (!di$type %in% c("bbd","lhs","Dopt","Dopt.blocked","Dopt.splitplot")){
              aus <- formula(paste(response, paste("(",paste(names(factor.names),collapse="+"),")^",degree,sep=""),sep="~"))
              }
              }
      }
      }
       aus   
}