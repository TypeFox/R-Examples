# Author: Babak Naimi, naimi.b@gmail.com
# Date :  May 2015
# Version 1.0
# Licence GPL v3

setMethod ('show' , 'sdmdata',
           function ( object ) {
             cat('class                                 :' , class(object), '\n')
             cat('===========================================================','\n')
             cat('number of species                     : ' , length(object@species.names) , '\n')
             cat('species names                         : ' , if (length(object@species.names) > 3) paste(c(object@species.names[1:3],'...'),collapse=', ') else paste(object@species.names,collapse=', ') , '\n')
             cat('number of fearures                    : ' , length(object@features.name), '\n')
             cat('fearure names                         : ' , if (length(object@features.name) > 0) {
               if (length(object@features.name) > 3) paste(c(object@features.name[1:3],'...'),collapse=', ') else paste(object@features.name,collapse=', ')
             } else NA, '\n')
             if (length(object@factors) > 0) {
               cat('which feature is categorical (factor) : ' , 
                 if (length(object@factors) > 3) paste(c(object@factors[1:3],'...'),collapse=', ') else paste(object@factors,collapse=', ')
               , '\n') 
             }
             
             typ <- c()
             for (n in object@species.names) typ <- c(typ,object@species[[n]]@type)
             if (length(unique(typ)) == 1) typ <- unique(typ)
             hasTest <- 'test' %in% .getGroupNames(object,levels=TRUE)
             hasTrain <- 'train' %in% .getGroupNames(object,levels=TRUE)
             cat('type                                  : ', if (length(typ) > 3) paste(c(typ[1:3],'...'),collapse=', ') else paste(typ,collapse=', '), '\n')
             cat('has independet test data?             : ' , hasTest, '\n')
             if (hasTrain) {
               cat('number of records                     : ', if (hasTest) paste('train-> ',length(.getGroupIndex(object,'train')),"; ",'test-> ',length(.getGroupIndex(object,'test')),sep='') else length(.getSpeciesIndex(object)),'\n')
             } else
             cat('number of records                     : ', if (is.null(object@features)) '0' else nrow(object@features),'\n')
             cat('has Coordinates?                      : ' , !is.null(object@info) && !is.null(object@info@coords)  , '\n')
           }
)
#-----------
setMethod ('show' , '.sdmCorSetting',
           function (object) {
             cv.n <- 1
             r.n <- 1
             r <- 1
             
             cat('class                                 :' , class(object), '\n')
             cat('========================================================','\n')
             cat('modelling methods                     : ' , paste(object@methods,collapse=', '), '\n')
             sn <- object@sdmFormula@species
             cat('species names                         : ' , if (length(sn) > 3) paste(length(sn),'species including:',paste(c(sn,'...'),collapse=', ')) else paste(sn,collapse=', '), '\n')
             cat('feature names                         : ' , if (length(object@sdmFormula@vars) > 3) paste(length(object@sdmFormula@vars),'variables including:',paste(c(object@sdmFormula@vars[1:3],'...'),collapse=', ')) else paste(object@sdmFormula@vars,collapse=', '), '\n')
             cat('feature types                         : ' , paste(unique(unlist(lapply(object@featuresFrame@feature.types,function(x) x@type))),collapse=', '), '\n')
             if (!is.null(object@replicate)) {
               r.n <- length(object@replicate)
               r <- object@n.replicates
               cat('replicate.methods (data partitioning) : ' , paste(object@replicate,collapse=','), '\n')
               if ('subsampleing' %in% object@replicate)
                 cat('test percentage                       : ' , object@test.percentage, '\n')
               cat('number of replicates (each method)    : ' , r, '\n')
               if ("cross-validation" %in% object@replicate) {
                 cat('n.folds in cross-validation           : ' ,object@cv.folds , '\n')
               }
             }
             n <- r * (r.n-1+object@cv.folds)
             cat('------------------------------------------\n')
             cat('number of runs                        : ' , paste(n,' for each model, ',n*length(object@methods),' in total... (per species)' ,sep=''), '\n')
             
           }
)
#---------------------

setMethod ('show' , 'sdmModels',
           function (object) {
             mi <- object@run.info
             for (i in c(2,3)) mi[,i] <- as.character(mi[,i])
             sp <- unique(mi$species)
             mo <- unique(mi$method)
             n <- 1
             
             cat('class                                 :' , class(object), '\n')
             cat('========================================================','\n')
             cat('number of species                     : ' , length(sp), '\n')
             cat('number of modelling methods           : ' , length(mo), '\n')
             cat('names of modelling methods            : ' , paste(mo,collapse=', '), '\n')
             
             if (!is.na(mi[1,4])) {
               r.n <- length(unique(mi$replication))
               n <- length(object@replicates[[1]])
               
               
               cat('replicate.methods (data partitioning) : ' , paste(unique(mi$replication),collapse=','), '\n')
               cat('number of replicates (each method)    : ' , object@setting@n.replicates, '\n')
               cat('toral number of replicates per model  : ' , paste(n,'(per species)'), '\n')
               if ('subsampling' %in% unique(mi$replication))
                 cat('test percentage (in subsampling)      : ' , object@setting@test.percentage, '\n')
               
             }
             cat('------------------------------------------\n')
             cat('model run success percentage (per species)  :\n')
             cat('------------------------------------------\n')
             
             p2 <- function(x) {
               o <- rep(NA,length(sp))
               for (i in seq_along(sp)) {
                 w1 <- mi$species == sp[i]
                 w2 <- mi$method == x
                 w <- mi$success[w1 & w2]
                 o[i] <- (length(which(w))/length(w))*100
               }
               o
             }
             a <- ''
             for (i in sp) a <- paste(a,i,paste(rep(' ',15 - length(unlist(strsplit(i,'')))),collapse=''),collapse='')
             cat(paste('method        ',a,collapse='|'),'\n')
             cat(paste(rep('-',length(unlist(strsplit(a,'')))+5),collapse=''),'\n')
             
             for (i in seq_along(mo)) {
               cat(paste(mo[i],paste(rep(' ',10 - length(unlist(strsplit(mo[i],'')))),collapse=''),' : ',paste(rep(' ',3),collapse=''),sep='')  , paste(p2(mo[i]),collapse='       |     '), '%\n')
             }
             wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
             
             
             cat('\n###################################################################\n')
             if (n > 1) {
               if (wtest == 'test.dep')
                 cat('model Mean performance (per species), using test dataset (generated using partitioning):\n')
               else if (wtest == 'test.indep')
                 cat('model Mean performance (per species), using independent test dataset:\n')
               else
                 cat('model Mean performance (per species), using training test dataset:\n')
             } else {
               if (wtest == 'test.dep')
                 cat('model performance (per species), using test dataset (generated using partitioning):\n')
               else if (wtest == 'test.indep')
                 cat('model performance (per species), using independent test dataset:\n')
               else
                 cat('model performance (per species), using training test dataset:\n')
              } 
               
             cat('--------------------------------------------------------------------\n')
             p <- function(x,sp) {
               a <- c()
               w1 <- mi$species == sp
               w2 <- mi$method == x
               id <- mi$modelID[w1 & w2 & mi$success]
               if (length(id) > 0) {
                 o <- ._getPerformance(object,id,wtest = NULL,s1=c('AUC','COR','Deviance'),s2='TSS',opt = 2)
                 a <- mean(sapply(o,function(x) x$AUC),na.rm=TRUE)
                 a <- c(a,mean(sapply(o,function(x) x$COR),na.rm=TRUE))
                 a <- c(a,mean(sapply(o,function(x) x$TSS),na.rm=TRUE))
                 a <- c(a,mean(sapply(o,function(x) x$Deviance),na.rm=TRUE))
                 a <- round(a,2)
                 a <- as.character(a)
               } else a <- as.character(c(NA,NA,NA,NA))
               b <- c()
               for (i in a) b <- c(b,paste(i,paste(rep(' ',7 - length(unlist(strsplit(i,'')))),collapse=''),collapse=''))
               paste(b,collapse='|     ')
             }
             
             for (spp in sp) {
               
               cat('\n ## species   : ',spp,'\n')
               cat('======================\n\n')
               cat(paste('methods',paste(rep(' ',3),collapse=''),' : ',paste(rep(' ',3),collapse=''),sep='')  , paste(c('AUC','COR','TSS','Deviance'),collapse='     |     '), '\n')
               cat('-------------------------------------------------------------------\n')
               for (i in seq_along(mo)) {
                 cat(paste(mo[i],paste(rep(' ',10 - length(unlist(strsplit(mo[i],'')))),collapse=''),' : ',paste(rep(' ',3),collapse=''),sep='')  , p(mo[i],spp), '\n')
               }
             }
           }
)
