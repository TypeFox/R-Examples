
prediction <- function(fit, data, data.Longi, t, window, group, MC.sample=0){

        if (missing(fit)) stop("Need a fit")
        if ((class(fit)!="frailtyPenal") & (class(fit)!="jointPenal") & class(fit)!='longiPenal' & class(fit)!='trivPenal') stop("The argument fit must be a frailtyPenal or jointPenal object")
        if (fit$istop != 1) stop("Attempting to do predictions with a wrong model")
        if (missing(data)) stop("Need data to do some predictions")
        if (missing(data.Longi) & (class(fit)=="longiPenal") & (class(fit)=="trivPenal")) stop("Need data.Longi to do predictions")
        if (missing(t) | missing(window)) stop("Need times and a window to do predictions")

        if (length(t)!=1 & length(window)!=1) stop("t and window can not be vector both at the same time")
        if (any(window <= 0)) stop("Window must be positive")
        if (is.unsorted(t)) stop("Last time of predictions must be greater than first one")
        if (any(t < 0)) stop("Be careful, negative time input")
        if ((MC.sample < 0) | (MC.sample > 1000))  stop("MC.sample needs to be positive integer up to 1000")

        if ((class(fit)=="jointPenal" | class(fit)=='longiPenal' | class(fit)=='trivPenal') & (!missing(group))) stop("No need for 'group' argument for predictions on a joint model")

  if(class(fit)=='jointPenal' | class(fit)=='trivPenal'){
        if (max(t+window) > max(fit$xR)) stop("Prediction times cannot exceed maximum time of observation")
        if (max(t+window) > max(fit$xD)) stop("Prediction times cannot exceed maximum time of observation")}

  if(class(fit)=='frailtyPenal'){
    if (max(t+window) > max(fit$x)) stop("Prediction times cannot exceed maximum time of observation")
  }

        if(class(fit)=='longiPenal'){
          if (max(t+window) > max(fit$xD)) stop("Prediction times cannot exceed maximum time of observation")
        }
        # seulement dans le cas du shared
        if (missing(group)) type <- "marginal"
        else type <- "conditional"

        # if (!(predTime >= min(fit$x1))) stop("predtime must be in the right range")
        # mettre un warning quand une variable est un factor dans le fit et pas dans le datapred => source d'erreur

        if (MC.sample==0) ICproba <- FALSE
        else ICproba <- TRUE

        np <- fit$npar
        b <- fit$b
        typeof <- fit$typeof
        nva1 <- fit$nvarRec
        nva2 <- fit$nvarEnd
  nva3 <- fit$nvarY
        ng <- fit$groups
        nst <- 2
        HIH <- fit$varHIHtotal

        # a definir meme si non utilise
        nz <- 1
        zi <- 0
        nbintervR <- 1
        nbintervDC <- 1
        time <- 0
        timedc <- 0

        if(typeof == 0){
                nz <- fit$n.knots.temp
                zi <- fit$zi
        }

        if(typeof == 1){
                nbintervR <- fit$nbintervR
                nbintervDC <- fit$nbintervDC
                time <- fit$time
                timedc <- fit$timedc
        }

        # nombre de predictions a faire pour chaque individu
        moving.window <- FALSE
        if (length(t)==1) moving.window <- TRUE

        if (moving.window){
                predTime <- t
                timeAll <- t+window #seq(predTime+window,predMax,by=window)
                if (class(fit) == "jointPenal" | class(fit)== "trivPenal") window <- 0
        }else{
                predTime <- t[1]
                timeAll <- t+window
        }
        ntimeAll <- length(timeAll)

        # recuperation des profils d'individus pour la prediction
        m <- fit$call
        m2 <- match.call()

        m$formula.terminalEvent <- m$formula.LongitudinalData <- m$data.Longi <- m$random <- m$id  <- m$link <- m$left.censoring <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$init.Random <- m$init.Eta <- m$Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$...  <- NULL



        m[[1]] <- as.name("model.frame")
        m3 <- m # pour recuperer les donnees du dataset initial en plus
        m[[3]] <- as.name(m2$data)

        
        if (class(fit) == "jointPenal" | class(fit)=="trivPenal"){

                temp <- as.character(m$formula[[2]])
                if (temp[1]=="Surv"){

                        if (length(temp) == 4) m$formula[[2]] <- paste(c("cbind(",temp[3],",",temp[4],")"),collapse=" ")
                        else if (length(temp) == 3) m$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],")"),collapse=" ")
                        else stop("Wrong Surv function")

                }else{ # SurvIC
                        if (length(temp) == 4) m$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],",",temp[4],")"),collapse=" ")
                        else if (length(temp) == 5) m$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],",",temp[4],",",temp[5],")"),collapse=" ")
                        else stop("Wrong SurvIC function")
                }
                m$formula <- unlist(strsplit(deparse(m$formula)," "))
                m$formula <- gsub("\"","",m$formula)

                ter <- grep("terminal",m$formula)
                if (ter==length(m$formula)) m$formula <- as.formula(paste(m$formula[-c(ter,max(which(m$formula=="+")))],collapse=""))
                else m$formula <- as.formula(paste(m$formula[-ter],collapse=""))

                if (fit$joint.clust==0){

                        m$formula <- unlist(strsplit(deparse(m$formula)," "))
                        clus <- grep("cluster",m$formula)
                        if (clus==length(m$formula)) m$formula <- as.formula(paste(m$formula[-c(clus,max(which(m$formula=="+")))],collapse=""))
                        else m$formula <- as.formula(paste(m$formula[-clus],collapse=""))
                }
        }else{
                if (fit$Frailty==TRUE ){
                        m$formula <- unlist(strsplit(deparse(m$formula)," "))
                        clus <- grep("cluster",m$formula)
      if (clus==length(m$formula)) m$formula <- as.formula(paste(m$formula[-c(clus,max(which(m$formula=="+")))],collapse=""))
                        else m$formula <- as.formula(paste(m$formula[-clus],collapse=""))
                }else if(class(fit)=='longiPenal'){
                  m$formula <- unlist(strsplit(deparse(m$formula)," "))
                  m$formula <- as.formula(paste(m$formula,collapse=""))

                }

                m$formula[[2]] <- NULL # pas besoin du Surv dans formula
        }

        dataset <- eval(m, sys.parent())

        dataset3 <- eval(m3, sys.parent())

        typeofY <- attr(model.extract(dataset3, "response"),"type")
        Y <- model.extract(dataset3, "response")

        if (typeofY=="right") tt1 <- Y[,1]
        else tt1 <- Y[,2]

        class(m$formula) <- "formula"
        special <- c("strata", "cluster", "subcluster", "terminal", "num.id", "timedep")

        Terms <- terms(m$formula, special, data = data)

        m$formula <- Terms

        dropx <- NULL

        if (class(fit) == "jointPenal" | class(fit) == 'trivPenal'){
                if (fit$joint.clust==1){ # joint classique
                        tempc <- untangle.specials(Terms, "cluster", 1:10)
                        dropx <- c(dropx,tempc$terms)
                        cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
                        uni.cluster <- unique(cluster)

                        ic <- model.extract(dataset, "response")[,2]
                        npred <- length(uni.cluster)
                        nrec <- max(table(cluster[ic==1]))
                       
                        if (temp[1]=="Surv"){
                                Y <- NULL
                               
                                for (i in uni.cluster) {
                                        temp <- model.extract(dataset, "response")[,1]
                                        temp <- temp[cluster==i & ic==1]
                                        Y <- c(Y,c(temp,rep(0,nrec-length(temp))))
                                }
                               
                                predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
                                trunctime <- rep(0,npred)
                                lowertime <- rep(0,npred)
                                uppertime <- rep(0,npred)
                                
                        }else{
                                stop("Predictions not allowed for interval-censored yet...") # a enlever plus tard
                                predtimerec <- matrix(0,nrow=npred)
                                if (length(temp) == 4){ # pas troncature
                                        temp <- model.extract(dataset, "response")
                                        trunctime <- rep(0,npred)
                                        lowertime <- temp[,1]
                                        uppertime <- temp[,2]
                                }
                                if (length(temp) == 5){ # troncature
                                        temp <- model.extract(dataset, "response")
                                        trunctime <- temp[,1]
                                        lowertime <- temp[,2]
                                        uppertime <- temp[,3]
                                }
                        }

                }else{ # joint cluster
                        tempnum <- untangle.specials(Terms, "num.id", 1:10)
                        dropx <- c(dropx,tempnum$terms)
                        num.id <- strata(dataset[, tempnum$vars], shortlabel = TRUE)
                        uni.num.id <- unique(num.id)

                        ic <- model.extract(dataset, "response")[,2]
                        npred <- length(uni.num.id)
                        nrec <- max(table(num.id[ic==1]))

                        if (temp[1]=="Surv"){
                                Y <- NULL
                                for (i in uni.num.id) {
                                        temp <- model.extract(dataset, "response")[,1]
                                        temp <- temp[num.id==i & ic==1]
                                        Y <- c(Y,c(temp,rep(0,nrec-length(temp))))
                                }
                                predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
                                trunctime <- rep(0,npred)
                                lowertime <- rep(0,npred)
                                uppertime <- rep(0,npred)
                        }else{
                                stop("Predictions not allowed for interval-censored yet...") # a enlever plus tard
                                predtimerec <- matrix(0,nrow=npred)
                                if (length(temp) == 4){ # pas troncature
                                        temp <- model.extract(dataset, "response")
                                        trunctime <- rep(0,npred)
                                        lowertime <- temp[,1]
                                        uppertime <- temp[,2]
                                }
                                if (length(temp) == 5){ # troncature
                                        temp <- model.extract(dataset, "response")
                                        trunctime <- temp[,1]
                                        lowertime <- temp[,2]
                                        uppertime <- temp[,3]
                                }
                        }
                }
        }else{
                if (fit$Frailty){
                        class(m3$formula) <- "formula"
                        Terms3 <- terms(m3$formula, special, data = data)
                        m3$formula <- Terms3

                        tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
                        # je recupere le cluster du dataframe de depart et non pas des predictions

                        cluster <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
                        uni.cluster <- unique(cluster)
                }
        }
 
        if (!is.null(dropx)) newTerms <- Terms[-dropx]
        else newTerms <- Terms

        X <- model.matrix(newTerms, dataset)

        if (ncol(X) > 1) X <- X[, -1, drop = FALSE]


        if (class(fit) == "jointPenal"){#### for joint frailty

                if (fit$joint.clust==1) vaxpred <- aggregate(X,by=list(cluster),FUN=function(x) {x[1]})[,-1]
                else vaxpred <- aggregate(X,by=list(num.id),FUN=function(x) {x[1]})[,-1]

                # recuperation des variables partie deces
                m3 <- fit$call
                m2 <- match.call()

   m3$formula.LongitudinalData <- m3$data.Longi <- m3$random <- m3$id <- m3$link <- m3$left.censoring <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa <- m3$maxit <- m3$hazard <- m3$nb.int <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$init.Theta <- m3$init.Alpha <- m3$Alpha <- m3$init.Random <- m3$init.Eta <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$... <- NULL

                m3$formula[[3]] <- m3$formula.terminalEvent[[2]]
                m3$formula.terminalEvent <- NULL
                m3[[1]] <- as.name("model.frame")
                m3[[3]] <- as.name(m2$data)

                temp <- as.character(m3$formula[[2]])
                if (temp[1]=="Surv"){
                        if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[3])
                        else if (length(temp) == 3) m3$formula[[2]] <- as.name(temp[2])
                        else stop("Wrong Surv function")
                }else{ # SurvIC
                        if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[2])
                        else if (length(temp) == 5) m3$formula[[2]] <- as.name(temp[3])
                        else stop("Wrong SurvIC function")
                }

                datasetdc <- eval(m3, sys.parent())

                class(m3$formula) <- "formula"
                special2 <- c("strata", "timedep")

                Terms2 <- terms(m3$formula, special2, data = data)

                X2 <- model.matrix(Terms2, datasetdc)
                if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]
                if (fit$joint.clust==1) vaxdcpred <- aggregate(X2,by=list(cluster),FUN=function(x) {x[1]})[,-1]
                else vaxdcpred <- aggregate(X2,by=list(num.id),FUN=function(x) {x[1]})[,-1]

                cat("\n")
                cat("Calculating the probabilities ... \n")
                if(fit$logNormal==0){
                ans <- .Fortran("predict",
                                as.integer(np),
                                as.double(b),
                                as.integer(nz),
                                as.integer(nbintervR),
                                as.integer(nbintervDC),
                                as.integer(nva1),
                                as.integer(nva2),
                                as.integer(nst),
                                as.integer(typeof),
                                as.double(zi),
                                as.double(HIH),
                                as.double(time),
                                as.double(timedc),
                                as.integer(ntimeAll),
                                as.integer(npred),
                                as.double(predTime),
                                as.double(window),
                                as.double(predtimerec),
                                as.integer(nrec),
                                as.double(as.matrix(vaxpred)),
                                as.double(as.matrix(vaxdcpred)),
                                pred1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                pred2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                pred3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                predlow1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                predhigh1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                predlow2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                predhigh2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                predlow3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                predhigh3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                                icproba=as.integer(ICproba),
                                as.integer(MC.sample),
                                as.integer(fit$intcens),
                                as.double(trunctime),
                                as.double(lowertime),
                                as.double(uppertime),
                                as.integer(moving.window),
                                as.double(timeAll),
                                PACKAGE = "frailtypack")
                }else{
     #AK: joint log-normal

      ans <- .Fortran("predict_LogN",
                      as.integer(np),
                      as.double(b),
                      as.integer(nz),
                      as.integer(nbintervR),
                      as.integer(nbintervDC),
                      as.integer(nva1),
                      as.integer(nva2),
                      as.integer(nst),
                      as.integer(typeof),
                      as.double(zi),
                      as.double(HIH),
                      as.double(time),
                      as.double(timedc),
                      as.integer(ntimeAll),
                      as.integer(npred),
                      as.double(predTime),
                      as.double(window),
                      as.double(predtimerec),
                      as.integer(nrec),
                      as.double(as.matrix(vaxpred)),
                      as.double(as.matrix(vaxdcpred)),
                      pred1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      pred2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      pred3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predlow1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predhigh1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predlow2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predhigh2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predlow3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predhigh3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      icproba=as.integer(ICproba),
                      as.integer(MC.sample),
                      as.integer(fit$intcens),
                      as.double(trunctime),
                      as.double(lowertime),
                      as.double(uppertime),
                      as.integer(moving.window),
                      as.double(timeAll),
                      PACKAGE = "frailtypack")
                }

                out <- NULL
                out$call <- match.call()
                out$name.fit <- match.call()[[2]]
                out$npred <- npred
                out$window <- window
                out$predtimerec <- predtimerec
                out$moving.window <- moving.window
                if (moving.window){
                        out$x.time <- timeAll
                        out$t <- predTime
                }else{
                        out$x.time <- timeAll - window
                }
                if (fit$joint.clust==1) out$group <- uni.cluster
                else out$group <- uni.num.id

                if (!fit$intcens){
                        out$pred1 <- matrix(ans$pred1,nrow=npred,ncol=ntimeAll)
                        rownames(out$pred1) <- paste("ind",out$group)
                        colnames(out$pred1) <- c("times",rep(" ",ntimeAll-1))
                        out$pred3 <- matrix(ans$pred3,nrow=npred,ncol=ntimeAll)
                        rownames(out$pred3) <- paste("ind",out$group)
                        colnames(out$pred3) <- c("times",rep(" ",ntimeAll-1))
                }
                out$pred2 <- matrix(ans$pred2,nrow=npred,ncol=ntimeAll)
                rownames(out$pred2) <- paste("ind",out$group)
                colnames(out$pred2) <- c("times",rep(" ",ntimeAll-1))
                out$icproba <- ICproba
                if (ICproba){
                        if (!fit$intcens){
                                out$predlow1 <- matrix(ans$predlow1,nrow=npred,ncol=ntimeAll)
                                out$predhigh1 <- matrix(ans$predhigh1,nrow=npred,ncol=ntimeAll)
                                rownames(out$predlow1) <- paste("ind",out$group)
                                colnames(out$predlow1) <- c("times",rep(" ",ntimeAll-1))
                                rownames(out$predhigh1) <- paste("ind",out$group)
                                colnames(out$predhigh1) <- c("times",rep(" ",ntimeAll-1))
                                out$predlow3 <- matrix(ans$predlow3,nrow=npred,ncol=ntimeAll)
                                out$predhigh3 <- matrix(ans$predhigh3,nrow=npred,ncol=ntimeAll)
                                rownames(out$predlow3) <- paste("ind",out$group)
                                colnames(out$predlow3) <- c("times",rep(" ",ntimeAll-1))
                                rownames(out$predhigh3) <- paste("ind",out$group)
                                colnames(out$predhigh3) <- c("times",rep(" ",ntimeAll-1))
                        }
                        out$predlow2 <- matrix(ans$predlow2,nrow=npred,ncol=ntimeAll)
                        out$predhigh2 <- matrix(ans$predhigh2,nrow=npred,ncol=ntimeAll)
                        rownames(out$predlow2) <- paste("ind",out$group)
                        colnames(out$predlow2) <- c("times",rep(" ",ntimeAll-1))
                        rownames(out$predhigh2) <- paste("ind",out$group)
                        colnames(out$predhigh2) <- c("times",rep(" ",ntimeAll-1))
                }
                out$joint.clust <- fit$joint.clust
                out$intcens <- fit$intcens

                cat("Predictions done for",npred,"subjects and",ntimeAll,"times \n")

                class(out) <- c("predJoint")
        }else if(class(fit)=="longiPenal" | class(fit)=="trivPenal"){## predictions joint for longitudinal data and terminal event or longitudinal dtaa, recurrent events and a terminal event
          cat("\n")
          cat("Calculating the probabilities ... \n")

          if(class(fit)=="longiPenal"){
      expBX <- exp(X %*% fit$coef[1:fit$nvarEnd])
          }else{

           vaxpred <- aggregate(X,by=list(cluster),FUN=function(x) {x[1]})[,-1]

            # recuperation des variables partie deces
            m3 <- fit$call
            m2 <- match.call()


            m3$formula.LongitudinalData <- m3$data.Longi <- m3$random <- m3$id <- m3$link <- m3$left.censoring <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa <- m3$maxit <- m3$hazard <- m3$nb.int <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$init.Theta <- m3$init.Alpha <- m3$Alpha <- m3$init.Random <- m3$init.Eta <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$... <- NULL

            m3$formula[[3]] <- m3$formula.terminalEvent[[2]]
            m3$formula.terminalEvent <- NULL
            m3[[1]] <- as.name("model.frame")
            m3[[3]] <- as.name(m2$data)

            temp <- as.character(m3$formula[[2]])

              if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[3])
              else if (length(temp) == 3) m3$formula[[2]] <- as.name(temp[2])
              else stop("Wrong Surv function")

           
            datasetdc <- eval(m3, sys.parent())

            class(m3$formula) <- "formula"
            special2 <- c("strata", "timedep")

            Terms2 <- terms(m3$formula, special2, data = data)

            X2 <- model.matrix(Terms2, datasetdc)
            if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]
             vaxdcpred <- aggregate(X2,by=list(cluster),FUN=function(x) {x[1]})[,-1]


          }

          # nombre de predictions a faire pour chaque individu
          if (moving.window){ # 2 facons differentes de faire des predictions, soit h evolue, soit t evolue
            sequence2 <- t+window #seq(predTime+window,predMax,by=window)
            sequence <- rep(predTime,times=length(sequence2))
          }else{
            sequence <- t #seq(predTime,predMax,length=50)
            sequence2 <- t+window #sequence+window
          }

          predMat <- NULL

    m2 <- fit$call
          m2$formula <-  m2$data <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$method.GH <- m2$... <- NULL

          special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")

          #========= Longitudinal Data preparation =========================
          class(m2$formula.LongitudinalData) <- "formula"

          TermsY <- terms(m2$formula.LongitudinalData, special, data = data.Longi)

          llY <- attr(TermsY, "term.labels")#liste des variables explicatives
          ord <- attr(TermsY, "order")

          #=========================================================>

          name.Y <- as.character(attr(TermsY, "variables")[[2]])
          yy <- data.Longi[,which(names(data.Longi)==name.Y)]

          # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
          ind.placeY <- which(llY%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llY)],function(x) length(levels(x)))>2)))
          vec.factorY <- NULL
          vec.factorY <- c(vec.factorY,llY[ind.placeY])


          mat.factorY <- matrix(vec.factorY,ncol=1,nrow=length(vec.factorY))

          # Fonction servant a prendre les termes entre "as.factor"
          vec.factorY <-apply(mat.factorY,MARGIN=1,FUN=function(x){
            if (length(grep("as.factor",x))>0){
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- length(unlist(strsplit(x,split="")))-1
              return(substr(x,start=pos1,stop=pos2))
            }else{
              return(x)
            }})

          ind.placeY <- grep(paste(vec.factorY,collapse="|"),llY)

          if(is.factor(data.Longi[,names(data.Longi)==llY[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llY[1]])-1
          else X_L<- data.Longi[,names(data.Longi)==llY[1]]

          if(length(llY)>1){
            for(i in 2:length(llY)){
              if(is.factor(data.Longi[,names(data.Longi)==llY[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llY[i]])-1)
              else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llY[i]])
            }}
          #X_L<- data.Longi[,names(data.Longi)%in%(llY)]

          if(sum(ord)>length(ord)){
            for(i in 1:length(ord)){
              if(ord[i]>1){
                v1 <- strsplit(as.character(llY[i]),":")[[1]][1]
                v2 <- strsplit(as.character(llY[i]),":")[[1]][2]
                if(is.factor(data.Longi[,names(data.Longi)==v1]) && length(levels(data.Longi[,names(data.Longi)==v1]))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
                if(is.factor(data.Longi[,names(data.Longi)==v2]) && length(levels(data.Longi[,names(data.Longi)==v2]))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
                if(is.factor(data.Longi[,names(data.Longi)==v1]) || !is.factor(data.Longi[,names(data.Longi)==v2])){
                  X_L <- cbind(X_L,(as.numeric(data.Longi[,names(data.Longi)==v1])-1)*data.Longi[,names(data.Longi)==v2])
                  llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v1])[2],sep="")
                }else if(!is.factor(data.Longi[,names(data.Longi)==v1]) || is.factor(data.Longi[,names(data.Longi)==v2])){
                  X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*(as.numeric(data.Longi[,names(data.Longi)==v2])-1))
                  llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v2])[2],sep="")
                }else{
                  X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*data.Longi[,names(data.Longi)==v2])
                }

              }
            }
          }

          if(dim(X_L)[2]!=length(llY))stop("The variables in the longitudinal part must be in the data.Longi")
          X_L <- as.data.frame(X_L)
          names(X_L) <- llY

          Intercept <- rep(1,dim(X_L)[1])

          if(fit$intercept)X_L <- cbind(Intercept,X_L)

          X_Lall<- X_L
          "%+%"<- function(x,y) paste(x,y,sep="")

          if(length(vec.factorY) > 0){
            for(i in 1:length(vec.factorY)){
              X_L <- cbind(X_L[,-(which(names(X_L)==vec.factorY[i]))],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1])
            }


            vect.factY<-names(X_L)[which(!(names(X_L)%in%llY))]

            occurY <- rep(0,length(vec.factorY))

            for(i in 1:length(vec.factorY)){
              #occur[i] <- sum(vec.factor[i] == vect.fact)
              occurY[i] <- length(grep(vec.factorY[i],vect.factY))
            }
          }

          if (ncol(X_L) == 0){
            noVarY <- 1
          }else{
            noVarY <- 0
          }
          #=========================================================>

          clusterY <- data.Longi$id
          maxy_rep <- max(table(clusterY))
          uni.cluster<-as.factor(unique(clusterY))
          npred <- length(uni.cluster)


          nvarY<-ncol(X_L) #nvar==1 correspond a 2 situations:


          varY <- as.matrix(sapply(X_L, as.numeric))


          #=======================================>
          #======= Construction du vecteur des indicatrice
          if(length(vec.factorY) > 0){
            #           ind.place <- ind.place -1
            k <- 0
            for(i in 1:length(vec.factorY)){
              ind.placeY[i] <- ind.placeY[i]+k
              k <- k + occurY[i]-1
            }
          }

   if(fit$link=="Random-effects")link <- 1
   if(fit$link=="Current-level") link <- 2

         if(fit$leftCensoring==FALSE){s_cag_id = 0
       s_cag = 0}else{
       s_cag_id = 1
       s_cag = fit$leftCensoring.threshold
    }

   if(class(fit)=="longiPenal"){
         ans <- .Fortran("predict_biv",
           as.integer(np),
                 as.double(b),
                 as.integer(nz),
                 as.integer(nva2),
                 as.integer(nva3),
           as.integer(fit$ne_re),
           as.integer(fit$netadc),
           as.integer(link),
                 as.integer(nst),
                 as.integer(typeof),
                 as.double(zi),
                 as.double(HIH),
                 as.integer(ntimeAll),
                 as.integer(npred),
                 as.double(predTime),
                 as.double(window),
                 as.integer(fit$max_rep),
           as.double(yy),
                 as.double(as.matrix(X)),
                 as.double(as.matrix(varY)),
                 as.integer(clusterY),
           as.integer(length(clusterY)),
            pred=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                  predlow=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                 predhigh=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                 icproba=as.integer(ICproba),
                 as.integer(MC.sample),
                 as.integer(moving.window),
                 as.double(timeAll),
           as.integer(s_cag_id),
           as.double(s_cag),
                 PACKAGE = "frailtypack")

         predMat <- matrix(ans$pred,nrow=nrow(data),ncol=ntimeAll)
         predMatLow <- matrix(ans$predlow,nrow=nrow(data),ncol=ntimeAll)
         predMatHigh <- matrix(ans$predhigh,nrow=nrow(data),ncol=ntimeAll)


         out <- NULL
         out$call <- match.call()
         out$name.fit <- match.call()[[2]]
         out$npred <- npred
         out$moving.window <- moving.window
         if (moving.window){
           out$x.time <- sequence2
           out$t <- predTime
         }else{
           out$x.time <- sequence
         }
         out$group <- uni.cluster
       
         out$pred <- predMat

         colnames(out$pred) <- c("times",rep(" ",dim(out$pred)[2]-1))
         rownames(out$pred) <- paste("ind",1:out$npred)

         out$icproba <- ICproba
         if (ICproba){
           out$predLow <- predMatLow
           out$predHigh <- predMatHigh
           colnames(out$predLow) <- c("times",rep(" ",dim(out$predLow)[2]-1))
           rownames(out$predLow) <- paste("ind",1:out$npred)
           colnames(out$predHigh) <- c("times",rep(" ",dim(out$predHigh)[2]-1))
           rownames(out$predHigh) <- paste("ind",1:out$npred)
         }
         out$window <- window
         out$trivariate <- FALSE
   }else if(class(fit)=="trivPenal"){
  
     ans <- .Fortran("predict_tri",
                     as.integer(np),
                     as.double(b),
                     as.integer(nz),
                     as.integer(nva1),
                     as.integer(nva2),
                     as.integer(nva3),
                     as.integer(fit$ne_re),
                     as.integer(fit$netar),
                     as.integer(fit$netadc),
                     as.integer(link),
                     as.integer(nst),
                     as.integer(typeof),
                     as.double(zi),
                     as.double(HIH),
                     as.integer(ntimeAll),
                     as.integer(npred),
                     as.double(predTime),
                     as.double(window),
                     as.double(predtimerec),
                     as.integer(nrec),
                     as.integer(fit$max_rep),
                     as.double(yy),
                     as.double(as.matrix(vaxpred)),
                     as.double(as.matrix(X)),
                     as.double(as.matrix(varY)),
                     as.integer(clusterY),
                     as.integer(length(clusterY)),
                     as.integer(npred),
                     pred=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                     predlow=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                     predhigh=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                     icproba=as.integer(ICproba),
                     as.integer(MC.sample),
                     as.integer(moving.window),
                     as.double(timeAll),
                     as.integer(s_cag_id),
                     as.double(s_cag),
                     PACKAGE = "frailtypack")


     out <- NULL
     out$call <- match.call()
     out$name.fit <- match.call()[[2]]
     out$npred <- npred
     out$window <- window
     out$predtimerec <- predtimerec
     out$moving.window <- moving.window
     if (moving.window){
       out$x.time <- timeAll
       out$t <- predTime
     }else{
       out$x.time <- timeAll - window
     }

   out$group <- uni.cluster


       out$pred <- matrix(ans$pred,nrow=npred,ncol=ntimeAll)
       rownames(out$pred) <- paste("ind",out$group)
       colnames(out$pred) <- c("times",rep(" ",ntimeAll-1))


     out$icproba <- ICproba
     if (ICproba){

         out$predLow <- matrix(ans$predlow,nrow=npred,ncol=ntimeAll)
         out$predHigh <- matrix(ans$predhigh,nrow=npred,ncol=ntimeAll)
         rownames(out$predLow) <- paste("ind",out$group)
         colnames(out$predLow) <- c("times",rep(" ",ntimeAll-1))
         rownames(out$predHigh) <- paste("ind",out$group)
         colnames(out$predHigh) <- c("times",rep(" ",ntimeAll-1))

     }
    out$trivariate <- TRUE
   }



   cat("Predictions done for",npred,"subjects and",ntimeAll,"times \n")


 class(out) <- "predLongi"



  }else if(class(fit)=="frailtyPenal"){## predictions shared
                cat("\n")
                cat("Calculating the probabilities ... \n")

                expBX <- exp(X %*% fit$coef)

                # nombre de predictions a faire pour chaque individu
                if (moving.window){ # 2 facons differentes de faire des predictions, soit h evolue, soit t evolue
                        sequence2 <- t+window #seq(predTime+window,predMax,by=window)
                        sequence <- rep(predTime,times=length(sequence2))
                }else{
                        sequence <- t #seq(predTime,predMax,length=50)
                        sequence2 <- t+window #sequence+window
                }

                predMat <- NULL

                if (fit$Frailty){
                        if (type=="marginal"){ ## marginal ##
                          if(fit$logNormal==0){ ## Gamma distribution
                            for (k in 1:nrow(data)){
                                          vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
                                          vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
                        pred <- 1-((1+fit$theta*(-log(vect.survival.X)))/(1+fit$theta*(-log(vect.survival.X.horizon))))**(1/fit$theta)
                                          predMat <- rbind(predMat,pred)
                                  }
        }else{ #AK: Normal distribution
          mat.survival.X <- NULL
          mat.survival.X.horizon <- NULL
          for (k in 1:nrow(data)){
            vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
            vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]

            mat.survival.X <- rbind( mat.survival.X,vect.survival.X)
            mat.survival.X.horizon <- rbind( mat.survival.X.horizon,vect.survival.X.horizon)
          }
          ans <- .Fortran("predict_LogN_sha",
                            as.integer(nrow(data)),
                            as.double(mat.survival.X),
                            as.double(mat.survival.X.horizon),
                            as.double(fit$sigma2),
                            pred=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                            as.integer(0),
                            as.integer(ntimeAll),
                            as.integer(MC.sample),
                            as.double(rep(0,MC.sample)),
                            as.double(matrix(0,nrow=nrow(data)*MC.sample,ncol=ntimeAll)),
                            as.double(matrix(0,nrow=nrow(data)*MC.sample,ncol=ntimeAll)),
                            predlow1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                            predhigh1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                            PACKAGE = "frailtypack")
            predMat <- matrix(ans$pred,nrow=nrow(data),ncol=ntimeAll)


          }

                        }else{ ## conditional ##
                                if (!(group %in% uni.cluster)) stop("Are you sure that the group is present in your cluster variable ?")

                                for (k in 1:nrow(data)){
                                        vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
                                        vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
          if (fit$logNormal==0){ # Gamma distribution
            pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[uni.cluster==group]}
          else{ #AK: Normal distribution
            pred <- 1-(vect.survival.X.horizon/vect.survival.X)**exp(fit$frailty.pred[uni.cluster==group])}

                                        predMat <- rbind(predMat,pred)
                                }
                        }
                }else{ ## for Cox model
                        if (!missing(group)) stop("No need for a group to predict on a proportionnal hazard model")

                        for (k in 1:nrow(data)){
                                vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
                                vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
                                pred <- 1-(vect.survival.X.horizon/vect.survival.X)
                                predMat <- rbind(predMat,pred)
                        }
                }

                # -------------------------------------------------------- #
                # calcul des bornes de confiances (methode de Monte Carlo) #
                # -------------------------------------------------------- #
                        if (ICproba){
                        balea <- mvrnorm(MC.sample,fit$b,fit$varHtotal)
                        #print(fit$b)
                        #print(apply(balea,2,mean))
                        #print(fit$b-apply(balea,2,mean))

                        if (fit$Frailty) { #AK: For Gamma we have variance theta and for Normal we have variance sigma2
        if(fit$logNormal==0)theta.mc <- balea[,fit$np-fit$nvar]^2
        if(fit$logNormal==1)sigma2.mc <- balea[,fit$np-fit$nvar]^2
                        }

                        aleaCoef <- balea[,(fit$np-fit$nvar+1):(fit$np)]
                        expBX.mc <- exp(X %*% t(aleaCoef))

                        # recuperation parametres de la fonction de risque/survie (splines,piecewise,weibull)
                        if (fit$typeof == 0){
                                para.mc <- balea[,1:(fit$n.knots+2)]^2
                                if(fit$n.strat == 2) para.mc2 <- balea[,(fit$n.knots+3):(2*(fit$n.knots+2))]^2
                                else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$n.knots+2)
                        }else if (fit$typeof == 1){
                                para.mc <- balea[,1:(fit$nbintervR)] # attention de ne pas elever au carre
                                if(fit$n.strat == 2) para.mc2 <- balea[,(fit$nbintervR+1):(2*fit$nbintervR)]
                                else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$nbintervR)
                        }else{
                                para.mc <- balea[,1:2]^2
                                if(fit$n.strat == 2) para.mc2 <- balea[,2:4]^2
                                else para.mc2 <- matrix(0,nrow=MC.sample,ncol=2)
                        }

                        survival.mc <- function(t,ObjFrailty,para1,para2){ # dans les trois cas para1 et para2 seront traites differemment
                                if (ObjFrailty$typeof == 0){ # splines
                                        nz <- ObjFrailty$n.knots
                                        zi <- ObjFrailty$zi
                                        res <- NULL
                                        nst <- ObjFrailty$n.strat
                                        out <- .Fortran("survival",as.double(t),as.double(para1),as.double(para2),as.integer(nz+2),
                                        as.double(zi),survival=as.double(c(0,0)),lam=as.double(c(0,0)),as.integer(nst),PACKAGE = "frailtypack") # lam ajoute suite aux modif de survival
                                        if(ObjFrailty$n.strat == 2){
                                                res <- c(res,out$survival)
                                        }else{
                                                res <- c(res,out$survival[1])
                                        }
                                        return(res)
                                }
                                if (ObjFrailty$typeof == 1){ # piecewise
                                        res <- NULL
                                        if (ObjFrailty$n.strat == 2) b <- c(para1,para2)
                                        else b <- para1
                                        time <- ObjFrailty$time
                                        out <- .Fortran("survival_cpm",as.double(t),as.double(b),
                                        as.integer(ObjFrailty$n.strat),as.integer(ObjFrailty$nbintervR),
                                        as.double(time),survival=as.double(c(0,0)),PACKAGE = "frailtypack")
                                        if(ObjFrailty$n.strat == 2){
                                                res <- c(res,out$survival)
                                        }else{
                                                res <- c(res,out$survival[1])
                                        }
                                        return(res)
                                }
                                if (ObjFrailty$typeof == 2){ # weibull
                                        res <- NULL
                                        sh1 <- para1[1]
                                        sc1 <- para1[2]
                                        res <- c(res,exp(-(t/sc1)^sh1))
                                        if(ObjFrailty$n.strat == 2){
                                                sh1 <- para2[1]
                                                sc1 <- para2[2]
                                                res <- c(res,exp(-(t/sc1)^sh1))
                                        }
                                        return(res)
                                }
                        }

                        # calcul de la somme des risques cumules juste pour le groupe defini
                        X3 <- model.matrix(newTerms, dataset3)
                        if (ncol(X3) > 1) X3 <- X3[, -1, drop = FALSE]

                        expBX3 <- exp(X3 %*% fit$coef)
                        if ((fit$Frailty) & (type=="conditional")) res1 <- sum((-log(sapply(tt1[which(cluster==group)],survival,ObjFrailty=fit))) %*% expBX3[which(cluster==group)])

                        predMatLow <- NULL
                        predMatHigh <- NULL

                        frailty.mc <- NULL

                        if (fit$Frailty){
                                if (type=="marginal"){ ## marginal ##
                                  if(fit$logNormal==0){ ## Gamma distribution
                                          for (k in 1:nrow(data)){
                                                  realisations <- NULL
                  for (i in 1:MC.sample){
                                                        vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                                                        vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                                                    pred <- 1-((1+theta.mc[i]*(-log(vect.survival.X)))/(1+theta.mc[i]*(-log(vect.survival.X.horizon))))**(1/theta.mc[i])
                realisations <- cbind(realisations,pred)
                                                }
              predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
                                                predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
                                        }
                                  }else{ #AK: Normal distribution

                                    mat.survival.X.mc <- NULL
                                    mat.survival.X.horizon.mc <- NULL

              for(i in 1:MC.sample){
                mat.survival.X.samp <- NULL
                  mat.survival.X.horizon.samp <- NULL
                for(k in 1:nrow(data)){
                   vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                   vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                   mat.survival.X.samp <- rbind(mat.survival.X.samp,vect.survival.X.samp)
                   mat.survival.X.horizon.samp <- rbind(mat.survival.X.horizon.samp,vect.survival.X.horizon.samp)

                }
                mat.survival.X.mc <- rbind(mat.survival.X.mc,mat.survival.X.samp)
                mat.survival.X.horizon.mc <- rbind(mat.survival.X.horizon.mc,mat.survival.X.horizon.samp)

              }
              ans <- .Fortran("predict_LogN_sha",
                                                    as.integer(nrow(data)),
                                                    as.double(mat.survival.X),
                                                    as.double(mat.survival.X.horizon),
                                                    as.double(fit$sigma2),
                                                    pred=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                                                    as.integer(1),
                                                    as.integer(ntimeAll),
                                                    as.integer(MC.sample),
                            as.double(sigma2.mc),
                                                    as.double(mat.survival.X.mc),
                            as.double(mat.survival.X.horizon.mc),
                                                    predlow1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                                                    predhigh1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                                                    PACKAGE = "frailtypack")

            predMatLow <- matrix(ans$predlow1,nrow=nrow(data),ncol=ntimeAll)
            predMatHigh <- matrix(ans$predhigh1,nrow=nrow(data),ncol=ntimeAll)

                                  }
                                }else{ ## conditional ##
                                  mi <- fit$n.eventsbygrp[uni.cluster==group]
                                        for (k in 1:nrow(data)){
                                                realisations <- NULL
                                                for (i in 1:MC.sample){

                                                        vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                                                        vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]

              if(fit$logNormal==0){## Gamma distribution

                if (k == 1) frailty.mc <- c(frailty.mc,rgamma(1,shape=fit$n.eventsbygrp[uni.cluster==group]+1/theta.mc[i],scale=1/(res1+1/theta.mc[i])))
                   pred <- 1-(vect.survival.X.horizon/vect.survival.X)**frailty.mc[i]
              }else{ ## AK: Normal distribution
                if (k==1){
                res<-.Fortran("frailpred_sha_nor_mc",
                            as.integer(fit$npar),
                             frail.out=as.double(0),
                             as.double(sigma2.mc[i]),
                             as.double(res1),
                             as.integer(mi),
                             PACKAGE = "frailtypack" )
                frailty.mc[i] <- res$frail.out
                }
                pred <- 1-(vect.survival.X.horizon/vect.survival.X)**exp(frailty.mc[i])
              }

                                                        realisations <- cbind(realisations,pred)
                                                }
                                                predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
                                                predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
                                        }
                                }
                        }else{ ## for a Cox model
                                for (k in 1:nrow(data)){
                                        realisations <- NULL
                                        for (i in 1:MC.sample){
                                                vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                                                vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                                                pred <- 1-(vect.survival.X.horizon/vect.survival.X)
                                                realisations <- rbind(realisations,pred)
                                        }
                                        predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
                                        predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
                                }
                        }
                }

                out <- NULL
                out$call <- match.call()
    out$name.fit <- match.call()[[2]]

                out$npred <- nrow(data)
                out$moving.window <- moving.window
                if (moving.window){
                        out$x.time <- sequence2
                        out$t <- predTime
                }else{
                        out$x.time <- sequence
                }

                out$pred <- predMat

    colnames(out$pred) <- c("times",rep(" ",dim(out$pred)[2]-1))
                rownames(out$pred) <- paste("ind",1:out$npred)

                out$icproba <- ICproba
                if (ICproba){
                        out$predLow <- predMatLow
                        out$predHigh <- predMatHigh
                        colnames(out$predLow) <- c("times",rep(" ",dim(out$predLow)[2]-1))
                        rownames(out$predLow) <- paste("ind",1:out$npred)
                        colnames(out$predHigh) <- c("times",rep(" ",dim(out$predHigh)[2]-1))
                        rownames(out$predHigh) <- paste("ind",1:out$npred)
                }
                if (fit$Frailty) out$type <- type
                out$window <- window
                if (type == "conditional") out$group <- group

                cat("Predictions done for",nrow(data),"subjects \n")

                class(out) <- "predFrailty"
        }
        out
} 