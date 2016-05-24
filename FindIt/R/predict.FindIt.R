predict.FindIt<-
    function(object,newdata,sort=TRUE,decreasing =TRUE,
             wts=1,unique=FALSE,...){
        treat.type   <- object$treat.type
        type <- object$type
        threshold <- object$threshold
        make.twoway <- object$make.twoway
        make.allway  <- object$make.allway
        main <- object$main
	nway <- object$nway
        model.treat <- object$model.treat
        if(main){
            model.main <- object$model.main
            model.int  <- object$model.int 
        }

        if(missing(newdata)){
            different <- 0
        }
        if(!missing(newdata)){
            different <- 1
            newdata <- newdata
        }
        
        if(different==0){
            coefs2       <- object$coefs.orig
            y            <- object$y
            treat        <- object$treat
            y.orig     <- object$y.orig
            treat.orig   <- object$treat.orig
            X.t          <- object$X.t
            if(main){
                X.c          <- object$X.c        
                X.int        <- object$X.int
                X.int.orig <- object$X.int.orig
                scale.c      <- object$scale.c
                scale.int    <- object$scale.int
                X.c.orig   <- object$X.c.orig
            }                      
        }
        if(different==1){
            terms.treat     <- terms(model.treat)
            terms.treat2 <- delete.response(terms.treat)
            treat.frame     <- model.frame(terms.treat2,data=newdata, na.action=NULL) 
            if(main){
                terms.main  <- terms(model.main)
                terms.int   <- terms(model.int)
                main.frame      <- model.frame(terms.main,data=newdata,na.action=NULL)
                int.frame       <- model.frame(terms.int,data=newdata,na.action=NULL)
            }
            ## NA manipulation: List-wise deletion
            data.which <- c()
            if(main){
                for(i in 1:nrow(treat.frame)){
                    if(any(is.na(main.frame[i,])) | any(is.na(int.frame[i,])) |
                       any(is.na(treat.frame[i,])) )
                        {data.which[i] <- 0}
                    else{data.which[i] <- 1}
                }
            }
            if(main==FALSE){
                for(i in 1:nrow(treat.frame)){
                    if(any(is.na(treat.frame[i,])))
                        {data.which[i] <- 0}
                    else{data.which[i] <- 1}
                }
            }

            treat.frame <- treat.frame[data.which==1,]
            treat <- treat.orig <- treat.frame
            
            ## y     <- y.orig      <- treat.frame[,1]
            ## if(type=="binary"){
            ##     y <- 2*y-1
            ## }
            ## y <- y*wts^.5
            
            if(main==TRUE){
                X.c   <- main.frame[data.which==1,]
                X.int <- int.frame[data.which==1,]
                X.int.orig <- X.int
                
                if(make.twoway==TRUE){
                    frame.meanPreC <- object$frame.meanC
                    XC      <- maketwoway(X.c,wts=wts,center=TRUE,frame.meanPre=frame.meanPreC,predict=TRUE)
                    X.c     <- XC$X
                    scale.c <- XC$scale.X
                    frame.meanPreInt <- object$frame.meanInt
                    XInt    <- maketwoway(X.int,wts=wts,center=TRUE,frame.meanPre=frame.meanPreInt,predict=TRUE)
                    X.int   <- XInt$X
                    scale.int <- XInt$scale.X
                    X.c <- as.matrix(X.c)
                    X.int <- as.matrix(X.int)

                    ## print("check after fix")
                    ## print(dim(X.c))
                    ## print(colnames(X.c))
                    ## print(dim(X.int))
                    ## print(colnames(X.c))
                }
                
                
                ## Remove all elements with no variance and add an intercept.
                ## X.c<-X.c[,apply(X.c,2,sd)>0]
                ## in the prediction model, it is allowed.
                X.c.m <- X.c
                ## X.c<-unique(X.c,MARGIN=2)
                scale.c <- scale.c[colnames(X.c.m)==colnames(X.c)]
                if(all(X.c[,1]==1)==FALSE){
                    X.c<-cbind(1,X.c)
                    colnames(X.c)[1]<-"Intercept"
                }
                X.c <- X.c*wts^.5
                ## X.int <- X.int*wts^.5  

                if(treat.type=="single"){
                    ## Remove all elements with no variance and add an intercept.
                    ## X.int<-X.int[,apply(X.int,2,sd)>0]
                    ## in the prediction, it is OK.
                    X.int.m <- X.int
                    ## X.int<-unique(X.int,MARGIN=2)
                    scale.int <- scale.int[colnames(X.int.m)==colnames(X.int)] 
                    if(all(X.int[,1]==1)==FALSE){
                        X.int<-cbind(1,X.int)
                        colnames(X.int)[1]<-"Intercept"
                    }       
                    X.t<-cbind(treat*1e4,(treat>0)*X.int[,-1])
                    X.t[treat!=0,-1]<-
                        apply(X.t[treat!=0,-1],2,
                              FUN=function(x) x-mean(x))
                    colnames(X.t)<-c("treat",paste("treat",
                                                   colnames(X.int)[-1],
                                                   sep=":"))
                }
                
            }
            
            if(treat.type=="multiple"){
                if(make.allway==TRUE){
                    treat <- as.matrix(treat)
                    Allway <- makeallway(treat,threshold,
                                         make.reference=FALSE,
                                         nway=nway)
                    X.t    <- Allway$FinalData
                    X.t    <- as.matrix(X.t)
                    reference.main <- Allway$reference
                }
                if(make.allway==FALSE){
                    X.t    <- as.matrix(treat)
                    ## reference.main <- reference.main
                    reference.main <- "No Reference"
                }
            }
            X.t <- X.t*wts^.5
        }

        
        ## Adjust the coefficients.
        ## When we have smaller variation in the new data.
        if(different==1){
            if(treat.type=="single"){
                coefs.orig         <- object$coefs.orig
                names(coefs.orig)  <- c(object$name.c,"treat",object$name.int[-1])      
                coefs.orig.c   <- coefs.orig[1:ncol(object$X.c)]
                coefs.orig.int <- coefs.orig[-c(1:ncol(object$X.c))]                  
                ## X.c include intercept.
                ## I use the same formula, so the order is the same.

                ## above is the original code
                ## X.c include intercept.
                ## I use the same formula, so the order is the same.
                coefs.orig.c    <- coefs.orig.c[is.element(object$name.c,
                                                           colnames(X.c))]
                coefs.orig.int    <- coefs.orig.int[is.element(object$name.int,
                                                               colnames(X.int))]

                
                ## For coefs.orig.int
                col2.int <- colnames(X.int)
                col2.int <- gsub(":",".",col2.int)
                col2.int <- gsub(" ",".",col2.int)
                col2.int[1] <- "treat"
                coef.int.name <- names(coefs.orig.int)
                coef.int.name <- gsub(":",".",coef.int.name)
                coef.int.name <- gsub(" ",".",coef.int.name)
                ##Coeffcicients Assignment
                equal.int <- (col2.int==coef.int.name)
                search.int <- seq(1:length(coefs.orig.int))[-which(equal.int)]
                coefs.orig.int.new <- c()

                for(i in 1:length(coefs.orig.int)){
                    if(equal.int[i]==TRUE){
                        coefs.orig.int.new[i] <- coefs.orig.int[i] 
                    }else{
                        for(j in search.int){
                            if(col2.int[i]==coef.int.name[j]){
                                coefs.orig.int.new[i] <- coefs.orig.int[j]
                            }
                        }
                    }
                }
                names(coefs.orig.int.new) <- colnames(X.int)
                names(coefs.orig.int.new)[1] <- "treat" 

                ## For coefs.orig.c
                col2.c <- colnames(X.c)
                col2.c <- gsub(":",".",col2.c)
                col2.c <- gsub(" ",".",col2.c)
                coef.c.name <- names(coefs.orig.c)
                coef.c.name <- gsub(":",".",coef.c.name)
                coef.c.name <- gsub(" ",".",coef.c.name)
                ##Coeffcicients Assignment
                equal.c <- (col2.c==coef.c.name)
                search.c <- seq(1:length(coefs.orig.c))[-which(equal.c)]
                coefs.orig.c.new <- c()                   

                for(i in 1:length(coefs.orig.c)){
                    if(equal.c[i]==TRUE){
                        coefs.orig.c.new[i] <- coefs.orig.c[i] 
                    }else{
                        for(j in search.c){
                            if(col2.c[i]==coef.c.name[j]){
                                coefs.orig.c.new[i] <- coefs.orig.c[j]
                            }
                        }
                    }
                }
                names(coefs.orig.c.new) <- colnames(X.c)
                
                coefs2    <- c(coefs.orig.c.new, coefs.orig.int.new)

                ## print("Check int")
                ## coefs.orig.int.check <- coefs.orig.int[order(coef.int.name)]
                ## coefs.int2.check <- coefs.orig.int.new[order(names(coefs.orig.int.new))]
                ## print(coefs.orig.int.check==coefs.int2.check)
                ## break

                ## print("Check c")
                ## print(coefs.orig.c)
                ## print(coefs.orig.c.new)
                ## coefs.orig.c.check <- coefs.orig.c[order(coef.c.name)]
                ## coefs.c2.check <- coefs.orig.c.new[order(names(coefs.orig.c.new))]
                ## print(coefs.orig.c.check==coefs.c2.check)
                ## break

            }
            if(treat.type=="multiple"){
                coefs.orig <- object$coefs.orig
                
                if(main==FALSE){
                    names(coefs.orig)  <- c("Intercept", object$name.t)
                    coefs.orig.t <- coefs.orig[-1]
                    coefs.orig.t <- coefs.orig.t[is.element(object$name.t,
                                                            colnames(X.t))]
                    col2 <- colnames(X.t)
                    col2 <- gsub(":",".",col2)
                    col2 <- gsub(" ",".",col2)
                    coef.t.name <- names(coefs.orig.t)
                    coef.t.name <- gsub(":",".",coef.t.name)
                    coef.t.name <- gsub(" ",".",coef.t.name)
                    ##Coefficients Assignment
                    equal <- (col2==coef.t.name)
                    search <- seq(1:length(coefs.orig.t))[-which(equal)]
                    coefs.orig.new <- c()

                    for(i in 1:length(coefs.orig.t)){
                        if(equal[i]==TRUE){
                            coefs.orig.new[i] <- coefs.orig.t[i] 
                        }else{
                            for(j in search){
                                if(col2[i]==coef.t.name[j]){
                                    coefs.orig.new[i] <- coefs.orig.t[j]
                                }
                            }
                        }
                    }
                    
                    coefs2s  <- c(coefs.orig[1], coefs.orig.new)
                    coefs.t2 <- coefs.orig.new
                    names(coefs.t2) <- colnames(X.t)

                    scale.out <- c(rep(1,length(coefs2s)))
                    coefs2 <- coefs2s * scale.out

                    ## coefs.orig.t.check <- coefs.orig.t[order(coef.t.name)]
                    ## coefs.t2.check <- coefs.t2[order(names(coefs.t2))]
                    ## print(coefs.orig.t.check==coefs.t2.check)
                    
                    names(coefs2) <- c("Intercept", colnames(X.t))
                }                
                if(main){
                    names(coefs.orig)  <- c(object$name.c,object$name.t)
                    coefs.orig.c   <- coefs.orig[1:ncol(object$X.c)]
                    coefs.orig.t   <- coefs.orig[-c(1:ncol(object$X.c))]              
                    ## X.c include intercept.
                    ## I use the same formula, so the order is the same.
                    coefs.orig.c    <- coefs.orig.c[is.element(object$name.c,
                                                               colnames(X.c))]
                    coefs.orig.t    <- coefs.orig.t[is.element(object$name.t,
                                                               colnames(X.t))]
                    
                    ## For coefs.orig.t
                    col2.t <- colnames(X.t)
                    col2.t <- gsub(":",".",col2.t)
                    col2.t <- gsub(" ",".",col2.t)
                    coef.t.name <- names(coefs.orig.t)
                    coef.t.name <- gsub(":",".",coef.t.name)
                    coef.t.name <- gsub(" ",".",coef.t.name)
                    ##Coeffcicients Assignment
                    equal.t <- (col2.t==coef.t.name)
                    search.t <- seq(1:length(coefs.orig.t))[-which(equal.t)]
                    coefs.orig.t.new <- c()

                    for(i in 1:length(coefs.orig.t)){
                        if(equal.t[i]==TRUE){
                            coefs.orig.t.new[i] <- coefs.orig.t[i] 
                        }else{
                            for(j in search.t){
                                if(col2.t[i]==coef.t.name[j]){
                                    coefs.orig.t.new[i] <- coefs.orig.t[j]
                                }
                            }
                        }
                    }
                    names(coefs.orig.t.new) <- colnames(X.t)

                    ## For coefs.orig.c
                    col2.c <- colnames(X.c)
                    col2.c <- gsub(":",".",col2.c)
                    col2.c <- gsub(" ",".",col2.c)
                    coef.c.name <- names(coefs.orig.c)
                    coef.c.name <- gsub(":",".",coef.c.name)
                    coef.c.name <- gsub(" ",".",coef.c.name)
                    ##Coeffcicients Assignment
                    equal.c <- (col2.c==coef.c.name)
                    search.c <- seq(1:length(coefs.orig.c))[-which(equal.c)]
                    coefs.orig.c.new <- c()                   

                    for(i in 1:length(coefs.orig.c)){
                        if(equal.c[i]==TRUE){
                            coefs.orig.c.new[i] <- coefs.orig.c[i] 
                        }else{
                            for(j in search.c){
                                if(col2.c[i]==coef.c.name[j]){
                                    coefs.orig.c.new[i] <- coefs.orig.c[j]
                                }
                            }
                        }
                    }
                    names(coefs.orig.c.new) <- colnames(X.c)
                                     
                    coefs2s    <- c(coefs.orig.c.new, coefs.orig.t.new)
                    scale.out <-c(1,scale.c,rep(1,dim(X.t)[2]))
                    coefs2 <- coefs2s * scale.out
                    
                    coefs.t2 <- coefs.orig.t.new
                    names(coefs.t2) <- colnames(X.t)

                    ## print("Check t")
                    ## coefs.orig.t.check <- coefs.orig.t[order(coef.t.name)]
                    ## coefs.t2.check <- coefs.orig.t.new[order(names(coefs.orig.t.new))]
                    ## print(coefs.orig.t.check==coefs.t2.check)


                    ## print("Check c")
                    ## print(coefs.orig.c)
                    ## print(coefs.orig.c.new)
                    ## coefs.orig.c.check <- coefs.orig.c[order(coef.c.name)]
                    ## coefs.c2.check <- coefs.orig.c.new[order(names(coefs.orig.c.new))]
                    ## print(coefs.orig.c.check==coefs.c2.check)
                    
                }
                X.t <- X.t[,is.element(colnames(X.t), names(coefs.t2))]
            }
        }

        if(main){
            if(any(is.element(colnames(X.c),object$name.c)==FALSE)){
                warning("The new data set has more variation
                         than the original model.")
            }
            
            if(any(is.element(colnames(X.int), object$name.int)==FALSE)){
                warning(
                    "The new data set has more variation
                     than the original model."
                    )
            }            
        }

        if(any(is.element(colnames(X.t), object$name.t)==FALSE)){
            warning(
                "The new data set has more variation than the original model."
                )
        }

        
        if(treat.type=="single"){        
            scale.out<-c(1,scale.c,1e-4,scale.int)
            coefs <- coefs2 * scale.out            
            ## Generate outcome treatment and no treatment.
            preds.treat<-cbind(X.c,cbind(1e4,X.int[,-1]))%*%coefs
            ## preds.treat<-cbind(X.c,cbind(1e4,X.c[,-1]))%*%coefs
            preds.control<-cbind(X.c,0*cbind(1e4,X.int[,-1]))%*%coefs
            ## preds.control<-cbind(X.c,0*cbind(1e4,X.c[,-1]))%*%coefs
            preds <- cbind(X.c,X.t)%*%coefs
            
            if(type=="binary"){
                preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                preds.control <- sign(preds.control)*pmin(abs(preds.control),1)
                preds         <- sign(preds)*pmin(abs(preds),1)
                ATE <- mean(preds.treat-preds.control)/2
                preds.diff<- (preds.treat-preds.control)/2
            }
            if(type=="continuous"){
                preds     <- preds
                preds.diff<-preds.treat-preds.control
                ATE<-mean(preds.treat-preds.control)
            }

            if(different==0){
                pred.data <- cbind(preds.diff,y.orig,treat.orig,X.int.orig)
                pred.data <- as.data.frame(pred.data)
                colnames(pred.data) <- c("Treatment.effect",
                                         "outcome",
                                         "treatment",
                                         colnames(X.int.orig))
            }
            if(different==1){
                pred.data <- cbind(preds.diff,treat.orig,X.int.orig)
                pred.data <- as.data.frame(pred.data)
                colnames(pred.data) <- c("Treatment.effect",
                                         "treatment",
                                         colnames(X.int.orig))
            }
            
            if(sort==TRUE){
                pred.data.out <- pred.data[order(pred.data$Treatment.effect,
                                                 decreasing=decreasing),]
            }
            else{pred.data.out <- pred.data}
            
        }
        if(treat.type=="multiple"){
            if(different==0){
                if(main){
                    scale.out<-c(1,scale.c,rep(1,dim(X.t)[2]))
                    coefs <- coefs2 * scale.out
                }
                if(main==FALSE){
                    scale.out <- c(1,rep(1,ncol(X.t)))
                    coefs <- coefs2 * scale.out
                }
            }else{
                coefs <- coefs2
            }
            
            if(type=="binary"){
                if(main){
                    X.t1 <- as.matrix(X.t)
                    preds.treat    <- cbind(X.c,X.t1) %*% coefs
                    preds.control  <- X.c %*% coefs[c(1:ncol(X.c))]
                    preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                    preds.control <- sign(preds.control)*pmin(abs(preds.control),1)
                    preds.diff <- (preds.treat-preds.control)/2
                }
                if(main==FALSE){
                    preds.treat <- X.t %*% coefs[-1]
                    preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                    preds.diff <- preds.treat/2
                }
                ATE <- mean(preds.diff)
            }
            if(type=="continuous"){
                if(main){
                    preds.diff    <- X.t %*% coefs[-c(1:(dim(X.c)[2]))]
                }
                if(main==FALSE){
                    preds.diff <- X.t %*% coefs[-1]
                }
                ATE <- mean(preds.diff)
            }

            ## Make the SATE for Multiple treatment version.            
            if(unique==FALSE){
                ##preds.control <- X.c %*% coefs[1:ncol(X.c)]
                ##preds <- cbind(X.c,X.t)%*%coefs
                ##preds.diff <- preds - preds.control
                preds.diff <- preds.diff
                pred.data <- cbind(preds.diff,treat.orig)
                colnames(pred.data) <- c("Treatment.effect",                          
                                         colnames(treat.orig)
                                         )
                pred.data <- as.data.frame(pred.data)
                
                if(sort==TRUE){
                    pred.data.out <- pred.data[order(pred.data$Treatment.effect,
                                                     decreasing=decreasing),]
                }
                else{pred.data.out <- pred.data}
                
            }
            if(unique==TRUE){
                X.t <- as.matrix(X.t)
                X.t.u <- X.t
                preds.diff.u <- as.data.frame(preds.diff)
                rownames(preds.diff.u) <- rownames(X.t.u) <- seq(1:nrow(X.t.u))
                treat.unique  <- unique(treat.orig,MARGIN=1)
                X.t.unique    <- unique(X.t.u, MARGIN=1)
                X.t.unique2   <- unique(X.t, MARGIN=1)
                preds.diff   <- preds.diff.u[rownames(preds.diff.u) %in%
                                             rownames(X.t.unique),]
                
                preds.diff <- as.data.frame(preds.diff)
                ##name back
                rownames(X.t.unique) <- rownames(X.t.unique2)
                rownames(preds.diff) <- rownames(X.t.unique2)

                pred.data <- cbind(preds.diff,treat.unique)
                colnames(pred.data) <- c("Treatment.effect",
                                         colnames(treat.orig)
                                         )
                pred.data <- as.data.frame(pred.data)
                pred.data.out <-
                    pred.data[order(pred.data$Treatment.effect,
                                    decreasing=decreasing
                                    ),]
            }
            
        }
        pred.data.out <- as.data.frame(pred.data.out)
        out <- list("treat.type"=treat.type, "ATE"=ATE, "data" =pred.data.out,
                    "coefs"=coefs,"orig.coef"=coefs2)
        class(out) <- "PredictFindIt"
        invisible(out)
    }
