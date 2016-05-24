## This file contains the functions used for analyzing causal heterogeneity
## using support vector machines under a double LASSO constraint.

## The functions consist of 
## scale.func: demeans and standardizes variables.
## scale.func.2: gives variables a standard deviation of one
## FindIt: a wrapper function for SVM.func
## maketwoway<-rescales covariates and makes all twoway interactions.
## SVM.func: the workhorse function.  Given values of y, control, treatment, and lambdas, 
## it returns betas, the loss (GCV), fitted values, and size of the margin.
## find.lambda: Conducts the alternating line search.
## A large value of lambda.c is selected (10), and alternating
## line searches of decreasing radius are conducted.
## lars: An adaptation of the LARS algorithm of Efron, et al. 2004.
## Due to the large numbers being considered in the fitting,
## the code is adapted to not crash.  The original code is commented out.
## The remainder of the code is untouched.

##################################################################
##################################################################
##################################################################

scale.func<-function(x) {
    x2<-x-mean(x)
    x2/sd(x2)}

scale.func.2<-function(x) {x2<-x
                           x2*sd(x2)}

#####################################################

##change here. 
FindIt <- function(model.treat, model.main,model.int,data=NULL,
                   type="binary",
                   treat.type="multiple",
                   nway,
                   search.lambdas=TRUE,lambdas=NULL,
                   make.twoway=TRUE,make.allway=TRUE,
                   wts=1,scale.c=1,scale.int=1,
                   fit.glmnet=TRUE,
                   make.reference=TRUE,
                   reference.main=NULL,                   
                   threshold=0.999999){

    ## if(make.allway==FALSE & is.null(reference.main)){
    ##     warning("Need to specify the corresponding reference of treatment matrix")
    ##     break
    ## }
    unique <- FALSE

    ## Check nway
    if(missing(nway) & treat.type=="multiple") {
        print("Need to specify nway")
        break
    }
    if(missing(nway) & treat.type=="single") {
        nway <- 2
    }

    if(sum(wts<0)>0) {
        print("Only non-negative weights allowed")
        break
    }

    ## Check whether the model has the main effect covariates or not. 
    if(!missing(model.main)|!missing(model.int)){
        main <- TRUE
    }    
    if(missing(model.int) & !missing(model.main)){
        model.int <- model.main
        main <- TRUE
    }    
    if(!missing(model.int) & missing(model.main)){
        model.main <- model.int
        main <- TRUE
    }    
    if(missing(model.main) & missing(model.int)){
        main <- FALSE
    }

    ## Extract the data
    treat.frame     <- model.frame(model.treat,data, na.action=NULL)
    if(main){
        main.frame  <- model.frame(model.main, data, na.action=NULL)
        int.frame   <- model.frame(model.int,  data, na.action=NULL)
    }

    ## NA : List-wise deletion.
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
    treat.orig <- treat.frame[,-1]
    treat <- treat.frame[,-1]
    treat <- as.matrix(treat)
    y     <- y.orig      <- treat.frame[,1]
    if(type=="binary"){
        y.orig <- y
        y <-(2*y-1)        
    }
    
    if(main){
        X.c   <- X.c.orig    <- main.frame[data.which==1,]
        X.int <- X.int.orig  <- int.frame[data.which==1,]
        X.c   <- X.c.orig <- as.data.frame(X.c)
        X.int <- X.int.orig <- as.data.frame(X.int)

        if(ncol(X.c)>1){
            if(any(apply(X.c,2,FUN=function(x) length(unique(x)) == 1))){
                warning("There is a column with no variation
                         in the main effect matrix.")
            }
        }
        if(ncol(X.int)>1){
            if(any(apply(X.int,2,FUN=function(x) length(unique(x)) == 1))){
                warning("There is a column with no variation
                         in the interaction matrix.")
            }   
        }
        if(treat.type=="single" & main==FALSE){
            warning("Need to specify the covariates to interact with treatment.")
        }
    }

    ## Make multiple-treatment matrix
    if(treat.type == "multiple"){
        if(make.allway==TRUE){
            Allway <- makeallway(treat,threshold,make.reference=make.reference,nway=nway)
            X.t    <- Allway$FinalData
            reference.main <- Allway$reference
        }
        if(make.allway==FALSE){
            X.t    <- as.matrix(treat)
            ## reference.main <- reference.main
            reference.main <- "No reference"
        }
    }

    ## print("X.t")
    ## print(head(X.t))
    ## print(head(X.c))
    ## print(head(X.int))

    ## Maketwoway    
    if(main){
        X.c.old <- X.c
        X.int.old <- X.int
        if(make.twoway==TRUE){
            if(ncol(X.c.old)>1){
                XC         <- maketwoway(X.c,center=TRUE,wts=wts)
                X.c        <- XC$X
                scale.c    <- XC$scale.X
                XC.deleted <- XC$D
                frame.meanC <- XC$frame.mean
            }
            if(ncol(X.c.old)==1){
                form <- ~ .
                X.c.old <- as.data.frame(X.c.old)
                X.c1 <- model.matrix(form,data=X.c.old)
                X.c1 <- as.data.frame(X.c1)
                X.c <- X.c1[,-1]
                if(is.vector(X.c)){
                    scale.c <- sd(X.c)
                    frame.meanC <- mean(X.c)
                    X.c <- (X.c-mean(wts^.5*X.c))/sd(X.c)    
                    X.c <- as.data.frame(X.c)
                    colnames(X.c) <- colnames(X.c1)[-1]                   
                }else{
                    scale.c <- apply(X.c,2,sd)
                    frame.meanC <- apply(X.c,2,function(x) mean(wts^.5*x))
                    X.c <- apply(X.c,2,FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                }
            }
            if(ncol(X.int.old)>1){
                XInt       <- maketwoway(X.int,center=TRUE,wts=wts)
                X.int      <- XInt$X
                scale.int  <- XInt$scale.X
                XInt.deleted <- XInt$D
                frame.meanInt <- XInt$frame.mean
            }
            if(ncol(X.int.old)==1){
                form <- ~ .
                X.int.old<- as.data.frame(X.int.old)
                X.int1 <- model.matrix(form,data=X.int.old)
                X.int1 <- as.data.frame(X.int1)
                X.int <- X.int1[,-1]
                if(is.vector(X.int)){
                    scale.int <- sd(X.int)
                    frame.meanInt <- mean(X.int)
                    X.int <- (X.int-mean(wts^.5*X.int))/sd(X.int) 
                    X.int <- as.data.frame(X.int)
                    colnames(X.int) <- colnames(X.int1)[-1]
                }
                else{
                    scale.int <- apply(X.int,2,sd)
                    frame.meanInt <- apply(X.int,2,function(x) mean(wts^.5*x))
                    X.int <- apply(X.int,2,FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                }
            }
        }
        
        ## Remove all elements with no variance and add an intercept.
        if(is.matrix(X.c)){
            scale.c <- scale.c[apply(X.c,2,sd)>0]
            X.c<-X.c.1 <- X.c[,apply(X.c,2,sd)>0]
            ## X.c<-unique(X.c,MARGIN=2)
            scale.c <- scale.c[colnames(X.c.1)==colnames(X.c)]
            frame.meanC <- frame.meanC[colnames(X.c.1)==colnames(X.c)]
        }
        if(is.vector(X.c)){
            X.c<-cbind(1,X.c)
            colnames(X.c)[1]<-"Intercept"
        }
        else{
            if(sd(X.c[,1])>0){
                X.c<-cbind(1,X.c)
                colnames(X.c)[1]<-"Intercept"
            }
        }
        name.c <- colnames(X.c)
        ## This X.c is used in the SVM.
        
        ## Remove all elements with no variance and add an intercept.
        if(is.matrix(X.int)){
            scale.int <- scale.int[apply(X.int,2,sd)>0]
            frame.meanInt <- frame.meanInt[apply(X.int,2,sd)>0]
            X.int<- X.int.1<- X.int[,apply(X.int,2,sd)>0]
            scale.int <- scale.int[colnames(X.int.1)==colnames(X.int)]
            frame.meanInt <- frame.meanInt[colnames(X.int.1)==colnames(X.int)]
        }
        if(is.vector(X.int)){
            X.int<-cbind(1,X.int)
            colnames(X.int)[1]<-"Intercept"
        }
        else{
            if(sd(X.int[,1])>0){
                X.int<-cbind(1,X.int)
                colnames(X.int)[1]<-"Intercept"
            }
        }
        name.int <- colnames(X.int)
        ## This X.int is used in the SVM

        ## Make X.t
        ## I changed X.c to X.int        
        if(treat.type == "single") {
            X.t<-cbind(treat*1e4,(treat>0)*X.int[,-1])
            X.t[treat!=0,-1]<-apply(X.t[treat!=0,-1],2, FUN=function(x) x-mean(x))
            ## This X.t is used.            
            colnames(X.t)<-c("treat",paste("treat",colnames(X.int)[-1],sep=":"))
        }
    }

###################################
    ## Create the FindIt.fit function.
################################### 

    FindIt.fit<-function(y,X.c,X.t,X.int,treat,treat.type,                     
                         search.lambdas,lambdas,
                         wts,scale.c,scale.int,
                         reference.main,
                         fit.glmnet,
                         type,
                         X.c.orig,X.int.orig,treat.orig,
                         threshold,
                         unique,
                         main){

        ## Give wts mean 1.
	wts<-wts/mean(wts)
        
	y <- y*wts^.5
        X.t <- X.t*wts^.5
        name.t <- colnames(X.t)
        if(main){
            X.c <- X.c*wts^.5
            ## X.int <- X.int*wts^.5
            name.c <- colnames(X.c)            
            name.int <- colnames(X.int)           
        }
        ## treat<-treat*wts^.5

        

        X.t <- as.matrix(X.t)
        if(main){
            X.c <- as.matrix(X.c)
            X.int <- as.matrix(X.int)
        }
	
	if(search.lambdas==TRUE){
            if(main){
                lambdas <-search.lambda(y,X.c,X.t,
                                        fit.glmnet=fit.glmnet,
                                        main=main,
                                        type=type)
            }
            if(main==FALSE){
                lambdas <- search.lambda.nocov(y=y,X.t=X.t,
                                               fit.glmnet=fit.glmnet,
                                               main=main,
                                               type=type)
            }
        }
        

        if(main){
            A<-SVM.func(y,X.c[,-1],X.t,
                        lambdas[1],lambdas[2],
                        fit.glmnet=fit.glmnet,
                        type=type,
                        main=main)
        }
        if(main==FALSE){
            A<-SVM.func(y=y,X.c[,-1],X.t=X.t,
                        lambda.c=1,lambda.t=lambdas,
                        fit.glmnet=fit.glmnet,
                        type=type,
                        main=main)
        }
        ## When we input X.c doesnt have intercept. 
        
	if(treat.type=="single"){
            ## Keep the names of coefficients.
            reference.coef <- as.data.frame(A$beta[-c(1:(dim(X.c)[2]))])
            rownames(reference.coef) <- c("treat",colnames(X.int[,-1]))
            if(type=="binary"){
                reference.coef <- reference.coef/2
            }
            ## Generate outcome treatment and no treatment.
            ## I changed the code.
            preds.treat<-cbind(X.c,cbind(1e4,X.int[,-1]))%*%A$beta
            preds.control<-cbind(X.c,0*cbind(1e4,X.int[,-1]))%*%A$beta
            ## preds <- cbind(X.c,X.t) %*% A$beta
            ## diff  <- mean(y.inter)- mean(preds)
            
            ## cbind(1e4,X.int[,-1]) is X.t when treat =1.
            ## preds.treat<-cbind(X.c,cbind(1e4,X.c[,-1]))%*%A$beta
            ## in here, X.c includes the intercept.
            ## 1e4 means treated group         
            ## preds.control<-cbind(X.c,0*cbind(1e4,X.c[,-1]))%*%A$beta

            if(type=="binary"){
                preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                preds.control <- sign(preds.control)*pmin(abs(preds.control),1)
                ATE <- mean(preds.treat-preds.control)/2
                preds.diff<- (preds.treat-preds.control)/2
            }
            if(type=="continuous"){
                preds.diff<-preds.treat-preds.control
                ATE<-mean(preds.treat-preds.control)
            }

            scale.out<-c(1,scale.c,1e-4,scale.int)
            Treatment.version <- NULL
            
	}
        if(treat.type == "multiple"){
            ## make the reference matrix.
            if(main){
                ## preds <- cbind(X.c,X.t) %*% A$beta
                ## diff <- mean(y.inter) - mean(preds)
                scale.out<-c(1,scale.c,rep(1,dim(X.t)[2]))
                coef.orig.ref <- A$beta[-c(1:ncol(X.c))]/scale.out[-c(1:ncol(X.c))]
            }
            if(main==FALSE){
                ## preds <- cbind(1,X.t) %*% A$beta
                ## diff <- mean(y.inter) - mean(preds)
                scale.out <- c(1,rep(1,ncol(X.t)))
                ## I need one more for intercept.
                coef.orig.ref <- A$beta[-1]/scale.out[-1]
            }
            
            if(type=="binary"){
                reference.coef <- as.data.frame(signif(coef.orig.ref/2,6))
            }
            if(type=="continuous"){
                reference.coef <- as.data.frame(signif(coef.orig.ref,6))
            }
            
            rownames(reference.coef) <- colnames(X.t)

###############################
            ## create the Treatment.version
###############################
            coefs <- A$beta
            
            if(type=="binary"){
                if(main){
                    preds.treat    <- cbind(X.c,X.t) %*% coefs
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
            
            if(unique==FALSE){
                ##preds.control <- X.c %*% coefs[1:ncol(X.c)]
                ##preds <- cbind(X.c,X.t)%*%coefs
                ##preds.diff <- preds - preds.control
                preds.diff <- preds.diff
                Treatment.version <- cbind(preds.diff,treat.orig)
                colnames(Treatment.version) <- c("Treatment.effect",
                                                 colnames(treat.orig))
                
                Treatment.version <- as.data.frame(Treatment.version)
                
            }
            if(unique==TRUE){
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

                Treatment.version <- cbind(preds.diff,treat.unique)
                colnames(Treatment.version) <- c("Treatment.effect",
                                                 colnames(treat.orig))
                
                Treatment.version <- as.data.frame(Treatment.version)
                Treatment.version<-
                    Treatment.version[order(Treatment.version$Treatment.effect,
                                            decreasing=TRUE),]
            }            
	}
        
        if(main){
            names.out <- c(colnames(X.c),colnames(X.t))
        }
        else{
            names.out <- c("Intercept",colnames(X.t))
        }

        if(main){
            output<-(list("coefs"=A$beta,"coefs.orig"=A$beta/scale.out,
                          "fit"=A$fit,
                          "scale.out"=scale.out,
                          "scale.c" = scale.c,
                          "scale.int"=scale.int,
                          "frame.meanC" = frame.meanC,
                          "frame.meanInt" = frame.meanInt,
                          "names.out"=names.out,
                          "y"=y,"X.c"=X.c,"X.t"=X.t,"X.int"=X.int,
                          "GCV"=exp(A$loss),
                          "name.c" =name.c,
                          "name.int" = name.int,
                          "name.t" =name.t,
                          "ATE" = ATE,
                          "Treatment.version" = Treatment.version,
                          "lambdas"=lambdas,
                          "reference.main"=reference.main,
                          "reference.coef"=reference.coef)
                     )
        }
        if(main==FALSE){
            output<-(list("coefs"=coefs,"coefs.orig"=A$beta/scale.out,
                          "fit"=A$fit,
                          "scale.out"=scale.out,
                          "names.out"=names.out,
                          "y"=y,"X.t"=X.t,
                          "GCV"=exp(A$loss),
                          "name.t" =name.t,
                          "ATE" = ATE,
                          "Treatment.version" = Treatment.version,
                          "lambdas"=lambdas,
                          "reference.main"=reference.main,
                          "reference.coef"=reference.coef)
                     )
        }
	class(output)<-c("FindIt","list")
	output
    }

###########################################
    ## The FindIt.fit function is finished
##########################################
    
    if(main==FALSE){
        X.c <- X.int <- scale.c <- scale.int <- frame.meanC <- frame.meanInt <- NULL
        X.c.orig <- X.int.orig <- NULL
    }

    
    Fit.output <- FindIt.fit(y=y,X.c=X.c,X.t=X.t,X.int=X.int,treat=treat,
                             treat.type=treat.type,
                             search.lambdas = search.lambdas,
                             lambdas = lambdas,
                             wts = wts,
                             scale.c = scale.c,
                             scale.int = scale.int,
                             reference.main = reference.main,
                             fit.glmnet = fit.glmnet,
                             type = type,
                             X.c.orig = X.c.orig,
                             X.int.orig = X.int.orig,
                             treat.orig = treat.orig,
                             threshold = threshold,
                             unique = unique,
                             main=main)

#########################################################
#### Make the reference matrix when treat.type is multiple.
#########################################################
    if(make.reference==TRUE & make.allway==TRUE){
        if(treat.type=="multiple"){
            reference.main <- Fit.output$reference.main
            reference.coef <- Fit.output$reference.coef
            coef   <- c()
            
            for(i in 1:nrow(reference.main)){
                if(rownames(reference.main)[i] %in% rownames(reference.coef)){
                    coef[i] <-
                        reference.coef[rownames(reference.coef)==
                                       rownames(reference.main)[i],1]
                }
                else{coef[i] <- NA
                 }
            }
            
            reference.Final <- as.data.frame(cbind(reference.main,coef))
            colnames(reference.Final) <- c(colnames(reference.main),
                                           "Coefficients"
                                           )
        }
    }else{
        reference.Final <- "No reference"
    }
    if(treat.type=="single"){
        reference.Final <- "No reference"
    }

    if(main){
        output<-(list("coefs"=Fit.output$coefs,"coefs.orig"=Fit.output$coefs.orig,
                      "fit"=Fit.output$fit,                     
                      "scale.out"=Fit.output$scale.out,
                      "scale.c" = Fit.output$scale.c,
                      "scale.int" = Fit.output$scale.int,
                      "frame.meanC" = frame.meanC,
                      "frame.meanInt" = frame.meanInt,
                      "names.out"=Fit.output$names.out,
                      "y"=Fit.output$y,"X.c"=Fit.output$X.c,
                      "X.t"=Fit.output$X.t, "X.int"=Fit.output$X.int,
                      "treat" = Fit.output$treat,
                      "GCV"=Fit.output$GCV,
                      "y.orig"=y.orig, "X.c.orig" = X.c.orig,
                      "X.int.orig" = X.int.orig, "treat.orig" = treat.orig,
                      "name.c" = Fit.output$name.c,
                      "name.int" = Fit.output$name.int,
                      "name.t" = Fit.output$name.t,
                      "treat.type" = treat.type,
                      "model.main"=model.main,
                      "model.int" = model.int,
                      "model.treat"=model.treat,
                      "ATE"=Fit.output$ATE,
                      "data"=data,
                      "Treatment.version" = Fit.output$Treatment.version,
                      "lambdas"=Fit.output$lambdas,
                      "threshold" = threshold,
                      "make.twoway" = make.twoway,
                      "make.allway" = make.allway,
                      "model.main"=model.main,
                      "model.int"=model.int, "model.treat"=model.treat,
                      "reference"= reference.Final,"type" = type,
                      "main" =main,
		      "nway" =nway
                      ))
    }
    if(main==FALSE){
        output<-(list("coefs"=Fit.output$coefs,"coefs.orig"=Fit.output$coefs.orig,
                      "fit"=Fit.output$fit,
                      "scale.out"=Fit.output$scale.out,
                      "names.out"=Fit.output$names.out,
                      "y"=Fit.output$y,                    
                      "X.t"=Fit.output$X.t,                     
                      "treat" = Fit.output$treat,
                      "GCV"=Fit.output$GCV,
                      "y.orig"=y.orig,                     
                      "treat.orig" = treat.orig,
                      "name.t" = Fit.output$name.t,
                      "treat.type" = treat.type,
                      "model.treat"=model.treat,
                      "ATE"=Fit.output$ATE,
                      "data"=data,
                      "Treatment.version" = Fit.output$Treatment.version,
                      "lambdas"=Fit.output$lambdas,                      
                      "threshold" = threshold,
                      "make.twoway" = make.twoway,
                      "make.allway" = make.allway,
                      "reference"= reference.Final,
                      "type" = type,
                      "main" =main,
		      "nway" =nway
                      ))
    }
    
    class(output)<-c("FindIt","list")
    output

}

#####################################################
SVM.func<-function(y,X.c,X.t,lambda.c,lambda.t,wts=1,fit.glmnet,type,main){

    ### From the construction, this X.c doesnt include the intercept. 

    n<-length(y)
    if(main){
        X<-X2<-cbind(X.c,X.t)
    }
    else{
        X <- X2 <- X.t 
    }

    ##Originally written to handle a string of lambda.c and lambda.t
    ##loss.run creates a matrix of losses, across values of lambda.
    ## loss.run<-matrix(0,nrow=length(lambda.c),ncol=length(lambda.t))
    ## min.params<-c(lambda.c[1],lambda.t[1])
    ## for(i in 1:length(lambda.c)){
    ## for(j in 1:length(lambda.t)){

    ##Scale the X matrix by lambdas.
    X3<-X2
    if(main){
        if(is.matrix(X.c)){
            X3[,1:dim(X.c)[2]]<-1/exp(lambda.c)*X2[,1:dim(X.c)[2]]
            X3[,-c(1:dim(X.c)[2])]<-1/exp(lambda.t)*X2[,-c(1:dim(X.c)[2])]
        }
        if(is.vector(X.c)){
            X3[,1]<-1/exp(lambda.c)*X2[,1]
            X3[,-1]<-1/exp(lambda.t)*X2[,-1]
        }
    }

    if(main==FALSE){
        X3 <- 1/exp(lambda.t)*X2
    }
    ##X3.c <- 1/exp(lambda.c[i])*X.c
    ##X3.t <- 1/exp(lambda.t[j])*X.t
    ##X3   <- cbind(X.c,X.t)
    

    ##Declare original values.
    X.new<-X3
    which.use<-rep(T,n)
    beta.new<-beta.curr<-rep(0,dim(X2)[2]+1)
    ## this +1 is for the intercept.

    ## The loop to fit the LASSO.
    ## Center on the margin, fit the lasso, update coefficients,
    ##then drop all outside the margin (y*yhat>1)
    if(type=="binary"){
        for(beta.run in 1:100){
            X.new.2<-apply(X.new[which.use,],MARGIN=2,FUN=function(x) x-mean(x))
            X.new.2<-apply(X.new.2,MARGIN=2,
                           FUN=function(x)
                           if(length(unique(x))==1) rnorm(sum(which.use))
                           else x)
            if(fit.glmnet==TRUE){
                glmnet.1<-glmnet(X.new.2,y[which.use]-mean(y[which.use]),
                                 family="gaussian",
                                 lambda=c(5,4,3,2,seq(1.5,1,-.1)),
                                 standardize=F)
                beta.new[-1]<-as.vector(glmnet.1$beta[,10])
            }else{
                ##This was the old L1 optimizer.
                lasso1<-lars(X.new.2,(y[which.use]-mean(y[which.use])),
                             max.steps=dim(X)[2]+2,
                             normalize=F,
                             type="lasso",
                             eps=0)
                beta.new[-1]<-as.vector(predict(lasso1,s=1,type="coef",
                                                mode="lambda")$coef)
                if(log(sum(which.use))/2*sum(beta.new[-1]!=0)>.9*sum(which.use)){
                    beta.new[-1]<-
                        as.vector(predict(lasso1,
                                          s=min(c(floor(.9*sum(which.use)*2/
                                              log(n)),dim(X.new)[2]*.8)),
                                          type="coef",mode="step")$coef)
                }
            }
            beta.curr<-.5*beta.new+.5*beta.curr
            
            beta.new[1]<-mean(y[which.use])-mean(X.new[which.use,]%*%beta.new[-1])
            ## This is for the intercept. This corrects the demean tricks.
            run.diff<-(mean((beta.new[-1]-beta.curr[-1])^2)/sum(beta.new[-1]^2+1e-10))
            if(run.diff<1e-6) break
            which.use<-(y*cbind(1,X.new)%*%beta.new)<=y^2
        }
        
        ##Find fitted values.
        fits<-pmin(abs(cbind(1,X.new)%*%beta.new),y^2)*sign((cbind(1,X.new)%*%beta.new))
        fits2<-cbind(1,X.new)%*%beta.new
        fits2[sign(fits2)==sign(y)]<-pmin(abs(fits2[sign(fits2)==sign(y)]),
                      y[sign(fits2)==sign(y)]^2)*sign(fits2[sign(fits2)==sign(y)])
        
    }
    if(type=="continuous"){
        X.new.2<-apply(X.new,MARGIN=2,FUN=function(x) x-mean(x))
        X.new.2<-apply(X.new.2,MARGIN=2,
                       FUN=function(x) if(length(unique(x))==1) rnorm(sum(which.use)) else x)
        if(fit.glmnet==TRUE){
            glmnet.1<-glmnet(X.new.2,y-mean(y),
                             family="gaussian",
                             lambda=c(5,4,3,2,seq(1.5,1,-.1)),
                             standardize=F)
            beta.new[-1]<-as.vector(glmnet.1$beta[,10])
        }else{
            ##This was the old L1 optimizer.
            lasso1<-lars(X.new.2,(y-mean(y)),
                         max.steps=dim(X)[2]+2,normalize=F,type="lasso",eps=0)
            beta.new[-1]<-as.vector(predict(lasso1,s=1,type="coef",mode="lambda")$coef)
            if(log(length(y))/2*sum(beta.new[-1]!=0)>.9*length(y)) 
                beta.new[-1]<-
                    as.vector(predict(
                        lasso1,
                        s=min(c(floor(.9*length(y)*2/log(n)),dim(X.new)[2]*.8)),
                        type="coef",mode="step")$coef)
        }
        beta.curr<-.5*beta.new+.5*beta.curr
        
        beta.new[1]<-mean(y)-mean(X.new%*%beta.new[-1])
        ## This is for the intercept
        ##run.diff<-(mean((beta.new[-1]-beta.curr[-1])^2)/sum(beta.new[-1]^2+1e-10))
        ##if(run.diff<1e-6) break
        ##which.use<-(y*cbind(1,X.new)%*%beta.new)<=y^2
        ##}
        ## We dont need loop. Because this is not SVM.
        
        ##Find fitted values.
        fits2 <- cbind(1,X.new)%*%beta.new
        fits  <- fits2
        ##fits<-pmin(abs(cbind(1,X.new)%*%beta.new),y^2)*sign((cbind(1,X.new)%*%beta.new))
        ##fits2<-cbind(1,X.new)%*%beta.new
        ##fits2[sign(fits2)==sign(y)]<-pmin(abs(fits2[sign(fits2)==sign(y)]),
        ##              y[sign(fits2)==sign(y)]^2)*sign(fits2[sign(fits2)==sign(y)])
    }
    
    ##Calculate degrees of freedom
    edf<-  1+sum(beta.new[-1]!=0)
    ##*log(n)/2#sum(diag(hat.mat))#+1-sum(beta.new[-c(1:(dim(X.c)[2]+1))]!=0)
    ##Left in for estimates of standard deviation.  	
    stdev<-NA
    
    ##GCV statistic
    if(type=="binary"){
        loss.run<-sum((y^2-y*fits2)^2)/(sum(which.use)-edf)^2*(mean(which.use))^2*n
    }
    if(type=="continuous"){
        loss.run<-sum((y-fits2)^2)/(n*((1-edf/n)^2))
    }
    
    ## Gather minimum loss function.
    ## if(i*j>1) if(loss.run[i,j]==min(loss.run[loss.run!=0]))
    ## min.params<-c(lambda.c[i],lambda.t[j])

    ##     }
    ## }
    
    ##Scale betas back.
    if(main){
        if(is.matrix(X.c)){
            beta.new[c(1:dim(X.c)[2]+1)]<-beta.new[c(1:dim(X.c)[2]+1)]/exp(lambda.c)
            ##for the main covariates and intercepts
            beta.new[-c(1:(dim(X.c)[2]+1))]<-beta.new[-c(1:(dim(X.c)[2]+1))]/exp(lambda.t)
            ##for the interaction covariates
        }
        if(is.vector(X.c)){
            beta.new[c(1:2)]<-beta.new[c(1:2)]/exp(lambda.c)
            ##for the main covariates and intercepts
            beta.new[-c(1:2)]<-beta.new[-c(1:2)]/exp(lambda.t)
            ##for the interaction covariates
        }
    }
    if(main==FALSE){
        beta.new <- beta.new/exp(lambda.t)
    }
    
    beta.new<-as.vector(beta.new)
    X<-as.matrix(X)

    ##Calculate intercept
    if(type=="binary"){
        beta.new[1]<-mean(y[which.use])-mean((X%*%beta.new[-1])[which.use])
    }
    if(type=="continuous"){
        beta.new[1]<-mean(y)-mean((X%*%beta.new[-1]))
    }

    output<-list(##"lambdas"=min.params,
                 "beta"=beta.new,"fits"=fits,
                 "loss"=log(loss.run),##"marg"=mean(which.use),
                 "edf"=edf,"sd"=stdev)
    invisible(output)
}


#####################################################

#####################################################


search.lambda<-function(y=y,X.c=X.c,X.t=X.t,fit.glmnet,main,type){
    lambda.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,
                                                      lambda[1],lambda[2],
                                                      fit.glmnet=fit.glmnet,
                                                      main=main,
                                                      type=type)$loss)
    lambda.c.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,
                                                        lambda,lambda.t,
                                                        fit.glmnet=fit.glmnet,
                                                        main=main,
                                                        type=type)$loss)
    lambda.t.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,
                                                        lambda.c,lambda,
                                                        fit.glmnet=fit.glmnet,
                                                        main=main,
                                                        type=type)$loss)
    lambda.c.old<-lambda.t.old<-999

    lambda.c.seek<-seq(-15,10,1)
    lambda.t.seek<-seq(-15,10,1)
    lambda.t<-25
    lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
    lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
    lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
    lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
    print(c(lambda.c,lambda.t))
    print(range(lambda.c.seek))
    print(range(lambda.t.seek))

    ##lambda.c<- -15
    ##lambda.t<-5
    lambda.c.old<-lambda.c
    lambda.t.old<-lambda.t

    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-5,lambda.c+5,1)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-5,lambda.t+5,1)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }



    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-2.5,lambda.c+2.5,.5)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-2.5,lambda.t+2.5,.5)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-1.5,lambda.c+1.5,.25)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-1.5,lambda.t+1.5,.25)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.5,lambda.c+.5,.1)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.5,lambda.t+.5,.1)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }



    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.25,lambda.c+.25,.05)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.25,lambda.t+.25,.05)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.15,lambda.c+.15,.025)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.15,lambda.t+.15,.025)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.1,lambda.c+.1,.01)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.1,lambda.t+.1,.01)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.05,lambda.c+.05,.005)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.05,lambda.t+.05,.005)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.025,lambda.c+.025,.0025)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.025,lambda.t+.025,.0025)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.01,lambda.c+.01,.001)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.01,lambda.t+.01,.001)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.0025,lambda.c+.0025,.0005)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.0025,lambda.t+.0025,.0005)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }




    
    optim.lambda<-c(lambda.c,lambda.t)


    invisible(optim.lambda)
}

#####################################################
search.lambda.nocov<-function(y=y,X.t=X.t,fit.glmnet,main,type){
    lambda.t.find<-function(lambda) as.numeric(SVM.func(y,X.c=NULL,X.t,
                                                        lambda.c=1,lambda,
                                                        fit.glmnet=fit.glmnet,
                                                        main=FALSE,
                                                        type=type)$loss)
    ## This lambda.c will not be used. 
    ## lambda.c.old<-lambda.t.old<-999
    lambda.t.old<-999

##############################
    ## The Initial Try
############################## 
    ## lambda.c.seek<-seq(-15,10,1)
    lambda.t.seek<-seq(-15,10,1)
    ## lambda.t<-25
    ## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
    ## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
    lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
    lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
    ## print(c(lambda.c,lambda.t))
    ## print(range(lambda.c.seek))
    print(range(lambda.t.seek))

    ##lambda.c<- -15
    ##lambda.t<-5
    ## lambda.c.old<-lambda.c
    lambda.t.old<-lambda.t

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-5,lambda.c+5,1)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-5,lambda.t+5,1)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-2.5,lambda.c+2.5,.5)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-2.5,lambda.t+2.5,.5)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
        ## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-1.5,lambda.c+1.5,.25)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-1.5,lambda.t+1.5,.25)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.5,lambda.c+.5,.1)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.5,lambda.t+.5,.1)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }



    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.25,lambda.c+.25,.05)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.25,lambda.t+.25,.05)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.15,lambda.c+.15,.025)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.15,lambda.t+.15,.025)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.1,lambda.c+.1,.01)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.1,lambda.t+.1,.01)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.05,lambda.c+.05,.005)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.05,lambda.t+.05,.005)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.025,lambda.c+.025,.0025)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.025,lambda.t+.025,.0025)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.01,lambda.c+.01,.001)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.01,lambda.t+.01,.001)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.0025,lambda.c+.0025,.0005)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.0025,lambda.t+.0025,.0005)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }
    
    ## optim.lambda<-c(lambda.c,lambda.t)
    optim.lambda<- lambda.t


    invisible(optim.lambda)
    
}


