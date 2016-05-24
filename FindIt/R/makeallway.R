makeallway<-function(X,threshold=0.999999,deletion=TRUE,
                     make.reference=TRUE,
                     sparse.use=FALSE,
		     nway){
    if(!is.vector(X) & !is.matrix(X)){
        warning("X should be a vector or matrix.")
    }
    X2<-NULL
    ## Make indicator variables for all columns.
    if(is.matrix(X)){
        matrix <- 1
        X <- data.frame(X)
        for(i in 1:dim(X)[2]) {
            X2[[i]]<-cbind(sapply(sort(unique(X[,i])),FUN=function(j) 1*(X[,i]==j)))
            colnames(X2[[i]])<-c(paste(names(X[i]), sort(unique(X[,i])),sep="_"))
        }
        n <- ncol(X)
    }

    ## if(is.vector(X)){
    ##     matrix <- 0
    ##     n <- length(unique(X))
    ##     X2 <- matrix(NA, ncol=n, nrow=length(X))
    ##     name <- c()
    ##     for(z in 1:n){
    ##         for(i in 1:length(X)){
    ##             X2[i,z]<- as.numeric(X[i]==sort(unique(X))[z])
    ##         }
    ##         name[z]<-c(paste(names(X),sort(unique(X))[z],sep="_"))              
    ##     }
    ##     X2 <- as.data.frame(X2)
    ##     colnames(X2) <- name
    ## }

    ## print("makeallway")
    ## print(head(X2[[1]]))

    
    ## NEW data sets
    ##One-way 
    one.way.data <- NULL
    one.way.data <- as.data.frame(X2)
    if(sparse.use==TRUE){
        one.way.data <- as.matrix(one.way.data)
        one.way.data <- Matrix(one.way.data,sparse=TRUE)
    }
    ## print(head(one.way.data))
    
    
    if(matrix == 1){
        ##Two-way
        if(ncol(X)>=2){
            two.way.data <- NULL
            two.way.name.w <- NULL
            for(j in 1:choose(n,2)){
                formula  <- ~ X2[[combn(n,2)[1,j]]]:X2[[combn(n,2)[2,j]]]
                two.way.data[[j]] <- model.matrix(formula)[,-1]
                if(sparse.use==TRUE){
                    two.way.data[[j]] <- Matrix(two.way.data[[j]],sparse=TRUE)
                }
                two.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,2)[1,j]]]),
                                                   colnames(X2[[combn(n,2)[2,j]]]))
                name <- c()
                for(i in 1:nrow(two.way.name.w[[j]])){
                    name[i] <- paste(two.way.name.w[[j]][i,1],
                                     two.way.name.w[[j]][i,2],sep=":")
                }
                colnames(two.way.data[[j]]) <- name
            }
            if(sparse.use==TRUE){
                two.way.data <- do.call(cBind,two.way.data)
            }else{
                two.way.data <- as.data.frame(two.way.data)
            }
            ## two.way.data <- as.matrix(two.way.data)
            ## two.way.data <- Matrix(two.way.data,sparse=TRUE)
        }
        
        ## Three-way
         if(nway>=3){
             three.way.data <- NULL
             three.way.name.w <- NULL
             for(j in 1:choose(n,3)){
                 formula  <- ~ X2[[combn(n,3)[1,j]]]:X2[[combn(n,3)[2,j]]]:
                     X2[[combn(n,3)[3,j]]]
                 three.way.data[[j]] <- model.matrix(formula)[,-1]
                 three.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,3)[1,j]]]),
                                                      colnames(X2[[combn(n,3)[2,j]]]),
                                                      colnames(X2[[combn(n,3)[3,j]]]))
                 name <- c()
                 for(i in 1:nrow(three.way.name.w[[j]])){
                     name[i] <- paste(three.way.name.w[[j]][i,1],
                                      three.way.name.w[[j]][i,2],
                                      three.way.name.w[[j]][i,3],
                                      sep=":")
                 }
                 colnames(three.way.data[[j]]) <- name
             }
             three.way.data <- as.data.frame(three.way.data)
             ## three.way.data <- as.matrix(three.way.data)
             ## three.way.data <- Matrix(three.way.data,sparse=TRUE)
         }
        
        
        ## Four-way
        if(nway>=4){
            four.way.data <- NULL
            four.way.name.w <- NULL
            for(j in 1:choose(n,4)){
                formula  <- ~ X2[[combn(n,4)[1,j]]]:X2[[combn(n,4)[2,j]]]:
                    X2[[combn(n,4)[3,j]]]:X2[[combn(n,4)[4,j]]]
                four.way.data[[j]] <- model.matrix(formula)[,-1]
                four.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,4)[1,j]]]),
                                                    colnames(X2[[combn(n,4)[2,j]]]),
                                                    colnames(X2[[combn(n,4)[3,j]]]),
                                                    colnames(X2[[combn(n,4)[4,j]]])
                                                    )
                name <- c()
                for(i in 1:nrow(four.way.name.w[[j]])){
                    name[i] <- paste(four.way.name.w[[j]][i,1],
                                     four.way.name.w[[j]][i,2],
                                     four.way.name.w[[j]][i,3],
                                     four.way.name.w[[j]][i,4],
                                     sep=":")
                }
                colnames(four.way.data[[j]]) <- name
            }
            four.way.data <- as.data.frame(four.way.data)
            ## four.way.data <- as.matrix(four.way.data)
            ## four.way.data <- Matrix(four.way.data,sparse=TRUE)
        }
        
        
        ## ## Five-way
        ## if(ncol(X)>=5){
        ##     five.way.data <- NULL
        ##     five.way.name.w <- NULL
        ##     for(j in 1:choose(n,5)){
        ##         formula  <- ~ X2[[combn(n,5)[1,j]]]:X2[[combn(n,5)[2,j]]]:
        ##             X2[[combn(n,5)[3,j]]]:X2[[combn(n,5)[4,j]]]:X2[[combn(n,5)[5,j]]]
        
        ##         five.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         five.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,5)[1,j]]]),
        ##                                             colnames(X2[[combn(n,5)[2,j]]]),
        ##                                             colnames(X2[[combn(n,5)[3,j]]]),
        ##                                             colnames(X2[[combn(n,5)[4,j]]]),
        ##                                             colnames(X2[[combn(n,5)[5,j]]])
        ##                                             )
        ##         name <- c()
        ##         for(i in 1:nrow(five.way.name.w[[j]])){
        ##             name[i] <- paste(five.way.name.w[[j]][i,1],
        ##                              five.way.name.w[[j]][i,2],
        ##                              five.way.name.w[[j]][i,3],
        ##                              five.way.name.w[[j]][i,4],
        ##                              five.way.name.w[[j]][i,5],
        ##                              sep=":")
        ##         }
        ##         colnames(five.way.data[[j]]) <- name
        ##     }
        ##     five.way.data <- as.data.frame(five.way.data)
        ##     five.way.data <- as.matrix(five.way.data)
        ##     five.way.data <- Matrix(five.way.data,sparse=TRUE)
        ## }
        

        ## ## Six-way
        ## if(ncol(X)>=6){
        ##     six.way.data <- NULL
        ##     six.way.name.w <- NULL
        ##     for(j in 1:choose(n,6)){
        ##         formula  <- ~ X2[[combn(n,6)[1,j]]]:X2[[combn(n,6)[2,j]]]:
        ##             X2[[combn(n,6)[3,j]]]:X2[[combn(n,6)[4,j]]]:
        ##                 X2[[combn(n,6)[5,j]]]:X2[[combn(n,6)[6,j]]]
        ##         six.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         six.way.name.w[[j]] <- expand.grid( colnames(X2[[combn(n,6)[1,j]]]),
        ##                                            colnames(X2[[combn(n,6)[2,j]]]),
        ##                                            colnames(X2[[combn(n,6)[3,j]]]),
        ##                                            colnames(X2[[combn(n,6)[4,j]]]),
        ##                                            colnames(X2[[combn(n,6)[5,j]]]),
        ##                                            colnames(X2[[combn(n,6)[6,j]]])
        ##                                            )
        ##         name <- c()
        ##         for(i in 1:nrow(six.way.name.w[[j]])){
        ##             name[i] <- paste(six.way.name.w[[j]][i,1],
        ##                              six.way.name.w[[j]][i,2],
        ##                              six.way.name.w[[j]][i,3],
        ##                              six.way.name.w[[j]][i,4],
        ##                              six.way.name.w[[j]][i,5],
        ##                              six.way.name.w[[j]][i,6],
        ##                              sep=":")
        ##         }
        ##         colnames(six.way.data[[j]]) <- name
        ##     }
        ##     six.way.data <- as.data.frame(six.way.data)
        ##     six.way.data <- as.matrix(six.way.data)
        ##     six.way.data <- Matrix(six.way.data,sparse=TRUE)
        ## }
        

        ## ## Seven-way
        ## if(ncol(X)>=7){
        ##     seven.way.data <- NULL
        ##     seven.way.name.w <- NULL
        ##     for(j in 1:choose(n,7)){
        ##         formula  <- ~ X2[[combn(n,7)[1,j]]]:X2[[combn(n,7)[2,j]]]:
        ##             X2[[combn(n,7)[3,j]]]:X2[[combn(n,7)[4,j]]]:
        ##                 X2[[combn(n,7)[5,j]]]:X2[[combn(n,7)[6,j]]]:
        ##                     X2[[combn(n,7)[7,j]]]
        ##         seven.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         seven.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,7)[1,j]]]),
        ##                                              colnames(X2[[combn(n,7)[2,j]]]),
        ##                                              colnames(X2[[combn(n,7)[3,j]]]),
        ##                                              colnames(X2[[combn(n,7)[4,j]]]),
        ##                                              colnames(X2[[combn(n,7)[5,j]]]),
        ##                                              colnames(X2[[combn(n,7)[6,j]]]),
        ##                                              colnames(X2[[combn(n,7)[7,j]]])
        ##                                              )
        ##         name <- c()
        ##         for(i in 1:nrow(seven.way.name.w[[j]])){
        ##             name[i] <- paste(seven.way.name.w[[j]][i,1],
        ##                              seven.way.name.w[[j]][i,2],
        ##                              seven.way.name.w[[j]][i,3],
        ##                              seven.way.name.w[[j]][i,4],
        ##                              seven.way.name.w[[j]][i,5],
        ##                              seven.way.name.w[[j]][i,6],
        ##                              seven.way.name.w[[j]][i,7],
        ##                              sep=":")
        ##         }
        ##         colnames(seven.way.data[[j]]) <- name
        ##     }
        ##     seven.way.data <- as.data.frame(seven.way.data)
        ##     seven.way.data <- as.matrix(seven.way.data)
        ##     seven.way.data <- Matrix(seven.way.data,sparse=TRUE)
        ## }
        

        ## ## Eight-way
        ## if(ncol(X)>=8){
        ##     eight.way.data <- NULL
        ##     eight.way.name.w <- NULL
        ##     for(j in 1:choose(n,8)){
        ##         formula  <- ~ X2[[combn(n,8)[1,j]]]:X2[[combn(n,8)[2,j]]]:
        ##             X2[[combn(n,8)[3,j]]]:X2[[combn(n,8)[4,j]]]:
        ##                 X2[[combn(n,8)[5,j]]]:X2[[combn(n,8)[6,j]]]:
        ##                     X2[[combn(n,8)[7,j]]]:X2[[combn(n,8)[8,j]]]
        ##         eight.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         eight.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,8)[1,j]]]),
        ##                                              colnames(X2[[combn(n,8)[2,j]]]),
        ##                                              colnames(X2[[combn(n,8)[3,j]]]),
        ##                                              colnames(X2[[combn(n,8)[4,j]]]),
        ##                                              colnames(X2[[combn(n,8)[5,j]]]),
        ##                                              colnames(X2[[combn(n,8)[6,j]]]),
        ##                                              colnames(X2[[combn(n,8)[7,j]]]),
        ##                                              colnames(X2[[combn(n,8)[8,j]]])
        ##                                              )
        ##         name <- c()
        ##         for(i in 1:nrow(eight.way.name.w[[j]])){
        ##             name[i] <- paste(eight.way.name.w[[j]][i,1],
        ##                              eight.way.name.w[[j]][i,2],
        ##                              eight.way.name.w[[j]][i,3],
        ##                              eight.way.name.w[[j]][i,4],
        ##                              eight.way.name.w[[j]][i,5],
        ##                              eight.way.name.w[[j]][i,6],
        ##                              eight.way.name.w[[j]][i,7],
        ##                              eight.way.name.w[[j]][i,8],
        ##                              sep=":")
        ##         }
        ##         colnames(eight.way.data[[j]]) <- name
        ##     }
        ##     eight.way.data <- as.data.frame(eight.way.data)
        ##     eight.way.data <- as.matrix(eight.way.data)
        ##     eight.way.data <- Matrix(eight.way.data,sparse=TRUE)
        ## }
    }
    
    ## if(is.vector(X)){FinalData <- as.data.frame(one.way.data)}
    if(ncol(X)==1){FinalData <- one.way.data}
    if(matrix==1){
        if(sparse.use==TRUE){
            if(ncol(X)>=2){FinalData <- cBind(one.way.data,two.way.data)}
        }else{
            if(ncol(X)>=2){FinalData <- cbind(one.way.data,two.way.data)}
        }
        
        if(nway==3){FinalData <- cbind(one.way.data,two.way.data,three.way.data)}
        if(nway==4){FinalData <- cbind(one.way.data,two.way.data,three.way.data,
                   four.way.data)}
        ## if(ncol(X)==5){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data)}
        ## if(ncol(X)==6){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data,six.way.data)}
        ## if(ncol(X)==7){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data,six.way.data,
        ##            seven.way.data)}
        ## if(ncol(X)==8){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data,six.way.data,
        ##            seven.way.data,eight.way.data)}
    }
    
    ## FinalData <- as.data.frame(FinalData)
    ## print(head(FinalData))
    
    if(deletion==TRUE){
        ## Discard the columns with no variation or the same variation
        ## with other columns,
        ## while keeping the column names
        No.variation <- colnames(FinalData[,apply(FinalData,2,sd)==0])
        T2 <- FinalData[,apply(FinalData,2,sd)>0]
        T2.sparse <- T2
        
        T2 <- as.matrix(T2)
        ## T2 <- t(T2)
        coldrop      <- sapply(1:ncol(T2),
                               FUN=function(x)
                               which(cor(T2[,x] ,T2)^2 > threshold))
        
        col.names.w  <- sapply(coldrop, FUN=function(x) colnames(T2)[x])
        col.names    <- list(Perfect.Correlated.Columns=col.names.w,
                             no.variation.columns=No.variation)
        
        ## it could be three columns are the same. 
        coldrop.cols <- sapply(coldrop,FUN=function(x) x[1])
        Keep.cols <- unique(coldrop.cols)
        ## if there are more than two, the first one is gotten.
        ## so the second or the third one is the discarded one.
        ## if there are only one, then it is the same.
        
        ## colnames(T2)[coldrop[[Keep.cols]][-1]]
        ## This gives the corresponding dropped coefficients.
        T3 <- T2.sparse[,Keep.cols]
        
        if(make.reference==TRUE){
            Corresponding <- list()
            ## Keep Discarded Data.
            Discarded.cols <- seq(1:ncol(T2))[-Keep.cols]
            
            for(i in Keep.cols){
                if(length(coldrop[[i]][-1])==0){
                    Corresponding[[i]] <- "No Match"
                }
                if(length(coldrop[[i]][-1])!=0){
                    Corresponding[[i]] <- colnames(T2)[coldrop[[i]][-1]]
                }
            }
            for(i in Discarded.cols){
                Corresponding[[i]] <- "Discarded"
            }
            
            C <- max(sapply(Corresponding, FUN=function(x) length(x)))
            
            CorrespondingM <- matrix(NA,nrow=ncol(T2), ncol=C)
            Variation <- colnames(T2)
            rownames(CorrespondingM) <- colnames(T2)
            for(j in 1:C){
                for(i in 1:ncol(T2)){
                    if(length(Corresponding[[i]])>=j){
                        CorrespondingM[i,j] <- Corresponding[[i]][j]
                    }
                }
            }
            
            Reference <- matrix(NA,nrow=ncol(FinalData), ncol=C)
            rownames(Reference) <- colnames(FinalData)
            
            for(i in 1:nrow(Reference)){
                if(rownames(Reference)[i] %in% No.variation){
                    Reference[i,1] <- "No Variation"
                }
                if(rownames(Reference)[i] %in% Variation){
                    Reference[i,]  <- CorrespondingM[rownames(CorrespondingM)==
                                                     rownames(Reference)[i]]
                }
            }
            
            names <- c("Matched Variables", rep("Other matched variables",C-1))
            colnames(Reference) <- names
            Reference <- as.data.frame(Reference)
        }
        else{
            Reference <- NULL}
    }
    if(deletion==FALSE){
        T3 <- FinalData
        Reference <- NULL
    }
    
    return(list(FinalData=T3, reference=Reference))
}
