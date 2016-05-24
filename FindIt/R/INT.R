INT <- function(object,target.data,column,dist="target", base, sort=TRUE,
                  compare=FALSE,order=2){

    if(all(class(object)!="FindIt")){
        warning("the class of object needs to be FindIt.")
    }

    Restriction <- NULL
    measure <- "TIE"
    
    if(compare==TRUE){
        if(missing(order)){
            warning("Need to specify the order")
        }
    }
    
    coefs <- object$coefs.orig
    names(coefs) <- c("Intercept",object$name.t)
    outcome.type <- object$type

    
    if(dist=="unique"){
        data <- cbind(object$y.orig,object$treat.orig)
        colnames(data)[-1] <- colnames(object$treat.orig)
        data <- unique(data,MARGIN=1)
    }
    if(dist=="sample"){
        data <- cbind(object$y.orig,object$treat.orig)
        colnames(data)[-1] <- colnames(object$treat.orig)
    }

    if(dist=="target"){
        data <- target.data
        print("Using Data representing the target population.")
    }else{
        print("Using the Data used to fit the model.")
    }

    ##Relevel
    if(!missing(column) & !missing(base)){
        for(i in 1:length(column)){
            ind.column <- which(colnames(data)==column[i])
            data[,ind.column] <- relevel(data[,ind.column],base[i])
        }
    }

    ## break

    if(!missing(column)){
        ind <- c()
        for(i in 1:length(column)){
            ind[i] <- which(colnames(data)==column[i])
        }
        column <- column[order(ind,decreasing=FALSE)]
    }

    base <- 1
    
    ##Marginal Effect
    Marginal <- function(data,base,Restriction){
    Effect.list <- list()
    Name.list <- list()
    range.mar <- c()
    for(i in 1:(ncol(data)-1)){
        columnM <- colnames(data)[i+1]
        A <- tapply(data[,1],data[,columnM],mean,simplify=FALSE)
        Effect1 <- unlist(A)
        if(base=="min"){
            base <- which(Effect1==min(Effect1))
        }
        TEffect <- Effect1 - Effect1[base]

        ## if(!missing(Restriction)| !is.null(Restriction)){
        ##     for(res in 1:length(Restriction)){
        ##         restriction <- Restriction[[res]]
        ##         if(is.element(columnM,restriction$col.restricted)){
        ##             ## base.ind<- which(is.element(names(Effect1),restriction$fac.restricted))
        ##             ## base <-min(seq(1:length(Effect1))[-base.ind])
                    
        ##             ## if(base=="min"){
        ##             ##     base <- which(Effect1==min(Effect1))
        ##             ## }
                    
        ##             TEffect <- Effect1 - Effect1[base]            
        ##             ind.AC <- which(is.element(data[,restriction$col.restricting],
        ##                                        restriction$fac.restricting))
        ##             dataAC <- data[ind.AC,]
        ##             AC <- tapply(dataAC[,1],dataAC[,columnM],mean,simplify=FALSE)
        ##             Effect1AC <- unlist(AC)
        ##             TEffectAC <- Effect1AC - Effect1AC[base]
        ##             for(z in 1:length(Effect1)){
        ##                 if(is.element(names(Effect1)[z],restriction$fac.restricted)){
        ##                     ## print(names(Effect1)[z])
        ##                     TEffect[z] <- TEffectAC[z]
        ##                 }
        ##             }            
        ##         }
        ##     }
        ## }

        range.mar[i] <- max(TEffect) - min(TEffect)
        Effect.list[[i]] <- TEffect
        Name.list[[i]] <- paste(colnames(data)[(i+1)],names(TEffect),sep="_")
    }
    names(range.mar) <- colnames(data)[2:ncol(data)]
    name <- unlist(Name.list)
    Treatment.Effect1 <- unlist(Effect.list)
    Treatment.Effect1 <- as.data.frame(Treatment.Effect1)
    colnames(Treatment.Effect1) <- "AMTE"
    rownames(Treatment.Effect1) <- name
    ## print(rownames(Treatment.Effect1))
    ## print(names(Treatment.Effect1))
    output <- list("Treatment.Effect1"=Treatment.Effect1,"range"=range.mar)
    return(output)
}

    if(compare==FALSE & missing(column)){
        Treatment.Effect1 <- Marginal(data=data,
                                      base=base,
                                      Restriction=Restriction)$Treatment.Effect1
        Treatment.Effect <- Treatment.Effect1
    }
    if(compare==TRUE & order==1){
        Marginal.range <- Marginal(data=data,
                                   base=base,
                                   Restriction=Restriction)$range
        print("Range of Marginal Effects")
        return(Marginal.range)
    }
    
    
    
    comb.prem <- function(data,column,Restriction,order){
        ## if(!missing(Restriction)|!is.null(Restriction)){
        ##     ##Restriction
        ##     ## Col.Restricted
        ##     Col.Restricted <- c()    
        ##     for(i in 1:length(Res)){
        ##         Col.Restricted <- c(Col.Restricted,Res[[i]]$col.restricted)
        ##     }
        ##     ## Col.Restricting
        ##     Col.Restricting <- c()    
        ##     for(i in 1:length(Res)){
        ##         Col.Restricting <- c(Col.Restricting,Res[[i]]$col.restricting)
        ##     }

        ##     ##Pattern
        ##     ## Now we only support until order 2.
        ##     Res.type <- "Normal"
        ##     ## One 
        ##     if(sum(is.element(column,Col.Restricted))==1 &
        ##        sum(is.element(column,Col.Restricting))==0){Res.type <- "Restricted1"}
        ##     if(sum(is.element(column,Col.Restricted))==0 &
        ##        sum(is.element(column,Col.Restricting))==1){Res.type <- "Restrict1"}
        ##     ## two
        ##     if(sum(is.element(column,Col.Restricted))==2){Res.type <- "Restricted2"}
        ##     if(sum(is.element(column,Col.Restricting))==2){Res.type <- "Restrict2"}
        ##     ## Mix
        ##     if(sum(is.element(column,Col.Restricted))==1 &
        ##        sum(is.element(column,Col.Restricting))==1){
        ##         col.restricted.simple <- column[is.element(column,Col.Restricted)==TRUE]
        ##         col.restricting.simple <- column[is.element(column,Col.Restricting)==TRUE]
        ##         ## Simple Case
        ##         for(resC in 1:length(Restriction)){
        ##             restriction.check <- Restriction[[resC]]
        ##             if((restriction.check$col.restricted==col.restricted.simple)&
        ##                (restriction.check$col.restricting==col.restricting.simple)){
        ##                 Res.type <- "Restriction.mix.simple"
        ##             }else{
        ##                 Res.type <- "Restriction.mix.hard"
        ##             }
        ##         }
                
        ##     }
            
        ##     ## Make the data for Full common Support
        ##     if(Res.type=="Restrict1"| Res.type=="Restriction.mix.simple"){
        ##         column.restrict1 <- column[is.element(column,Col.Restricting)==TRUE]
        ##         for(z in 1:length(Restriction)){
        ##             restrictionR1 <- Restriction[[z]]
        ##             if(restrictionR1$col.restricting==column.restrict1){
        ##                 use.restriction <- restrictionR1
        ##                 break
        ##             }
        ##         }
        ##         ind.res1 <- which(is.element(data[,use.restriction$col.restricted],
        ##                                      use.restriction$fac.restricted))
        ##         data.use <- data[-ind.res1,]
        ##     }
        ##     if(Res.type=="Restrict2"){
        ##         ind.use <- c()
        ##         column.restrict2 <- column
        ##         for(z in 1:length(Restriction)){
        ##             restrictionR2 <- Restriction[[z]]
        ##             if(is.element(restrictionR2$col.restricting,column.restrict2)){
        ##                 use.restriction <- restrictionR2
        ##                 ind.res2 <- which(is.element(data[,use.restriction$col.restricted],
        ##                                              use.restriction$fac.restricted))
        ##                 ind.use <- union(ind.use,ind.res2)
        ##             }
        ##         }       
        ##         data.use <- data[-ind.use,]
        ##     }
        ##     if(Res.type=="Restricted1"){
        ##         column.restricted1 <- column[is.element(column,Col.Restricted)==TRUE]
        ##         for(z in 1:length(Restriction)){
        ##             restrictionRt1 <- Restriction[[z]]
        ##             if(restrictionRt1$col.restricted==column.restricted1){
        ##                 use.restriction <- restrictionRt1
        ##                 break
        ##             }
        ##         }
        ##         ind.rest1 <- which(is.element(data[,use.restriction$col.restricting],
        ##                                       use.restriction$fac.restricting))
        ##         data.use <- data[ind.rest1,]
        ##     }
        ##     if(Res.type=="Restricted2"){
        ##         ind.use <- seq(1:nrow(data))
        ##         column.restricted2 <- column
        ##         for(z in 1:length(Restriction)){
        ##             restrictionRt2 <- Restriction[[z]]
        ##             if(is.element(restrictionRt2$col.restricted,column.restricted2)){
        ##                 use.restriction <- restrictionRt2
        ##                 ind.rest2 <- which(is.element(data[,use.restriction$col.restricting],
        ##                                               use.restriction$fac.restricting))
        ##                 ind.use <- intersect(ind.use,ind.rest2)
        ##             }
        ##         }    
        ##         data.use <- data[ind.use,]
        ##     }
        ##     if(Res.type=="Normal"){
        ##         data.use <- data
        ##     }
        ##     if(Res.type=="Restriction.mix.hard"){
        ##         column.restrictedM <- column[is.element(column,Col.Restricted)==TRUE]
        ##         column.restrictM <- column[is.element(column,Col.Restricting)==TRUE]
        ##         for(z in 1:length(Restriction)){
        ##             restrictionRtM <- Restriction[[z]]
        ##             if(restrictionRtM$col.restricted==column.restrictedM){
        ##                 use.restricted <- restrictionRtM               
        ##             }
        ##             if(restrictionRtM$col.restricting==column.restrictM){
        ##                 use.restricting <- restrictionRtM               
        ##             }
        ##         }
        ##         ## Keep ind.restrictedM
        ##         ind.restrictedM <- which(is.element(data[,use.restricted$col.restricting],
        ##                                             use.restricted$fac.restricting))
        ##         data.use1 <- data[ind.restrictedM,]
        ##         ## Remove ind.restrictM
        ##         ind.restrictM <- which(is.element(data.use1[,use.restricting$col.restricted],
        ##                                           use.restricting$fac.restricted))
        ##         data.use <- data.use1[-ind.restrictM,]                
        ##     }

        ##     data <- data.use
        ## }

        ## Column name and Restriction Type
        ## print("Column name")
        ## print(column)
        ## print("Restriction Type")
        ## print(Res.type)
        
        ##Treatment Effect
        A <- tapply(data[,1],data[,column],mean,simplify=FALSE)
        A2 <- tapply(data[,1],data[,column],mean,simplify=TRUE)
        A[is.na(A2)] <- NA
        Effect1 <- unlist(A)
        if(base=="min"){
            base <- which(Effect1==min(Effect1))
        }        
        
        ## Combination Effect
        Comb.Effect <- Effect1 - Effect1[base]
        Treatment.Effect <- cbind(Comb.Effect,
                                  expand.grid(dimnames(A)))
        Treatment.Effect <- na.omit(Treatment.Effect)
        Comb.Effect <- Treatment.Effect[,1]
        Combination.name <- Treatment.Effect[,-1]
        Treatment.Effect.print <- Treatment.Effect
        
        for(j in 2:ncol(Treatment.Effect)){
            Treatment.Effect[,j] <- paste(colnames(Treatment.Effect)[j],
                                          Treatment.Effect[,j],
                                          sep="_")
        }

        ## Sum of Marginals
        Treatment.Effect1 <- Marginal(data=data,
                                      base=base,
                                      Restriction=Restriction)$Treatment.Effect1
        Main.comb <- c()
        for(i in 1:nrow(Treatment.Effect)){
            main <- c()
            for(j in 1:(ncol(Treatment.Effect)-1)){
                main[j] <- Treatment.Effect1[rownames(Treatment.Effect1)==
                                             Treatment.Effect[i,(j+1)],1]
            }
            Main.comb[i] <- sum(main)
        }
        Sum.Mar.Effect <- Main.comb - Main.comb[base]

        if(order==2){
            TIE <- round(Comb.Effect - Sum.Mar.Effect, digits=8)

        }
        if(order==3){
            ##Two combination Effect List
            data.column.three <- data[,column]
            Two.way.comb.list <- list()
            for(j in 1:3){
                column.com.three <- c()
                Two.way.comb <- list()
                column.com.three <-
                    colnames(data.column.three)[c(combn(length(data.column.three),2)[,j])]
                A.Three <- tapply(data[,1],data[,column.com.three],mean,simplify=FALSE)
                A2.Three <- tapply(data[,1],data[,column.com.three],mean,simplify=TRUE)
                A.Three[is.na(A2.Three)] <- NA
                EffectTwo <- unlist(A.Three)
                Two.way.comb <- EffectTwo - EffectTwo[base]
                Two.way.comb2 <- cbind(Two.way.comb,
                                          expand.grid(dimnames(A.Three)))
                for(k in 2:3){
                    Two.way.comb2[,k] <- paste(colnames(Two.way.comb2)[k],
                                               Two.way.comb2[,k],sep="_")
                }
                colnames(Two.way.comb2) <- c("Two.way.comb","first","second")
                Two.way.comb.list[[j]]  <- na.omit(Two.way.comb2)
                ## print(j)
                ## print("th")
                ## print(Two.way.comb.list[[j]])
            }
            Two.way.comb.Final <- do.call(rbind,Two.way.comb.list)


            Two.Way.sum <- c()
            for(i in 1:nrow(Treatment.Effect)){
                two.sum <- c()
                for(j in 1:nrow(Two.way.comb.Final)){
                    ind.two <- sum(is.element(Two.way.comb.Final[j,2:3],
                                         Treatment.Effect[i,2:4]))
                    if(ind.two==2){
                        two.sum[j] <- Two.way.comb.Final[j,1]
                    }else{
                        two.sum[j] <- 0
                    }
                }
                Two.Way.sum[i] <- sum(two.sum)
            }
            Sum.Two.Effect <- Two.Way.sum - Two.Way.sum[base]
            TIE <- round(Comb.Effect - Sum.Two.Effect + Sum.Mar.Effect, digits=8)
        }
       
        ## print("The range of TIE")
        range.change <- max(TIE) - min(TIE)
        ## Order diff
        Order.data <- cbind(Sum.Mar.Effect,Comb.Effect)
        Order.data <- Order.data[order(Order.data[,1],decreasing=TRUE),]
        Order.data <- as.data.frame(Order.data)
        Order.data$order.main <- seq(1:nrow(Order.data))
        Order.comb <- Order.data[order(Order.data[,2],decreasing=TRUE),"order.main"]
        Order.diff <- Order.comb - seq(1:nrow(Order.data))
        max.order.diff <- max(abs(Order.diff)) 
        ## print(range.change)
        return(list("A"=A,
                    "Treatment.Effect"=Treatment.Effect,
                    "Comb.Effect"=Comb.Effect,
                    "Combination.name"=Combination.name,
                    "Sum.Mar.Effect"=Sum.Mar.Effect,
                    "TIE"=TIE,
                    "range.change"=range.change,
                    "max.order.diff"=max.order.diff))
    }

    if(compare==FALSE){
        if(!missing(column)){
            X <- comb.prem(data=data,column=column,Restriction=Restriction,order=order)
            A <- X$A
            Treatment.Effect <- X$Treatment.Effect
            Comb.Effect <- X$Comb.Effect
            Combination.name <- X$Combination.name
            Sum.Mar.Effect <- X$Sum.Mar.Effect
            TIE <- X$TIE
            range.change <- X$range.change     
            
        }
    }else{        
        data.column <- data[,-1]
        column.com <- list()
        Range.com <- list()
        for(j in 1:choose(length(data.column),order)){
            column.com[[j]] <-
                colnames(data.column)[c(combn(length(data.column),order)[,j])]
            if(measure=="TIE"){
                Range.com[[j]] <- comb.prem(data=data,
                                            column=column.com[[j]],
                                            Restriction=Restriction,
                                            order=order)$range.change
            }
            if(measure=="Order.Diff"){
                Range.com[[j]] <- comb.prem(data=data,
                                            column=column.com[[j]],
                                            Restriction=Restriction,
                                            order=order)$max.order.diff
            }
            ## print(Range.com[[j]])
            ## print(column.com[[j]])
        }
        column.com <- as.data.frame(matrix(unlist(column.com),ncol=order,byrow=TRUE))
        Compare <- as.data.frame(cbind(column.com,unlist(Range.com)))
        Compare <- Compare[order(Compare[,(order+1)],decreasing=TRUE),]
        if(measure=="TIE"){
            colnames(Compare) <- c(colnames(Compare)[1:order],"range.TIE")
        }
        if(measure=="Order.Diff"){
            colnames(Compare) <- c(colnames(Compare)[1:order],"max.Order.diff")
        }
        return(Compare)
    }


    ##Put the coefficients together
    if(!missing(column)){
        if(length(column)==2){
            Treatment.coefs.name <- Treatment.Effect[,2]
            if(ncol(Treatment.Effect)>=3){
                for(i in 1:nrow(Treatment.Effect)){
                    for(j in 3:ncol(Treatment.Effect)){
                        Treatment.coefs.name[i] <-
                            paste(Treatment.coefs.name[i],
                                  Treatment.Effect[i,j],
                                  sep=".")
                    }
                }
            }
            Treatment.coefs.name <- gsub(" ", ".", Treatment.coefs.name)
        }
    }

    if(missing(column)){
        Treatment.coefs.name <- rownames(Treatment.Effect1)
        Treatment.coefs.name <- gsub(" ", ".", Treatment.coefs.name)
    }

    ##print("This is point")
    ##print(Treatment.coefs.name)
    ## print(Treatment.coefs.name2)
    ##print(names(coefs))
    ##print(Treatment.coefs.name[1])
    ##print(dim(Treatment.Effect))
    ##print(is.element(names(coefs), Treatment.coefs.name))
    ##print(coefs)
    ##print(coefs[names(coefs)==Treatment.coefs.name[1]])

    if(missing(column)){
        Coefficients <- c()
        for(i in 1:length(Treatment.coefs.name)){
            if(Treatment.coefs.name[i] %in% names(coefs)){
                Coefficients[i] <- coefs[names(coefs)==Treatment.coefs.name[i]]                
            }else{
                Coefficients[i] <- NA
            }
        }
        if(outcome.type=="binary"){
            Coefficients <- Coefficients/2
        }
    }
    if(!missing(column)){
        if(length(column)==2){
            Coefficients <- c()
            for(i in 1:length(Treatment.coefs.name)){
                if(Treatment.coefs.name[i] %in% names(coefs)){
                    Coefficients[i] <- coefs[names(coefs)==Treatment.coefs.name[i]]
                }else{
                    Coefficients[i] <- NA
                }
            }
        
            if(outcome.type=="binary"){
                Coefficients <- Coefficients/2
            }
        }
    }



    ## print(Treatment.coefs)

    if(!missing(column)){
        if(length(column)==2){
            Treatment.Effect.print <- cbind(Comb.Effect,
                                            TIE,
                                            ## Coefficients,
                                            Combination.name)
            colnames(Treatment.Effect.print) <- c("ATCE","AMTIE",colnames(Combination.name))
            Main.matrix <- cbind(Sum.Mar.Effect,
                                 TIE,
                                 ## Coefficients,
                                 Combination.name)
            colnames(Main.matrix) <- c("Sum of AMTEs","AMTIE",colnames(Combination.name))
            Inequality.Cor <- cor(Sum.Mar.Effect,TIE)
            Relative.change <- cbind(TIE,Combination.name)
            colnames(Relative.change) <- c("AMTIE",colnames(Combination.name))
            Relative.change <- Relative.change[order(Relative.change[,1],decreasing=TRUE),]
        }
        if(length(column)>=3){
            Treatment.Effect.print <- cbind(Comb.Effect,
                                            TIE,
                                            Combination.name)
            colnames(Treatment.Effect.print) <- c("ATCE","AMTIE",colnames(Combination.name))
            Main.matrix <- cbind(Sum.Mar.Effect,
                                 TIE,
                                 Combination.name)
            colnames(Main.matrix) <- c("Sum of AMTEs","AMTIE",colnames(Combination.name))
            Relative.change <- cbind(TIE,Combination.name)
            colnames(Relative.change) <- c("AMTIE",colnames(Combination.name))
            Relative.change <- Relative.change[order(Relative.change[,1],decreasing=TRUE),]
        }
    }else{
        ## Treatment.Effect.print <- cbind(Treatment.Effect1,
        ##                                 Coefficients
        ##                                 )
        Treatment.Effect.print <- Treatment.Effect1                                       
        Man.matrix <- NULL
    }

    Treatment.Effect.print <- as.data.frame(Treatment.Effect.print)
    if(!missing(column)){
        Main.matrix <- as.data.frame(Main.matrix)
    }

    if(sort==TRUE){
        Treatment.Effect.print1 <- Treatment.Effect.print
        Treatment.Effect.print <-
            Treatment.Effect.print[order(Treatment.Effect.print[,1],
                                         decreasing=TRUE),]
        
        if(!missing(column)){
            ## Order.data <- cbind(Main.matrix$Sum,Treatment.Effect.print1[,1])
            ## Order.data <- Order.data[order(Order.data[,1],decreasing=TRUE),]
            ## Order.data <- as.data.frame(Order.data)
            ## Order.data$order.main <- seq(1:nrow(Order.data))
            ## Order.comb <- Order.data[order(Order.data[,2],decreasing=TRUE),"order.main"]
            ## Order.diff <- Order.comb - seq(1:nrow(Order.data)) 
            
            Main.matrix <-
                Main.matrix[order(Main.matrix$Sum,
                                  decreasing=TRUE),]
            
            Treatment.Effect.print <- cbind(
                Treatment.Effect.print[,1:2],
                ## Order.diff,
                Treatment.Effect.print[,3:ncol(Treatment.Effect.print)])
        }
    }

    ## if(!missing(column)){
    ##     if(range.change>0){
    ##         print("Interaction matters")
    ##     }else{
    ##         print("No Interaction Effect")
    ##     }
    ## }

    ## if(!missing(column)){
    ##     if(all(Treatment.Effect.print[,"Coefficients"]==Main.matrix[,"Coefficients"])){
    ##         print("No Significant Interaction Effect")
    ##     }else{
    ##         print("Interaction Effect changes the order")
    ##     }
    ## }
    if(!missing(column)){
        if(sort){
            return(list("Range of AMTIE"=range.change,
                        "AMTIE"=Relative.change,
                        ## "Inequality Correlation"=Inequality.Cor,
                        ## "Absolute Order change"=max(abs(Order.diff)),
                        "ATCE"=Treatment.Effect.print,
                        "Sum of AMTEs"=Main.matrix
                        ))
        }else{
            return(list("Range of AMTIE"=range.change,
                        "AMTIE"=Relative.change,
                        ## "Inequality Correlation"=Inequality.Cor,
                        "ATCE"=Treatment.Effect.print,
                        "Sum of AMTEs"=Main.matrix
                        ))
        }
    }else{
        return(list("Treatment.Effect"=Treatment.Effect.print))
    }
}

## Restriction should look like this.
## restriction1 <- list("col.restricted"="FeatReason",
##                     "fac.restricted"="Escape",
##                     "col.restricting"="FeatCountry",
##                     "fac.restricting"=c("China","Sudan","Somalia","iraq"))

## restriction2 <- list("col.restricted"="FeatJob",
##                     "fac.restricted"=c(
##                         "Doctor","Research scientist",
##                         "Computer programmer","Financial analyst"),
##                     "col.restricting"="FeatEd",
##                     "fac.restricting"=c("Twoyears.college","College","Graduate"))

## Res <- list(restriction1,restriction2)
