reptolong <- function(design){
     if (!"design" %in% class(design)) stop("This function is applicable for class design objects only.")
     if (!design.info(design)$format == "repeatedMeasuresWide")
        stop("design is not of wide repeated measures format")
     di <- design.info(design)
     rl <- di$responselist
     restl <- di$restlist
     ro <- run.order(design)
     desnum <- desnum(design)
     ## use product here in the hope that it works with messed-up crossed designs
     if (di$type=="FrF2.blocked") di$wbreps<-di$wbreps*nrow(rl)
        else di$replications <- di$replications*nrow(rl)
     di$repeat.only <- TRUE
     di$response.names <- colnames(rl)
     di <- di[-which(names(di) %in% c("format","responselist"))]
     ro <- data.frame(run.no.in.std.order=rep(ro$run.no.in.std.order,each=nrow(rl)),
                          run.no = 1:(nrow(design)*nrow(rl)),
                          run.no.std.rp=paste(rep(ro$run.no.in.std.order,each=nrow(rl)),rep(1:nrow(rl),nrow(design)),sep="."))
     pasteform <- function(name){
           aus <- paste("interleave(",name)
             for (i in 2:nrow(rl)) aus <- paste(aus, name, sep=",")
           aus <- paste(aus, ")")
           aus
           }
     hilf <- undesign(design)
     designlong <- eval(parse(text=pasteform("hilf")))
     desnumlong <- eval(parse(text=pasteform("desnum")))
     rownames(designlong) <- rownames(desnumlong) <- rownames(ro) <- ro$run.no
     ## treat all columns that vary between repeated measurements
     ## these are all columns execpt blocks and factors 
     ##     that have not been specifically marked as constant in reptowide
     ## all in rl and restl
     if (!is.null(restl)) rl <- cbind(rl, restl)
     for (i in 1:ncol(rl)){
         hilf <- NULL
         cn <- colnames(rl)[i]
           hilf <- c(t(as.matrix(design[,rl[,i]])))
             ##stacks rows
         ## change first column name to long version  
         colnames(designlong)[which(colnames(designlong)==rl[1,i])] <- cn
         ## assign correct response values to it
         designlong[,colnames(rl)[i]] <- hilf
         ## remove other columns for the variable
         hilf <- which(colnames(designlong) %in% rl[,i])
         designlong <- designlong[,-hilf]
     }
     for (i in 1:ncol(rl)){
         hilf <- NULL
         cn <- colnames(rl)[i]
           hilf <- c(t(as.matrix(desnum[,rl[,i]])))
         colnames(desnumlong)[which(colnames(desnumlong)==rl[1,i])] <- cn
         desnumlong[,colnames(rl)[i]] <- hilf
         hilf <- which(colnames(desnumlong) %in% rl[,i])
         desnumlong <- desnumlong[,-hilf]
     }
     class(designlong) <- c("design","data.frame")
     desnum(designlong) <- desnumlong
     run.order(designlong) <- ro
     design.info(designlong) <- di
   designlong
}