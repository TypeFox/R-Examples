fold.design <- function(design, columns="full", ...){
## currently does not cover splitplot, block and parameter designs
## splitplot should be easiest to include as well

## full foldover:
## standard order for mirror runs is actually mirrored, i.e. e.g. for 8 to 16, mirror of 1 is 16, mirror of 2 is 15 etc.
## base factors remain base factors
## some generators are previous generators plus nruns, others remain unchanged, depending on the number of 
##    factors involved

## any other foldover
## analogous

   ## error checks
   if (!"design" %in% class(design))
      stop("design must be of class design")
   di <- design.info(design)
    if (!(substr(di$type,1,4)=="FrF2" | substr(di$type,1,2)=="pb" | 
             (length(grep("full factorial",di$type))>0 & all(di$nlevels==2))))
      stop("fold.design is applicable for FrF2 or pb designs only")
#   if (di$type %in% c("FrF2.blocked", "FrF2.param", "FrF2.splitplot"))
#      stop("blocked, split-plot and parameter designs cannot be treated with fold.design")
   if (di$type %in% c("FrF2.blocked", "FrF2.param"))
      stop("blocked designs and parameter designs in long format cannot be treated with fold.design")
   if (length(grep("center", di$type)) > 0)
      stop("currently, designs with center points cannot be treated with fold.design")
   if (length(grep("full factorial",di$type))>0)
      stop("currently, full factorial designs cannot be treated with fold.design")
   if (!(identical(columns, "full") | all(columns %in% names(di$factor.names)) | is.numeric(columns)))
      stop("columns must be full or a vector with factor position numbers or factor names")
   methcall <- columns
   if (is.character(columns) & !identical(columns, "full"))
      columns <- which(di$factor.names %in% columns)
   if (identical(columns, "full")) columns <- 1:di$nfactors
   if (is.numeric(columns))
      if (any(columns < 1 | columns>di$nfactors)) stop("elements of columns must between 1 and nfactors")
   ## now, columns is a vector of column numbers for columns to be mirrored
   fn <- names(factor.names(design))
   on <- setdiff(colnames(design),fn)  ## columns other than experimental factors
   k <- round(log2(di$nruns))  ## meaningful for FrF2 only, but for convenience generally available
   hilf <- undesign(design)[fn]
   hilf.mirror <- hilf
   for (i in 1:length(columns))
      hilf.mirror[[columns[i]]][1:nrow(design)] <- as.character(factor(3-as.numeric(hilf.mirror[[columns[i]]]), 
           labels=levels(hilf.mirror[[columns[i]]])))
   ## create new design matrix
   ## place fold factor at the end of the factors (pb), 
   ##       at the beginning of the factors (FrF2.splitplot), or after the base factors (other FrF2)
     if (substr(di$type,1,2)=="pb" | di$nfactors <= k){
        if (length(on)>0){
           hilf.mirror <- cbind(hilf.mirror,fold=rep("mirror",nrow(design)),matrix(NA,ncol=length(on),nrow=nrow(design)))
           hilf <- cbind(hilf, fold=rep("original",nrow(design)), design[,on])
           }
        else{
           hilf.mirror <- cbind(hilf.mirror,fold=rep("mirror",nrow(design)))
           hilf <- cbind(hilf, fold=rep("original",nrow(design)))
         }
        colnames(hilf.mirror) <- colnames(hilf)
        di$factor.names <- c(di$factor.names,list(fold=c("original","mirror")))
     }
     else{ 
       ## FrF2 type designs with more factors than base ones
       
       ## splitplot special, because WP-factors in front even if not base factors
       ## Fold factor after the WP factors because it is WP factor!!!
       if (di$type=="FrF2.splitplot"){
         if (length(on)>0){
           hilf.mirror <- cbind(fold=rep("mirror",nrow(design)),hilf.mirror,matrix(NA,ncol=length(on),nrow=nrow(design)))
           hilf <- cbind(fold=rep("original",nrow(design)), hilf, design[,on])
         }
         else{
           hilf.mirror <- cbind(fold=rep("mirror",nrow(design)),hilf.mirror)
           hilf <- cbind(fold=rep("original",nrow(design)),hilf)
           }
          di$factor.names <- c(list(fold=c("original","mirror")),di$factor.names)
           } 
       else {
         if (length(on)>0){
             ## further variables to be accomodated
             ## more than k factors
             hilf.mirror <- cbind(hilf.mirror[,1:k],
                                  fold=rep("mirror",nrow(design)),
                                  hilf.mirror[,(k+1):ncol(hilf.mirror),drop=FALSE],
                                  matrix(NA,ncol=length(on),nrow=nrow(design)))
             hilf <- cbind(hilf[,1:k], fold=rep("original",nrow(design)), hilf[,(k+1):ncol(hilf),drop=FALSE], design[,on])
         }
         else{ hilf.mirror <- cbind(hilf.mirror[,1:k],fold=rep("mirror",nrow(design)),
                             hilf.mirror[,(k+1):ncol(hilf.mirror),drop=FALSE])
             hilf <- cbind(hilf[,1:k], fold=rep("original",nrow(design)), hilf[,(k+1):ncol(hilf),drop=FALSE])
         }
         if (!di$nfactors>k)
            di$factor.names <- c(di$factor.names,list(fold=c("original","mirror")))
         else di$factor.names <- c(di$factor.names[1:k],list(fold=c("original","mirror")),di$factor.names[(k+1):di$nfactors])

         }
         colnames(hilf.mirror) <- colnames(hilf)
     }
#   if (length(on)>0)
#     hilf <- cbind(hilf, fold=rep("original",nrow(design)), design[,on])
#   else
#     hilf <- cbind(hilf, fold=rep("original",nrow(design)))
     hilf <- rbind(hilf, hilf.mirror)
     for (nam in fn){ 
        hilf[[nam]] <- factor(hilf[[nam]])
        contrasts(hilf[[nam]]) <- contr.FrF2(2)
        }
     contrasts(hilf$fold) <- contr.FrF2(2)
   
   gen <- NULL

   ## extract and remove old generator information
   ## also in form catlg.entry (or base.design)
   if (!is.null(di$catlg.entry)){
      gen <- di$catlg.entry[[1]]$gen
      if (is.null(di$base.design) & !is.null(di$catlg.entry)) di$base.design <- names(di$catlg.entry)
      di$catlg.entry <- NULL
   }
   else if (!is.null(di$generators)) 
        gen <- sapply(strsplit(di$generators,"="), 
        function(obj) if (substr(obj[2],1,1)=="-") -1*which(names(Yates)==substring(obj[2],2)) 
                      else which(names(Yates) == obj[2]))

  ## now gen is a vector of signed columns 
  ## as it was valid for the original design
  
  ## walk through all generated columns to determine whether they are moved to the right or not
  if (!is.null(gen)){
  for (i in 1:length(gen)){
     ## decide about signs of generators
     ## decide is 0 for words of even length
     ##   and 1 for words of odd length
     decide <- length(intersect(
                   c(round( log2( Yates[[abs(gen[i])]] )) + 1 , k+i), columns
                   ))%%2
     if (decide == 1){ 
            ## odd word
            gen[i] <- -sign(gen[i])*(abs(gen[i]) + di$nruns)  ## nruns columns to the right, reverse sign
            }
  }
  }
  
   ## below, numbering includes fold column
   di$nruns <- 2*di$nruns
   di$nfactors <- 1+di$nfactors
   di$nfac.WP <- 1+di$nfac.WP
   di$nWPs <- 2*di$nWPs
   di$res.WP <- "unknown"  ### ??? have to work on this
   ## new k would be k+1, not changed

  if (!is.null(gen)){
   ## create new generator information for display
   minus <- sign(gen)
   cminus <- rep("", length(gen))
   cminus[minus<0] <- "-"
   if (di$nfactors <=50)
      di$generators <- paste(Letters[(k+2):di$nfactors],paste(cminus,names(Yates)[abs(gen)],sep=""),sep="=")
   else di$generators <- paste(paste("F",(k+2):di$nfactors, sep=""),
               paste(cminus, sapply(gen,function(obj) paste(paste("F",abs(obj),sep=""),collapse=":")),sep=""), sep="=")
   }
   
   ## process alias information
   alorder <- length(di$aliased) - 1
   if (!alorder %in% c(2,3)) alorder <- 2 ## set for the case that it was unconfounded already
    if (di$nfactors<=50) legend <- paste(Letters[1:di$nfactors],names(di$factor.names),sep="=")
       else legend <- paste(paste("F",1:di$nfactors,sep=""),names(di$factor.names),sep="=")
   
   ### easy for generators, other scenarios more difficult
   if (!is.null(gen)) di$aliased <- c(legend=list(legend), alias3fi(k+1, gen, order=alorder))
   else{ 
      ## estimable is handled here
          ## split-plot designs might also be treatable here later
          ## blocked designs with generators should also be possible later
          ## (current priority for extending is low)
      dat <- hilf[,names(di$factor.names)]
      stupid.workaround <- rnorm(nrow(dat))
      if (alorder==2) fit <- lm(stupid.workaround~(.)^2,data=dat)
      if (alorder==3) fit <- lm(stupid.workaround~(.)^3,data=dat)
      if (!substr(di$type,1,2)=="pb") di$aliased <- aliases(fit, code=TRUE, condense=TRUE)
   } 
   

   desnum <- model.matrix(as.formula(paste("~",paste(names(di$factor.names),collapse="+"))), hilf)[,-1]
    if (length(on)>0){
      desnumo <- desnum(design)[,on,drop=FALSE]
      desnumo2 <- desnumo*NA
      desnum <- cbind(desnum, rbind(desnumo, desnumo2))
      }

   ro <- run.order(design)
   run.no <- c(ro$run.no, ro$run.no+max(ro$run.no))
   ## this should be adjusted for regular fractional factorials
   if (substr(di$type,1,2)=="pb"){ 
     run.no.in.std.order <- c(paste(ro$run.no.in.std.order,"orig",sep="."),paste(ro$run.no.in.std.order,"mirror",sep="."))
     run.no.std.rp <- c(paste(ro$run.no.std.rp,"orig",sep="."),paste(ro$run.no.std.rp,"mirror",sep="."))}
   else{
     if (any(1:k %in% columns)){
        ll <- rep(list(c(-1,1)), k)
        for (i in 1:k) if (i %in% columns) ll[[i]] <- -ll[[i]]
        v <- ord(expand.grid(ll)[,k:1]) 
     }
     else v <- 1:k
     ### now the reordering of standard order runs for the mirror portion has been calculated
     ### must be accomodated including replicates
     ## is easy to do as v is a permutation that solely swaps pairs of runs 
     
     ## so complicated because of split-plot designs
     run.no.in.std.order <- sapply(strsplit(as.character(ro$run.no.in.std.order),".",fixed=TRUE),function(obj) as.numeric(obj[1]))
     reordered <- v[run.no.in.std.order]
     
     ## create the replaced rp version for run order
     run.no.std.rp2 <- mapply("sub",run.no.in.std.order,reordered+round(di$nruns/2),ro$run.no.std.rp)
     run.no.in.std.order <- c(run.no.in.std.order, reordered+round(di$nruns/2))
     run.no.std.rp <- c(as.character(ro$run.no.std.rp), run.no.std.rp2)
   }
     
   ro <- data.frame(run.no.in.std.order,run.no, run.no.std.rp)
   di$type <- paste(di$type,"folded",sep=".")
   aus <- hilf
   class(aus) <- c("design","data.frame")
   desnum(aus) <- desnum
   run.order(aus) <- ro
   di$creator <- list(di$creator, fold=methcall)
   if (length(grep("splitplot",di$type))>0){
         hilfresp <- 1:di$nruns
         hilf <- aliases(lm(hilfresp~(.)^2,aus[,names(di$factor.names)[1:di$nfac.WP]]),code=TRUE,condense=TRUE)
         if(!is.null(hilf$main)) di$res.WP <- 3
         else if (!is.null(hilf$fi2)) di$res.WP <- 4
         if (di$res.WP=="unknown"){
             hilf <- aliases(lm(hilfresp~(.)^3,aus[,names(di$factor.names)[1:di$nfac.WP]]),code=TRUE,condense=TRUE)
             if(!is.null(hilf$fi2)) di$res.WP <- 5
             else if (!is.null(hilf$fi3)) di$res.WP <- 6
                  else {
                     di$res.WP <- 7
                     if (!di$nfac.WP < di$res.WP) 
                     message("resolution of whole plot portion of design is at least 7,
                     \nres.WP element of design.info attribute has been set to Inf")
                     di$res.WP <- Inf
                  }
         } 
      }
   design.info(aus) <- di
   aus
}