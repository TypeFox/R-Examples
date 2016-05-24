### required:    ID
###          OR  factor.names (as list) 
###          OR  nlevels (as vector)
###          OR  nfactors with nlevels as number
### min.residual.df extra runs over and above number of main effects parameters

### accomodate oa.design within the wrapper package

oa.design <- function(ID=NULL, nruns=NULL, nfactors=NULL, nlevels=NULL, 
      factor.names = if (!is.null(nfactors)) {
        if (nfactors <= 50) Letters[1:nfactors]
           else paste("F", 1:nfactors, sep = "")} 
        else NULL, columns="order",  
        replications=1, repeat.only=FALSE,
        randomize=TRUE, seed=NULL, min.residual.df=0, levordold=FALSE){
        ## ID identifies the design
        ## nruns, nfactors, factor.names self-explanatory
        ## nlevels is a numeric vector with the numbers of levels of the factors in the experiment
        ##         or a single number of levels (e.g. 3)
        ## factor.types is a character vector with "qual" or "quant" for each factor
            ## default is that all factors are qualitative
        ## columns is a vector with columns to be assigned to the factors (default is from left to right within each number of levels

        ## nruns is completely superfluous, it is just a check in case a user is unsure about the number of runs
        ## nfactors takes precedence over other ways to specify number of factors
        ## factor.names, nlevels or columns can also be used to specify number of factors
        ## no specification --> maximum number of factors for the array

        ## if more than one of the entries are given:
        ## compatibility checks necessary
      creator <- sys.call()

      generating.oa <- deparse(substitute(ID))  ## document selected OA; NULL if NULL
      if (generating.oa=="NULL") generating.oa <- NULL
      if ("try-error" %in% class(try(is.null(ID), silent = TRUE))){ 
          if (!generating.oa %in% oacat$name) stop("invalid ID")
          ID <- "NULL"   ## made NULL only after the check for nlevels etc.
          }
      
      if (is.null(ID) & is.null(factor.names) & is.null(nlevels)) 
         stop("ID or factor.names or nlevels must be specified!")
      if (identical(ID, "NULL")) ID <- NULL
      
      if (is.null(ID) & is.null(generating.oa)){
          ## determine array, if not explicitly given or non-existing in the data base
          ## array not given: both ID and generating.oa are NULL
          ## array non-existing: ID is NULL, generating.oa is a text string of the array
          if (is.numeric(columns)) 
                stop("column numbers must not be specified, if ID is omitted")
          if (!is.null(nlevels)) {
                if (!is.numeric(nlevels)) stop("nlevels must be numeric")
                if (!all(floor(nlevels)==nlevels)) stop("nlevels must be integer")
                if (length(nlevels)==1 & is.null(nfactors) & is.null(factor.names))
                    stop("designs for one factor only are not implemented")
                if (length(nlevels)==1 & is.null(nfactors)) nlevels <- rep(nlevels,length(factor.names))
                if (length(nlevels)==1 & !is.null(nfactors)) nlevels <- rep(nlevels,nfactors)
                }
          if(is.null(nlevels) & !is.list(factor.names))
             stop("number of levels for each factor must be specified via ID or nlevels or level specifications in factor.names!")
          if (is.null(nlevels) & !is.null(factor.names)) {
             nlevels <- sapply(factor.names, "length")
             if (any(nlevels < 2)) 
                stop("If ID and nlevels are not given, factor.names must contain level entries for at least two levels for each factors.")
          }
          ## now nlevels is known and is a vector longer than 1

          if (!is.null(nfactors)){ 
               if (!nfactors==length(nlevels)) 
                    stop("mismatch between nfactors and nlevels or factor.names")
          }
          else nfactors <- length(nlevels)

          if (is.list(factor.names)){ 
                    if (!(all(nlevels==sapply(factor.names,length) | sapply(factor.names,length)==1)))
                          stop("Entries in nlevels do not match entries in factor.names.") 
            if (is.null(names(factor.names))){ 
                    if (nfactors<=50) names(factor.names) <- Letters[1:nfactors] 
                       else names(factor.names) <- paste("F",1:nfactors,sep="")
                           }}
          ## make names valid under all circumstances 
          ## needed here already because of some return() instances below
          if (!is.null(names(factor.names)))
          names(factor.names) <- make.names(names(factor.names), unique=TRUE)

          ## determine ID from nlevels and other entries!
          hilf <- table(nlevels)
            names(hilf) <- paste("n",names(hilf),sep="")
          minnrun <- sum(nlevels) - nfactors + 1 + min.residual.df
          ffnrun <- prod(nlevels)
          if (!is.null(nruns)) {
                  if (nruns < minnrun) stop("Your request requires at least ", minnrun, ">", nruns, " runs.")
                  if (nruns > ffnrun) warning("You are requesting more runs than needed for a full factorial!")
                  ## full factorial, if this number of runs is requested
                  if (nruns==ffnrun) return(fac.design(nfactors=nfactors, nlevels=nlevels, factor.names=factor.names,
                         replications=replications, repeat.only=repeat.only, randomize=randomize, seed=seed))
                  cand <- oacat[oacat$nruns >= nruns,]
               }
          else cand <- oacat[oacat$nruns >= minnrun,]
          if (nrow(cand)==0) 
              stop("Currently, DoE.base does not contain an orthogonal array that can cover your request.")
          else {
             for (i in 1:length(hilf))
                cand <- eval(parse(text=paste("cand[cand$",names(hilf)[i],">=",hilf[i],",]",sep="")))
## bug fix: cand <=ffnrun was stupid, because it returns a worse oa in case a full factorial is possible
             cand <- cand[cand$nruns<ffnrun,]
             if (nrow(cand)==0){ 
                  if ((replications > 1 & repeat.only)){
                  if (minnrun <= ffnrun)
                      return(fac.design(nfactors=nfactors, nlevels=nlevels, factor.names=factor.names,
                         replications=replications, repeat.only=repeat.only, randomize=randomize, seed=seed))
                      else {warning("A full factorial without real replications does not fulfill your request for degrees of freedom. repeat.only has been set to FALSE.")
                           repeat.only <- FALSE}
                       }
                  return(fac.design(nfactors=nfactors, nlevels=nlevels, factor.names=factor.names,
                      replications=ceiling(minnrun/ffnrun), randomize=randomize, seed=seed))
                  }
             else {
                 if (cand[1,]$lineage == "") 
                       ID <- get(as.character(cand[1,1]))
                   else ID <- genChild(parseArrayLine(cand[1,]))
                       generating.oa <- as.character(cand[1,1])
             }
          }
      }  ## end if (is.null(ID) & is.null(generating.oa)), now both are non-null, if both were NULL
      
      if (is.null(ID)) ID <- genChild(parseArrayLine(oacat[which(oacat$name==generating.oa),]))

          des <- ID
          if (!("oa" %in% class(des))) stop("ID does not specify an orthogonal array.")
          if (!is.null(columns)) {
              if (is.numeric(columns)){
## bugfix: columns of length 1 threw an uncaptured error
                 if (length(columns)==1) stop("oa.design only works for at least two design columns")
                 if (!is.null(nfactors)){ 
                     if (!length(columns)==nfactors) stop("mismatch between columns and nfactors")}
                 else nfactors <- length(columns)
                 if (!is.null(nlevels)){
                     if (length(nlevels)==1) nlevels <- rep(nlevels,nfactors)
                     else if (!length(nlevels)==length(columns)) stop("mismatch between nlevels and columns") 
                         else if (any(!apply(des[,columns,drop=FALSE],2,function(obj) length(table(obj)))==nlevels))
                             stop("mismatch between nlevels and columns")}
                 else nlevels <- apply(des[,columns,drop=FALSE],2,function(obj) length(table(obj)))
                 }  ## end of numeric columns
              if (is.character(columns))
                 if (!columns %in% c("order","min3","min34","min3.rela","min34.rela","minPFT","minRPFT","minRelProjAberr")) 
                      stop("invalid choice for columns")
             }
          if (!(is.null(nruns))){
             if (!nrow(des)==nruns) 
                  stop("The design ", generating.oa, " has ", nrow(des), " runs, mismatch to specified nruns!")
             if (!is.numeric(nruns)) stop("nruns must be a number.")
             if (!nruns==floor(nruns)) stop("nruns must be an integer number.")
          }
          else (nruns <- nrow(des))
          if (!(is.null(nfactors))){
             if (!is.numeric(nfactors)) stop("nfactors must be a number.")
             if (!nfactors==floor(nfactors)) stop("nfactors must be an integer number.")
             if (nfactors > ncol(des)) stop("too many factors for design ", ID)
          }
          ## default: all columns are used, order is left to right
          if (is.null(nfactors) & is.null(factor.names) & is.null(nlevels) & !is.numeric(columns)){
               if (!(identical(columns,"order") | is.null(columns))) 
                  warning("All columns of the design are used, columns option was ignored")
               nfactors <- ncol(des)
               nlevels <- apply(des, 2, max)  ## assuming coding as 1:number of levels
               columns <- 1:ncol(des)
               factor.names <- as.list(rep("",nfactors))
               if (nfactors <= 50) names(factor.names) <- Letters[1:ncol(des)]
                    else names(factor.names) <- paste("F",1:nfactors,sep="")
               ## now, nfactors, nlevelsand factor.names are all non-null,
               ## and columns is numeric
               }  
         ## factor number is nfactors, otherwise given by (in this order of precedence) factor.names, columns, nlevels
         ## factor types is not permitted here, because at least one of these others also needs to be available
         if (is.null(nfactors)){ 
             if (!is.null(factor.names)) nfactors <- length(factor.names) else 
             if (is.numeric(columns)) nfactors <- length(columns) else 
             if (is.null(nfactors) & !is.null(nlevels)) 
                 if (length(nlevels) > 1) nfactors <- length(nlevels) else 
                 stop("A length 1 entry for nlevels can only be used,\nif at least factor.names, nfactors or a numeric vector for columns is specified as well.")
             if (is.null(nfactors))
                 stop("There is not enough information regarding the number of factors.")
             if (nfactors==1) stop("Designs with only one factor are not supported")
             }
               ## at least one of the above has an entry, i.e. nfactors is now non-null
               ## or an error has been thrown
## bug fix: ID and a length 1 nlevels threw an uncaptured error
         if (length(nlevels)==1) nlevels <- rep(nlevels, nfactors)
         
          

      ## correct number of valid factor types ?
#      if (!is.null(factor.types)) {
#           if (!length(factor.types)==nfactors) 
#               stop("If present, factor types must have an entry for each factor.")
#           if (any(!factor.types %in% c("qual","quant"))) 
#               stop("invalid entries in factor.types; only strings qual or quant are valid")
#           }
#      else factor.types <- rep("qual",nfactors)  ## default: qualitative

      ## check compatibility issues between different ways of specifying factor numbers
      ## and levels
      ## and determine appropriate values for null parameters
      ## nfactors is now non-null, since it has been adjusted otherwise  
         if (is.null(factor.names)) if (nfactors <= 50) factor.names <- Letters[1:nfactors]
             else factor.names <- paste("F", 1:nfactors, sep = "")
         ## factor.names is now also non-null

        if (!is.null(ID)) if (nfactors > ncol(des)) 
              stop("The design ", ID, " accomodates at most ", ncol(des), 
                    " factors, mismatch to specified nfactors!")
        if (!is.null(nlevels)) 
              if (!length(nlevels)==nfactors) 
                 stop("nlevels must have exactly one entry for each factor")
        if (!length(factor.names)==nfactors) 
                 stop("factor.names must have exactly one entry for each factor")
        if (is.numeric(columns)){
            if (any(table(columns) > 1)) stop("columns contains duplicates!")
            if (any(!columns %in% 1:ncol(des))) stop("invalid entry in columns")
            if (!is.null(nlevels)) 
              if (!all(apply(des[,columns],2,max)==nlevels)) 
                   stop("Mismatch between columns and nlevels!")
        }
      if (!(is.character(factor.names) | is.null(nlevels))){
          hilf <- sapply(factor.names,"length")
          if (!all(nlevels[hilf>1]==hilf[hilf>1])) 
               stop("mismatch between nlevels and level entries in factor.names")
      }  
      if (is.list(factor.names) & is.null(nlevels))
          nlevels <- sapply(factor.names,"length")  

      if (is.character(factor.names)) {
          hilf <- factor.names
          factor.names <- as.list(rep("",length(factor.names)))
          names(factor.names) <- hilf
      }   ## factor.names is now a list, but perhaps without level entries
      ## make factor names valid
      names(factor.names) <- make.names(names(factor.names), unique=TRUE)


      ### check whether the requested nlevels can be accomodated in design ID
           ### compare two tables, make sure that apples are compared to apples!
      if (!is.numeric(columns)){
         ## with columns has already been checked otherwise 
         hilf <- table(apply(des,2,max))
         if (is.null(nlevels)) nlevels <- sapply(factor.names,length)
         if (any(nlevels==1)) stop("If nlevels is not given, levels must be specified for all factors.")
         hilf2 <- table(nlevels)
         if (!all(names(hilf2) %in% names(hilf))) 
               stop("The chosen array does not offer factors with ", setdiff(names(hilf2),names(hilf)), " levels.")
         else if (length(hilf <- which(hilf[names(hilf2)] < hilf2))>0) 
               stop("The chosen design does not offer enough columns with ", hilf2[hilf], "-level factors.")
      }         

      ## nlevels should now be known
      ## make list factor.names always have level entries 
      for (i in 1:length(nlevels)) 
          if (identical(factor.names[[i]],"")) factor.names[[i]] <- 1:nlevels[i]
      
      ### arrange columns of oa in order needed for design
      if (is.numeric(columns)) {
         des <- des[,columns]
         desorigcode <- des
         origorder <- 1:nfactors
         names(nlevels) <- NULL
         names(columns) <- NULL}
      else{
          origorder <- (1:nfactors)[order(nlevels)]
          factor.names <- factor.names[order(nlevels)]
          nlevels <- nlevels[order(nlevels)]
          if ((columns=="order" | is.null(columns)) & !nfactors==ncol(ID))
             message("The columns of the array have been used in order of appearance. \n", 
             "For designs with relatively few columns, \n", 
             "the properties can sometimes be substantially improved \n", 
             "using option columns with min3 or even min34.\n")
          nutze <- NULL
          if (columns=="min3") nutze <- try(oa.min3(ID, nlevels)$column.variants[1,], silent=TRUE)
          else if (columns=="min3.rela") nutze <- try(oa.min34(ID, nlevels, rela=TRUE)$column.variants[1,], silent=TRUE)
          else if (columns=="min34") nutze <- try(oa.min34(ID, nlevels)$column.variants[1,], silent=TRUE)
          else if (columns=="min34.rela") nutze <- try(oa.min34(ID, nlevels, rela=TRUE)$column.variants[1,], silent=TRUE)
          else if (columns=="minRPFT") nutze <- try({
              hilf <- oa.min3(ID, nlevels, crit="worst", rela=TRUE)
              if (hilf[[1]]==4) hilf <- oa.maxGR(ID, nlevels, variants=hilf$column.variants)
              hilf2 <- which(bestPFT(matrix.fromPFTs(PFTs.from.variants(ID, hilf$column.variants))))[1]
              hilf$column.variants[hilf2,]
              }, silent=TRUE)
          else if (columns=="minRelProjAberr") nutze <- try(oa.minRelProjAberr(ID, nlevels)$column.variants[1,], silent=TRUE)
          if ("try-error" %in% class(nutze)) {
               columns <- "order"
               warning("resources were not sufficient for optimizing column selection")
               nutze <- NULL
          }
          ## select optimum columns (in the order of nlevels)
          if (!is.null(nutze)) des <- des[,nutze]
          ## assign design in requested order and column numbers
          columns <- numeric(0)
          hilf <- table(nlevels)
          for (i in as.numeric(names(hilf)))
              columns <- c(columns,which(apply(des,2,max)==i)[1:hilf[paste(i)]])
          des <- des[,columns]
          if (!is.null(nutze)) columns <- nutze[columns]
          names(columns) <- NULL
          names(nlevels) <- NULL
      desorigcode <- des[,order(origorder)]
      nlevels <- nlevels[order(origorder)]
      factor.names <- factor.names[order(origorder)]
#      factor.types <- factor.types[order(origorder)]
      columns <- columns[order(origorder)]
      }
      design <- as.data.frame(desorigcode)
      colnames(design) <- names(factor.names)
      quant <- sapply(factor.names, "is.numeric")
      for (i in 1:ncol(design)){
              ## recode 
              ## could later make (here and elsewhere)
              ## all factors if as.factor.result=TRUE
              ## all non-factors if FALSE
              ## character factors and numeric not if NULL
              ## in that case, quantitative factors as numeric contrasts
              #if (!is.numeric(factor.names[[i]])) 
              #    design[,i] <- des.recode(desorigcode[,i], paste(1:nlevels[i],"=",
              #       factor.names[[i]],sep="",collapse=";"),
              #       TRUE,TRUE)
              #else design[,i] <- des.recode(desorigcode[,i], paste(1:nlevels[i],"=",
              #       factor.names[[i]],sep="",collapse=";"),
              #       TRUE,FALSE)
              #    
              if (levordold){
              if (!is.numeric(factor.names[[i]])) 
                  design[,i] <- des.recode.old(desorigcode[,i], paste(1:nlevels[i],"=",
                     factor.names[[i]],sep="",collapse=";"),
                     TRUE,TRUE)
              else design[,i] <- des.recode.old(desorigcode[,i], paste(1:nlevels[i],"=",
                     factor.names[[i]],sep="",collapse=";"),
                     FALSE,FALSE)
              }
              else{
              if (!is.numeric(factor.names[[i]])) 
                  design[,i] <- des.recode(desorigcode[,i], paste(1:nlevels[i],"=",
                     factor.names[[i]],sep="",collapse=";"),
                     TRUE,TRUE)
              else design[,i] <- des.recode(desorigcode[,i], paste(1:nlevels[i],"=",
                     factor.names[[i]],sep="",collapse=";"),
                     FALSE,FALSE)
              }
              if (!is.factor(design[,i]))
                 design[,i] <- factor(design[,i],levels=factor.names[[i]]) 
              if (nlevels[i]==2) contrasts(design[,i]) <- contr.FrF2(2)
                    else if (quant[i]) contrasts(design[,i]) <- contr.poly(nlevels[i],scores=factor.names[[i]])
              ##  if (as.factor.result & factor.types[i]=="quant") 
              ##    if (!factor.names[i]=="") if (is.numeric(factor.names[i])) 
              ##            contrasts(design[,i]) <- contr.poly(nlevels[i], scores=factor.names[i])
              ##                      else contrasts(design[,i]) <- contr.poly(nlevels[i])
              }

      ## run.order is in line with the strategy for other designs
      rand.ord <- rep(1:nrow(design),replications)
      if (replications > 1 & repeat.only) rand.ord <- rep(1:nrow(design),each=replications)
      if (randomize & !is.null(seed)) set.seed(seed)
      if (randomize & !repeat.only) for (i in 1:replications) 
                  rand.ord[((i-1)*nrow(design)+1):(i*nrow(design))] <- sample(nrow(design))
      if (randomize & repeat.only) rand.ord <- rep(sample(1:nrow(design)), each=replications)

      aus <- design[rand.ord,]
      # GWP <- c("3"=length3(aus), "4"=length4(aus))  ## can be very resource intensive, stopped for the moment
      # resolution <- "III"
      
      ## extract run number in standard order
      ## remove uniqueness appendix
      orig.no <- orig.no.rp <- sapply(strsplit(rownames(aus),".",fixed=TRUE),function(obj) obj[1])
      if (replications>1) {
           if (repeat.only) orig.no.rp <- paste(orig.no.rp,rep(1:replications,nruns),sep=".")
           else orig.no.rp <- paste(orig.no.rp,rep(1:replications,each=nruns),sep=".")
       }
    ## row added 27 01 2011 (for proper ordering of design)
    orig.no.levord <- sort(as.numeric(orig.no),index=TRUE)$ix
      desmat <- model.matrix(~.,data=aus)[,-1,drop=FALSE]
      rownames(aus) <- rownames(desmat) <- 1:nrow(aus)

      attr(aus,"desnum") <- desmat
      ## change 27 Jan 2011: leave orig.no as a factor, but with better-ordered levels
      orig.no <- factor(orig.no, levels=unique(orig.no[orig.no.levord]))
      attr(aus,"run.order") <- data.frame("run.no.in.std.order"=orig.no,"run.no"=1:nrow(aus),"run.no.std.rp"=orig.no.rp)
      attr(aus,"design.info") <- list(type="oa",
              nruns=nruns,
              nfactors=nfactors,
              nlevels=nlevels,
              generating.oa=generating.oa,
              selected.columns=columns,
              origin=attr(ID,"origin"),
              comment=attr(ID,"comment"),
              residual.df = nruns - 1 - sum(nlevels) + length(nlevels),
              factor.names = factor.names, 
              replications=replications, 
              repeat.only=repeat.only, 
              randomize=randomize, 
              seed=seed, creator=creator)
      class(aus) <- c("design","data.frame")
      aus
}


