lhs.design <- function(nruns, nfactors, type="optimum", factor.names=NULL, seed=NULL, digits=NULL, nlevels=nruns, 
     default.levels = c(0,1), randomize = FALSE, ...){
     creator <- sys.call()
     if (!type %in% c("genetic","improved","maximin","optimum","random","dmax","faure","strauss","fact"))
           stop("invalid type")

## ??? open: error control for specification of required arguments for dmaxDesign and straussDesign

     ## calls functions geneticLHS ... randomLHS from package lhs
     ## or              dmaxDesign, runif.faure, straussDesign or factDesign from package DiceDesign
     ## additional parameters:  pop, gen and pMut for genetic
     ##                         dup for improved and maximin
     ##                         maxSweeps and eps for optimum
     ##                         none for random and runif.faure
     ##                         range and niter_max for dmaxDesign
     ##                         RND for straussDesign
     ##                         randomize for factDesign
     if (!is.null(digits)) if (!(length(digits)==nfactors | length(digits)==1)){
         stop("if specified, digits must be a single number or must contain an entry for each factor")
         if (!is.numeric(digits)) stop("digits must be numeric")
         if (!all(floor(digits)==digits)) stop("digits must be integer")
         }
     ## nruns can be missing for type "fact", otherwise it must be identical to nlevels or already correctly calculated
     if (missing(nruns) & !type=="fact") stop("nruns must be specified for all types except fact")
     if (missing(nfactors)) stop("nfactors must be specified")
     if (missing(nruns)){ 
           if (!is.numeric(nlevels)) stop("nlevels must be numeric")
           if (!length(nlevels) %in% c(1,nfactors)) stop("nlevels must be a number or a vector of length nfactors.")
           if (!all(nlevels%%1==0)) stop("All entries of nlevels must be integer numbers.")
           if (length(nlevels)==1) nlevels <- rep(nlevels, nfactors)
           nruns <- prod(nlevels)
         }

     if (!(is.numeric(nruns) & is.numeric(nfactors))) stop("nruns and nfactors must be numeric")
     if (!length(nruns)==1) stop("nruns must be a single number")
     if (!(nruns%%1==0 & nfactors%%1==0)) stop("nruns and nfactors must be integer")
     
     ## treat the case for type "fact", where nlevels is given in nruns
     if (type=="fact" & identical(nruns,nlevels)){
            nlevels <- rep(nlevels, nfactors)
            nruns <- prod(nlevels)
         }

     if (is.null(factor.names)){ 
          factor.names <- rep(list(default.levels),nfactors)
          names(factor.names) <- paste("X",1:nfactors,sep="")
          }
     if (!length(factor.names)==nfactors) stop("factor.names must have one entry for each factor.")
     if (is.list(factor.names)){
        if (!all(sapply(factor.names,"is.numeric"))) 
               stop("The list factor.names must give the scale ends for the quantitative variables, if given.")
        if (!all(sapply(factor.names,"length")==2)) 
               stop("The list factor.names must give the scale ends for the quantitative variables, if given.")
        if (is.null(names(factor.names))) names(factor.names) <- paste("X",1:nfactors,sep="")
        }
     if (is.character(factor.names)){
        hilf <- rep(list(default.levels),nfactors)
        names(hilf) <- factor.names
        factor.names <- hilf
     }
     ## make all factor names valid R names
     names(factor.names) <- make.names(names(factor.names), unique=TRUE)

     if (!is.null(seed)) set.seed(seed)
     ## determine lhs design
     if (type %in% c("genetic","improved","maximin","optimum","random"))
        desnum <- eval(parse(text=paste(type,"LHS(n=nruns,k=nfactors, ...)",sep="")))
     else{ 
        ## determine DiceDesign design
        if (type %in% c("dmax","strauss"))
           DD <- eval(parse(text=paste(type,"Design(n=nruns,dimension=nfactors, ...)",sep="")))
        if (type =="fact")
           DD <- factDesign(levels=nlevels,dimension=nfactors)
        if (type =="faure")
           DD <- runif.faure(n=nruns,dimension=nfactors)
        desnum <- DD$design
     }
        ## colnames needed as V1, ..., for recoding
        colnames(desnum) <- paste("V",1:nfactors,sep="")
        rownames(desnum) <- 1:nrow(desnum)
        design <- as.data.frame(desnum)
        design <- decode.data(coded.data(2*design-1, formulas = make.formulas(paste("V",1:nfactors,sep=""),factor.names)))
        ## colnames back to standard colnames
        colnames(desnum) <- names(factor.names)
        if (!is.null(digits)) if (length(digits)==1) design <- round(design, digits) else
            for (i in 1:nfactors) design[,i] <- round(design[,i],digits=digits[i])
        class(design) <- c("design",class(design))
        desnum(design) <- desnum
     if (type %in% c("genetic","improved","maximin","optimum","random"))
        di <- list(type="lhs", subtype=paste(type,"LHS",sep=""),
           nruns=nruns, nfactors=nfactors, factor.names=factor.names, 
           quantitative=rep(TRUE,nfactors), randomize=TRUE, seed=seed, replications=1, repeat.only=FALSE, digits=digits,
           creator=creator)
     else di <- append(list(type="lhs", subtype=paste(type,"DiceDesign",sep="."),
           nruns=nruns, nfactors=nfactors, factor.names=factor.names, 
           quantitative=rep(TRUE,nfactors), randomize=TRUE, seed=seed, replications=1, repeat.only=FALSE, digits=digits,
           creator=creator),DD[-(1:3)])
     if (type %in% c("fact","faure")){
         if (randomize){
             if (!is.null(seed)) set.seed(seed)
             rand.ord <- sample(nrow(design))
             design <- design[rand.ord,]
         }
         di$randomize <- randomize 
         if (type=="fact"){ 
            di$nlevels <- di$levels
            di$levels <- NULL
            }
         }
     ## append optimality information
     di$optimality.criteria <- list(S=Scalc(desnum),
                          mindist=mindist(desnum),
                          coverage=coverage(desnum),
                          det.cor=det(cor(desnum)))
## omitted meshRatio because of excessive duration

     ## add coding information for potential use with function rsm from package rsm
     di$coding <- lapply(make.formulas(paste("x",1:nfactors,sep=""), factor.names),
                       "as.formula",env=NULL)
     names(di$coding) <- paste("x",1:nfactors,sep="")
     design.info(design) <- di

        run.no.in.std.order <- run.no <- run.no.std.rp <- 1:nruns
        run.order(design) <- data.frame(run.no.in.std.order, run.no, run.no.std.rp)
        
        design
}

lhs.augment <- function(lhs, m=1, type="optAugment", seed=NULL, ...){
        ## type is one of augment, optSeeded or optAugment
        ## additional parameters:   maxSweeps, eps for optSeeded
        ##                          Mult for optAugment
        ## lhs is a design created before with lhs.design
        ## m is the number of additional points
        creator <- sys.call()
        if (!type %in% c("augment", "optSeeded", "optAugment"))
           stop("invalid type")
        if (!"design" %in% class(lhs)) stop("lhs must be of class design")
        di <- design.info(lhs)
        if (!di$type=="lhs") stop("lhs must be a design of type lhs")
        if (type=="optSeeded" & !is.null(di$response.names)) 
                stop("lhs.augment with type optSeeded does not work for designs with response data")

        factor.names <- di$factor.names
        digits <- di$digits
        nfactors <- di$nfactors
        nruns <- di$nruns
        oldseed <- di$seed
        if (!is.list(oldseed)) oldseed <- list(oldseed)
        if (!is.null(seed)) set.seed(seed)
        ## apply the appropriate augment function to existing lhs design
        ## must apply to desnum, because requires coding in unit cube
        ## must restrict to factors
        desnum <- eval(parse(text=paste(type,"LHS(desnum(lhs)[,names(factor.names)],m=m, ...)",sep="")))
        ## colnames needed as V1, ..., for recoding
        colnames(desnum) <- paste("V",1:nfactors,sep="")
        rownames(desnum) <- c(rownames(lhs), (nrow(lhs)+1):nrow(desnum))
        ## append previously added data
        design <- coded.data(as.data.frame(desnum*2-1), formulas = make.formulas(paste("V",1:nfactors,sep=""),factor.names))
        design <- decode.data(design)
        ## colnames back to standard colnames
        colnames(desnum) <- names(factor.names)
        if (!is.null(digits)) if (length(digits)==1) design <- round(design, digits) else
            for (i in 1:nfactors) design[,i] <- round(design[,i],digits=digits[i])
        diffnam <- setdiff(colnames(lhs), names(factor.names))
        ladd <- length(diffnam)
        if (ladd > 0){
            if (type=="optSeeded") 
                stop("lhs.augment with type optSeeded does not work for designs that have data additional to the experimental variables")
           appenddesnum <- rbind(desnum(lhs)[,diffnam,drop=FALSE], matrix(NA,nrow=m, ncol=ladd, 
             dimnames = list(NULL, diffnam)))
           desnum <- cbind(desnum, appenddesnum)
           appenddesign <- rbind(undesign(lhs)[,diffnam,drop=FALSE], matrix(NA,nrow=m, ncol=ladd, 
             dimnames = list(NULL, diffnam)))
           design <- cbind(design, appenddesign)
           }
        class(design) <- c("design",class(design))
        desnum(design) <- desnum
        di$creator <- list(di$creator, creator)
        di$seed <- append(oldseed, seed)
        di$nruns <- di$nruns+m
        di$subtype <- c(di$subtype, paste(type,"LHS",m,sep=""))
        design.info(design) <- di
        run.no.in.std.order <- run.no <- run.no.std.rp <- c(run.order(lhs)$run.no, (nruns+1):(nruns+m))
        run.order(design) <- data.frame(run.no.in.std.order, run.no, run.no.std.rp)
        design
}
