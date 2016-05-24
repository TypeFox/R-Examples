ccd.augment <- function(cube, ncenter = 4, columns="all", block.name="Block.ccd",
        alpha = "orthogonal", randomize=TRUE, seed=NULL, ...){
    ## design carries info on all factors etc.
    ## ccd can generate full design and then throw the cube portion away
    ## does not handle splitplot designs
    ## does not allow to expand the numeric factors in a design with mixed factor types
    
    ## would it make sense (be required) to allow bbreps to be different for star portion ?
    ## currently, a blocked design would always have one star block only
    creator <- sys.call()
    di <- design.info(cube)
   # if (length(grep("blocked",di$type))>0) stop("ccd.augment currently does not work for blocked designs")
    if (length(grep("splitplot",di$type))>0) stop("ccd.augment does not work for split-plot designs")
    if (!(substr(di$type,1,4)== "FrF2" | (length(grep("full factorial",di$type))>0 & all(di$nlevels==2))))
        stop("this is not a regular (fractional) factorial 2-level design")
    if (!columns=="all") stop("columns has not yet been implemented")
    if (!is.numeric(ncenter)) stop("ncenter must be numeric")
    if (!all(ncenter==floor(ncenter))) stop("ncenter must be integer")
    if (length(grep("blocked",di$type))> 0) 
       if (di$bbreps*di$wbreps > 1) stop("replicated blocked designs can not yet be treated with function ccd.augment")

    if (length(grep("center",di$type))==0){
      cube <- add.center(cube, ncenter[1])
      di <- design.info(cube)
    }

    more <- setdiff(colnames(cube),c(names(di$factor.names),di$block.name))
    moredn <- more
#    moredn <- setdiff(colnames(desnum(cube)),names(di$factor.names)) ## as block factor has different names!
    planvars <- colnames(cube)
    if (length(more)>0){
       addedvars <- cube[,more,drop=FALSE]
       planvars <- setdiff(planvars, more)
       }

    if (length(ncenter)==1) ncenter <- c(di$ncenter,ncenter)
    if (!length(ncenter)==2) stop("ncenter must have one or two elements")
       else if (!ncenter[1]==di$ncenter) stop("the first element of ncenter must correspond to number of center points in cube")
    nfactors <- di$nfactors
    factor.names <- di$factor.names
    bbreps <- di$bbreps
    if (is.null(bbreps)) if (!di$repeat.only) bbreps <- di$replications else
             stop("designs with repeat.only replications cannot be augmented to become ccd designs")
    wbreps <- di$wbreps
    if (is.null(wbreps)) wbreps <- 1
    #if (!di$repeat.only) wbreps <- di$replications else
    #         stop("designs with repeat.only replications cannot be augmented to become ccd designs")
    nlevels <- rep(2, nfactors)
    n.c <- di$ncube
    k <- round(log2(di$ncube))
    if (nfactors>k)
    generators <- generators.from.design(cube)
    else generators <- "full factorial" ## problem: if early vector generated from later one, ccd cannot handle this
                                                ## shuffle back and forth to solve this issue!!!

    if (randomize & !is.null(seed)) set.seed(seed)

    if (length(grep("estimable",di$type))>0)
        map <- di$map
        else map <- list(1:nfactors)

  ### cases without blocks
    if (is.null(di$block.gen)){
    ### treat case with generators
    if (!k >= nfactors){
        aus <- .ccd.1.41(k,
            generators=generators,
            blocks = block.name, n0 = ncenter, alpha = alpha,
            wbreps = wbreps, bbreps = bbreps, randomize = randomize,
            coding = make.formulas(paste("x",1:nfactors,sep=""),factor.names[map[[1]]])
        #      coding = make.formulas(paste("x",map[[1]],sep=""),factor.names[map[[1]]])
            )
        }
    else{
    ### treat case without generators
    if (k>=nfactors) wbreps <- 2^(k-nfactors)*wbreps
    aus <- .ccd.1.41(nfactors,
        blocks = block.name, n0 = ncenter, alpha = alpha,
        wbreps = wbreps, bbreps = bbreps, randomize = randomize,
        coding = make.formulas(paste("x",map[[1]],sep=""),factor.names[map[[1]]])
        )
        }
    }
    else{
    ### treat case with blocks
    ### LHS is the block name
    ### RHS is vector of block generator formulae (no warnings, if these are dependent)
    block.form <- di$block.gen
    if (is.vector(block.form)) block.form <- Yates[block.form]
    block.form <- paste("c(",paste(sapply(block.form, function(obj) paste(paste("x", obj, sep=""),collapse="*")),collapse=","),")")
    if (!k >= nfactors){
        aus <- .ccd.1.41(k,
            generators=generators,
            blocks = as.formula(paste(block.name,block.form,sep="~")), n0 = ncenter, alpha = alpha,
            wbreps = wbreps, bbreps = bbreps, randomize = randomize,
            coding = make.formulas(paste("x",map[[1]],sep=""),factor.names[map[[1]]])
            )
        }
    else{
    ### treat case without generators
    if (k > nfactors) wbreps <- 2^(k-nfactors)*wbreps
    aus <- .ccd.1.41(nfactors,
        blocks = as.formula(paste(block.name,block.form,sep="~")), n0 = ncenter, alpha = alpha,
        wbreps = wbreps, bbreps = bbreps, randomize = randomize,
        coding = make.formulas(paste("x",map[[1]],sep=""),factor.names[map[[1]]])
        )
        }
    }
  #  center.positions.cube <- which(apply(desnum[1:((n.c*wbreps+ncenter[1])*bbreps),],1,function(obj) sum(abs(obj))==0))
    if (is.null(di$blocks)) nblocks <- 1 else nblocks <- di$nblocks
    star.points <- ((n.c*wbreps+ncenter[1]*nblocks)*bbreps+1):nrow(aus)

    ### uses block numbers generated by ccd to make row names unique
    ### relates to run number in standard order as created by FrF2, which is different from that in ccd
    rostd <- as.character(run.order(cube)$run.no.in.std.order)   ## change 18 Feb 2013
    if (!is.numeric(rostd)) rostd <- sapply(strsplit(rostd, ".", fixed=TRUE), function(obj) as.numeric(obj[1]))
    rn <- c(paste(paste("C", aus[[block.name]][-star.points],sep=""), rostd,sep="."), rownames(aus)[star.points])
    ## error in R CMD check, but not in R itself; why?
    cubecenter <- which(!iscube(cube))
    rn[cubecenter] <- paste(paste("C", as.character(aus[[block.name]])[cubecenter],sep=""),(n.c+1):(n.c+ncenter[1]),sep=".")
    design <- decode.data(aus)[,-1]   #[,FrF2 : : : invperm(map[[1]])]
    if (length(more)>0) design <- cbind(design, matrix(NA, nrow=nrow(design), ncol=length(more), dimnames=list(rn, more)))
    design <- rbind(cube[,c(names(factor.names),more)],design[star.points,])
    design <- cbind(aus[[block.name]],design)
    colnames(design)[1] <- block.name
    rownames(design) <- rn
    desnum <- coded.data(design, formulas=attr(aus,"coding"))
    class(design) <- c("design","data.frame")

    #desnum <- aus[,-1]
    attr(desnum,"codings") <- NULL
    #desnum <- desnum[,FrF2 : : : invperm(map[[1]])]   ## rearrange columns to match factor order in case of estimable
    #colnames(desnum) <- names(factor.names)
    #if (length(more)>0) 
        #desnum <- cbind(desnum, matrix(NA, nrow=nrow(desnum), ncol=length(moredn), dimnames=list(rownames(desnum), moredn)))
    desnum <- model.matrix(~.,model.frame(~.,desnum, na.action=na.pass))[,-1]
    #hilf <- desnum(cube)[,c(names(factor.names),moredn)]
    #oldnames <- colnames(desnum)
    #colnames(hilf)[1:nfactors] <- oldnames 
    #desnum <- rbind(hilf,desnum[star.points,])
    #desnum <- cbind(aus[,1],desnum)
    colnames(desnum)[1] <- block.name
    colnames(desnum)[2:(1+nfactors)] <- names(factor.names)
    
    desnum(design) <- desnum
    ## bug fix Nov 15: aus replaced by design (mismatch of run.order with row.names)
    run.order(design) <- data.frame(run.no.in.std.order=row.names(design),run.no=1:nrow(design),run.no.std.rp =rownames(design))
    di$type <- "ccd"
    di$block.name <- block.name
    di$coding <- lapply(attr(aus,"coding"),"as.formula",env=NULL)[map[[1]]]
    di$cube.gen <- generators
    di$nstar <- (nrow(aus)/bbreps-sum(ncenter))/wbreps-n.c
    di$ncenter <- ncenter
    di$creator <- append(di$creator, creator)
    di$nruns <- sum(di$ncenter)+di$ncube+di$nstar
    design.info(design) <- di
    design
}
