oa.maxGR <- function (ID, nlevels, variants=NULL)
{
    tab.needed <- table(nlevels)
    GR <- 3

    ## retrieve child array or array identified by character string
          ## gsub for case where ID is character string
    IDname <- gsub("\"","",deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)){
    if (!exists(IDname))
          ID <- eval(parse(text=paste("oa.design(",IDname,")")))
    else if (is.character(ID))
          ID <- eval(parse(text=paste("oa.design(",IDname,")")))
    }

    ## identify match between available and requested levels
    nlevID <- apply(ID, 2, function(obj) length(table(obj)))
    tab.available <- table(nlevID)[names(tab.needed)]
    if (any(is.na(names(tab.available)))) stop("not all levels can be accomodated")
    col.lists <- lapply(names(tab.needed), function(obj) which(nlevID ==
        as.numeric(obj)))
    spielraum <- tab.available - tab.needed
         if (any(spielraum < 0))
             stop("design does not have enough factors with ",
                  paste(names(spielraum)[which(spielraum<0)], collapse=" and "), " levels")

    ## provide candidate column list to be looped through
    cand.lists <- mapply(nchoosek, tab.available, tab.needed, SIMPLIFY=FALSE)
    cand.lists <- mapply(function(obj1, obj2) matrix(obj1[obj2],
        nrow = nrow(obj2), ncol = ncol(obj2)), col.lists, cand.lists,
        SIMPLIFY = FALSE)

    ## provide full factorial for all combinations of subsets,
    ## e.g. combining each variant of 3 2-level factors with each variant of 4 3-level factors
    hilf <- lapply(cand.lists, function(obj) 1:ncol(obj))
    hilf <- expand.grid(hilf)
    if (!is.null(variants)) hilf <- variants

    ## initialize curMax
    curMax <- -Inf
    MaxVariants <- numeric(0)
    MaxProj <- vector("list",0)
    for (i in 1:nrow(hilf)) {
        if (is.null(variants)) spalten <- c(unlist(mapply(function(obj1, obj2) obj1[,
            obj2], cand.lists, hilf[i, ])))
            else spalten <- hilf[i,]
        cur3 <- GR(ID[, spalten], digits=4)
        if (cur3$GR == curMax){
            MaxVariants <- rbind(MaxVariants, spalten)
            MaxProj <- c(MaxProj, list(cur3$RPFT))
            }
        else if (cur3$GR > curMax) {
            curMax <- cur3$GR
            MaxVariants <- matrix(spalten, nrow = 1)
            MaxProj <- list(cur3$RPFT)
        }
    }
    rownames(MaxVariants) <- 1:nrow(MaxVariants)
    list(GR = c(GR=curMax), column.variants = MaxVariants, RPFTs = MaxProj)
}

oa.minRelProjAberr <- function(ID, nlevels, maxGR=NULL){
    ## retrieve child array or array identified by character string
          ## gsub for case where ID is character string
    IDname <- gsub("\"","",deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)){
    if (!exists(IDname))
          ID <- eval(parse(text=paste("oa.design(",IDname,")")))
    else if (is.character(ID))
          ID <- eval(parse(text=paste("oa.design(",IDname,")")))
    }
    ## determine maxGR, if not handed to the function from previous call
     if (is.null(maxGR)) maxGR <- oa.min3(ID, nlevels, crit="worst", rela=TRUE)
     
     if (!is.list(maxGR)) stop("maxGR must be a list")
     if (!all(c("GR","column.variants") %in% names(maxGR)))
         stop("maxGR is not of the appropriate form")
     ## oa.min3 with crit="worst" not appropriate for resolution IV or higher designs
     GR <- maxGR$GR
     if (GR==4 & !"RPFTs" %in% names(maxGR)) maxGR <- oa.maxGR(ID, nlevels, variants=maxGR$column.variants)
     GR <- maxGR$GR
     
     if (GR==5) {
         hilf <- c("3"=0,"4"=0)
         aus <- list(GWP=hilf, column.variants=maxGR$column.variants, complete=TRUE)
     }
     else{
     ## the more frequent case of a resolution 3 or 4 design
     reso <- floor(maxGR$GR)
     ## one single entry only
     if (nrow(maxGR$column.variants)==1) {
         if (reso==3)
         hilf <- c("3.relative"=length3(ID[,maxGR$column.variants], rela=TRUE),"4"=length4(ID[,maxGR$column.variants]))
         if (reso==4)
         hilf <- c("3"=0,"4.relative"=length4(ID[,maxGR$column.variants], rela=TRUE))
         aus <- list(GWP=hilf, column.variants=maxGR$column.variants, complete=TRUE)
     }
     else{
     ## reduce maxGR to best rA3/rA4 design
     if (reso==3) minrA <- oa.min3(ID, nlevels, variants=maxGR$column.variants, rela=TRUE)
     else if (reso==4) minrA <- oa.min34(ID, nlevels, variants=maxGR$column.variants, rela=TRUE)
     
     if (!"RPFTs" %in% names(maxGR)){
         maxGR$column.variants <- minrA$column.variants
         maxGR$RPFTs <- PFTs.from.variants(ID, maxGR$column.variants, R=reso, rela=TRUE)
         }
     else{
       zeile <- 1
       hilf <- maxGR$column.variants
       auswahl <- rep(FALSE, nrow(hilf))
       for (i in 1:nrow(minrA$column.variants)){
          found <- FALSE
          while (!found){
          if (all(minrA$column.variants[i,]==hilf[zeile,])) {
             auswahl[zeile] <- TRUE
             zeile <- zeile + 1
             found <- TRUE
             }
          }
       }
       maxGR$column.variants <- maxGR$column.variants[auswahl,,drop=FALSE]
       maxGR$RPFTs <- maxGR$RPFTs[auswahl]
     }
         
     ## optimizing RPFTs
     RPFTs <- maxGR$RPFTs
     best <- which(bestPFT(matrix.fromPFTs(RPFTs)))
     maxGR$column.variants <- maxGR$column.variants[best,,drop=FALSE]
     maxGR$RPFTs <- maxGR$RPFTs[best]
     
     ## resolving final ties with A4
     if (length(best)>1 & reso==3){ 
        maxGR[[1]] <- length3(ID[,maxGR$column.variants[1,]], rela=TRUE)
        names(maxGR)[1] <- "GWP3"; names(maxGR[[1]]) <- "rA3"
        maxGR$complete <- TRUE
        maxGR <- oa.min34(ID, nlevels=nlevels, min3=maxGR)
        }
     
    aus <- maxGR
     }
     }
   c(GR=GR, aus)
}