oa.min3 <- function (ID, nlevels, all = FALSE, rela = FALSE, variants=NULL, crit="total")
{
    if (!is.logical(rela)) stop("rela must be logical")
    if (!is.logical(all)) stop("all must be logical")
    if (! crit %in% c("total","worst")) stop("invalid crit")
  ## might make sense to have more versions here 
  ##  if (! crit %in% c("total","worst","wt","tw")) stop("invalid crit")
    tab.needed <- table(nlevels)
    ## retrieve child array or array identified by character string
          ## gsub for case where ID is character string
    IDname <- gsub("\"", "", deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)) {
        if (!exists(IDname))
            ID <- eval(parse(text = paste("oa.design(", IDname,
                ", randomize=FALSE)")))
        else if (is.character(ID))
            ID <- eval(parse(text = paste("oa.design(", IDname,
                ", randomize=FALSE)")))
    }
    ## identify match between available and requested levels
    nlevID <- apply(ID, 2, function(obj) length(table(obj)))
    if (!is.null(variants)){
       if (!is.numeric(variants)) stop("If given, variants must be numeric.")
       if (!is.matrix(variants)) stop("If given, variants must be a matrix.")
       nr <- nrow(variants); nc <- ncol(variants)
       if (!nc==length(nlevels)) stop("variants has the wrong length")
       if (!all(matrix(nlevID[variants],nr,nc)==matrix(nlevels,nr,nc,byrow=TRUE))) 
            stop("variants has invalid entries")
    }
    tab.available <- table(nlevID)[names(tab.needed)]
    if (any(is.na(names(tab.available))))
        stop("not all levels can be accomodated")
    col.lists <- lapply(names(tab.needed), function(obj) which(nlevID ==
        as.numeric(obj)))
    spielraum <- tab.available - tab.needed
    if (any(spielraum < 0))
        stop("design does not have enough factors with ", paste(names(spielraum)[which(spielraum <
            0)], collapse = " and "), " levels")
    triples <- nchoosek(length(nlevels), 3)  ## for later selection from spalten
    ## provide candidate column list to be looped through
    if (!is.null(variants)) hilf <- t(variants)
    else{
    cand.lists <- mapply(nchoosek, tab.available, tab.needed,
        SIMPLIFY = FALSE)
    cand.lists <- mapply(function(obj1, obj2) matrix(obj1[obj2],
        nrow = nrow(obj2), ncol = ncol(obj2)), col.lists, cand.lists,
        SIMPLIFY = FALSE)
    ## provide full factorial for all combinations of subsets,
    ## e.g. combining each variant of 3 2-level factors with each variant of 4 3-level factors
    hilf <- lapply(cand.lists, function(obj) 1:ncol(obj))
    hilf <- expand.grid(hilf)
    }
    curMin <- Inf
    ## would be for two-part crit
    ## if (crit %in% c("tw","wt")) curMin <- c(Inf,Inf)
    MinVariants <- numeric(0)
    ll <- length3(ID, J=TRUE)  # J characteristics
    ll <- split(ll, names(ll))
    ll <- sapply(ll, function(obj) sum(obj^2))  ## length 3 words for all triples
    div <- sapply(names(ll), function(obj) min((nlevID-1)[as.numeric(unlist(strsplit(obj,":")))]))
    if (rela) ll <- ll/div   ## relative length 3 words for all triples
    if (is.null(variants))
    hilf <- apply(hilf,1,function(obj) c(unlist(mapply(function(obj1, obj2) obj1[,
            obj2], cand.lists, obj))))
    ## hilf is a matrix, the columns of which contain column selections to be checked
        ## spalten contains all column numbers of the current choice
        ## ordered by number of levels (for unspecified variants)
    for (i in 1:ncol(hilf)) {
        spalten <- hilf[,i]
        nam.need <- apply(triples, 2, function(obj) paste(sort(spalten[obj]), collapse=":"))
        if (crit=="total") cur3 <- round(sum(ll[nam.need]),4)
        else cur3 <- round(max(ll[nam.need]),4)
        ## else if (crit=="wt") cur3 <- c(round(max(ll[nam.need]),4),round(sum(ll[nam.need]),4))
        ## else if (crit=="tw") cur3 <- c(round(sum(ll[nam.need]),4),round(max(ll[nam.need]),4))
        if (cur3 == curMin)
            MinVariants <- rbind(MinVariants, spalten)
        ## else if (cur3[1] < curMin[1] | (cur3[1] <= curMin[1] & cur3[2]<curMin[2])){
        else if (cur3 < curMin){
            curMin <- cur3
            MinVariants <- matrix(spalten, nrow = 1)
            if (curMin == 0 & !all){
                if (crit=="total")
                return(list(GWP3 = 0, column.variants = matrix(spalten,
                  nrow = 1), complete = FALSE))
            else
                return(list(GR = 4, column.variants = matrix(spalten,
                  nrow = 1), complete = FALSE))
            }
        }
    }
    rownames(MinVariants) <- 1:nrow(MinVariants)
    if (!rela){
    if (crit=="total")
        names(curMin) <- "A3"
    else names(curMin) <- "worst"
    ## else if (crit=="wt") names(curMin) <- c("worst", "A3")
    ## else if (crit=="tw") names(curMin) <- c("A3","worst")
    }
    else {
        if (crit=="total") names(curMin) <- "rA3"
        else {
          curMin <- 3+1-round(sqrt(curMin),2)
          names(curMin) <- "GR"
        }
        ##if (crit=="wt") {
        ##  curMin[1] <- 3+1-round(sqrt(curMin[1]),2)
        ##  names(curMin=c("GR","rA3"))
        ##}
        ##if (crit=="tw") {
        ##  curMin[2] <- 3+1-round(sqrt(curMin[2]),2)
        ##  names(curMin=c("rA3","GR"))
        ##}
    }
    aus <- list(optimum = curMin, column.variants = MinVariants, complete = TRUE)
    if (crit == "total" & rela) names(aus)[1] <- "GWP3"
    if (crit == "total" & !rela) names(aus)[1] <- "GWP3"
    if (crit == "worst" & !rela) names(aus)[1] <- "worst.a3"
    if (crit == "worst" & rela) names(aus)[1] <- "GR"
    aus
}

oa.max3 <- function (ID, nlevels, rela = FALSE)
{
    tab.needed <- table(nlevels)
    ## retrieve child array or array identified by character string
          ## gsub for case where ID is character string
    IDname <- gsub("\"", "", deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)) {
        if (!exists(IDname))
            ID <- eval(parse(text = paste("oa.design(", IDname,
                ")")))
        else if (is.character(ID))
            ID <- eval(parse(text = paste("oa.design(", IDname,
                ")")))
    }
    ## identify match between available and requested levels
    nlevID <- apply(ID, 2, function(obj) length(table(obj)))
    tab.available <- table(nlevID)[names(tab.needed)]
    if (any(is.na(names(tab.available))))
        stop("not all levels can be accomodated")
    col.lists <- lapply(names(tab.needed), function(obj) which(nlevID ==
        as.numeric(obj)))
    spielraum <- tab.available - tab.needed
    if (any(spielraum < 0))
        stop("design does not have enough factors with ", paste(names(spielraum)[which(spielraum <
            0)], collapse = " and "), " levels")
    triples <- nchoosek(length(nlevels), 3)  ## for later selection from spalten

    ## provide candidate column list to be looped through
    cand.lists <- mapply(nchoosek, tab.available, tab.needed,
        SIMPLIFY = FALSE)
    cand.lists <- mapply(function(obj1, obj2) matrix(obj1[obj2],
        nrow = nrow(obj2), ncol = ncol(obj2)), col.lists, cand.lists,
        SIMPLIFY = FALSE)

    ## provide full factorial for all combinations of subsets,
    ## e.g. combining each variant of 3 2-level factors with each variant of 4 3-level factors
    hilf <- lapply(cand.lists, function(obj) 1:ncol(obj))
    hilf <- expand.grid(hilf)

    ## initialize curMax
    curMax <- -Inf
    MaxVariants <- numeric(0)
    ll <- length3(ID, J=TRUE)  # J characteristics
    ll <- split(ll, names(ll))
    ll <- sapply(ll, function(obj) sum(obj^2))  ## length 3 words for all triples
    div <- sapply(names(ll), function(obj) min((nlevID-1)[as.numeric(unlist(strsplit(obj,":")))]))
    if (rela) ll <- ll/div   ## relative length 3 words for all triples
    hilf <- apply(hilf,1,function(obj) c(unlist(mapply(function(obj1, obj2) obj1[,
            obj2], cand.lists, obj))))
    ## hilf is a matrix, the columns of which
        ## spalten contain all column numbers of the current choice
        ## ordered by number of levels
    for (i in 1:ncol(hilf)) {
        spalten <- hilf[,i]
        nam.need <- apply(triples, 2, function(obj) paste(sort(spalten[obj]), collapse=":"))
        cur3 <- round(sum(ll[nam.need]),4)
        if (cur3 == curMax)
            MaxVariants <- rbind(MaxVariants, spalten)
        else if (cur3 > curMax) {
            curMax <- cur3
            MaxVariants <- matrix(spalten, nrow = 1)
        }
    }
    rownames(MaxVariants) <- 1:nrow(MaxVariants)
    if (rela)
        names(curMax) <- "3.relative"
    else names(curMax) <- "3"
    list(GWP3 = curMax, column.variants = MaxVariants, complete = TRUE)
}


oa.max4 <- function (ID, nlevels, rela=FALSE)
{
## oa.max4 only makes sense for selection from a resolution IV design
    tab.needed <- table(nlevels)

    ## retrieve child array or array identified by character string
          ## gsub for case where ID is character string
    IDname <- gsub("\"","",deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)){
    if (!exists(IDname))
          ID <- eval(parse(text=paste("oa.design(",IDname,", randomize=FALSE)")))
    else if (is.character(ID))
          ID <- eval(parse(text=paste("oa.design(",IDname,", randomize=FALSE)")))
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
    quadruples <- nchoosek(length(nlevels), 4)  ## for later selection from spalten

    ## provide candidate column list to be looped through
    cand.lists <- mapply(nchoosek, tab.available, tab.needed, SIMPLIFY=FALSE)
    cand.lists <- mapply(function(obj1, obj2) matrix(obj1[obj2],
        nrow = nrow(obj2), ncol = ncol(obj2)), col.lists, cand.lists,
        SIMPLIFY = FALSE)

    ## provide full factorial for all combinations of subsets,
    ## e.g. combining each variant of 3 2-level factors with each variant of 4 3-level factors
    hilf <- lapply(cand.lists, function(obj) 1:ncol(obj))
    hilf <- expand.grid(hilf)

    ## initialize curMax
    curMax <- -Inf
    MaxVariants <- numeric(0)
    ll <- length4(ID, J=TRUE)  # J characteristics
    ll <- split(ll, names(ll))
    ll <- sapply(ll, function(obj) sum(obj^2))  ## length 4 words for all quadruples
    div <- sapply(names(ll), function(obj) min((nlevID-1)[as.numeric(unlist(strsplit(obj,":")))]))
    if (rela) ll <- ll/div   ## relative length 4 words for all quadruples
    hilf <- apply(hilf,1,function(obj) c(unlist(mapply(function(obj1, obj2) obj1[,
            obj2], cand.lists, obj))))
    ## hilf is a matrix, the columns of which
        ## spalten contain all column numbers of the current choice
        ## ordered by number of levels
    for (i in 1:ncol(hilf)) {
        spalten <- hilf[,i]
        nam.need <- apply(quadruples, 2, function(obj) paste(sort(spalten[obj]), collapse=":"))
        cur4 <- round(sum(ll[nam.need]),4)
        if (cur4 == curMax)
            MaxVariants <- rbind(MaxVariants, spalten)
        else if (cur4 > curMax) {
            curMax <- cur4
            MaxVariants <- matrix(spalten, nrow = 1)
        }
    }
    rownames(MaxVariants) <- 1:nrow(MaxVariants)
    if (rela) names(curMax) <- "4.relative" else names(curMax) <- "4"
    list(GWP4 = curMax, column.variants = MaxVariants, complete = TRUE)
}