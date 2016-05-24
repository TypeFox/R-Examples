nextItem<-function (itemBank, model = NULL, theta = 0, out = NULL, x = NULL, 
    criterion = "MFI", method = "BM", priorDist = "norm", priorPar = c(0, 
        1), D = 1, range = c(-4, 4), parInt = c(-4, 4, 33), infoType = "observed", 
    randomesque = 1, rule = "length", thr = 20, SETH = NULL, 
    AP = 1, nAvailable = NULL, maxItems = 50, cbControl = NULL, 
    cbGroup = NULL) 
{
    crit <- switch(criterion, MFI = "MFI", bOpt = "bOpt", MLWI = "MLWI", 
        MPWI = "MPWI", MEI = "MEI", MEPV = "MEPV", random = "random", 
        progressive = "progressive", proportional = "proportional", 
        KL = "KL", KLP = "KLP", thOpt = "thOpt", GDI = "GDI", GDIP = "GDIP")
    if (is.null(crit)) 
        stop("invalid 'criterion' name", call. = FALSE)
    if (!is.null(model)) {
        mod <- switch(model, GRM = 1, MGRM = 2, PCM = 3, GPCM = 4, 
            RSM = 5, NRM = 6)
        if (is.null(mod)) 
            stop("invalid 'model' type!", call. = FALSE)
    }
    if (is.null(cbControl)) 
        OUT <- out
    else {
        if (is.null(cbGroup)) 
            stop("'cbGroup' argument must be provided for content balancing!", 
                call. = FALSE)
        if (sum(cbControl$props) != 1) 
            cbControl$props <- cbControl$props/sum(cbControl$props)
        nrGroup <- length(cbControl$names)
        if (is.null(out)) 
            empProp <- rep(0, nrGroup)
        else {
            empProp <- NULL
            for (i in 1:nrGroup) empProp[i] <- length(out[cbGroup[out] == 
                cbControl$names[i]])
            empProp <- empProp/sum(empProp)
        }
        thProp <- cbControl$props
        if (min(empProp) == 0) {
            indGroup <- (1:nrGroup)[empProp == 0]
            selGroup <- ifelse(length(indGroup) == 1, indGroup, 
                sample(indGroup, 1))
        }
        else {
            indGroup <- (1:nrGroup)[(thProp - empProp) == max(thProp - 
                empProp)]
            selGroup <- ifelse(length(indGroup) == 1, indGroup, 
                sample(indGroup, 1))
        }
        OUT <- unique(c(out, (1:length(cbGroup))[cbGroup != cbControl$names[selGroup]]))
    }
    if (!is.null(nAvailable)) {
        discard <- unique(c(OUT, which(nAvailable == 0)))
        OUT <- discard
    }
    if (crit == "MFI") {
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        info <- Ii(theta, itemBank, model = model, D = D)$Ii
        ranks <- rank(info)
        nrIt <- min(c(randomesque, sum(items)))
        keepRank <- sort(ranks[items == 1], decreasing = TRUE)[1:nrIt]
        keep <- NULL
        for (i in 1:length(keepRank)) keep <- c(keep, which(ranks == keepRank[i] & items == 1))
        select <- ifelse(length(keep) == 1, keep, sample(c(keep), 1))
        res <- list(item = select, par = itemBank[select, ], 
            info = info[select], criterion = criterion, randomesque = randomesque)
    }
    if (crit == "bOpt") {
        if (!is.null(model)) 
            stop("bOpt's rule cannot be considered with polytomous items", 
                call. = FALSE)
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        distance <- abs(itemBank[, 2] - theta)
        ranks <- rank(distance)
        ranks[OUT] <- -1
        nrIt <- min(c(randomesque, sum(items)))
        keepRank <- sort(ranks[items == 1], decreasing = FALSE)[1:nrIt]
        keepRank <- unique(keepRank)
        keep <- NULL
        for (i in 1:length(keepRank)) keep <- c(keep, which(ranks == keepRank[i] & items == 1))
        select <- ifelse(length(keep) == 1, keep, sample(keep, 
            1))
        res <- list(item = select, par = itemBank[select, ], 
            info = distance[select], criterion = criterion, randomesque = randomesque)
    }
    if (crit == "MLWI" | crit == "MPWI") {
        if (length(out) == 1) 
            par <- rbind(itemBank[out, ])
        else par <- itemBank[out, ]
        ITEMS <- rep(1, nrow(itemBank))
        ITEMS[OUT] <- 0
        likInfo <- rep(0, nrow(itemBank))
        for (i in 1:nrow(itemBank)) {
            if (ITEMS[i] == 1) 
                likInfo[i] <- MWI(itemBank, i, x, it.given = par, 
                  model = model, type = criterion, lower = parInt[1], 
                  upper = parInt[2], nqp = parInt[3], priorDist = priorDist, 
                  priorPar = priorPar, D = D)
        }
        likVal <- sort(likInfo, decreasing = TRUE)[min(c(randomesque,sum(ITEMS)))]
        keep <- (1:length(ITEMS))[likInfo >= likVal]
        select <- ifelse(length(keep) == 1, keep, sample(keep,1))
        res <- list(item = select, par = itemBank[select, ], 
            info = likInfo[select], criterion = criterion, randomesque = randomesque)
    }
    if (crit == "KL" | crit == "KLP") {
        if (length(out) == 1) 
            par <- rbind(itemBank[out, ])
        else par <- itemBank[out, ]
        ITEMS <- rep(1, nrow(itemBank))
        ITEMS[OUT] <- 0
        klvalue <- rep(0, nrow(itemBank))
        L <- function(th, r, param) prod(Pi(th, param, D = D)$Pi^r * 
            (1 - Pi(th, param, D = D)$Pi)^(1 - r))
        X <- seq(from = parInt[1], to = parInt[2], length = parInt[3])
        LL <- function(th, r, param, model, D = D) {
          if (dim(param)[1]==0)
            res<-1
          else {
              prob <- Pi(th, param, model = model, D = D)$Pi
              res <- 1
              for (i in 1:length(r)) res <- res * prob[i, r[i] + 1]
              }
          return(res)
        }
        if (is.null(model)) 
            LF <- sapply(X, L, x, par)
        else
            LF <- sapply(X, LL, x, par, model = model, D = D)
        for (i in 1:nrow(itemBank)) {
            if (ITEMS[i] == 1) 
                klvalue[i] <- KL(itemBank, i, x, it.given = par, 
                  model = model, theta = theta, type = criterion, 
                  lower = parInt[1], upper = parInt[2], nqp = parInt[3], 
                  priorDist = priorDist, priorPar = priorPar, 
                  lik = LF, X = X, D = D)
        }
        klVal <- sort(klvalue, decreasing = TRUE)[min(c(randomesque, 
            sum(ITEMS)))]
        keep <- (1:length(ITEMS))[klvalue >= klVal]
        select <- ifelse(length(keep) == 1, keep, sample(keep, 
            1))
        res <- list(item = select, par = itemBank[select, ], 
            info = klvalue[select], criterion = criterion, randomesque = randomesque)
    }
    
    
    
    if (crit == "GDI" | crit == "GDIP") {
      if (length(out) == 1) 
        par <- rbind(itemBank[out, ])
      else par <- itemBank[out, ]
      ITEMS <- rep(1, nrow(itemBank))
      ITEMS[OUT] <- 0
      gdivalue <- rep(0, nrow(itemBank))
      L <- function(th, r, param) prod(Pi(th, param, D = D)$Pi^r * (1 - Pi(th, param, D = D)$Pi)^(1 - r))
      X <- seq(from = parInt[1], to = parInt[2], length = parInt[3])
      LLL <- function(th, r, param, model, D = 1) {
        if (dim(param)[1]==0)
          res<-1
        else {
          prob <- Pi(th, param, model = model, D = D)$Pi
          res <- 1
          for (i in 1:length(r)) res <- res * prob[i, r[i] + 1]
        }
        return(res)
      }
      if (is.null(model)) 
        LF <- sapply(X, L, x, par)
      else
        LF <- sapply(X, LLL, x, par, model = model, D = D)
      for (i in 1:nrow(itemBank)) {
        if (ITEMS[i] == 1) 
          gdivalue[i] <- GDI(itemBank, i, x, it.given = par, 
                           model = model, type = criterion, 
                           lower = parInt[1], upper = parInt[2], nqp = parInt[3], 
                           priorDist = priorDist, priorPar = priorPar, 
                           lik = LF, X = X, D = D)
      }
      gdiVal <- sort(gdivalue, decreasing = TRUE)[min(c(randomesque,sum(ITEMS)))]
      keep <- (1:length(ITEMS))[gdivalue >= gdiVal]
      select <- ifelse(length(keep) == 1, keep, sample(keep,1))
      res <- list(item = select, par = itemBank[select, ], 
                  info = gdivalue[select], criterion = criterion, randomesque = randomesque)
    }

    if (crit == "MEI") {
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        infos <- rep(0, length(items))
        for (i in 1:length(items)) {
            if (items[i] > 0) 
                infos[i] <- MEI(itemBank, item = i, x = x, theta = theta, 
                  it.given = itemBank[out, ], model = model, 
                  method = method, priorDist = priorDist, priorPar = priorPar, 
                  D = D, range = range, parInt = parInt, infoType = infoType)
        }
        infoVal <- sort(infos, decreasing = TRUE)[min(c(randomesque, 
            sum(items)))]
        keep <- (1:nrow(itemBank))[infos >= infoVal]
        select <- ifelse(length(keep) == 1, keep, sample(keep, 
            1))
        res <- list(item = select, par = itemBank[select, ], 
            info = infos[select], criterion = criterion, randomesque = randomesque)
    }
    if (crit == "MEPV") {
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        epvs <- rep(1000, length(items))
        for (i in 1:length(items)) {
            if (items[i] > 0) 
                epvs[i] <- EPV(itemBank, item = i, x = x, theta = theta, 
                  it.given = itemBank[out, ], model = model, 
                  priorDist = priorDist, priorPar = priorPar, 
                  D = D, parInt = parInt)
        }
        epVal <- sort(epvs)[min(c(randomesque, sum(items)))]
        keep <- (1:nrow(itemBank))[epvs <= epVal]
        select <- ifelse(length(keep) == 1, keep, sample(keep, 
            1))
        res <- list(item = select, par = itemBank[select, ], 
            info = epvs[select], criterion = criterion, randomesque = randomesque)
    }
    if (crit == "random") {
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        gen <- as.integer(runif(1, 0, 1) * (sum(items))) + 1
        ind <- (1:nrow(itemBank))[items > 0][gen]
        res <- list(item = ind, par = itemBank[ind, ], info = NA, 
            criterion = criterion, randomesque = randomesque)
    }
    if (crit == "progressive") {
        items_administered <- length(out)
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        info <- Ii(theta, itemBank, model = model, D = D)$Ii
        itemMaxInfo <- max(info[items == 1])
        randomValues <- runif(length(items), 0, itemMaxInfo)
        wq <- 0
        if (rule == "precision") {
            infostop <- (1/thr)^2
            cuminfo <- (1/SETH)^2
            if (items_administered > 0) 
                wq <- max(cuminfo/infostop, items_administered/(maxItems - 
                  1))^AP
        }
        if (rule == "length") {
            if (items_administered > 0) {
                numerador <- sum((1:items_administered)^AP)
                denominador <- sum((1:(thr - 1))^AP)
                wq <- numerador/denominador
            }
        }
        funcPR <- info * wq + randomValues * (1 - wq)
        funcPR[OUT] <- 0
        keep <- which(funcPR == max(funcPR))
        select <- ifelse(length(keep) == 1, keep, sample(keep, 
            1))
        res <- list(item = select, par = itemBank[select, ], 
            info = info[select], criterion = criterion, randomesque = randomesque)
    }
    if (crit == "proportional") {
        items_administered <- length(out)
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        wq <- 0
        if (rule == "precision") {
            infostop <- (1/thr)^2
            cuminfo <- (1/SETH)^2
            if (items_administered > 0) 
                wq <- infostop * max(cuminfo/infostop, items_administered/(maxItems - 
                  1))^AP
        }
        if (rule == "length") 
            if (items_administered > 0) {
                numerador <- sum((1:items_administered)^AP)
                denominador <- sum((1:(thr - 1))^AP)
                wq <- thr * numerador/denominador
            }
        info <- Ii(theta, itemBank, model = model, D = D)$Ii
        infoPR <- info^wq
        infoPR[OUT] <- 0
        totalInfoPR <- sum(infoPR[items == 1])
        probSelect <- infoPR/totalInfoPR
        select <- sample(1:length(items), size = 1, prob = probSelect)
        res <- list(item = select, par = itemBank[select, ], 
            info = info[select], criterion = criterion, randomesque = randomesque)
    }
    if (crit == "thOpt") {
        if (!is.null(model)) 
            stop("'thOpt' rule cannot be considered with polytomous items", 
                call. = FALSE)
        items <- rep(1, nrow(itemBank))
        items[OUT] <- 0
        u <- -3/4 + (itemBank[, 3] + itemBank[, 4] + -2 * itemBank[, 
            3] * itemBank[, 4])/2
        v <- (itemBank[, 3] + itemBank[, 4] - 1)/4
        xstar <- 2 * sqrt(-u/3) * cos(acos(-v * sqrt(-27/u^3)/2)/3 + 
            4 * pi/3) + 1/2
        thstar <- itemBank[, 2] + log((xstar - itemBank[, 3])/(itemBank[, 
            4] - xstar))/(D * itemBank[, 1])
        distance <- abs(thstar - theta)
        ranks <- rank(distance)
        ranks[OUT] <- -1
        nrIt <- min(c(randomesque, sum(items)))
        keepRank <- sort(ranks[items == 1], decreasing = FALSE)[1:nrIt]
        keepRank <- unique(keepRank)
        keep <- NULL
        for (i in 1:length(keepRank)) {
            keep <- c(keep, which(ranks == keepRank[i]))
        }
        select <- ifelse(length(keep) == 1, keep, sample(keep, 
            1))
        res <- list(item = select, par = itemBank[select, ], 
            info = distance[select], criterion = criterion, randomesque = randomesque)
    }
    if (is.null(cbControl)) 
        res[[6]] <- res[[7]] <- res[[8]] <- NA
    else {
        res[[6]] <- empProp
        postProp <- NULL
        for (i in 1:nrGroup) postProp[i] <- length(c(res$item, 
            out)[cbGroup[c(res$item, out)] == cbControl$names[i]])
        res[[7]] <- postProp/sum(postProp)
        res[[8]] <- thProp
    }
    names(res)[6:8] <- c("prior.prop", "post.prop", "cb.prop")
    return(res)
}
