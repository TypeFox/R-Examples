# ----------------------------------------------------
# Function to produce the "strata" dataframe
# starting from the available sampling frame
# Author: Giulio Barcaroli
# Date: 23 September 2015
# ----------------------------------------------------
buildStrataDF <- function(dataset) {
    # stdev1 is for sampling data
    stdev1 <- function(x, w) {
        mx <- sum(x * w)/sum(w)
        sqrt(sum(w * (x - mx)^2)/(sum(w) - 1))
    }
    # stdev2 is for population data
    stdev2 <- function(x, w) {
        mx <- sum(x * w)/sum(w)
        sqrt(sum(w * (x - mx)^2)/(sum(w)))
    }
    colnames(dataset) <- toupper(colnames(dataset))
    nvarX <- length(grep("X", names(dataset)))
    nvarY <- length(grep("Y", names(dataset)))
    if (length(grep("WEIGHT", names(dataset))) == 1) {
        cat("\nComputations have been done on sampling data\n")
        stdev <- "stdev1"
    }
    if (length(grep("WEIGHT", names(dataset))) == 0) {
        dataset$WEIGHT <- rep(1, nrow(dataset))
        stdev <- "stdev2"
        cat("\nComputations have been done on population data\n")
    }
    numdom <- length(levels(as.factor(dataset$DOMAINVALUE)))
    stratatot <- NULL
    # begin domains cycle
    for (d in (1:numdom)) {
		dom <- levels(as.factor(dataset$DOMAINVALUE))[d]
		domain <- dataset[dataset$DOMAINVALUE == dom, ]
        listX <- NULL
        namesX <- NULL
        for (i in 1:nvarX) {
            name <- paste("X", i, sep = "")
            namesX <- cbind(namesX, name)
            if (i < nvarX) 
                listX <- paste(listX, "domain$X", i, ",", sep = "") else listX <- paste(listX, "domain$X", i, sep = "")
        }
        listM <- NULL
        listS <- NULL
        for (i in 1:nvarY) {
            listM <- paste(listM, "M", i, ",", sep = "")
            listS <- paste(listS, "S", i, ",", sep = "")
        }
        stmt <- paste("domain$STRATO <- as.factor(paste(", listX, 
            ",sep='*'))", sep = "")
        eval(parse(text = stmt))
        for (i in 1:nvarY) {
            WEIGHT <- NULL
            STRATO <- NULL
            Y <- NULL
            stmt <- paste("Y <- domain$Y", i, "[!is.na(domain$Y", 
                i, ")]", sep = "")
            eval(parse(text = stmt))
            stmt <- paste("WEIGHT <- domain$WEIGHT[!is.na(domain$Y", 
                i, ")]", sep = "")
            eval(parse(text = stmt))
            stmt <- paste("STRATO <- domain$STRATO[!is.na(domain$Y", 
                i, ")]", sep = "")
            eval(parse(text = stmt))
            STRATO <- factor(STRATO)
            stmt <- paste("M", i, " <- tapply(WEIGHT * Y,STRATO,sum) / tapply(WEIGHT,STRATO,sum)", 
                sep = "")
            eval(parse(text = stmt))
            samp <- NULL
            stmt <- paste("samp <- domain[!is.na(domain$Y", i, 
                "),]", sep = "")
            eval(parse(text = stmt))
            l.split <- split(samp, samp$STRATO, drop = TRUE)
            stmt <- paste("S", i, " <- sapply(l.split, function(df,x,w) ", 
                stdev, "(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')", 
                sep = "")
            eval(parse(text = stmt))
            stmt <- paste("stratirid <- unlist(attr(M", i, ",'dimnames'))", 
                sep = "")
            eval(parse(text = stmt))
            strati <- data.frame(X1 = levels(domain$STRATO))
            stmt <- paste("m <- data.frame(cbind(X1=stratirid,X2=M", 
                i, "))", sep = "")
            eval(parse(text = stmt))
            m <- merge(strati, m, by = c("X1"), all = TRUE)
            m$X2 <- as.character(m$X2)
            m$X2 <- as.numeric(m$X2)
            m$X2 <- ifelse(is.na(m$X2), 0, m$X2)
            stmt <- paste("M", i, " <- m$X2", sep = "")
            eval(parse(text = stmt))
            stmt <- paste("s <- data.frame(cbind(X1=stratirid,X2=S", 
                i, "))", sep = "")
            eval(parse(text = stmt))
            s <- merge(strati, s, by = c("X1"), all = TRUE)
            s$X2 <- as.character(s$X2)
            s$X2 <- as.numeric(s$X2)
            s$X2 <- ifelse(is.na(s$X2), 0, s$X2)
            stmt <- paste("S", i, " <- s$X2", sep = "")
            eval(parse(text = stmt))
        }
        N <- tapply(domain$WEIGHT, domain$STRATO, sum)
        STRATO <- domain$STRATO
        COST <- rep(1, length(levels(domain$STRATO)))
        CENS <- rep(0, length(levels(domain$STRATO)))
        DOM1 <- rep(as.character(dom), length(levels(domain$STRATO)))
        stmt <- paste("strata <- as.data.frame(cbind(STRATO=levels(STRATO),N,", 
            listM, listS, "COST,CENS,DOM1))")
        eval(parse(text = stmt))
        for (i in 1:nvarX) {
            stmt <- paste("strata$X", i, " <- rep(0, length(levels(domain$STRATO)))", 
                sep = "")
            eval(parse(text = stmt))
        }
        strata$STRATO <- as.character(strata$STRATO)
        for (i in 1:nrow(strata)) {
            strata[i, c(namesX)] <- unlist(strsplit(strata$STRATO[i], 
                "\\*"))
        }
        stratatot <- rbind(stratatot, strata)
    }  # end domain cycle
    colnames(stratatot) <- toupper(colnames(stratatot))
    stratatot$DOM1 <- as.factor(stratatot$DOM1)
    write.table(stratatot, "strata.txt", quote = FALSE, sep = "\t", 
        dec = ".", row.names = FALSE)
    stratatot <- read.delim("strata.txt")
    return(stratatot)
}
