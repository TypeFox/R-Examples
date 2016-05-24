`pertables` <-
function (data, index = NULL, nsim = 100) 
{
    index0 <- c("Indet", "indet", "", " ", as.character(c(1:100)), 
        "sp", paste("sp", as.character(c(1:100)), sep = ""), 
        paste("sp", as.character(c(1:100)), sep = " "))
    index <- c(index0, index)
    colnames(data)[1:3] <- c("Family", "Genus", "Specific")
    data <- cbind(Species = paste(data$Genus, data$Specific), 
        data)
    cond <- rep(0, length(data$Species))
    cond <- as.numeric(data$Specific %in% index)
    cond <- cond + as.numeric(data$Genus %in% index)
    cond <- cond + as.numeric(data$Family %in% index)
    results <- list()
    for (p in 1:nsim) {
        print(paste("This is simulation number", p))
        data.A <- data
        for (i in 1:length(data.A[cond == 3, 1])) {
            pr <- as.numeric(c(rep(0, 4), data.A[cond == 3, -c(1:4)][i, 
                ]))
            pr <- ifelse(is.na(pr) == "TRUE", 0, pr)
            temp1 <- apply(as.data.frame(data.A[, pr > 0]), 1, 
                sum)
            data.A$Species[cond == 3][i] <- sample(rep(c(as.character(data.A$Species[cond == 
                3][i]), as.character(data.A$Species[temp1 == 
                0])), 2), 1)
        }
        for (i in 1:length(data.A[cond == 2, 1])) {
            pr2 <- as.numeric(c(rep(0, 4), data.A[cond == 2, 
                -c(1:4)][i, ]))
            pr2 <- ifelse(is.na(pr2) == "TRUE", 0, pr2)
            data.f <- data.A[data.A$Family == data.A$Family[cond == 
                2][i], ]
            temp1 <- apply(as.data.frame(data.f[, pr2 > 0]), 
                1, sum)
            data.A$Species[cond == 2][i] <- ifelse(dim(data.f)[1] == 
                0, NA, sample(rep(c(as.character(data.A$Species[cond == 
                2][i]), as.character(data.f$Species[temp1 == 
                0])), 2), 1))
        }
        for (i in 1:length(data.A[cond == 1, 1])) {
            pr3 <- as.numeric(c(rep(0, 4), data.A[cond == 1, 
                -c(1:4)][i, ]))
            pr3 <- ifelse(is.na(pr3) == "TRUE", 0, pr3)
            data.g <- data.A[data.A$Genus == data.A$Genus[cond == 
                1][i], ]
            temp1 <- apply(as.data.frame(data.g[, pr3 > 0]), 
                1, sum)
            data.A$Species[cond == 1][i] <- ifelse(dim(data.g)[1] == 
                0, NA, sample(rep(c(as.character(data.A$Species[cond == 
                1][i]), as.character(data.g$Species[temp1 == 
                0])), 2), 1))
        }
        fun1 <- function(arg1) {
            apply(arg1, 2, sum)
        }
        data.A2 <- do.call(rbind, by(data.A[, -c(1:4)], data.A$Species, 
            fun1))
        results[[p]] <- t(data.A2)
    }
    n <- table(cond)
    names(n) <- ifelse(names(n) == "0", "Fully identified", names(n))
    names(n) <- ifelse(names(n) == "1", "Identified to genus", 
        names(n))
    names(n) <- ifelse(names(n) == "2", "Identified to family", 
        names(n))
    names(n) <- ifelse(names(n) == "3", "Fully undetermined", 
        names(n))
    attributes(dimnames(n))$names <- "Taxonomic uncertainty"
    raw <- t(data[, -c(1:4)])
    raw <- raw[, cond == 0]
    final <- list(taxunc = n, pertables = results, raw = raw)
    class(final) <- c("pertables", class(final))
    return(final)
}

