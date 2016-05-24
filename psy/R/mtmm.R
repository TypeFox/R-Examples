mtmm <-
function (datafile, x, color = FALSE, itemTot = FALSE, graphItem = FALSE, stripChart=FALSE, namesDim=NULL) 
{
 
 ###############################################################################
 #          Graph of Correlation within and between Item of subscale           #
 ###############################################################################
   par.old <- par("mfrow")
    k <- 1
    name <- nam <- c()
    r2 <- c()
    nn <- c()
    nbdim <- length(x)
    X <- list()
    n <- as.data.frame(matrix(nrow = 1, ncol = nbdim))
    datafile <- na.omit(datafile[,unlist(x)])
    for (i in 1:nbdim) {
        X[[i]] <- datafile[, x[[i]]]
    }
# Names of variables and matrix of correlations
    for (i in 1:nbdim) {
        n[i] <- length(X[[i]])
    }
    nn <- c(nn, n)
    k <- 1
    for (i in 1:nbdim) {
        nameD <- as.data.frame(matrix(nrow = n[[i]], ncol = nbdim))
        nameD[i] <- names(X[[i]])
        for (j in 1:nbdim) {
            r <- as.data.frame(matrix(nrow = n[[i]], ncol = n[[j]]))
            r <- cor(X[i][[1]], X[j][[1]])
            r2[k] <- list(r)
            k <- k + 1
        }
        name <- c(name, nameD[i])
        nam <- c(nam, name[[i]])
    }
	# Separation of both types of matrix
    a <- seq(from = 1, to = length(r2), by = (nbdim + 1))
    X1 <- r2[a]
    X2 <- r2[-a]

    p1 <- length(X1)
    p2 <- length(X2)

# Correlation of item with its own dimension
 # Initialisation
    correlation <- c()
    dimension <- c()
    mat <- c()
    V1 <- c()
    V2 <- c()

# Number of item

    k1 <- c()
    for (q in 1:nbdim) {
        w <- seq(from = 1, to = nn[[q]])
        k = 0
        for (i in 1:(nn[[q]] - 1)) {
            k <- w[i] + k
        }
        k1 <- c(k1, k)
    }
    for (l in 1:p1) {
        mat1 <- matrix(ncol = 1, nrow = k1[l])
        mat2 <- matrix(ncol = 1, nrow = k1[l])
        mat3 <- matrix(ncol = 1, nrow = k1[l])
        Var1 <- matrix(ncol = 1, nrow = k1[l])
        Var2 <- matrix(ncol = 1, nrow = k1[l])
        k = 1
        for (i in 1:(nn[[l]] - 1)) {
            for (j in (i + 1):nn[[l]]) {
                mat1[k, 1] <- X1[[l]][j, i]
                mat2[k, 1] <- l
                mat3[k, 1] <- l
                Var1[k, 1] <- attributes(X1[[l]])$dimnames[[1]][i]
                Var2[k, 1] <- attributes(X1[[l]])$dimnames[[1]][j]
                k = k + 1
            }
        }
        correlation <- append(correlation, mat1)
        dimension <- append(dimension, c(mat2))
        mat <- append(mat, c(mat3))
        V1 <- append(V1, c(Var1))
        V2 <- append(V2, c(Var2))
    }

# Format
    er <- cbind(dimension, mat, V1, V2)
    er <- as.data.frame(er)
    er$correlation <- correlation

# Correlation of item with other dimensions

    er1 <- c()
    V1 <- c()
    V2 <- c()
    for (l in 1:p2) {
        k <- 1
        mat1 <- matrix(ncol = 1, nrow = ((dim(X2[[l]])[2]) * 
            (dim(X2[[l]])[1])))
        Var1 <- matrix(ncol = 1, nrow = ((dim(X2[[l]])[2]) * 
            (dim(X2[[l]])[1])))
        Var2 <- matrix(ncol = 1, nrow = ((dim(X2[[l]])[2]) * 
            (dim(X2[[l]])[1])))
        for (i in 1:(dim(X2[[l]])[1])) {
            for (j in 1:(dim(X2[[l]])[2])) {
                mat1[k, 1] <- X2[[l]][i, j]
                Var1[k, 1] <- rep(attributes(X2[[l]])$dimnames[[1]][i], 
                  each = (dim(X2[[l]])[2]))[j]
                Var2[k, 1] <- attributes(X2[[l]])$dimnames[[2]][j]
                k <- k + 1
            }
        }
        er1 <- append(er1, c(mat1))
        V1 <- append(V1, c(Var1))
        V2 <- append(V2, c(Var2))
    }
  # Format
    err <- cbind(V1, V2)
    err <- as.data.frame(err)
    err$correlation <- er1

# Names of matrix
    i2 <- c()
    i1 <- c()
    for (k in 1:(nrow(err))) {
        for (l in 1:p1) {
            if (err$V2[k] %in% name[[l]]) 
                i2[k] <- l
            if (err$V1[k] %in% name[[l]]) 
                i1[k] <- l
        }
    }
    err$mat <- i2
    err$dimension <- i1

# grouping of matrix

    corre <- rbind(er, err)

#****************************************************************************
#********************************* Colors ***********************************
#****************************************************************************
    col <- Y <- c()
    for (i in 1:nbdim) {
        if (color == FALSE) {
            col <- rep("white", times = nbdim)
            col[i] <- c("grey")
        }
        else col <- 2:nbdim
        Y[i] <- list(col)
    }
#****************************************************************************
#********************************  plots ***********************************
#****************************************************************************
    if (graphItem == FALSE) {
        if (stripChart==FALSE) {  
            for (i in 1:p1) {
                ifelse(is.null(namesDim), maintitle <- paste("Scale", i, sep = " "),
                maintitle <- paste(i, namesDim[i], sep = " "))
                Dim <- subset(corre, (corre$dimension == i), drop = TRUE)
                boxplot(Dim$correlation ~ Dim$mat, col = Y[[i]], 
                    main = maintitle, xlab = "i", 
                    ylab = paste("Correlation of Items of Scale i with Items of Scale",
                    i), ylim = c(min(corre$correlation), max(corre$correlation)))
            }
        }
        if (stripChart==TRUE) {  
            for (i in 1:p1) {
                ifelse(is.null(namesDim), maintitle <- paste("Scale", i, sep = " "),
                maintitle <- paste(i, namesDim[i], sep = " "))
                Dim <- subset(corre, (corre$dimension == i), drop = TRUE)
                stripchart(Dim$correlation ~ Dim$mat, vertical=TRUE, method="jitter", 
                    jitter=0.05, pch=1, cex=1.5, main = maintitle, xlab = "i", 
                    ylab = paste("Correlation of Items of Scale i with Items of Scale",
                    i), ylim = c(min(corre$correlation), max(corre$correlation)))
                for (j in 1:p1) points(j, median(Dim$correlation[Dim$mat==j]), pch="-", cex=4)
             }
        }
    }
###############################################################################
#          Value of Correlation of Item of Scale with Scale                  #
##############################################################################
    n2 <- nrow(X[[1]])
    Score <- matrix(nrow = n2, ncol = nbdim)
    for (i in 1:nbdim) {
        k <- 1
        for (j in 1:n2) {
            Score[k, i] <- sum(X[[i]][j, ])
            k <- k + 1
        }
    }

# score and correlation
    Scoreaj <- c()
    mat <- c()
    mat2 <- c()
    for (i in 1:nbdim) {
        Scorea <- matrix(nrow = n2, ncol = nn[[i]])
        for (j in 1:(nn[[i]])) {
            for (l in 1:n2) {
                Scorea[l, j] <- Score[l, i] - X[[i]][l, j]
            }
        }
        Scoreaj[i] <- list(Scorea)
        mat <- rep(Score[, i], each = nn[[i]])
        dim(mat) <- c(nn[[i]], n2)
        mat2[i] <- list(t(mat))
    }

# correlations
    Scale <- ScaleI <- Item <- co <- c <- c()
    t <- 0
    for (i in 1:nbdim) {
        t <- t + as.numeric(nn[i])
        for (j in 1:nbdim) {
            for (k in 1:nn[[i]]) {
                if (i == j) {
                  Score <- Scoreaj[[i]][, k]
                  c <- c(c, cor(Score, X[[i]][, k]))
                  Item <- c(Item, names(X[[i]][k]))
                  ScaleI <- c(ScaleI, i)
                  Scale <- c(Scale, j)
                }
                else {
                  c <- c(c, cor(X[[i]], mat2[[j]][, 1])[k])
                  Item <- c(Item, names(X[[i]][k]))
                  ScaleI <- c(ScaleI, i)
                  Scale <- c(Scale, j)
                }
            }
        }
    }
    co <- c(co, c)
    err <- as.data.frame(cbind(Item, Scale))
    err$ScaleI <- ScaleI
    err$correlation <- co
    tabl <- reshape(err, direction = "wide", timevar = "Scale", 
        idvar = "Item", v.names = c("correlation"))
    for (i in 1:nbdim) {
        colnames(tabl)[i + 2] <- ifelse(is.null(namesDim), paste("Scale", i, sep = " "),
          paste(i, namesDim[i], sep = " "))
    }
    if (itemTot == TRUE & graphItem == FALSE) {
        x11()
        par(mfrow = par.old)
        if (stripChart==FALSE) {  
            for (i in 1:nbdim) {
              ifelse(is.null(namesDim), maintitle <- paste("Scale", i, sep = " "),
              maintitle <- paste(i, namesDim[i], sep = " "))
              Dim <- subset(err, (err$Scale == i), drop = TRUE)
                boxplot(Dim$correlation ~ Dim$ScaleI, col = Y[[i]], 
                    main = maintitle, ylim = c(min(err$correlation), 
                      max(err$correlation)), cex.axis = 0.75, xlab = "i", ylab = 
                      paste("Corrected Correlation of Item i with Total Score",i))
            }
        }
        if (stripChart==TRUE) {  
            for (i in 1:nbdim) {
                ifelse(is.null(namesDim), maintitle <- paste("Scale", i, sep = " "),
                maintitle <- paste(i, namesDim[i], sep = " "))
                Dim <- subset(err, (err$Scale == i), drop = TRUE)
                stripchart(Dim$correlation ~ Dim$ScaleI, vertical=TRUE, method="jitter", 
                    jitter=0.05, pch=1, cex=1.5, main = maintitle, xlab = "i", 
                    ylab = paste("Corrected Correlation of Item i with Total Score",
                    i), ylim = c(min(err$correlation), max(err$correlation)))
                for (j in 1:nbdim) points(j, median(Dim$correlation[Dim$ScaleI==j]), pch="-", cex=4)
             }
        }
    }
    if (graphItem == TRUE) {
        for (i in 1:nbdim) {
            Dim <- subset(err, (err$Scale == i), drop = TRUE)
            Dim$Item <- factor(Dim$Item, levels = unique(as.character(Dim$Item)))
            with(Dim, stripchart(correlation ~ Item, vertical = TRUE, 
                add = FALSE, main = paste("Score of Scale", i, 
                  sep = " "), xlab = "Item", ylab = "Correlation", 
                ylim = c(min(corre$correlation), max(corre$correlation)), cex.axis = 1, las = 3, 
                pch = 20, cex = 0.75))
            abline(h = median(Dim$correlation), col = "red")
        }
    }
    return(tabl)
}
