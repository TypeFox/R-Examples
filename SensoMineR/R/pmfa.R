pmfa<-function (matrice, matrice.illu = NULL, mean.conf = NULL, dilat = TRUE,
    graph.ind = TRUE, graph.mfa = TRUE, lim = c(60, 40), coord = c(1, 2), cex = 0.8)
{
    procrustes <- function(amat, target, orthogonal = F, translate = F,
        magnify = FALSE) {
        for (i in nrow(amat):1) {
            if (any(is.na(amat)[i, ]) | any(is.na(target)[i,
                ])) {
                amat <- amat[-i, ]
                target <- target[-i, ]
            }
        }
        dA <- dim(amat)
        dX <- dim(target)
        if (length(dA) != 2 || length(dX) != 2)
            stop("arguments amat and target must be matrices")
        if (any(dA != dX))
            stop("dimensions of amat and target must match")
        if (length(attr(amat, "tmat")))
            stop("oblique loadings matrix not allowed for amat")
        if (orthogonal) {
            if (translate) {
                p <- dX[1]
                target.m <- (rep(1/p, p) %*% target)[, ]
                amat.m <- (rep(1/p, p) %*% amat)[, ]
                target.c <- scale(target, center = target.m,
                  scale = FALSE)
                amat.c <- scale(amat, center = amat.m, scale = FALSE)
                j <- svd(crossprod(target.c, amat.c))
            }
            else {
                amat.c <- amat
                j <- svd(crossprod(target, amat))
            }
            rot <- j$v %*% t(j$u)
            if (magnify)
                beta <- sum(j$d)/sum(amat.c^2)
            else beta <- 1
            B <- beta * amat.c %*% rot
            if (translate)
                B <- B + rep(as.vector(target.m), rep.int(p,
                  dX[2]))
            value <- list(rmat = B, tmat = rot, magnify = beta)
            if (translate)
                value$translate <- target.m - (rot %*% amat.m)[,
                  ]
        }
        else {
            b <- solve(amat, target)
            gamma <- sqrt(diag(solve(crossprod(b))))
            rot <- b * rep(gamma, rep.int(dim(b)[1], length(gamma)))
            B <- amat %*% rot
            fcor <- solve(crossprod(rot))
            value <- list(rmat = B, tmat = rot, correlation = fcor)
        }
        value
    }
 #------------------------------------------------------------------------------
 
    nbjuge <- ncol(matrice)/2
    matricemoyenne<-colMeans(matrice)
    matrice <- scale(matrice, center = TRUE, scale = FALSE)
    if (!is.null(matrice.illu)) matrice.illu = matrice.illu[rownames(matrice),] 
    do.mfa = FALSE
    if (is.null(mean.conf)) {
        do.mfa = TRUE
        if (is.null(matrice.illu)) res.afm <- MFA(as.data.frame(matrice), group = rep(2, nbjuge), ncp = max(coord),type=rep("c",nbjuge),graph=FALSE)
        else res.afm <- MFA(cbind.data.frame(matrice,matrice.illu), group = c(rep(2, nbjuge),ncol(matrice.illu)), ncp = max(coord),type=c(rep("c",nbjuge),"s"),num.group.sup = nbjuge+1,graph=FALSE)
        mean.conf <- res.afm$ind$coord
    }
    mean.conf <- as.matrix(mean.conf[, coord])
    res <- matrix(0, nbjuge, 1)
 #-----------------------------------------------------------------------------#
    for (j in 1:nbjuge) {
        atourner <- as.matrix(matrice[, (2 * (j - 1) + 1):(2 *
            j)])
        if ((dilat == TRUE) & (do.mfa == TRUE)) {
            eig <- eigen(1/nrow(atourner) * t(scale(atourner,scale=FALSE)) %*% scale(atourner,scale=FALSE))
            res.procrustes <- procrustes(atourner, mean.conf,
                orthogonal = TRUE, translate = TRUE, magnify = FALSE)
            magnify <- sqrt(res.afm$eig[1,1])/sqrt(eig$values[1])
        }
        else {
            res.procrustes <- procrustes(atourner, mean.conf,
                orthogonal = TRUE, translate = TRUE, magnify = TRUE)
            magnify <- res.procrustes$magnify
        }
        tourne <- scale(atourner,scale=FALSE) %*% res.procrustes$tmat * magnify
        res[j] <- coeffRV(mean.conf, tourne)$rv
        if (graph.ind == TRUE) {
            dd = cbind(mean.conf, tourne)
            nappe <- rbind((matrix(c(0, 0, 0, lim[2], lim[1], lim[2], lim[1], 0),ncol=2,byrow = TRUE)
                        -cbind(rep( matricemoyenne[(2 * (j - 1) + 1)],4),rep( matricemoyenne[(2 * j)],4))) %*% res.procrustes$tmat * magnify)
                
               
            if (j != 1)
            dev.new()
            plot(rbind(tourne, mean.conf, nappe), type = "n",
                xlab = paste("Dim", coord[1]), ylab = paste("Dim",
                  coord[2]), asp = 1, main = colnames(matrice)[2 *
                  j], sub = paste("RV between the mean representation and the representation of",
                  colnames(matrice)[2 * j], ": ", signif(res[j], 4)), cex.sub = cex)
            for (i in 1:nrow(mean.conf)) points(mean.conf[i, 1], mean.conf[i, 2], cex = cex, pch = 20)
            for (i in 1:nrow(mean.conf)) text(mean.conf[i, 1], mean.conf[i, 2], rownames(matrice)[i], cex = cex,
                pos = 1, offset = 0.5)
            lines(nappe[c(1:4, 1), ], col = 3, lty = 2)
            for (i in 1:nrow(mean.conf)) points(tourne[i, 1],  tourne[i, 2], cex = cex, pch = 20, col = 3)
            for (i in 1:nrow(mean.conf)) text(tourne[i, 1], tourne[i, 2], rownames(matrice)[i], col = 3, font = 3,
                cex = cex, pos = 2, offset = 0.2)
        }
    }
    if (do.mfa&graph.mfa){
      plot(res.afm,choix="var",habillage="var",axes=coord,new.plot=TRUE)
      if (!is.null(matrice.illu)){
        plot(res.afm,choix="var",invisible="sup",habillage="var",axes=coord,new.plot=TRUE)
        plot(res.afm,choix="var",invisible="actif",axes=coord,new.plot=TRUE)
      }
      plot(res.afm,choix="ind",partial="all",habillage="group",axes=coord,new.plot=TRUE)
      plot(res.afm,choix="ind",axes=coord,new.plot=TRUE)
      plot(res.afm,choix="group",axes=coord,new.plot=TRUE)
    }
    dimnames(res) <- list(colnames(matrice)[(1:(ncol(matrice)/2)) * 2], "RV coeff")
    return(res)
}
