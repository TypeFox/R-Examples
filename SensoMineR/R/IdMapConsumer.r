IdMapConsumer = function (dataset.id, dataset.signa, col.p, col.j, col.lik, num.col.var.signa, conf.level=0.95, id.recogn, nbchoix = NULL,
    nbsimul = 500, alpha = 0.05, coord = c(1, 2), precision = 0.1,
    levels.contour = NULL, color = FALSE, simusigni = 500) 
{
#rownames
aux<-rownames(dataset.signa)
dataset.signa<-cbind.data.frame(aux, dataset.signa)#ajout de la colonne contenant les juges
for(f in 1:ncol(dataset.signa)) dataset.signa[, f] <- as.factor(as.character(dataset.signa[, f]))
col.j.signa = 1
num.col.var.signa<-num.col.var.signa+1# +1 car on a ajouté la colonne des juges
var.signa.test<-num.col.var.signa #on fera le test global (inertie) sur toutes les variables
dataset <- as.data.frame(merge(dataset.id, dataset.signa, by.x = colnames(dataset.id)[col.j], by.y = colnames(dataset.signa)[col.j.signa]))
num.col.var <- ncol(dataset.id) + num.col.var.signa-1
	
	for(f in num.col.var) dataset[, f] <- as.factor(as.character(dataset[, f]))
   forprefmap <- function(dataset, col.p, col.j, col.lik, id.recogn, 
        rm.j = TRUE) {
        dataset[, col.p] <- as.factor(as.character(dataset[, col.p]))
        product <- levels(dataset[, col.p])
        nbprod <- length(product)
        dataset[, col.j] <- as.factor(as.character(dataset[, col.j]))
        juge <- levels(dataset[, col.j])
        nbjuge <- length(juge)
        desc <- dataset[, c(col.j, col.p)]
        id.pos <- grep(id.recogn, colnames(dataset))
        intensity <- dataset[, id.pos - 1]
        int.data <- cbind(desc, intensity)
        liking <- as.matrix(dataset[, col.lik])
        colnames(liking) <- "liking"
        lik.data <- cbind(desc, liking)
        int.avg <- averagetable(int.data, formul = paste("~", 
            colnames(dataset)[col.p], "+", colnames(dataset)[col.j], 
            sep = ""), firstvar = 3, method = "coeff")
        lik.j <- matrix(0, nbprod, 0)
        rownames(lik.j) <- product

        for (j in 1:nbjuge) {
            lik.j.temp <- lik.data[lik.data[, 1] == juge[j], 
                ]
            lik.j.temp.rn <- lik.j.temp[, 2]
            lik.j.temp <- as.matrix(lik.j.temp[, 3])
            rownames(lik.j.temp) <- lik.j.temp.rn
            lik.j <- merge(lik.j, lik.j.temp, all = T, by = 0, 
                sort = F)
            lik.j.rn <- lik.j[, 1]
            lik.j <- as.matrix(lik.j[, -1])
            rownames(lik.j) <- lik.j.rn
            colnames(lik.j)[ncol(lik.j)] <- juge[j]
        }
        if (rm.j) {
            rm.c <- NULL
            for (j in 1:ncol(lik.j)) if (var(na.omit(lik.j[, 
                j])) == 0) 
                rm.c <- c(rm.c, j)
            if (!is.null(rm.c)) {
                lik.j <- lik.j[, -rm.c]
                juge.rm <- juge[rm.c]
            }
            else {
                juge.rm <- NULL
            }
        }
        res <- list()
        res$senso <- int.avg #fonction qui retourne l'average table
        res$hedo <- as.data.frame(lik.j)#et les notes hédoniques
        if (rm.j) 
            juge.rm
        return(res)
    }
    procrustes <- function(amat, target, orthogonal = FALSE, 
        translate = FALSE, magnify = FALSE) {
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
        return(value)
    }
    oo <- order(dataset[, col.p])
    dataset <- dataset[oo, ]
    oo <- order(dataset[, col.j])
    dataset <- dataset[oo, ]
    dataset[, col.p] <- as.factor(dataset[, col.p])
    product <- levels(dataset[, col.p])
    nbprod <- length(product)
    dataset[, col.j] <- as.factor(dataset[, col.j])
    juge <- levels(dataset[, col.j])
    nbjuge <- length(juge)
    info <- dataset[, c(col.j, col.p)]
    id.pos <- grep(id.recogn, colnames(dataset))
    intensity <- dataset[, id.pos - 1]
    attribut <- colnames(intensity)
    nbatt <- length(attribut)
    ideal <- dataset[, id.pos]
    int.data <- cbind(info, intensity)
    id.data <- cbind(info, ideal)
    data.cut <- forprefmap(dataset, col.p = col.p, col.j = col.j, 
        id.recogn = id.recogn, col.lik = col.lik)#découpage du jdd en average table et notes hédoniques (liste)
    int.p.avg <- scale(data.cut$senso, scale = F)
    int.j.avg <- averagetable(int.data, formul = paste("~", colnames(dataset)[col.j], 
        "+", colnames(dataset)[col.p]), firstvar = 3, method = "coeff")
    id.j.avg <- averagetable(id.data, formul = paste("~", colnames(dataset)[col.j], 
        "+", colnames(dataset)[col.p]), firstvar = 3, method = "coeff")
    ideal.juge <- vector("list", nbjuge)
    names(ideal.juge) <- juge
    for (j in 1:nbjuge) {
        ideal.j <- id.data[id.data[, 1] == juge[j], ]
        rownames(ideal.j) <- ideal.j[, 2]
        ideal.j <- ideal.j[, -c(1, 2)]
        temp <- as.matrix(int.j.avg[j, ])
        ideal.juge[[j]] <- sweep(ideal.j, 2, as.vector(as.matrix(int.j.avg[j, 
            ])), FUN = "-")
    }
    data.j.cplt <- matrix(0, 0, nbatt)
    colnames(data.j.cplt) <- attribut

    l = 0
	
    for (j in 1:nbjuge) {
        data.j.cplt <- rbind(data.j.cplt, ideal.juge[[j]])
        rownames(data.j.cplt)[c((l + 1):nrow(data.j.cplt))] <- paste(rownames(ideal.juge[[j]]), 
            "_", juge[j], sep = "")
        l = nrow(data.j.cplt)
    }
    colnames(data.j.cplt) <- colnames(int.p.avg)
    data.pca <- rbind.data.frame(int.p.avg, data.j.cplt)
    res.pca <- PCA(data.pca, ind.sup = c((nbprod + 1):nrow(data.pca)), 
        graph = F, ncp = Inf)#acp pour construire l'espace produit
    id.j.avg.cor <- id.j.avg - int.j.avg
    colnames(id.j.avg.cor) <- colnames(int.p.avg)
    data.pcab <- rbind.data.frame(int.p.avg, id.j.avg.cor)
    res.pcab <- PCA(data.pcab, ind.sup = c((nbprod + 1):nrow(data.pcab)), 
        graph = F, ncp = Inf)#acp pour projeter en supplémentaire les profils idéaux moyens de chaque juge
    plot.PCA(res.pcab, choix = "ind", cex = 0.8, label = "ind.sup", 
        new.plot = T, title = "Projection of the individual averaged ideal profiles")
    eig <- res.pca$eig
    data.pca2 <- merge(data.cut$senso, data.cut$hedo, all = T, 
        by = 0, sort = F)
    rownames(data.pca2) <- data.pca2[, 1]
    data.pca2 <- data.pca2[, -1]
    res.pca2 <- PCA(data.pca2, quanti.sup = (ncol(data.cut$senso) + 
        1):ncol(data.pca2), graph = F)#acp pour projeter en supplémentaire les notes hédoniques de chaque juge
    plot.PCA(res.pca2, choix = "var", invisible = "var", label = "quanti.sup", 
        cex = 0.9, new.plot = T, title = "Projection of the individual hedonic scores")
    ideal.j.dim <- vector("list", nbjuge)
    names(ideal.j.dim) <- juge
    l = 0
    for (j in 1:nbjuge) {
        ideal.j.dim[[j]] <- res.pca$ind.sup$coord[c((l + 1):(l + 
            nrow(ideal.juge[[j]]))), ]
        rownames(ideal.j.dim[[j]]) <- paste(product, "_", juge[j], 
            sep = "")
        l <- l + nrow(ideal.juge[[j]])
    }
    if (is.null(nbchoix)) 
        nbchoix = nbprod
    simul <- matrix(0, nbchoix, 0)
    rownames(simul) <- paste("prod", 1:nbchoix, sep = "")
    for (sim in 1:nbsimul) simul <- cbind(simul, as.matrix(sample(nbprod, 
        nbchoix, replace = T)))
    colnames(simul) <- paste("Simul.", 1:nbsimul, sep = "")
    ponder = res.pcab$call$col.w
    estim.ncp <- max(max(coord), estim_ncp(sweep(int.p.avg, 2, 
        sqrt(ponder), FUN = "*"), scale = FALSE, ncp.min = 0, 
        ncp.max = min(10, ncol(int.p.avg)))$ncp)
    target.pca <- res.pcab$ind.sup$coord[, 1:estim.ncp]
    jdd <- target.pca
    for (sim in 1:nbsimul) {
        juge.sample <- matrix(0, nbjuge, nbatt)
        rownames(juge.sample) <- juge
        colnames(juge.sample) <- colnames(int.p.avg)
        for (j in 1:nbjuge) juge.sample[j, ] <- apply(ideal.juge[[j]][as.vector(simul[, 
            sim]), ], 2, mean)
        data.pca.temp <- rbind.data.frame(int.p.avg, juge.sample)
        res.pca.temp <- PCA(data.pca.temp, ind.sup = c((nbprod + 
            1):nrow(data.pca.temp)), graph = F, ncp = estim.ncp)$ind.sup$coord
        aux <- procrustes(res.pca.temp, target.pca, orthogonal = TRUE, 
            translate = TRUE, magnify = FALSE)$rmat
        colnames(aux) <- colnames(jdd)
        jdd = rbind.data.frame(jdd, aux)
    }
    truc = cbind.data.frame(jdd[-(1:nrow(target.pca)), ], rep(rownames(target.pca), 
        nbsimul))
    res.simul = list()
    res.simul$moy$simul = truc[order(truc[, ncol(truc)]), ]
    res.simul$moy$J = cbind.data.frame(target.pca, rownames(target.pca))
    res.simul$moy$J = res.simul$moy$J[order(res.simul$moy$J[, 
        ncol(res.simul$moy$J)]), ]
    plotellipse2 <- function(mat, alpha = 0.05, coord = c(1, 
        2), eig, cex = 1, color = NULL) {
        res <- plotellipseinter2(mat, alpha = alpha, coord = coord, 
            nbgroup = 1, moy = T, eig = eig, color = color, cex = cex)
        if (length(mat$partiel) != 0) {
            dev.new()
            nbgroup = length(levels(mat$partiel$simul[, ncol(mat$partiel$simul)]))/length(levels(mat$moy$simul[, 
                ncol(mat$moy$simul)]))
            plotellipseinter2(mat, alpha = alpha, coord = coord, 
                nbgroup = nbgroup, moy = F, eig = eig, color = color, 
                cex = cex)
        }
        return(res)
    }
    plotellipseinter2 <- function(mat, alpha = 0.05, coord = c(1, 
        2), nbgroup = 1, moy = TRUE, eig, cex = 1, color = NULL) {
        if (moy == T) {
            matJ = cbind.data.frame(mat$moy$J[, 1:2], mat$moy$J[, 
                ncol(mat$moy$J)])
            matJP = cbind.data.frame(mat$moy$JP[, 1:2], mat$moy$JP[, 
                ncol(mat$moy$JP)])
            matsimul = cbind.data.frame(mat$moy$simul[, 1:2], 
                mat$moy$simul[, ncol(mat$moy$simul)])
        }
        if (moy == F) {
            matmoyJ = cbind.data.frame(mat$moy$J[, 1:2], mat$moy$J[, 
                ncol(mat$moy$J)])
            matmoyJP = cbind.data.frame(mat$moy$JP[, 1:2], mat$moy$JP[, 
                ncol(mat$moy$JP)])
            matmoysimul = cbind.data.frame(mat$moy$simul[, 1:2], 
                mat$moy$simul[, ncol(mat$moy$simul)])
            matJ = cbind.data.frame(mat$partiel$J[, 1:2], mat$partiel$J[, 
                ncol(mat$partiel$J)])
            matJP = cbind.data.frame(mat$partiel$JP[, 1:2], mat$partiel$JP[, 
                ncol(mat$partiel$JP)])
            matsimul = cbind.data.frame(mat$partiel$simul[, 1:2], 
                mat$partiel$simul[, ncol(mat$partiel$simul)])
        }
        nbp <- nrow(matJ)
        nbjuge <- nbp/nbgroup
        coord.ellipse.a.tracer <- matrix(0, 402, 2 * nbp)
        p <- 2
        nbprod <- nrow(matJP)/nrow(matJ)
        nbsimul <- nrow(matsimul)/nrow(matJ)
        for (i in 1:nbp) {
            VX <- var(matsimul[((i - 1) * nbsimul + 1):(i * nbsimul), 
                1:2])
            coord.ellipse.a.tracer[, (1 + 2 * (i - 1)):(2 * i)] <- ellipse2(as.numeric(t(matJ[i, 
                1:2])), VX, alpha)
        }
        minx <- min(coord.ellipse.a.tracer[, 1 + 2 * (0:(nbp - 
            1))], na.rm = T)
        maxx <- max(coord.ellipse.a.tracer[, 1 + 2 * (0:(nbp - 
            1))], na.rm = T)
        miny <- min(coord.ellipse.a.tracer[, 2 * (1:nbp)], na.rm = T)
        maxy <- max(coord.ellipse.a.tracer[, 2 * (1:nbp)], na.rm = T)
        plot(0, 0, xlab = paste("Dim ", coord[1], " (", round(eig[coord[1], 
            2], 2), "%)", sep = ""), ylab = paste("Dim ", coord[2], 
            " (", round(eig[coord[2], 2], 2), "%)", sep = ""), 
            xlim = c(minx * 1.05, maxx * 1.05), ylim = c(1.05 * 
                miny, 1.05 * maxy), col = "white", asp = 1)
        if (moy == T) 
            title(main = "Individual ideal confidence ellipses")
        abline(v = 0, lty = 2)
        abline(h = 0, lty = 2)
        if (moy == F) {
            points(matmoyJ[, 1], matmoyJ[, 2], cex = 0.8 * cex, 
                col = "blue3", pch = 15)
            text(matmoyJ[, 1], matmoyJ[, 2], matmoyJ[, ncol(matmoyJ)], 
                cex = 0.8 * cex, pos = 4, offset = 0.2, col = "blue3")
        }
        if (moy == T) 
            text(matJ[, 1], matJ[, 2], matJ[, ncol(matJ)], cex = 0.8 * 
                cex, pos = 4, offset = 0.2, col = "blue3")
        for (j in 1:nbgroup) {
            for (i in 1:nbjuge) {
                points(matJ[(j - 1) * nbjuge + i, 1], matJ[(j - 
                  1) * nbjuge + i, 2], cex = 0.8 * cex, col = "blue3", 
                  pch = 20)
                if (moy == F) 
                  lines(c(matJ[(j - 1) * nbjuge + i, 1], matmoyJ[i, 
                    1]), c(matJ[(j - 1) * nbjuge + i, 2], matmoyJ[i, 
                    2]), col = "lightblue", lty = j)
                lines(coord.ellipse.a.tracer[, (1 + 2 * ((i + 
                  (j - 1) * nbjuge) - 1)):(2 * (i + (j - 1) * 
                  nbjuge))], col = "lightblue", lty = j)
            }
        }
        return(coord.ellipse.a.tracer)
    }
    "ellipse2" <- function(loc, cov, alpha) {
        A <- cov
        detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
        dist <- sqrt(qchisq(1 - alpha/2, 2))
        ylimit <- sqrt(A[2, 2]) * dist
        y <- seq(-ylimit, ylimit, 0.01 * ylimit)
        sqrt.discr <- sqrt(detA/A[2, 2]^2 * abs(A[2, 2] * dist^2 - 
            y^2))
        sqrt.discr[c(1, length(sqrt.discr))] <- 0
        b <- loc[1] + A[1, 2]/A[2, 2] * y
        x1 <- b - sqrt.discr
        x2 <- b + sqrt.discr
        y <- loc[2] + y
        return(rbind(cbind(x1, y), cbind(rev(x2), rev(y))))
    }
    dev.new()
    res.ellipse <- round(plotellipse2(res.simul, alpha = 0.05, 
        coord = coord, eig = res.pca$eig, cex = 1, color = NULL), 
        1)
    minx <- min(min(res.ellipse[, 1 + 2 * (0:(nbjuge - 1))], 
        na.rm = T), floor(min(res.pca$ind$coord[, coord[1]])))
    maxx <- max(max(res.ellipse[, 1 + 2 * (0:(nbjuge - 1))], 
        na.rm = T), ceiling(max(res.pca$ind$coord[, coord[1]])))
    miny <- min(min(res.ellipse[, 2 * (1:nbjuge)], na.rm = T), 
        floor(min(res.pca$ind$coord[, coord[2]])))
    maxy <- max(max(res.ellipse[, 2 * (1:nbjuge)], na.rm = T), 
        ceiling(max(res.pca$ind$coord[, coord[2]])))
    juge.mat <- vector("list", nbjuge)#liste contenant autant d'éléments que de juges
    names(juge.mat) <- juge
    lim.minx <- floor(minx)#bornes pour le quadrillage de l'idmap
    lim.maxx <- ceiling(maxx)
    lim.miny <- floor(miny)
    lim.maxy <- ceiling(maxy)
    nbrow <- length(seq(lim.minx, lim.maxx, precision))#nombre de cases en abcsisse du quadrillage
    nbcol <- length(seq(lim.miny, lim.maxy, precision))#nombre de cases en ordonnée du quadrillage
    juge.tot <- matrix(0, nbcol, nbrow)#matrice aux dimensions du quadrillage de la carte idéale, en effet "~tracé en chiffres"
    rownames(juge.tot) <- round(seq(lim.miny, lim.maxy, precision), 
        1)
    colnames(juge.tot) <- round(seq(lim.minx, lim.maxx, precision), 
        1)
    cons.wgt <- matrix(0, 1, nbjuge)#poids de  chaque juge
    rownames(cons.wgt) <- "weight"
    colnames(cons.wgt) <- juge
    for (j in 1:nbjuge) {
        juge.mat[[j]] <- matrix(0, nbcol, nbrow)
        rownames(juge.mat[[j]]) <- round(seq(lim.miny, lim.maxy, 
            precision), 1)
        colnames(juge.mat[[j]]) <- round(seq(lim.minx, lim.maxx, 
            precision), 1)
        ellipse.x <- res.ellipse[, (1 + 2 * (j - 1))]
        ellipse.y <- res.ellipse[, (2 * j)]
        for (i in 1:length(ellipse.x)) {
            posx <- grep(ellipse.x[i], colnames(juge.mat[[j]]))
            if (length(posx) > 1) {
                posx.temp = NULL
                for (l in 1:length(posx)) if (colnames(juge.mat[[j]])[posx[l]] == 
                  ellipse.x[i]) 
                  posx.temp = posx[l]
                if (!is.null(posx.temp)) {
                  posx = posx.temp
                }
                else {
                  stop("Not convenient posx definition")
                }
            }
            posy <- grep(ellipse.y[i], rownames(juge.mat[[j]]))
            if (length(posy) > 1) {
                posy.temp = NULL
                for (l in 1:length(posy)) if (rownames(juge.mat[[j]])[posy[l]] == 
                  ellipse.y[i]) 
                  posy.temp = posy[l]
                if (!is.null(posy.temp)) {
                  posy = posy.temp
                }
                else {
                  stop("Not convenient posy definition")
                }
            }
            juge.mat[[j]][posy, posx] = 1
        }
        for (i in 1:nrow(juge.mat[[j]])) {
            pos1 <- grep(1, juge.mat[[j]][i, ])
            if (length(pos1) >= 2) 
                juge.mat[[j]][i, c(pos1[1]:pos1[length(pos1)])] = 1
        }
       
       
            cons.wgt[1, j] <- 1#poids: tous les individus sont pris en compte dans le tracé de la carte globale
        
        juge.mat[[j]] <- juge.mat[[j]] * cons.wgt[1, j]#pondération
		juge.tot = juge.tot+ juge.mat[[j]]#matrice pour tracer le fond de carte (idmap globale sur tous les individus)
		
		juge.tot.list<-vector("list",length(num.col.var)) #liste stockant les différents juge.tot pour chaque modalite de toutes les variables dans le cas où les cartes"partielles" sont tracée avec la méthode ellipse 
    }
    juge.tot.rn <- paste("Y_", rownames(juge.tot), sep = "")
    juge.tot.cn <- paste("X_", colnames(juge.tot), sep = "")
    rownames(juge.tot) <- juge.tot.rn
    colnames(juge.tot) <- juge.tot.cn
    juge.tot <- 100 * round(juge.tot/sum(cons.wgt[1, ]), 3)
#remplissage de juge.tot.list : on ne somme que les indiv appartenenant à la modalité	
			# for(v in 1:length(num.col.var))#boucle sur les variables
		# {
			# for(m in 1:nlevels(dataset[,num.col.var[v]]))#boucle sur les modalites
			# {
			# juge.tot.list[[v]][[m]] <- matrix(0, nbcol, nbrow)
			# rownames(juge.tot.list[[v]][[m]]) <- round(seq(lim.miny, lim.maxy, precision), 1)
			# colnames(juge.tot.list[[v]][[m]]) <- round(seq(lim.minx, lim.maxx, precision), 1)
			# somme = 0
			# for (j in 1:nbjuge){#boucle sur les juges
			# if(dataset[nbprod*j,num.col.var[v]] == levels(dataset[,num.col.var[v]])[m]){
			# juge.tot.list[[v]][[m]] = juge.tot.list[[v]][[m]] + juge.mat[[j]]#le juge "participe" à la construction de juge.tot de la modalite ssi il la possede
			# somme=somme+1#compteur de nb de juges ayant la modalite
			# }			
			# }
	
			# juge.tot.list[[v]][[m]] <- 100 * round(juge.tot.list[[v]][[m]]/somme, 3)#matrice d'altitide pour l'ideal des indiv ayant la modalite m de la variable v
			# }
		# }
	
    f1 <- seq(lim.minx, lim.maxx, precision)#x : quadrillage des abscisses de l'idmap
    f2 <- seq(lim.miny, lim.maxy, precision)#y : quadrillage des ordonnées de l'idmap
    if (!is.null(levels.contour)) {
        if (min(levels.contour) < 0 || max(levels.contour) > 
            100 || length(levels.contour) < 2) {
            warning("Not convenient 'levels.contour' definition: the default value will be used")
            levels.contour = NULL
        }
        else {
            oo <- order(levels.contour)
            levels.contour <- levels.contour[oo]
        }
    }
    if (is.null(levels.contour)) 
        levels.contour <- seq(10, 5 * floor(max(juge.tot)/5), 
            5)
    # dev.new()
    
       
    
   
    if (color) {#si couleur zone d'idéal par modalités tracées avec la méthode ellipse sur fond de carte de l'idmap global
         titrecolor <- "Ideal Mapping"
		image(f1, f2, t(juge.tot), asp = 1,col = terrain.colors(200), 
            xlab = paste("Dim ", coord[1], "(", round(res.pca$eig[coord[1], 
                2], 2), "%)", sep = ""), ylab = paste("Dim ", 
                coord[2], "(", round(res.pca$eig[coord[2], 2], 
                  2), "%)", sep = ""), main = titrecolor)
        contour(f1, f2, t(juge.tot), nlevels = length(levels.contour), 
            levels = levels.contour, add = T, labex = 0, asp=1)
        for (i in 1:nrow(res.pca$ind$coord)) {
            points(res.pca$ind$coord[i, coord[1]], res.pca$ind$coord[i, 
                coord[2]], pch = 15)
            text(res.pca$ind$coord[i, coord[1]], res.pca$ind$coord[i, 
                coord[2]], rownames(res.pca$ind$coord)[i], pos = 4, 
                offset = 0.2, cex = 0.7)
		}#fin de tracé idmap globale en couleur
		# coul<-0
		# modalite<-list()
		# for(v in 1:length(num.col.var))
		# {
			# for(m in 1:nlevels(dataset[,num.col.var[v]]))
			# {
			# coul<-coul+1
			# contour(f1, f2, t(juge.tot.list[[v]][[m]]),levels=max( t(juge.tot.list[[v]][[m]])),add=TRUE,col=coul, lwd = 3) #tracé des contours de la zone d'ideal max de chaque modalite
			# effectif<- table(dataset.signa[,colnames(dataset)[num.col.var[v]]] == levels(dataset[,num.col.var[v]])[m])[2]
			# modalite[[coul]]<-paste(levels(dataset[,num.col.var[v]])[m], "(", effectif, ")")
			# }
		
		# }
		
		# legend("topright", title = paste("Ideal areas by ", colnames(dataset)[num.col.var[v]]), legend = modalite, fill = 1:length(modalite), text.col = "black")
        abline(v = 0, lty = 2)
        abline(h = 0, lty = 2)
    }
    else {
		
###########################################################################trace ac ellipses des zones ideales des modalites		
	# if (ellipse)
	# {
		# coul<-1
		# modalite<-list()
		# for(v in 1:length(num.col.var))
		# {
			# for(m in 1:nlevels(dataset[,num.col.var[v]]))
			# {
			# coul<-coul+1
			# contour(f1, f2, t(juge.tot.list[[v]][[m]]),levels=max( t(juge.tot.list[[v]][[m]])),add=TRUE,col=coul, lwd = 3)
			# effectif<- table(dataset.signa[,colnames(dataset)[num.col.var[v]]] == levels(dataset[,num.col.var[v]])[m])[2]
			# modalite[[coul-1]]<-paste(levels(dataset[,num.col.var[v]])[m], "(", effectif, ")")
			# }
		
		# }
		
		# legend("topright", title = paste("Ideal areas by ", colnames(dataset)[num.col.var[v]]), legend = modalite, fill = 2:(length(modalite)+1), text.col = "black")
        # abline(v = 0, lty = 2, col = "white")
        # abline(h = 0, lty = 2, col = "white")
    # }
	# if (densite)############################################################### trace ac densite pr les modalites
	# {debut du if densite
##		require(KernSmooth)#package contenant la fonction bkde2D pour estimer la densité
		coul<-1
		modalite<-list()
		listpointsvrai<-vector("list",length(num.col.var.signa))#liste contenant les coordonnées sur le 1er plan factoriel des idéaux de toutes les variables/toutes les modalités
		listpointssimu<-vector("list",length(num.col.var.signa))#idem pour les idéaux simulés (affectation aléatoire des individus aux modalités)
		cont=1
	    coord.ideaux = res.pca$ind.sup$coord #en ligne les indiv, en col les coord de leur ideal par pdt sur l'espace produit
		nom.juge<-rep(0,nrow(coord.ideaux))
		coord.ideaux<-cbind.data.frame(nom.juge,coord.ideaux)#ajout de la colonne du nom du juge
		for (i in 1:nrow(coord.ideaux)) coord.ideaux[i,1]<-strsplit(rownames(coord.ideaux),"_")[[i]][2]#remplissage de la colonne par le nom du juge
		for(v in num.col.var.signa)#boucle sur les variables
		{
			for(m in 1:nlevels(dataset.signa[,v]))#boucle sur les modalités
			{
			coul<-coul+1
			lignes<-which(dataset.signa[,v]== levels(dataset.signa[,v])[m])#numero des lignes prenant la modalite m pour la variable v
			jugemod<-dataset.signa[lignes,col.j.signa]#dataframe avec les juges ayant la modalité m
			larecup<-list()
			for (i in 1:length(jugemod)) larecup[[i]]<-which(coord.ideaux[,"nom.juge"]== jugemod[i]) #recupération des lignes correspondant aux juges ayant la modalite m
			larecup<-unlist(larecup)
			coord.ideaux.moda<-as.data.frame(coord.ideaux[larecup,c(2,3)])#coordonnées des idéaux SUR DIM1 ET DIM2 des juges ayant la modalité m
			effectif<- table(dataset.signa[,v] == levels(dataset.signa[,v])[m])[2]#effectif de la modalite m
			est<-bkde2D(coord.ideaux.moda, gridsize = c(100, 100), bandwidth=c(0.5, 0.5))#estimation de la densite, idée envisagée : largeur de bande = f(effectif)
			# contour(est$x1, est$x2, est$fhat, add=TRUE, levels=max(est$fhat),col=coul, lwd = 3)
			# modalite[[coul-1]]<-paste(levels(dataset.signa[,v])[m], "(", effectif, ")")
			#recuperation des coord du max
			est$fhat = as.data.frame(est$fhat)
			rownames(est$fhat) = est$x2#ordonnée
			colnames(est$fhat) = est$x1#abscisse
			mat.pos.maxI = which(est$fhat == max(est$fhat), arr.ind = TRUE)#abscisse et ordonnée de la position du maximum de densité donc de l'idéal
			abs.maxI = as.numeric(colnames(est$fhat)[mat.pos.maxI[,2]])#abs des points de densité max
			ord.maxI = as.numeric(rownames(est$fhat)[mat.pos.maxI[,1]])#ord des points de densité max
			# points(abs.maxI, ord.maxI, pch = 17, col = coul)
			listpointsvrai[[cont]][[m]] = cbind(abs.maxI, ord.maxI)#coordonnées de l'ideal de la var v, modalite m
			names(listpointsvrai)[cont] = colnames(dataset.signa)[v]
			names(listpointsvrai[[cont]])[m] = levels(dataset.signa[,v])[m]
			}
		cont=cont+1
		# legend("topright", title = paste("Ideal areas by ", colnames(dataset)[v]), legend = modalite, fill = 2:(length(modalite)+1), text.col = "black")
        # abline(v = 0, lty = 2, col = "white")
        # abline(h = 0, lty = 2, col = "white")
		}
		

###simulations sans remise, en respectant l'effectif des classes
listpointssimu<-vector("list",length(num.col.var.signa))#idem que liste point vrai mais contient un niveau supplémentaire pour les simulation
contb=1
for(v in num.col.var.signa){#boucle sur les variables
	for(s in 1:simusigni){#boucle sur les simulations
		dataset.signa.simu<-dataset.signa
		pointssimu = list()
		for(m in 1:nlevels(dataset.signa.simu[,v])){#boucle sur les modalités
			aenlever<-list()#liste qui contiendra les num des lignes à enlever car les juges ont déjà été affecté à une modalité de v
			nblignes<-length(which(dataset.signa[,v]== levels(dataset.signa[,v])[m]))#nb de personnes ayant la modalite m pour la variable v
			jugemodsim<-sample(dataset.signa.simu[,col.j.signa],nblignes)#tirage aléatoire de nblignes juges:lignes tirées, les juges sont affectés à la modalité m
			larecupsim<-list()
			for (i in 1:length(jugemodsim)){ 
			 larecupsim[[i]]<-which(coord.ideaux[,"nom.juge"]== jugemodsim[i])
			 aenlever[[i]]<-which(dataset.signa.simu[,col.j.signa] == jugemodsim[i])
			}
			larecupsim<-unlist(larecupsim)
			aoter<-unlist(aenlever)
			coord.ideaux.moda.sim<-as.data.frame(coord.ideaux[larecupsim,c(2,3)])
			# effectif<- table(dataset.signa[,v] == levels(dataset.signa[,v])[m])[2]
			est<-bkde2D(coord.ideaux.moda.sim, gridsize = c(100, 100), bandwidth=c(0.5, 0.5)) ##########################largeur de bande = f(effectif)
			est$fhat = as.data.frame(est$fhat)
			rownames(est$fhat) = est$x2
			colnames(est$fhat) = est$x1
			mat.pos.max = which(est$fhat == max(est$fhat), arr.ind = TRUE)#matrice avec 2 colonnes : numero de lignes et num de colonnes de toutes le cases contenant le max
			abs.max = as.numeric(colnames(est$fhat)[mat.pos.max[,2]])#abs des points de densité max
			ord.max = as.numeric(rownames(est$fhat)[mat.pos.max[,1]])#ord des points de densité max
			abs.moy = mean(abs.max)
			ord.moy = mean(ord.max)
			pointssimu[[m]] = cbind(abs.moy, ord.moy) #liste stockant les resultats de la simu : cad coordonnées des modalités
			listpointssimu[[contb]][[s]] = pointssimu #liste dans laquelle on accumule les résultats de toutes les simulations. Liste à 3 niveaux : la variable, la simulation, et toutes les moda correspondant a la simulation
			# names(pointssimu)[s] = paste(s,"_", colnames(dataset.signa.simu)[v],"_",levels(dataset.signa.simu[,v])[m])
			dataset.signa.simu = dataset.signa.simu[-aoter,] #on enleve les juges affectés à la modalité m du dataframe pour refaire la prochaine affectation(pour la meme simulation, meme variable)
			}
		}
	contb = contb+1
	}

###########TEST GLOBAL (calcul d'inertie entre les modalités d'une même variable)ET TEST PAR PAIRES (calcul de distances entre les modalités 2 à 2 d'une variable)
cv=1
resultat.test = list()
inertie.obs<-as.data.frame(rep(0,length(var.signa.test)))
inertie.critique<-as.data.frame(rep(0,length(var.signa.test)))
inertie.pvalue<-as.data.frame(rep(0,length(var.signa.test)))
for(v in var.signa.test)
	{
		paires<-combn(1:nlevels(dataset.signa[,v]),2)#en colonnes toutes les paires possibles
		resultat.paires = as.data.frame(matrix(0,ncol(paires),4))
		colnames(resultat.paires) = c("Categories","Dc", "Dobs","p-value")
		#test global
		X=list()
		Y = list()
		poids = list()
		for (m in 1:nlevels(dataset.signa.simu[,v]))
		{
			X[[m]] = listpointsvrai[[cv]][[m]][1]#abscisse de l'idéal de la modalite m
			Y[[m]] = listpointsvrai[[cv]][[m]][2]#ordonnée de l'idéal de la modalite m
			poids[[m]] = (table(dataset.signa[,v] == levels(dataset.signa[,v])[m])[2])/nbjuge#pourcentage d'individus ayant la modalité m
		}
		X=do.call(rbind,X)#vecteur contenant les abscisses de l'idéal de toutes les modalités
		Y=do.call(rbind,Y)#idem ordonnées
		poids=do.call(rbind,poids)#idem poids
		rownames(inertie.obs)[cv]<-colnames(dataset.signa)[v]
		colnames(inertie.obs) = "inertie.obs"
		mX = sum(poids*X)#pondération
		mY = sum(poids*Y)
		inertie.obs[cv,1]<-sum(poids*(X-mX)^2)+ sum(poids*(Y-mY)^2)#calcul de l'inertie observée
		##calcul de inertie intra variables pour chaque simulation
		vectinertiesimu<-rep(0,simusigni)
		for(s in 1:simusigni)
			{
			X.simu=list()
			Y.simu=list()
				for (m in 1:nlevels(dataset.signa.simu[,v]))
				{
				X.simu[[m]] = listpointssimu[[cv]][[s]][[m]][1]
				Y.simu[[m]] = listpointssimu[[cv]][[s]][[m]][2]
				}
				X.simu = do.call(rbind, X.simu)
				Y.simu = do.call(rbind, Y.simu)
				mX.simu = sum(poids*X.simu)
				mY.simu = sum(poids*Y.simu)
				vectinertiesimu[s]<-sum(poids*(X.simu-mX.simu)^2)+ sum(poids*(Y.simu-mY.simu)^2)#vecteur contenant les simmusigni(500) inerties simulées
			}
		inertie.critique[cv,1]<-quantile(vectinertiesimu, probs = conf.level)
		rownames(inertie.critique)[cv]<-colnames(dataset.signa)[v]
		colnames(inertie.critique) = "inertie.critique"
		inertie.pvalue[cv,1]<-sum(vectinertiesimu>inertie.obs[cv,1])/simusigni
		rownames(inertie.pvalue)[cv]<-colnames(dataset.signa)[v]
		colnames(inertie.pvalue) = "p-value"
		#tracé de l'idmap si test global sur inertie intra var significatif
		if(inertie.obs[cv,1]>inertie.critique[cv,1])
		{
		dev.new()
    
        titre <- "Ideal Mapping"
			#####################trace du fond de carte
	image(f1, f2, t(juge.tot), asp=1, col = grey(1:max(juge.tot)/100), 
            xlab = paste("Dim ", coord[1], "(", round(res.pca$eig[coord[1], 
                2], 2), "%)", sep = ""), ylab = paste("Dim ", 
                coord[2], "(", round(res.pca$eig[coord[2], 2], 
                  2), "%)", sep = ""), main = titre)
        contour(f1, f2, t(juge.tot), nlevels = length(levels.contour), 
            levels = levels.contour, add = T, labex = 0, col = "white", asp=1)
        for (i in 1:nrow(res.pca$ind$coord)) {
            points(res.pca$ind$coord[i, coord[1]], res.pca$ind$coord[i, 
                coord[2]], pch = 15, col = "white")
            text(res.pca$ind$coord[i, coord[1]], res.pca$ind$coord[i, 
                coord[2]], rownames(res.pca$ind$coord)[i], pos = 4, 
                offset = 0.2, cex = 0.7, col = "white")
				
		}		
	#####################trace des modalites
	for (m in 1:nlevels(dataset.signa[,v])) 
	{points(listpointsvrai[[cv]][[m]][1], listpointsvrai[[cv]][[m]][2], pch = 17, col = m+1)
	effectif<- table(dataset.signa[,v] == levels(dataset.signa[,v])[m])[2]
	modalite[[m]]<-paste(levels(dataset.signa[,v])[m], "(", effectif, ")")
	}
	legend("topright", title = paste("Ideal areas by ", colnames(dataset.signa)[v]), legend = modalite, fill = 2:(length(modalite)+1), text.col = "black")
        abline(v = 0, lty = 2, col = "white")
        abline(h = 0, lty = 2, col = "white")
		}#fin du if trace

		
	##test par paires de modalites intra variable(distance)
	for(colo in 1:ncol(paires))#pour la paire de modalité correspondant à la 1ere colonne de toutes les combinaisons possibles (rangées dans paires)
		{
		resultat.paires[colo,1] = paste(colnames(dataset.signa)[v],"_", levels(dataset.signa[,v])[paires[1,colo]], ":",colnames(dataset.signa)[v],"_", levels(dataset.signa[,v])[paires[2,colo]])
		vectdsimu = rep(0,simusigni)
		for(s in 1:simusigni)
			{
			vectdsimu[s] = dist(rbind(listpointssimu[[cv]][[s]][[paires[1,colo]]], listpointssimu[[cv]][[s]][[paires[2,colo]]]))
			}
		Dc<-quantile(vectdsimu, probs = conf.level)#H1: la dist est superieure au cas où elle est due au hasard
		Dobs<-dist(rbind(listpointsvrai[[cv]][[paires[1,colo]]], listpointsvrai[[cv]][[paires[2,colo]]]))
		resultat.paires[colo,2] = Dc
		resultat.paires[colo,3] = Dobs
		resultat.paires[colo,4] = sum(vectdsimu>Dobs)/simusigni
		# if(Dobs>Dc) {resultat.paires[colo,4] = 1-conf.level} else {resultat.paires[colo,4] = "NS"}
		}
	resultat.test[[cv]] = resultat.paires
	names(resultat.test)[cv] = colnames(dataset.signa)[v]
	cv=cv+1
	}
	}	

    maxval <- max(juge.tot)
    res.id <- matrix(0, 0, 2)
    colnames(res.id) <- c("X", "Y")
    for (i in 1:nrow(juge.tot)) for (j in 1:ncol(juge.tot)) if (juge.tot[i, 
        j] == maxval) 
        res.id <- rbind(res.id, matrix(c(i, j), 1, 2))
    rownames(res.id) <- paste("Ideal_", LETTERS[1:nrow(res.id)], 
        sep = "")
    id.profile <- matrix(0, 0, nbatt)
    colnames(id.profile) <- attribut
    juge.max <- vector("list", nrow(res.id))
    names(juge.max) <- rownames(res.id)
    for (i in 1:nrow(res.id)) for (j in 1:nbjuge) if (juge.mat[[j]][res.id[i, 
        1], res.id[i, 2]] == 1) 
        juge.max[[i]] <- c(juge.max[[i]], j)
    id.j.avg <- averagetable(id.data, formul = paste("~", colnames(dataset)[col.j], 
        "+", colnames(dataset)[col.p]), firstvar = 3, method = "mean")
    for (i in 1:nrow(res.id)) id.profile <- rbind(id.profile, 
        t(as.matrix(apply(id.j.avg[juge.max[[i]], ], 2, mean))))
    id.profile <- t(as.matrix(apply(id.profile, 2, mean)))
    rownames(id.profile) <- "Ideal"
    res <- vector("list")
    res$PCA <- res.pca
    res$PCA$data <- data.pca
    res$PCA$dim <- coord
    res$idmap$data <- juge.tot
    res$idmap$j.weight <- cons.wgt
    res$idmap$precision <- precision
    res$ideal$profiles <- id.profile
    res$ideal$pct.conso <- maxval
	# res$exdensiteinertieH0<-vectinertiesimu #exemple permettant de tracer la densité de l'inertie sous H0 (cas de la derniere variable)
	# res$coord.ideaux<-res.pca$ind.sup.coord#peut être enlevée car contenue dans res$PCA, mais au moins y a pas a chercher les coordonnées des p idéaux de tous les juges (si p produits) qui sont en ind supp
	if(color == "FALSE") {
	  res$test.paires<-resultat.test#ne se fait que si on a choisi l'option color=false, si color = true, test.paires n'existe pas
	  res$coordobs<-listpointsvrai#coorddonnées des idéaux observés
	  res$test.global<-cbind.data.frame(inertie.obs, inertie.critique, inertie.pvalue)
    }
	class(res) <- c("IdMap", "list")
    return(res)
}
