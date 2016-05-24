MS.test.clust <-
function (data_tot, nclust) 
{
attach(what=NULL, name="e1")
if(missing(data_tot)){
data_tot<-tk_choose.files(default=file.path(getwd(),"*.txt"),caption="Please, select a file for MS.test.clust", multi=FALSE, filters=matrix(c("your file","*.txt"),ncol=2))
data_tot<-read.table(data_tot, header=TRUE, check.names=FALSE)
}
    library(fpc)
    Dunn.test <- matrix(nrow = 7, ncol = 5)
    rownames(Dunn.test) <- c("average", "single", "complete", 
        "ward", "centroid", "diana", "PAM")
    colnames(Dunn.test) <- c("cor", "mink_1/2", "mink_1/3", "manhattan", 
        "euclidean")
    silhouette.width <- matrix(nrow = 7, ncol = 5)
    rownames(silhouette.width) <- c("average", "single", "complete", 
        "ward", "centroid", "diana", "PAM")
    colnames(silhouette.width) <- c("cor", "mink_1/2", "mink_1/3", 
        "manhattan", "euclidean")
    for (i in c("average", "single", "complete", "ward", "centroid")) {
        Dis <- as.dist((1 - cor(t(data_tot[, 4:ncol(data_tot)]), 
            method = "spearman")))
        Hc <- hclust(Dis, method = i)
        cl <- cutree(Hc, k = nclust)
        Test_valid <- cluster.stats(Dis, cl)
        Dunn.test[rownames(Dunn.test) == i, colnames(Dunn.test) == 
            "cor"] <- Test_valid$dunn
        silhouette.width[rownames(silhouette.width) == i, colnames(silhouette.width) == 
            "cor"] <- Test_valid$avg.silwidth
        Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
            p = 1/2)
        Hc <- hclust(Dis, method = i)
        cl <- cutree(Hc, k = nclust)
        Test_valid <- cluster.stats(Dis, cl)
        Dunn.test[rownames(Dunn.test) == i, colnames(Dunn.test) == 
            "mink_1/2"] <- Test_valid$dunn
        silhouette.width[rownames(silhouette.width) == i, colnames(silhouette.width) == 
            "mink_1/2"] <- Test_valid$avg.silwidth
        Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
            p = 1/3)
        Hc <- hclust(Dis, method = i)
        cl <- cutree(Hc, k = nclust)
        Test_valid <- cluster.stats(Dis, cl)
        Dunn.test[rownames(Dunn.test) == i, colnames(Dunn.test) == 
            "mink_1/3"] <- Test_valid$dunn
        silhouette.width[rownames(silhouette.width) == i, colnames(silhouette.width) == 
            "mink_1/3"] <- Test_valid$avg.silwidth
        Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
            p = 1)
        Hc <- hclust(Dis, method = i)
        cl <- cutree(Hc, k = nclust)
        Test_valid <- cluster.stats(Dis, cl)
        Dunn.test[rownames(Dunn.test) == i, colnames(Dunn.test) == 
            "manhattan"] <- Test_valid$dunn
        silhouette.width[rownames(silhouette.width) == i, colnames(silhouette.width) == 
            "manhattan"] <- Test_valid$avg.silwidth
        Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
            p = 2)
        Hc <- hclust(Dis, method = i)
        cl <- cutree(Hc, k = nclust)
        Test_valid <- cluster.stats(Dis, cl)
        Dunn.test[rownames(Dunn.test) == i, colnames(Dunn.test) == 
            "euclidean"] <- Test_valid$dunn
        silhouette.width[rownames(silhouette.width) == i, colnames(silhouette.width) == 
            "euclidean"] <- Test_valid$avg.silwidth
    }
    Dis <- as.dist((1 - cor(t(data_tot[, 4:ncol(data_tot)]), 
        method = "spearman")))
    Hc <- diana(Dis, diss = TRUE)
    cl <- cutree(Hc, k = nclust)
    Test_valid <- cluster.stats(Dis, cl)
    Dunn.test[rownames(Dunn.test) == "diana", colnames(Dunn.test) == 
        "cor"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "diana", colnames(silhouette.width) == 
        "cor"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 1/2)
    Hc <- diana(Dis, diss = TRUE)
    cl <- cutree(Hc, k = nclust)
    Test_valid <- cluster.stats(Dis, cl)
    Dunn.test[rownames(Dunn.test) == "diana", colnames(Dunn.test) == 
        "mink_1/2"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "diana", colnames(silhouette.width) == 
        "mink_1/2"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 1/3)
    Hc <- diana(Dis, diss = TRUE)
    cl <- cutree(Hc, k = nclust)
    Test_valid <- cluster.stats(Dis, cl)
    Dunn.test[rownames(Dunn.test) == "diana", colnames(Dunn.test) == 
        "mink_1/3"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "diana", colnames(silhouette.width) == 
        "mink_1/3"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 1)
    Hc <- diana(Dis, diss = TRUE)
    cl <- cutree(Hc, k = nclust)
    Test_valid <- cluster.stats(Dis, cl)
    Dunn.test[rownames(Dunn.test) == "diana", colnames(Dunn.test) == 
        "manhattan"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "diana", colnames(silhouette.width) == 
        "manhattan"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 2)
    Hc <- diana(Dis, diss = TRUE)
    cl <- cutree(Hc, k = nclust)
    Test_valid <- cluster.stats(Dis, cl)
    Dunn.test[rownames(Dunn.test) == "diana", colnames(Dunn.test) == 
        "euclidean"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "diana", colnames(silhouette.width) == 
        "euclidean"] <- Test_valid$avg.silwidth
    Dis <- as.dist((1 - cor(t(data_tot[, 4:ncol(data_tot)]), 
        method = "spearman")))
    Hc <- pamk(Dis, krange = nclust)
    Test_valid <- cluster.stats(Dis, Hc$pamobject$clustering)
    Dunn.test[rownames(Dunn.test) == "PAM", colnames(Dunn.test) == 
        "cor"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "PAM", colnames(silhouette.width) == 
        "cor"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 1/2)
    Hc <- pamk(Dis, krange = nclust)
    Test_valid <- cluster.stats(Dis, Hc$pamobject$clustering)
    Dunn.test[rownames(Dunn.test) == "PAM", colnames(Dunn.test) == 
        "mink_1/2"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "PAM", colnames(silhouette.width) == 
        "mink_1/2"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 1/3)
    Hc <- pamk(Dis, krange = nclust)
    Test_valid <- cluster.stats(Dis, Hc$pamobject$clustering)
    Dunn.test[rownames(Dunn.test) == "PAM", colnames(Dunn.test) == 
        "mink_1/3"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "PAM", colnames(silhouette.width) == 
        "mink_1/3"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 1)
    Hc <- pamk(Dis, krange = nclust)
    Test_valid <- cluster.stats(Dis, Hc$pamobject$clustering)
    Dunn.test[rownames(Dunn.test) == "PAM", colnames(Dunn.test) == 
        "manhattan"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "PAM", colnames(silhouette.width) == 
        "manhattan"] <- Test_valid$avg.silwidth
    Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
        p = 2)
    Hc <- pamk(Dis, krange = nclust)
    Test_valid <- cluster.stats(Dis, Hc$pamobject$clustering)
    Dunn.test[rownames(Dunn.test) == "PAM", colnames(Dunn.test) == 
        "euclidean"] <- Test_valid$dunn
    silhouette.width[rownames(silhouette.width) == "PAM", colnames(silhouette.width) == 
        "euclidean"] <- Test_valid$avg.silwidth
    homogeneite <- matrix(nrow = 7, ncol = 5)
    rownames(homogeneite) <- c("average", "single", "complete", 
        "ward", "centroid", "diana", "PAM")
    colnames(homogeneite) <- c("cor", "mink_1/2", "mink_1/3", 
        "manhattan", "euclidean")
    k <- 1
    lis <- c("average", "single", "complete", "ward", "centroid")
    lisd <- c(1/2, 1/3, 1, 2)
    for (i in 1:length(lis)) for (j in 1:length(lisd)) {
        nb <- 0
        Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
            p = lisd[j])
        Hc <- hclust(Dis, method = lis[i])
        cl <- cutree(Hc, k = nclust)
        res <- as.data.frame(cbind(as.vector(data_tot[, 1]), 
            cl))
        clust <- levels(res$cl)
        homog_clust <- list()
        for (l in 1:length(clust)) {
            res_temp <- subset(res, res$cl == clust[l])
            temp_mol <- levels(as.factor(as.vector(res_temp$V1)))
            homog_clust[[l]] <- temp_mol
        }
        for (g in (1:length(homog_clust))) {
            if (length(homog_clust[[g]]) == 1) 
                nb <- nb + 1
        }
        homogeneite[i, j + 1] <- nb/nclust
    }
    for (i in 1:length(lis)) {
        nb <- 0
        Dis <- as.dist((1 - cor(t(data_tot[, 4:ncol(data_tot)]), 
            method = "spearman")))
        Hc <- hclust(Dis, method = lis[i])
        cl <- cutree(Hc, k = nclust)
        res <- as.data.frame(cbind(as.vector(data_tot[, 1]), 
            cl))
        clust <- levels(res$cl)
        homog_clust <- list()
        for (l in 1:length(clust)) {
            res_temp <- subset(res, res$cl == clust[l])
            temp_mol <- levels(as.factor(as.vector(res_temp$V1)))
            homog_clust[[l]] <- temp_mol
        }
        for (g in (1:length(homog_clust))) {
            if (length(homog_clust[[g]]) == 1) 
                nb <- nb + 1
        }
        homogeneite[i, 1] <- nb/nclust
    }
    for (j in 1:length(lisd)) {
        nb <- 0
        Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
            p = lisd[j])
        Hc <- diana(Dis, diss = TRUE)
        cl <- cutree(Hc, k = nclust)
        res <- as.data.frame(cbind(as.vector(data_tot[, 1]), 
            cl))
        clust <- levels(res$cl)
        homog_clust <- list()
        for (l in 1:length(clust)) {
            res_temp <- subset(res, res$cl == clust[l])
            temp_mol <- levels(as.factor(as.vector(res_temp$V1)))
            homog_clust[[l]] <- temp_mol
        }
        for (g in (1:length(homog_clust))) {
            if (length(homog_clust[[g]]) == 1) 
                nb <- nb + 1
        }
        homogeneite[6, j + 1] <- nb/nclust
    }
    for (j in 1:length(lisd)) {
        nb <- 0
        Dis <- dist(data_tot[, 4:ncol(data_tot)], method = "minkowski", 
            p = lisd[j])
        Hc <- pamk(Dis, krange = nclust)
        cl <- Hc$pamobject$clustering
        res <- as.data.frame(cbind(as.vector(data_tot[, 1]), 
            cl))
        clust <- levels(res$cl)
        homog_clust <- list()
        for (l in 1:length(clust)) {
            res_temp <- subset(res, res$cl == clust[l])
            temp_mol <- levels(as.factor(as.vector(res_temp$V1)))
            homog_clust[[l]] <- temp_mol
        }
        for (g in (1:length(homog_clust))) {
            if (length(homog_clust[[g]]) == 1) 
                nb <- nb + 1
        }
        homogeneite[7, j + 1] <- nb/nclust
    }
    nb <- 0
    Dis <- as.dist((1 - cor(t(data_tot[, 4:ncol(data_tot)]), 
        method = "spearman")))
    Hc <- diana(Dis, diss = TRUE)
    cl <- cutree(Hc, k = nclust)
    res <- as.data.frame(cbind(as.vector(data_tot[, 1]), cl))
    clust <- levels(res$cl)
    homog_clust <- list()
    for (l in 1:length(clust)) {
        res_temp <- subset(res, res$cl == clust[l])
        temp_mol <- levels(as.factor(as.vector(res_temp$V1)))
        homog_clust[[l]] <- temp_mol
    }
    for (g in (1:length(homog_clust))) {
        if (length(homog_clust[[g]]) == 1) 
            nb <- nb + 1
    }
    homogeneite[6, 1] <- nb/nclust
    nb <- 0
    Dis <- as.dist((1 - cor(t(data_tot[, 4:ncol(data_tot)]), 
        method = "spearman")))
    Hc <- pamk(Dis, krange = nclust)
    cl <- Hc$pamobject$clustering
    res <- as.data.frame(cbind(as.vector(data_tot[, 1]), cl))
    clust <- levels(res$cl)
    homog_clust <- list()
    for (l in 1:length(clust)) {
        res_temp <- subset(res, res$cl == clust[l])
        temp_mol <- levels(as.factor(as.vector(res_temp$V1)))
        homog_clust[[l]] <- temp_mol
    }
    for (g in (1:length(homog_clust))) {
        if (length(homog_clust[[g]]) == 1) 
            nb <- nb + 1
    }
    homogeneite[7, 1] <- nb/nclust
	
	#creation of a directory 
	st <- strsplit(date(), " ")[[1]]
    stBis <- strsplit(st[4], ":")[[1]]
    Hour <- paste(stBis[1], stBis[2], stBis[3], sep = "-")
    Date <- paste(st[1], st[2], st[3], Hour, sep = "_")
    Mypath <- paste("output_MStestclust", "_", "result", Date, sep = "")
    dir.create(Mypath)
	
	
		par(mfrow=c(2,1))
		test<-matrix(homogeneite,nrow=7)
		ar<-names(homogeneite[,1])
		ar2<-names(homogeneite[1,])
		matplot( test, type = "b", lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), cex = 0.8, col = c("red", "blue","black","green","orange"), xlab = NULL, ylab = "Matching coefficient", main = "Matching coefficient", axes=FALSE)
		axis(1,seq(7), labels=ar, las=2)
		axis(2)
		legend(cex=0.5,"bottomright",legend=ar2, col=c("red", "blue","black","green","orange"), lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), bg="white")
		test<-matrix(silhouette.width,nrow=7)
		ar<-names(silhouette.width[,1])
		ar2<-names(silhouette.width[1,])
		matplot( test, type = "b", lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), cex = 0.8, col = c("red", "blue","black","green","orange"), xlab = NULL, ylab = "Silhouette width", main = "Silhouette Width", axes=FALSE)
		axis(1,seq(7), labels=ar, las=2)
		axis(2)
		legend(cex=0.5,"bottomright",legend=ar2, col=c("red", "blue","black","green","orange"), lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), bg="white")
	
	#creation of a pdf file with graphics
	pdf(paste(Mypath, "/", "Graph_MStestClust.pdf", sep = ""))
		par(mfrow=c(2,1))
		test<-matrix(homogeneite,nrow=7)
		ar<-names(homogeneite[,1])
		ar2<-names(homogeneite[1,])
		matplot( test, type = "b", lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), cex = 0.8, col = c("red", "blue","black","green","orange"), xlab = NULL, ylab = "Matching coefficient", main = "Matching coefficient", axes=FALSE)
		axis(1,seq(7), labels=ar, las=2)
		axis(2)
		legend(cex=0.5,"bottomright",legend=ar2, col=c("red", "blue","black","green","orange"), lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), bg="white")
		test<-matrix(silhouette.width,nrow=7)
		ar<-names(silhouette.width[,1])
		ar2<-names(silhouette.width[1,])
		matplot( test, type = "b", lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), cex = 0.8, col = c("red", "blue","black","green","orange"), xlab = NULL, ylab = "Silhouette width", main = "Silhouette Width", axes=FALSE)
		axis(1,seq(7), labels=ar, las=2)
		axis(2)
		legend(cex=0.5,"bottomright",legend=ar2, col=c("red", "blue","black","green","orange"), lty = c("dotted", "solid"), lwd = 2, pch = c(1, 2,3), bg="white")
	dev.off()
	par(mfrow=c(1,1))
    assign("Dunn.test",Dunn.test ,envir=as.environment("e1"))
	assign("silhouette.width", silhouette.width,envir=as.environment("e1"))
	assign("matching.coef",homogeneite,envir=as.environment("e1"))	
	print(paste("A pdf file have been generated in the path:", Mypath, cat("\n")))

}

