`ensemble.centroids` <- function(
    presence.raster=NULL, x=NULL, categories.raster=NULL,
    an=10000, ext=NULL, name="Species001", 
    pca.var=0.95, centers=0, use.silhouette=TRUE,
    plotit=FALSE, dev.new.width=7, dev.new.height=7 
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(presence.raster) == T) {stop("value for parameter presence.raster is missing (RasterLayer object)")}
    if(inherits(presence.raster, "RasterLayer") == F) {stop("x is not a RasterLayer object")}
    if(raster::minValue(presence.raster) != 0) {cat(paste("\n", "WARNING: minValue of presence.raster not 0, so this raster layer possibly does not indicate presence-absence", "\n", sep = ""))}
    if(raster::maxValue(presence.raster) != 1) {cat(paste("\n", "WARNING: maxValue of presence.raster not 1, so this raster layer possibly does not indicate presence-absence", "\n", sep = ""))}
    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x, "RasterStack") == F) {stop("x is not a RasterStack object")}

# only for plotting
predict.zone <- function(object=centroid.model, newdata=newdata) {
    centroids <- object$centroids
    cov.mahal <- object$cov.mahal
    nc <- nrow(centroids)
    result <- data.frame(array(0, dim=c(nrow(newdata), nc)))
    for (i in 1:nc) {
        result[,i] <- mahalanobis(newdata, center=as.numeric(centroids[i,]), cov=cov.mahal)
    }
    p <- apply(result[, 1:nc], 1, which.min)
    p <- as.numeric(p)
    return(p)
}

# same extent for predictors and presence map
    if (is.null(ext) == F) {
        if(length(x@title) == 0) {x@title <- "stack1"}
        title.old <- x@title
        x <- raster::crop(x, y=ext, snap="in")
        x@title <- title.old
        presence.raster <- raster::crop(presence.raster, y=ext, snap="in")
        if (is.null(categories.raster) == F) {categories.raster <- raster::crop(categories.raster, y=ext, snap="in")}
    }

# mask the presence area
    presence.raster <- raster::mask(presence.raster, presence.raster, inverse=T, maskvalue=1)

# create background data
    a <- dismo::randomPoints(presence.raster, n=an, p=NULL, excludep=F)
    a <- data.frame(a)
    background.data <- raster::extract(x, a)
    background.data <- data.frame(background.data)
    TrainValid <- complete.cases(background.data)
    a <- a[TrainValid,]
    background.data <- background.data[TrainValid,]

# PCA of scaled variables
    rda.result <- vegan::rda(X=background.data, scale=T)
# select number of axes
    ax <- 2
    while ( (sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result))) < pca.var ) {ax <- ax+1}
    cat(paste("\n", "Percentage of variance of the selected axes (1 to ", ax, ") of principal components analysis", "\n", sep = ""))
    print(100*sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result)))
    rda.scores <- vegan::scores(rda.result, display="sites", scaling=1, choices=c(1:ax))

# K-means analysis
    if (is.null(categories.raster) == T) {
        if (centers < 1) {    
            cat(paste("\n", "Number of centroids determined with cascadeKM (2 to 16 clusters, Calinski-Harabasz criterion)", "\n", sep = ""))
            cascaderesult <- vegan::cascadeKM(rda.scores, inf.gr=2, sup.gr=16, iter=200, criterion="calinski")
            print(cascaderesult$results)
            groupnumbers <- (as.numeric(gsub(" groups", "", colnames(cascaderesult$results))))
            w <- cascaderesult$results[2,]
            maxx = which.max(w[])
            centers <- groupnumbers[maxx]
            cat(paste("\n", "Selected number of centroids: ", centers, "\n", sep = ""))
        }
        kmeans.result <- kmeans(rda.scores, centers=centers, iter.max=1000)
        clusters <- kmeans.result$cluster
        clusters.remember <- clusters
        clusters.names <- c(1:centers)

# predefined categories
    }else{
        categories.data <- raster::extract(categories.raster, a)        
        categories.data <- data.frame(categories.data)
        categories <- levels(as.factor(categories.data[,1]))
        categories <- categories[is.na(categories) == F]
        clusters.names <- as.numeric(categories)
        clusters <- match(categories.data[,1], categories)
        clusters.remember <- clusters
        centers <- length(categories)
    }

# centroids
    centroid.data <- background.data[1:centers,]
    centroid.rda <- rda.scores[1:centers,]
    centroid.data[] <- NA
    centroid.rda[] <- NA
    if (use.silhouette == T) {
        background.silhouette <- cluster::silhouette(clusters, vegan::vegdist(rda.scores, method="euc"))
        clusters[background.silhouette[,"sil_width"] <= 0] <- -1
    }
    for (i in c(1:centers)) {
        centroid.data[i,] <- colMeans(background.data[as.numeric(clusters)==i,])
        centroid.rda[i,] <- colMeans(rda.scores[as.numeric(clusters)==i,])
    }

# calculate variance for Mahalanobis distance
    cov.mahal <- cov(background.data)

# find analog locations that are closest to the centroids in PCA space
    centroid.analogs <- cbind(a[1:centers,], background.data[1:centers,], pca.dist=as.numeric(rep(NA, centers)))
    centroid.analogs[] <- NA
    remember.closest <- numeric(centers)
    for (i in c(1:centers)) {
        centroid.rda1 <- rbind(centroid.rda[i,], rda.scores)
        euc.dist <- as.matrix(vegan::vegdist(centroid.rda1, "euc"))
        euc.dist <- euc.dist[1,]
        euc.dist <- euc.dist[-1]
	closest <- as.numeric(which.min(euc.dist))
        remember.closest[i] <- closest
        centroid.analogs[i, as.numeric(na.omit(match(names(a), names(centroid.analogs))))] <- a[closest,]
        centroid.analogs[i, as.numeric(na.omit(match(names(background.data), names(centroid.analogs))))] <- background.data[closest,]
        centroid.analogs[i, as.numeric(na.omit(match("pca.dist", names(centroid.analogs))))] <- euc.dist[closest]
    }

# output
    centroid.model <- list(centroids=centroid.data, centroid.analogs=centroid.analogs, cov.mahal=cov.mahal, name=name, clusters.names=clusters.names)  

# plotting
    if (plotit == T) {
        par.old <- graphics::par(no.readonly=T)
        if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
        graphics::par(mfrow=c(2,2))
        zone <- predict.zone(centroid.model, newdata=background.data)
        # simple plot in geographic space with K-means clustering
        graphics::plot(a[, 2] ~ a[, 1], pch=15, col=grDevices::rainbow(centers)[as.numeric(clusters.remember)], 
            main="zones (based on K-means clustering in PCA) and locations of centroid analogs", xlab=names(a)[1], ylab=names(a)[2])
        graphics::points(centroid.analogs[, 2] ~ centroid.analogs[, 1], pch=20)
        graphics::text(centroid.analogs[, 2] ~ centroid.analogs[, 1], labels=c(1:centers), pos=3)
        graphics::legend(x="topright", legend=c(1:centers), pch=rep(15, centers), col=grDevices::rainbow(centers))
        # simple plot in geographic space with Mahalanobis clustering
        graphics::plot(a[, 2] ~ a[, 1], pch=15, col=grDevices::rainbow(centers)[as.numeric(zone)], 
            main="zones (based on Mahalanobis distance from centroid) and locations of centroid analogs", xlab=names(a)[1], ylab=names(a)[2])
        graphics::points(centroid.analogs[, 2] ~ centroid.analogs[, 1], pch=20)
        graphics::text(centroid.analogs[, 2] ~ centroid.analogs[, 1], labels=c(1:centers), pos=3)
        graphics::legend(x="topright", legend=c(1:centers), pch=rep(15, centers), col=grDevices::rainbow(centers))
        # plot in PCA space
        graphics::plot(rda.scores[, 2] ~ rda.scores[, 1], pch=15, col=grDevices::rainbow(centers)[as.numeric(zone)], main="locations of centroids (bullets) and analogs (asterisks)")
        graphics::points(centroid.rda[, 2] ~ centroid.rda[, 1], pch=20)
        graphics::text(centroid.rda[, 2] ~ centroid.rda[, 1], labels=c(1:centers), pos=3)
        graphics::points(rda.scores[remember.closest, 2] ~ rda.scores[remember.closest, 1], pch=8)
        graphics::legend(x="topright", legend=c(1:centers), pch=rep(15, centers), col=grDevices::rainbow(centers))
        if (ax >= 4) {
            graphics::plot(rda.scores[, 4] ~ rda.scores[, 3], pch=15, col=grDevices::rainbow(centers)[as.numeric(zone)], main="locations of centroids (bullets) and analogs (asterisks)")
            graphics::points(centroid.rda[, 4] ~ centroid.rda[, 3], pch=20)
            graphics::text(centroid.rda[, 4] ~ centroid.rda[, 3], labels=c(1:centers), pos=3)
            graphics::points(rda.scores[remember.closest, 4] ~ rda.scores[remember.closest, 3], pch=8)
            graphics::legend(x="topright", legend=c(1:centers), pch=rep(15, centers), col=grDevices::rainbow(centers))
        }
        graphics::par(par.old)
    }

# output
    return(centroid.model)
}

