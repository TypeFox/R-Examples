BSL <-
function(band3, band4, method = "quantile", ulimit = .99, llimit = .005)
{
    # find Bare Soil Line and vegetation peak

    if(is.character(band3)) {
        band3 <- read.asciigrid(band3)
        band3 <- band3@data[,1]
    } else {
        if(class(band3) == "SpatialGridDataFrame") {
            band3 <- band3@data[,1]
        } else {
            band3 <- as.vector(as.matrix(band3))
        }
    } 
    
    if(is.character(band4)) {
        band4 <- read.asciigrid(band4)
        band4 <- band4@data[,1]
    } else {
        if(class(band4) == "SpatialGridDataFrame") {
            band4 <- band4@data[,1]
        } else {
            band4 <- as.vector(as.matrix(band4))
        }
    } 


    # find joint minimum and maximum
    bsl.joint <- cbind(band3, band4)
    bsl.joint <- bsl.joint[apply(bsl.joint, 1, function(x)!any(is.na(x))), ]
    bsl.joint <- bsl.joint[apply(bsl.joint, 1, function(x)all(x < 255)), ]

    ratio43 <- bsl.joint[,2]/bsl.joint[,1]

    if(method == "quantile") {
        bsl.lmodel2 <- lmodel2(bsl.joint[ratio43 < quantile(ratio43, llimit), 2] ~ bsl.joint[ratio43 < quantile(ratio43, llimit), 1])
    }
    else if(method == "minimum") {
        # want lowest band4 value for each band3 value (lowest NIR for each red)
        bsl.min <- factor(bsl.joint[,1], levels=1:254)
        bsl.min <- split(bsl.joint[,2], bsl.min, drop=TRUE)
        bsl.min <- lapply(bsl.min, min)

        bsl.lmodel2 <- lmodel2(as.numeric(bsl.min) ~ as.numeric(names(bsl.min)))
    }
    else {
        stop("Method not found.\n")
    }

    bsl.lm <- unlist(bsl.lmodel2$regression.results[2, 2:3])
    names(bsl.lm) <- c("Intercept", "Slope")

    ### next, find top vegetation point

    bsl.test <- bsl.joint
    bsl.test[,2] <- 255 - bsl.test[,2] # want high values of band 4
    bsl.test <- apply(bsl.test, 1, sum)

    # want high veg cover
    bsl.top <- bsl.joint[ratio43 > quantile(ratio43, ulimit, na.rm=TRUE), ]
    bsl.test <- bsl.test[ratio43 > quantile(ratio43, ulimit, na.rm=TRUE)]
    bsl.top <- bsl.top[bsl.test == min(bsl.test), ]
    if(!is.null(dim(bsl.top))) bsl.top <- bsl.top[sample(1:nrow(bsl.top), 1),]

    list(BSL=bsl.lm, top=bsl.top)
}

