

test.umarkedMultFrame.subset <- function() {

    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)

    dat <- as(umf1, "data.frame")

    umf1.obs1 <- umf1[,1]
    checkEquals(umf1.obs1@y, y[,1:3])
    checkEquals(umf1.obs1@siteCovs, sc)
    checkEqualsNumeric(unlist(umf1.obs1@obsCovs),
                       as.numeric(t(oc[[1]][,1:3])))
    checkEqualsNumeric(unlist(umf1.obs1@yearlySiteCovs), ysc[[1]][,1])
    checkEquals(umf1.obs1@numPrimary, 1)

    umf1.obs1and3 <- umf1[,c(1,3)]

    umf1.site1 <- umf1[1,]
    checkEquals(umf1.site1@y, y[1,, drop=FALSE])
    checkEquals(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    checkEqualsNumeric(unlist(umf1.site1@obsCovs), oc$x3[1,])
    checkEqualsNumeric(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    checkEquals(umf1.site1@numPrimary, 3)

    umf1.sites1and3 <- umf1[c(1,3),]

    }






test.umarkedFrameGMM.subset <- function() {

    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedFrameGMM(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3,
        type="removal")

    dat <- as(umf1, "data.frame")

    umf1.site1 <- umf1[1,]
    checkEquals(umf1.site1@y, y[1,, drop=FALSE])
    checkEquals(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    checkEqualsNumeric(unlist(umf1.site1@obsCovs), oc$x3[1,])
    checkEqualsNumeric(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    checkEquals(umf1.site1@numPrimary, 3)

    umf1.sites1and3 <- umf1[c(1,3),]

    checkEquals(class(umf1.site1)[1], "unmarkedFrameGMM")

    umf1.sites1and1 <- umf1[c(1,1),]

    umf1.obs1and2 <- umf1[,c(1,2)]

    checkEqualsNumeric(dim(getY(umf1.obs1and2)), c(3,6))
    checkEqualsNumeric(dim(siteCovs(umf1.obs1and2)), c(3,1))
    checkEqualsNumeric(dim(obsCovs(umf1.obs1and2)), c(18,1))

    umf1.sites1and2.obs1and2 <- umf1[c(1,2),c(1,2)]
    checkEqualsNumeric(dim(getY(umf1.sites1and2.obs1and2)), c(2,6))
    checkEqualsNumeric(dim(siteCovs(umf1.sites1and2.obs1and2)), c(2,1))
    checkEqualsNumeric(dim(obsCovs(umf1.sites1and2.obs1and2)), c(12,1))

    # THis doesn't work
    umf1.sites1and1.obs1and1 <- umf1[c(1,1),c(1,1)]


    }




test.umarkedFrameGDS.subset <- function() {

    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))

    umf1 <- unmarkedFrameGDS(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        numPrimary = 3,
        survey="point",
        dist.breaks=c(0, 10, 20, 30),
        unitsIn="m")

    dat <- as(umf1, "data.frame")
    checkEquals(nrow(dat), nrow(y))

    umf1.site1 <- umf1[1,]
    checkEquals(umf1.site1@y, y[1,, drop=FALSE])
    checkEquals(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    checkEqualsNumeric(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    checkEquals(umf1.site1@numPrimary, 3)
    checkEquals(umf1.site1@survey, "point")

    umf1.sites1and3 <- umf1[c(1,3),]


    umf2 <- unmarkedFrameGDS(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        numPrimary = 3,
        survey="line",
        dist.breaks=c(0, 10, 20, 30),
        tlength=rep(1,nrow(y)),
        unitsIn="m")

    dat <- as(umf2, "data.frame")

    umf2.site1 <- umf2[1,]
    checkEquals(umf2.site1@y, y[1,, drop=FALSE])
    checkEquals(umf2.site1@siteCovs, sc[1,, drop=FALSE])
    checkEqualsNumeric(unlist(umf2.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    checkEquals(umf2.site1@numPrimary, 3)
    checkEquals(umf2.site1@survey, "line")

    umf2.sites1and3 <- umf2[c(1,3),]

    }









test.umarkedFrameGPC.subset <- function() {

    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedFrameGPC(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)

    dat <- as(umf1, "data.frame")

    umf1.site1 <- umf1[1,]
    checkEquals(umf1.site1@y, y[1,, drop=FALSE])
    checkEquals(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    checkEqualsNumeric(unlist(umf1.site1@obsCovs), oc$x3[1,])
    checkEqualsNumeric(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    checkEquals(umf1.site1@numPrimary, 3)

    umf1.sites1and3 <- umf1[c(1,3),]

    checkEquals(class(umf1.site1)[1], "unmarkedFrameGPC")

    umf1.sites1and1 <- umf1[c(1,1),]

    umf1.obs1and2 <- umf1[,c(1,2)]

    checkEqualsNumeric(dim(getY(umf1.obs1and2)), c(3,6))
    checkEqualsNumeric(dim(siteCovs(umf1.obs1and2)), c(3,1))
    checkEqualsNumeric(dim(obsCovs(umf1.obs1and2)), c(18,1))

    umf1.sites1and2.obs1and2 <- umf1[c(1,2),c(1,2)]
    checkEquals(class(umf1.sites1and2.obs1and2)[1], "unmarkedFrameGPC")
    checkEqualsNumeric(dim(getY(umf1.sites1and2.obs1and2)), c(2,6))
    checkEqualsNumeric(dim(siteCovs(umf1.sites1and2.obs1and2)), c(2,1))
    checkEqualsNumeric(dim(obsCovs(umf1.sites1and2.obs1and2)), c(12,1))

    # THis doesn't work
    umf1.sites1and1.obs1and1 <- umf1[c(1,1),c(1,1)]


    }







test.umarkedFramePCO.subset <- function() {

    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedFramePCO(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)

    dat <- as(umf1, "data.frame")

    umf1.site1 <- umf1[1,]
    checkEquals(umf1.site1@y, y[1,, drop=FALSE])
    checkEquals(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    checkEqualsNumeric(unlist(umf1.site1@obsCovs), oc$x3[1,])
    checkEqualsNumeric(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    checkEquals(umf1.site1@numPrimary, 3)
    checkEquals(class(umf1.site1)[1], "unmarkedFramePCO")

    umf1.sites1and3 <- umf1[c(1,3),]

    }
