#!/usr/bin/env Rscript

test_that("test default  (simple)", {
    library("raster")
    env.one <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.one[] <- c(3,2,1,4,3,2,5,4,3)

    env.two <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.two[] <- c(3,3,5,5,3,4,9,5,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two))
    result.layer <- pickMean(env.stack)

    standard.result.layer <- raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    standard.result.layer[] <- c(1,1,0,1,1,1,0,1,1)

    extent.result <- compareRaster(result.layer, standard.result.layer)
    if (extent.result) {
        vector.one <- values(result.layer)
        vector.two <- as.logical(values(standard.result.layer))
        expect_that(vector.one, equals(vector.two))
    } else {
        expect_that(extent.result, is_true())
    }
})

test_that("test subset", {
    library("raster")
    env.one <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.one[] <- c(3,2,1,4,3,2,5,4,3)

    env.two <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.two[] <- c(3,3,5,5,3,4,9,5,3)

    env.three <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.three[] <- c(3,3,3,3,3,3,3,3,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two, env_three=env.three))
    subset <- c("env_one", "env_two")
    result.layer <- pickMean(env.stack, subset)

    standard.result.layer <- raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    standard.result.layer[] <- c(1,1,0,1,1,1,0,1,1)

    extent.result <- compareRaster(result.layer, standard.result.layer)
    if (extent.result) {
        vector.one <- values(result.layer)
        vector.two <- as.logical(values(standard.result.layer))
        expect_that(vector.one, equals(vector.two))
    } else {
        expect_that(extent.result, is_true())
    }
})

test_that("test stack (simple)", {
    library("raster")
    env.one <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.one[] <- c(3,2,1,4,3,2,5,4,3)

    env.two <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.two[] <- c(3,3,5,5,3,4,9,5,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two))

    result.stack <- pickMean(env.stack, stack=TRUE)

    layer.number <- nlayers(result.stack)

    expect_that(layer.number, equals(2))
})
