#!/usr/bin/env Rscript

test_that("test RasterLayer", {
    library("raster")
    raster.object <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    raster.object[] <- c(3,2,1,4,3,2,5,4,3)

    result.layer <- rescale(raster.object)

    standard.result.layer <- raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    standard.result.layer[] <- c(0.6,0.4,0.2,0.8,0.6,0.4,1.00,0.8,0.6)

    extent.result <- compareRaster(result.layer, standard.result.layer)
    if (extent.result) {
        vector.one <- values(result.layer)
        vector.two <- values(standard.result.layer)
        expect_that(vector.one, equals(vector.two))
    } else {
        expect_that(extent.result, is_true())
    }
})

test_that("test RasterStack (simple)", {
    library("raster")
    env.one <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.one[] <- c(3,2,1,4,3,2,5,4,3)

    env.two <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.two[] <- c(3,3,5,5,3,4,9,5,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two))
    result.stack <- rescale(env.stack)

    layer.number <- nlayers(result.stack)

    expect_that(layer.number, equals(2))
})
