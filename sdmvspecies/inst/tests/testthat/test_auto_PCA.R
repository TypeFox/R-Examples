#!/usr/bin/env Rscript

test_that("test default", {
    library("raster")
    env.one <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.one[] <- c(3,2,1,4,3,2,5,4,3)

    env.two <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.two[] <- c(3,3,5,5,3,4,9,5,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two))

    result.stack <- autoPCA(env.stack, nfactors=1)
    result.layer <- result.stack[[1]]

    standard.result.layer <- raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    standard.result.layer[] <- c(3.3778134876722929469,
                                 2.8148445730602440484,
                                 3.3778134876722929469,
                                 5.0667202315084391984,
                                 3.3778134876722929469,
                                 3.3778134876722929469,
                                 7.8815648045686836909,
                                 5.0667202315084391984,
                                 3.3778134876722929469)

    extent.result <- compareRaster(result.layer, standard.result.layer)
    if (extent.result) {
       vector.one <- values(result.layer)
       vector.two <- values(standard.result.layer)
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

    env.three <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.three[] <- c(3,8,3,4,1,3,2,8,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two, env_three=env.three))

    result.stack <- autoPCA(env.stack, nfactors=2)

    layer.number <- nlayers(result.stack)

    expect_that(layer.number, equals(2))
})
