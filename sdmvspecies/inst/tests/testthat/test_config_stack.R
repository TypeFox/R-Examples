#!/usr/bin/env Rscript

test_that("test default  (simple)", {
    library("raster")
    env.one <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.one[] <- c(1,1,1,1,1,1,1,1,1)

    env.two <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.two[] <- c(2,2,2,2,2,2,2,2,2)

    env.three <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.three[] <- c(3,3,3,3,3,3,3,3,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two, env_three=env.three))
    config <- list(c("env_one", 3, 2), c("env_two", 4, 2))
    result.layer <- configStack(env.stack, config)

    layer.number <- nlayers(result.layer)

    expect_that(layer.number, equals(2))
})
