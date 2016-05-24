#!/usr/bin/env Rscript

test_that("test default", {
    library("raster")
    env.one <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.one[] <- c(3,2,1,4,3,2,5,4,3)

    env.two <- raster::raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    env.two[] <- c(3,3,5,5,3,4,9,5,3)

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two))
    config <- list(c("env_one", "1", 2), c("env_two", "2", 4))
    result.layer <- nicheSynthese(env.stack, config)

    standard.result.layer <- raster(nrow=3, ncol=3, xmn=0, xmx=3, ymn=0, ymx=3, crs=NA)
    standard.result.layer[] <- c(0.3333333333333333148296,
                                 0.1082174891194499083413,
                                 0.2259252210683029560290,
                                 0.3304397113416721043500,
                                 0.3333333333333333148296,
                                 0.2193286002305609994067,
                                 0.6703696655127474590685,
                                 0.3304397113416721043500,
                                 0.3333333333333333148296)

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

    env.stack <- raster::stack(list(env_one=env.one, env_two=env.two))
    config <- list(c("env_one", 3, 2), c("env_two", 4, 2))
    result.stack <- artificialBellResponse(env.stack, config, stack=TRUE)

    layer.number <- nlayers(result.stack)

    expect_that(layer.number, equals(2))
})
