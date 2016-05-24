

context('Brownian bridge dyn')
test_that('dbbmm error handling',{
  data(leroy)
  expect_error(
    brownian.bridge.dyn(
      object = leroy2, location.error = 23.5, dimSize = 150, ext = .3, time.step =
        600
    ), "object 'leroy2' not found"
  )
  data2 <- spTransform(leroy[1:50,], center = T)
  expect_error(
    brownian.bridge.dyn(data2, dimSize = 150 ,location.error = '2346'), 'column indicated for location error probably does not exist'
  )# character loc
  expect_error(
    brownian.bridge.dyn(data2,  location.error = 1:5),'The location error vector has not the same length as the move object'
  )
})
test_that('dbbmm error handling in relation to projections',{
  data(leroy)
  data <- leroy
  expect_error(
    brownian.bridge.dyn(
      object = data, location.error = 23.5, dimSize = 150, ext = .3, time.step =
        600
    )
    , 'You can not use longitude latitude projection for this function. To transform your coordinates use the spTransform function.'
  )
  r <- raster(
    nrows = 100, ncols = 100, xmn = 0, xmx = 10
  )
  expect_error(
    brownian.bridge.dyn(
      data2 <-
        spTransform(data[1:50,], center = T), raster = r, location.error = 20
    ), 'The projection of the raster and the Move object are not equal'
  )# equal projection
})
test_that('brownian bridge dyn for bursted',{
  data(leroy)
  data <- leroy
  dataP <- spTransform(data[1:100,], center = T)
  dataPB <- move::burst(dataP, round((1:(n.locs(
    dataP
  ) - 1) / 50)))
  expect_error(brownian.bridge.dyn(
    dataPB, dimSize = 150, location.error = 23, ext = .3, time.step = 4, window.size =29, burstType="l"), 'none of the burstTypes is in the data')
  expect_equivalent(class(brownian.bridge.dyn(
    dataPB, dimSize = 150, location.error = 23, ext = .3, time.step = 4, window.size =29, burstType=c("l","X1"))), 'DBBMMBurstStack')
  
    expect_message(
    udS <-
      brownian.bridge.dyn(
        dataPB, dimSize = 150, location.error = 23, ext = .3, time.step = 4, window.size =
          29
      ), 'Computational size'
  )
  expect_is(udS,"DBBMMBurstStack")
  expect_true(validObject(udS))
  expect_message(
    ud <-
      brownian.bridge.dyn(
        dataP, dimSize = 150, location.error = 23, ext = .3, time.step = 4, window.size =
          29
      ), 'Computational size'
  )
  expect_equal(udS@DBMvar@means,
               ud@DBMvar@means)
  expect_true(validObject(ud))
  b <- as.numeric(getZ(udS))
  expect_equal(colSums(values(udS)),
               b / sum(b), check.attributes = F)
})
test_that("verbosity brownian bridge works",{
  skip("will only work with testthat .11")
  p <- seq(0, 2 * pi, length.out = 49)
  tmp <-
    move(sin(p), cos(p), as.POSIXct(1:length(p), origin = '1970-1-1'), proj =
           '+proj=aeqd +ellps=WGS84')
  t <- .05
  expect_message(
    u <-
      brownian.bridge.dyn(
        tmp, dimSize = 200, location.error = .1, time.step = t
        , verbose = T
      ), 'Computational size'
  )
  expect_silent(brownian.bridge.dyn(
    tmp, dimSize = 200, location.error = .1, time.step = t, verbose = F
  ))
})
test_that('brownian bridge dyn value comparison bursted',{
  p <- seq(0, 2 * pi, length.out = 49)
  tmp <-
    move(sin(p), cos(p), as.POSIXct(1:length(p), origin = '1970-1-1'), proj =
           '+proj=aeqd +ellps=WGS84')
  t <- .05
  expect_message(
    u <-
      brownian.bridge.dyn(
        tmp, dimSize = 200, location.error = .1, time.step = t
      ), 'Computational size'
  )

  expect_message(
    us <-
      brownian.bridge.dyn(
        move::burst(tmp, round(1:length(p[-1]) / 30)), dimSize = 200, location.error =
          .1, time.step = t
      ), 'Computational size'
  )
  expect_equal(values(u),
               c(values(u)))
  p <- seq(0, 2 * pi, length.out = 199)
  set.seed(3245)
  tmp <- move(
    sin(p) + rnorm(length(p)) * round(sin(p)) * .03,
    cos(p),
    as.POSIXct(1:length(p),
               origin = '1970-1-1'),
    proj = '+proj=aeqd +ellps=WGS84'
  )
  t <- .25
  expect_message(
    u <- brownian.bridge.dyn(
      tmp, dimSize = 500,
      location.error = .01, time.step = t,
      ext = .2, margin = 15
    ), 'Computational size'
  )
  expect_message(
    us <- brownian.bridge.dyn(
      move::burst(tmp,
                  round(1:length(p[-1]) / 30)),
      dimSize = 500, location.error = .01,
      time.step = t, ext = .2, margin = 15
    ), 'Computational size'
  )
  expect_true(all(abs(values(sum(
    us
  )) - values(u)) < .Machine$double.eps))
})
test_that('Brownian bridge, running with character and vector input to location error',{
  data(leroy)
  data <- spTransform(leroy[1:51,], center = T)
  data2 <- data
  lc <- 6.54
  data2$asdf <- lc
  expect_message(
    a <-
      brownian.bridge.dyn(
        data, location.error = lc, ext = 2, dimSize = 50
      ), 'Computational size'
  )
  a@DBMvar$asdf <- lc
  expect_equal(a,
               suppressMessages(
                 brownian.bridge.dyn(
                   data2,
                   location.error = "asdf",
                   ext = 2, dimSize = 50
                 )
               ))
  set.seed(3245)
  lc <- rweibull(n.locs(data),2,3)
  data2$asdf <- lc
  expect_message(
    a <-
      brownian.bridge.dyn(
        data, location.error = lc, ext = 2, dimSize = 50
      ), 'Computational size'
  )
  a@DBMvar$asdf <- lc
  expect_message(
    xx <- brownian.bridge.dyn(
      data2,
      location.error = "asdf",
      ext = 2, dimSize =
        50
    ), 'Computational size'
  )
  
  expect_equal(a, xx)
})

test_that('Brownian bridge, running with character and vector input to location error for stacks',{
  data(leroy)
  data(ricky)
  data <-
    spTransform(moveStack(list(leroy[1:75,], ricky[1:73,])), center = T)
  data2 <- data
  lc <- 6.54
  data2$asdf <- lc
  expect_message(
    a <-
      brownian.bridge.dyn(
        data, location.error = lc, ext = 2, dimSize = 50
      ), 'Computational size'
  )
  expect_true(validObject(a))
  
  a@DBMvar$asdf <- lc
  expect_equal(a,suppressMessages(
    brownian.bridge.dyn(
      data2,
      location.error = "asdf",
      ext = 2, dimSize = 50
    )
  ))
  set.seed(3245)
  lc <- rweibull(sum(n.locs(data)),2,3)
  data2$asdf <- lc
  expect_message(
    a <-
      brownian.bridge.dyn(
        data, location.error = lc, ext = 2, dimSize = 50
      ), 'Computational size'
  )
  a@DBMvar$asdf <- lc
  expect_equal(a,
               suppressMessages(
                 brownian.bridge.dyn(
                   data2,
                   location.error = "asdf",
                   ext = 2, dimSize = 50
                 )
               ))
  expect_error(
    suppressMessages(
      brownian.bridge.dyn(
        data,
        location.error = lc[-1], ext = 2, dimSize = 50
      )
    ),
    'Location error needs to be the same length as the number of locations'
  )
})
