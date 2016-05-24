context("Sequent Peak Algorithm")
# they handle leap years quite strange...
sp.storage <- read.csv2("tallaksen-sequent-peak-storage.csv")[425:1883, ]
sp.storage$time <- as.Date(sp.storage$time, format = "%d.%m.%y")

sp.summary <- read.csv2("tallaksen-sequent-peak-summary.csv")
sp.summary$start <- as.Date(sp.summary$start, format = "%d.%m.%Y")
sp.summary <- sp.summary[sp.summary$start >=min(sp.storage$time) &
                         sp.summary$start <= max(sp.storage$time), ]

ng <- xts(x = data.frame(discharge = sp.storage$streamflow),
          order.by = sp.storage$time)
xtsAttributes(ng)[["unit"]] <- "m^3/s"
ng <- .check_xts(ng)

deficit <- pool_sp(find_droughts(ng, threshold = 5.18))


test_that("internal storage is computed correctly", {
  # vector of base flows is as long as input
  expect_equal(nrow(deficit), nrow(sp.storage))

  deficit$storage <- 0
  for(i in setdiff(unique(deficit$event.no), 0)) {
    rng <- deficit$event.no == i
    deficit$storage[rng] <- cumsum(deficit$def.increase[rng])
  }


  # values according to Tallaksen and van Lanen.
  # fails, because Tallaksen and van Lanen exclude Feb 28th
  expect_equal2(as.vector(deficit$storage), sp.storage$storage * 86400,
                tolerance = 1e-2,
                label = "Deficit volumes given in Tallakesen is equal to computed deficit volumes")
})


test_that("deficit volume is computed correctly", {

  # values according to Tallaksen and van Lanen.
  expect_equal(summary(deficit, drop_minor = 0)$volume,
               sp.summary$def.volume * 86400, tolerance = 1e-3)
})

test_that("deficit duration is computed correctly", {

  # values according to Tallaksen and van Lanen.
  expect_equal(summary(deficit, drop_minor = 0)$duration,
               sp.summary$duration)
})
