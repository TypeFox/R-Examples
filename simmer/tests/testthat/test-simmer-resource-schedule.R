context("resource-schedule")

test_that("capacity & queue size change", {
  inf_sch <- schedule(c(8, 16, 24), c(1, 2, 3), Inf)
  fin_sch <- schedule(c(8, 16, 24), c(1, 2, 3), 24)
  
  limits <- simmer() %>%
    add_resource("dummy", inf_sch) %>%
    run(16) %>% reset() %>% run(48) %>%
    get_mon_resources("limits")
  
  expect_equal(limits$time, c(8, 16, 24))
  expect_equal(limits$server, c(1, 2, 3))
  
  limits <- simmer() %>%
    add_resource("dummy", fin_sch) %>%
    run(16) %>% reset() %>% run(48) %>%
    get_mon_resources("limits")
  
  expect_equal(limits$time, c(8, 16, 24, 32, 40, 48))
  expect_equal(limits$server, c(1, 2, 3, 1, 2, 3))
})

test_that("queue size changes", {
  inf_sch <- schedule(c(8, 16, 24), c(1, 2, 3), Inf)
  fin_sch <- schedule(c(8, 16, 24), c(1, 2, 3), 24)
  
  limits <- simmer() %>%
    add_resource("dummy", 1, inf_sch) %>%
    run(16) %>% reset() %>% run(48) %>%
    get_mon_resources("limits")
  
  expect_equal(limits$time, c(8, 16, 24))
  expect_equal(limits$queue, c(1, 2, 3))
  
  limits <- simmer() %>%
    add_resource("dummy", 1, fin_sch) %>%
    run(16) %>% reset() %>% run(48) %>%
    get_mon_resources("limits")
  
  expect_equal(limits$time, c(8, 16, 24, 32, 40, 48))
  expect_equal(limits$queue, c(1, 2, 3, 1, 2, 3))
})

test_that("arrivals 1) are dequeued when resource's capacity increases and 2) remain in server when it decreases", {
  t <- create_trajectory() %>% 
    seize("dummy", 1) %>%
    timeout(2) %>%
    release("dummy", 1)
  
  inf_sch <- schedule(c(0, 1, 2), c(1, 3, 1), Inf)
  
  arrivals <- simmer() %>%
    add_resource("dummy", inf_sch) %>%
    add_generator("asdf", t, at(0, 0, 0)) %>%
    run() %>%
    get_mon_arrivals()
  
  expect_equal(arrivals$end_time, c(2, 3, 3))
  expect_equal(arrivals$activity_time, c(2, 2, 2))
})

test_that("arrivals 1) are dequeued when resource's capacity increases and 2) remain in server when it decreases", {
  t <- create_trajectory() %>% 
    seize("dummy", 1) %>%
    timeout(2) %>%
    release("dummy", 1)
  
  inf_sch <- schedule(c(0, 1, 2), c(1, 3, 1), Inf)
  
  arrivals <- simmer() %>%
    add_resource("dummy", inf_sch) %>%
    add_generator("asdf", t, at(0, 0, 0)) %>%
    run() %>%
    get_mon_arrivals()
  
  expect_equal(arrivals$end_time, c(2, 3, 3))
  expect_equal(arrivals$activity_time, c(2, 2, 2))
})

test_that("arrivals are preempted when resource's capacity decreases", {
  t <- create_trajectory() %>% 
    seize("dummy", 1, restart=TRUE) %>%
    timeout(2) %>%
    release("dummy", 1)
  
  inf_sch <- schedule(c(0, 1, 2), c(1, 3, 1), Inf)
  
  arrivals <- simmer() %>%
    add_resource("dummy", inf_sch, preemptive=TRUE) %>%
    add_generator("asdf", t, at(0, 0, 0)) %>%
    run() %>%
    get_mon_arrivals()
  
  expect_equal(arrivals$end_time, c(2, 3, 5))
  expect_equal(arrivals$activity_time, c(2, 2, 3))
})
