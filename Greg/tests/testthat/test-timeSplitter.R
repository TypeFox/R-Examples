library(testthat)

context("Test timeSplitter")
library(magrittr)
library(dplyr)

test_that("Check the number of events and rows are correct",{
  test_data <- data.frame(
    id = 1:4,
    time = c(4, 3.5, 1, 5),
    event = c("censored", "dead", "alive", "dead"),
    age = c(62.2, 55.3, 73.7, 46.3),
    date = as.Date(
      c("2003-01-01", 
        "2010-04-01", 
        "2013-09-20",
        "2002-02-23"))
  )
  
  split_data <- 
    test_data %>% 
    select(id, event, time, age, date) %>% 
    timeSplitter(by = 2, # The time that we want to split by
                 event_var = "event",
                 time_var = "time",
                 event_start_status = "alive",
                 time_related_vars = c("age", "date"))  
  
  expect_equal(nrow(split_data), sum(ceil(test_data$time / 2)))
  expect_equal(sum(split_data$event == "dead"), 
               sum(test_data$event == "dead"))
  expect_equal(sum(split_data$event == "censored"), 
               sum(test_data$event == "censored"))
  expect_more_than(sum(split_data$event == "alive"), 
                   sum(test_data$event == "alive"))
})


test_that("Check that labels are preserved",{
  test_data <- data.frame(
    id = 1:4,
    time = c(4, 3.5, 1, 5),
    event = c("censored", "dead", "alive", "dead"),
    age = c(62.2, 55.3, 73.7, 46.3),
    date = as.Date(
      c("2003-01-01", 
        "2010-04-01", 
        "2013-09-20",
        "2002-02-23"))
  )
  
  library(Hmisc)
  label(test_data$age) <- "Age (years)"
  split_data <- 
    test_data %>% 
    select(id, event, time, age, date) %>% 
    timeSplitter(by = 2, # The time that we want to split by
                 event_var = "event",
                 time_var = "time",
                 event_start_status = "alive",
                 time_related_vars = c("age", "date"))  
  
  expect_equivalent(label(split_data$age), 
                    label(test_data$age))
})


test_that("Check that age and calendar time is updated",{
  test_data <- data.frame(
    id = 1:4,
    time = c(4, 3.5, 1, 5),
    event = c("censored", "dead", "alive", "dead"),
    age = c(62.2, 55.3, 73.7, 46.3),
    date = as.Date(
      c("2003-01-01", 
        "2010-04-01", 
        "2013-09-20",
        "2002-02-23"))
  )
  
  library(Hmisc)
  label(test_data$age) <- "Age (years)"
  split_data <- 
    test_data %>% 
    select(id, event, time, age, date) %>% 
    timeSplitter(by = 2, # The time that we want to split by
                 event_var = "event",
                 time_var = "time",
                 event_start_status = "alive",
                 time_related_vars = c("age", "date"))  
  test_data <- cal.yr(test_data)
  library(dplyr)
  for (i in 1:nrow(test_data)){
    row <- test_data[i,]
    last_age <- row$age + (ceil(row$time/2) - 1)*2
    age <- split_data[split_data$id == row$id, "age"]
    expect_true(max(age) - last_age < .Machine$double.eps,
                info = paste0("Could not identify the split max age (", max(age), ")",
                              " to be equal to the calculated (", last_age, ")",
                              " for row no ", i))
    last_cal <- row$date + (ceil(row$time/2) - 1)*2
    cal <- split_data[split_data$id == row$id, "date"]
    expect_true(max(age) - last_age < .Machine$double.eps,
                info = paste0("Could not identify the split max date (", max(cal), ")",
                              " to be equal to the calculated (", last_cal, ")",
                              " for row no ", i))
  }
})

  