# Assessment data

context("Test that lags works")

test_that("Object fails correctly", {
  test_data <- expand.grid(id = sample(letters, 10), 
                           time = 1:10)
  vals <- data.frame(value1 = rnorm(100), value2 = runif(100),
                     value3 = rpois(100, 4))
  test_data <- cbind(test_data, vals); rm(vals)
  
  expect_error(lag_data(test_data, group="id", time="time", 
                     values=c("value1", "value2"), periods = "A"), 
               "Periods must be a numeric or integer value")
  expect_error(lag_data(test_data, group=id, time="time", 
                        values=c("value1", "value2"), periods = 2))
  expect_error(lag_data(test_data, group="id", time="time", 
                        values=c("value1", "id"), periods = 2), 
               "Cannot lag a grouping or a time variable")
  expect_warning(lag_data(test_data, group="id", time="time", 
                          values=c("value1", "value2"), periods = 2.16),
                          "parameter periods has been forced to integer of value 2")
})

test_that("Object values are correct", {
  test_data <- expand.grid(id = sample(letters, 10), 
                           time = 1:10)
  vals <- data.frame(value1 = rnorm(100), value2 = runif(100),
                     value3 = rpois(100, 4))
  test_data <- cbind(test_data, vals); rm(vals)
  new_data <- lag_data(test_data, group="id", time="time", 
           values=c("value1", "value2"), periods = 3)
  
  for(i in unique(test_data$id)){
    for(j in unique(test_data$time)){
      if(j + 3 > max(test_data$time)){
        
      } else{
        expect_equal(test_data[test_data$id == i & test_data$time == j, ]$value1, 
                     new_data[new_data$id == i & new_data$time == j+3, ]$value1.lag3)
        expect_equal(test_data[test_data$id == i & test_data$time == j, ]$value2, 
                     new_data[new_data$id == i & new_data$time == j+3, ]$value2.lag3)
      }
    }
  }
  
  new_data2 <- lag_data(test_data, group="id", time="time", 
           values=c("value1", "value2"), periods = -2)
  
  for(i in unique(test_data$id)){
    for(j in unique(test_data$time)){
      if(j -2 < min(test_data$time)){
        
      } else{
        expect_equal(test_data[test_data$id == i & test_data$time == j, ]$value1, 
                     new_data2[new_data$id == i & new_data2$time == j-2, ]$value1.lead2)
        expect_equal(test_data[test_data$id == i & test_data$time == j, ]$value2, 
                     new_data2[new_data$id == i & new_data2$time == j-2, ]$value2.lead2)
      }
    }
  }  
})

context("Test crossplot functionality")

test_that("Crossplots work and get right answers", {
  set.seed(1213)
  sampDat <- data.frame(cbind(x=seq(1,3,by=1), y=sample(LETTERS[6:8], 60, 
                                                        replace=TRUE)),
                        fac=sample(LETTERS[1:4], 60, replace=TRUE))
  varnames<-c('Quality','Grade')
  
  out <- crosstabs(sampDat, "y", "fac", varnames = varnames)
  expect_is(out, "list")
  expect_equal(names(out), c("TABS", "PROPORTIONS", "TABSPROPORTIONS"))
  
  
})

context("Test mosaictabs with labels")

test_that("Plots work correctly", {
  sampDat <- data.frame(cbind(x=seq(1,3,by=1), y=sample(LETTERS[6:8], 60, 
                                                        replace=TRUE)),
                        fac=sample(LETTERS[1:4], 60, replace=TRUE))
  varnames<-c('Quality','Grade')
  test_plot_file <- "crosstabplot.png"
  png(test_plot_file)
  crosstabplot(sampDat, "y", "fac", varnames = varnames)
  dev.off()
  expect_true(file.exists(test_plot_file))
  unlink(test_plot_file)
  
  test_plot_file <- "crosstabplot.png"
  png(test_plot_file)
  crosstabplot(sampDat, "y", "fac", varnames = varnames, shade = TRUE)
  dev.off()
  expect_true(file.exists(test_plot_file))
  unlink(test_plot_file)
  
  
  test_plot_file <- "crosstabplot.png"
  png(test_plot_file)
  crosstabplot(sampDat, "y", "fac", varnames = varnames, shade = FALSE)
  dev.off()
  expect_true(file.exists(test_plot_file))
  unlink(test_plot_file)
  
  test_plot_file <- "crosstabplot.png"
  png(test_plot_file)
  crosstabplot(sampDat, "y", "fac", varnames = varnames, label = TRUE)
  dev.off()
  expect_true(file.exists(test_plot_file))
  unlink(test_plot_file)
  
  
})

context("Test proficiency polygons")

test_that("profpoly.data produces correct objects", {
  
  grades<-c(3,4,5,6,7,8)
  g <- length(grades)
  LOSS <- rep(200, g)
  HOSS <- rep(650, g)
  basic <- c(320,350,370,390,420,440)
  minimal <- basic-30
  prof <- c(380,410,430,450,480,500)
  adv <- c(480,510,530,550,580,600)
  
  z <- profpoly.data(grades, LOSS, minimal, basic, proficient = prof, 
                     advanced = adv, HOSS)
  expect_is(z, "data.frame")
  
})

test_that("profpoly makes valid ggplot objects", {
  grades<-c(3,4,5,6,7,8)
  g <- length(grades)
  LOSS <- rep(200, g)
  HOSS <- rep(650, g)
  basic <- c(320,350,370,390,420,440)
  minimal <- basic-30
  prof <- c(380,410,430,450,480,500)
  adv <- c(480,510,530,550,580,600)
  
  z <- profpoly.data(grades, LOSS, minimal, basic, proficient = prof, 
                     advanced = adv, HOSS)
  p1 <- profpoly(z)
  expect_is(p1, c("gg", "ggplot"))
  
#   test_plot_file <- "profpoly.png"
#   png(test_plot_file)
#   profpoly(z)
#   dev.off()
#   expect_true(file.exists(test_plot_file))
#   unlink(test_plot_file)
})

