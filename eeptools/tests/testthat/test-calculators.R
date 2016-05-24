# age, retention, and moves calculations

context("Test age calculator")

test_that("Leap year calculations work", {
  # from @larmarange
  expect_equal(age_calc(as.Date('2004-01-15'), as.Date('2004-02-16')), 1.034483, 
               tol = .00001)
  expect_equal(age_calc(as.Date('2005-01-15'), as.Date('2005-02-16')), 1.035714, 
               tol = .00001)
  expect_equal(age_calc(as.Date('1995-01-15'), as.Date('2003-02-16')), 
               age_calc(as.Date('1994-01-15'), as.Date('2002-02-16')))
  expect_false(age_calc(as.Date('1996-01-15'), as.Date('2004-02-16')) ==
                 age_calc(as.Date('1994-01-15'), as.Date('2002-02-16')))
})

test_that("All function parameters result in a numeric calculations with sane inputs", {
  tests <- expand.grid(precise = c(TRUE, FALSE), 
                       units = c("days", "months", "years"), 
                       dob = c("atomic", "vector"), 
                       enddate = c("atomic", "vector"))
  
  safe.ifelse <- function(cond, yes, no) structure(ifelse(cond, yes, no), class = class(yes))
  
  for(i in 1:nrow(tests)){
    atomDOB <- as.Date(as.POSIXct('1987-05-29 018:07:00'))
    vecDOB <- as.Date(seq(as.POSIXct('1987-05-29 018:07:00'), len=26, by="21 day"))
    vecED <- as.Date(seq(as.POSIXct('2017-05-29 018:07:00'), len=26, by="21 day"))
    atomED <- as.Date(as.POSIXct('2017-05-29 018:07:00'))
    
    dob <- safe.ifelse(tests[i, "dob"] == "atomic", atomDOB, vecDOB)
    enddate <- safe.ifelse(tests[i, "enddate"] == "atomic", atomED, vecED)
    
    out <- age_calc(dob = dob, enddate = enddate, units = tests[i, ]$units, 
                    precise = tests[i, ]$precise)
    expect_true(class(out) %in% c("difftime", "numeric"))
  }

})

test_that("Bad inputs yield correct errors", {
  expect_error(age_calc('2004-01-15', '2004-02-16'), 
               "Both dob and enddate must be Date class objects")
  expect_error(age_calc(as.Date('2004-01-15'), '2004-02-16'), 
               "Both dob and enddate must be Date class objects")
  expect_error(age_calc('2004-01-15', as.Date('2004-02-16')), 
               "Both dob and enddate must be Date class objects")
  expect_error(age_calc(as.Date('2004-02-16'), as.Date('2004-01-15')), 
               "End date must be a date after date of birth")
  
})


context("Test retention calculator")

test_that("standard cases work", {
  x <- data.frame(sid = c(101, 101, 102, 103, 103, 103, 104, 105, 105, 106, 106),
                 grade = c(9, 10, 9, 9, 9, 10, 10, 8, 9, 7, 7))
  expect_is(retained_calc(x), "data.frame")
  expect_equal(nrow(retained_calc(x)), 4)
  z <- data.frame(stuid = c(101, 101, 102, 103, 103, 103, 104, 105, 105, 106, 106),
                  grade_cd = c(9, 10, 9, 9, 9, 10, 10, 8, 9, 7, 7))
  expect_is(retained_calc(z, sid = "stuid", grade = "grade_cd"), "data.frame")
  expect_identical(retained_calc(z, sid = "stuid", grade = "grade_cd"), 
                   retained_calc(x))
  tests <- data.frame(grade_val = 1:12, expected_val = NA)
  
  test_dat <- data.frame(stuid = rep(101:130, each = 12), 
                         grade = rep(seq(1:12), 30))
  test_dat <- test_dat[order(test_dat$stuid, test_dat$grade),]
  test_dat$stuid <- as.character(test_dat$stuid)
  
  test_dat$grade[test_dat$stuid == "120"] <- c(1, 1, 2:11) 
  test_dat$grade[test_dat$stuid == "121"] <- c(1, 2, 2, 3:11) 
  test_dat$grade[test_dat$stuid == "122"] <- c(1:3, 3, 4:11) 
  test_dat$grade[test_dat$stuid == "123"] <- c(1:4, 4, 5:11) 
  test_dat$grade[test_dat$stuid == "124"] <- c(1:5, 5, 6:11) 
  test_dat$grade[test_dat$stuid == "125"] <- c(1:6, 6, 7:11) 
  test_dat$grade[test_dat$stuid == "126"] <- c(1:7, 7, 8:11) 
  test_dat$grade[test_dat$stuid == "127"] <- c(1:8, 8, 9:11) 
  test_dat$grade[test_dat$stuid == "128"] <- c(1:9, 9, 10:11) 
  test_dat$grade[test_dat$stuid == "129"] <- c(1:10, 10, 11) 
  test_dat$grade[test_dat$stuid == "130"] <- c(1:11, 11) 
  test_dat$grade[test_dat$stuid == "102"] <- c(1:4, 4, 3, 6:11) 
  test_dat$grade[test_dat$stuid == "103"] <- c(1:8, 6, 10, 11, 12)
  test_dat$grade[test_dat$stuid == "104"] <- c(1:11, 9)
  test_dat$grade[test_dat$stuid == "105"] <- c(1:5, 3, 4, 5:9)
  
  tests$expected_val[tests$grade_val == 1] <- list(c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "Y", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 2] <- list(c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "N", "Y", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 3] <- list(c("N", "Y", "N", "N", "Y", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "N", "N", "Y", "N", "N", "N", 
                                                "N", "N", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 4] <- list(c("N", "Y", "N", "N", "Y", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Y", "N", "N", 
                                                "N", "N", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 5] <- list(c("N", "N", "N", "Y", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Y", "N", "N", 
                                                "N", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 6] <- list(c("N", "N", "Y", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                               "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Y", 
                                               "N", "N", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 7] <- list(c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "Y", "N", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 8] <- list(c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "Y", "N", "N", "N"))
  tests$expected_val[tests$grade_val == 9] <- list(c("N", "N", "Y", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                "N", "Y", "N", "N"))
  tests$expected_val[tests$grade_val == 10] <- list(c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                 "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                 "N", "N", "Y", "N"))
  tests$expected_val[tests$grade_val == 11] <- list(c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                 "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                 "N", "N", "N", "Y"))
  tests$expected_val[tests$grade_val == 12] <- list(c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
                                                 "N", "N", "N", "N"))
  
  for(i in tests$grade_val){
    expect_identical(retained_calc(test_dat, sid = "stuid", grade = "grade", grade_val = i)$retained, 
              unlist(tests$expected_val[tests$grade_val == i]))
    
  }
})

test_that("Nonstandard cases fail", {
  x <- data.frame(sid = c(101, 101, 102, 103, 103, 103, 104, 105, 105, 106, 106),
                  grade = c(9, 10, 9, 9, 9, 10, 10, 8, 9, 7, 7))
  expect_is(retained_calc(x, grade_val = 13), "data.frame")
  expect_equal(nrow(retained_calc(x, grade_val = 13)), 0)
  expect_is(retained_calc(x, grade_val = -2), "data.frame")
  expect_equal(nrow(retained_calc(x, grade_val = -2)), 0)
})


context("Test moves calculator")

test_that("Zed", {
  df <- data.frame(sid = c(rep(1,3), rep(2,4), 3, rep(4,2)),
                   schid = c(1, 2, 2, 2, 3, 1, 1, 1, 3, 1),
                   enroll_date = as.Date(c('2004-08-26',
                                           '2004-10-01',
                                           '2005-05-01',
                                           '2004-09-01',
                                           '2004-11-03',
                                           '2005-01-11',
                                           '2005-04-02',
                                           '2004-09-26',
                                           '2004-09-01',
                                           '2005-02-02'),
                                         format='%Y-%m-%d'),
                   exit_date = as.Date(c('2004-08-26',
                                         '2005-04-10',
                                         '2005-06-15',
                                         '2004-11-02',
                                         '2005-01-10',
                                         '2005-03-01',
                                         '2005-06-15',
                                         '2005-05-30',
                                         NA,
                                         '2005-06-15'),
                                       format='%Y-%m-%d'))
  moves <- moves_calc(df)
  moves
  expect_is(moves, "data.frame")
  expect_equal(nrow(moves), 4)
  # moves <- moves_calc(df, enrollby='2004-10-15', gap=22)
  # moves
  # moves <- moves_calc(df, enrollby='2004-10-15', exitby='2005-05-29')
  # moves
})

