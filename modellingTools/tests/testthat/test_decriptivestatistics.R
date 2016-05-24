require(modellingTools, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(foreach, quietly = TRUE, warn.conflicts = FALSE)


context("Descriptive Statistics")

### proc_freq

x <- data_frame(v1 = c(rep(1,5),rep(2,5)),
                v2 = 1:10,
                v3 = letters[1:10],
                v4 = factor(letters[1:10])
)

pf1 <- proc_freq(x,"v1")
pf2 <- proc_freq(x,"v2")
pf3 <- proc_freq(x,"v3")
pf4 <- proc_freq(x,"v4")


test_that("proc_freq gives as many levels as there are unique variable values", {
  expect_equal(length(column_vector(pf1,"level")),
               length(unique(column_vector(x,"v1")))
  )
  expect_equal(length(column_vector(pf2,"level")),
               length(unique(column_vector(x,"v2")))
  )
  expect_equal(length(column_vector(pf3,"level")),
               length(unique(column_vector(x,"v3")))
  )
  expect_equal(length(column_vector(pf4,"level")),
               length(unique(column_vector(x,"v4")))
  )

})

test_that("proc_freq gives the class of levels the same as the class of the variable", {
  expect_equal(class(column_vector(pf1,"level")),
               class(column_vector(x,"v1"))
  )
  expect_equal(class(column_vector(pf2,"level")),
               class(column_vector(x,"v2"))
  )
  expect_equal(class(column_vector(pf3,"level")),
               class(column_vector(x,"v3"))
  )
  expect_equal(class(column_vector(pf4,"level")),
               class(column_vector(x,"v4"))
  )

})

pfb1 <- proc_freq(x,"v1",bins = 4)
pfb2 <- proc_freq(x,"v2",bins = 4)
pfb3 <- proc_freq(x,"v3",bins = 4)
pfb4 <- proc_freq(x,"v4",bins = 4)


test_that("proc freq gives correct number of bins", {
  expect_equal(length(column_vector(pfb1,"level")),2)
  expect_equal(length(column_vector(pfb2,"level")),4)
  expect_equal(length(column_vector(pfb3,"level")),10)
  expect_equal(length(column_vector(pfb4,"level")),10)
})

test_that("count sums to the sample size", {
  expect_equal(sum(column_vector(pf1,"count")),10)
  expect_equal(sum(column_vector(pf2,"count")),10)
  expect_equal(sum(column_vector(pf3,"count")),10)
  expect_equal(sum(column_vector(pf4,"count")),10)

  expect_equal(sum(column_vector(pfb1,"count")),10)
  expect_equal(sum(column_vector(pfb2,"count")),10)
  expect_equal(sum(column_vector(pfb3,"count")),10)
  expect_equal(sum(column_vector(pfb4,"count")),10)

})

### get_top_corrs

# Change datatypes and add NAs to mtcars
mt_dat <- mtcars %>%
            tbl_df() %>%
            mutate(vs = factor(vs),
                   am = factor(am),
                   gear = as.character(gear))

times(nrow(mt_dat)) %do%
  if (runif(1) < .1) {
    mt_dat[sample(1:nrow(mt_dat),1),"hp"] <- NA
  }

mtcor <- get_top_corrs(mt_dat,"mpg")

test_that("get_top_corrs considers only and all numeric variables", {
  expect_true(all(column_vector(mtcor,"var_name") %in%
                    c("wt","cyl","disp","hp","drat","carb","qsec"))
  )
  expect_true(all(c("wt","cyl","disp","hp","drat","carb","qsec") %in%
                    column_vector(mtcor,"var_name"))
  )

})

test_that("correlations are in proper range and not missing", {
  expect_false(any(is.na(column_vector(mtcor,"correlation"))))
  expect_true(all(column_vector(mtcor,"correlation") > -1 &&
                  column_vector(mtcor,"correlation") < 1))
})















