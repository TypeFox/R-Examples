context("utilities")

test_that("utilities loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))

    library(codyn)

    # Load our example data set
    data(knz_001d)
    expect_that(names(knz_001d)[4], equals("abundance"))

    #take a subset
    dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")

    # Does transpose_community correctly transform a long dataframe into a matrix?
    com_wide_df <- transpose_community(df = dat1, "year", "species", "abundance")
    expect_true(is.data.frame(com_wide_df))
    
    df2 = subset(knz_001d, subplot == "A_1")
    names(df2)=c("spp", "yr", "subplot", "abund")
    com_wide_df <- transpose_community(df = df2, "yr", "spp", "abund")
    expect_true(is.data.frame(com_wide_df))
    # rownames should be timevar
    expect_true(all(as.numeric(rownames(com_wide_df)) %in% dat1[["year"]]))
    # colnames should be species
    expect_true(all(colnames(com_wide_df) %in% dat1[["species"]]))
})

test_that("Name checking works", {
  data(knz_001d)
  #knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"),
  #                     sep=",", header=TRUE)
  expect_error(check_names(given = c("AAA", "year", "subplot", "abundance"),
                           data = knz_001d),
                 "data does not have name .*")
  expect_error(check_names(given = c("species", "BBB", "subplot", "abundance"),
                           data = knz_001d),
               "data does not have name .*")
  expect_error(check_names(given = c("species", "year", "CCC", "abundance"),
                           data = knz_001d),
               "data does not have name .*")
  expect_error(check_names(given = c("species", "year", "subplot", "DDD"),
                           data = knz_001d),
               "data does not have name .*")
})

test_that("Single-record checks for species work", {
  data(knz_001d)
 
  expect_null(check_single(knz_001d, time.var = "year", species.var = "species", replicate.var = "subplot"))

  knz_001d2<-knz_001d
  names(knz_001d2)=c("spp", "yr", "subp", "abund")
  expect_null(check_single(knz_001d2, time.var = "yr", species.var = "spp", replicate.var = "subp"))
  
  dftog = rbind(knz_001d, knz_001d[nrow(knz_001d),])
  expect_error(check_single(dftog, time.var = "year", species.var = "species", replicate.var = "subplot"))
  
  dftog2 = rbind(knz_001d2, knz_001d2[nrow(knz_001d2),])
  expect_error(check_single(dftog2, time.var = "yr", species.var = "spp", replicate.var = "subp"))
  
  
  dfa = subset(knz_001d, subplot == "A_1")
  expect_null(check_single_onerep(dfa, time.var = "year", species.var = "species"))
  dfab = subset(knz_001d, subplot == "A_1" | subplot == "A_2")
  expect_warning(check_single_onerep(dfab, time.var = "year", species.var = "species"))
  
  df2 = subset(knz_001d, subplot == "A_1")
  names(df2)=c("spp", "yr", "subplot", "abund")
  expect_null(check_single_onerep(df2, time.var = "yr", species.var = "spp"))
  
  #check check_sppvar
  comdf <- as.data.frame(matrix(1, 5, 5))
  expect_error(check_sppvar(comdf))
  
  comdf2 <-as.data.frame(matrix(c(1, 2, 3, 4, 5), 5, 5))
  check_sppvar(comdf2)
  
  })
