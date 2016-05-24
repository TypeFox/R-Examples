context("Testing aggregation by categories of exposure")

test_that("prepExpo works in the simple case", {
  skip_on_cran()
  library(Epi)
  
  df <- data.frame(id = "A", birth  = c(1952.4534), 
                   entry = c(1965.4746, 1972.42845, 1991.78643),
                   exit = c(1968.56346, 1979.32478, 1997.32432), fail = 0)
  
  # Define as Lexis object with timescales calendar time and age
  x <- Lexis( entry = list(work = 0, per=entry ),
                 exit = list( per=exit, age=exit-birth ),
                 exit.status = fail, id = id,
                 data = df )
  
  x2 <- prepExpo(x, freezeScales = "work", 
                 cutScale = "per", 
                 entry = 1964, 
                 exit = 2012, by = "lex.id")
  cd <- cumsum(x$lex.dur)
  exp_work <-  c(0, 0, cd[1], cd[1], cd[2], cd[2], cd[3])
  expect_equal(x2$work, exp_work)
  expect_equal(x2$per, 1964+c(0,cumsum(x2$lex.dur)[-nrow(x2)]))
  
  BL <- list(work = 0:50, age = c(0,18,Inf), per = 1963:2014)
  x2 <- prepExpo(x, freezeScales = "work", 
                 cutScale = "per", 
                 entry = 1964, 
                 # verbose = TRUE,
                 exit = 2012, by = "lex.id", 
                 breaks = BL)
  ag <- aggre(x2, by = list(lex.id, per, age))
  
  xx <- Lexis(entry = list(per = 1964, age = 1964-birth), exit = list(per=2012), data = df[1,])
  ag2 <- splitMulti(xx, breaks = BL[c("per","age")])
  ag2 <- aggre(ag2, by = list(lex.id, per, age))
  
  setkeyv(ag, c("lex.id","per"))
  setkeyv(ag2, c("lex.id","per"))
  
  expect_equal(ag$pyrs, ag2$pyrs)
  
  
  
  
})
