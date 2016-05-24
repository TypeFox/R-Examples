context("Compare splitMulti results with splitLexis results")

test_that("splitMulti and splitLexis are congruent", {
  skip_on_cran()
  skip_on_travis()
  library(Epi)
  
  sire2 <- copy(sire)
  sire2[, dg_yrs := get.yrs(dg_date, "actual")]
  sire2[, ex_yrs := get.yrs(ex_date, "actual")]
  sire2[, bi_yrs := get.yrs(bi_date, "actual")]
  sire2[, id := 1:.N]
  
  BL1 <- list(fot = 0:5, per = 1990:1995, age = c(0, 60, Inf))
  
  BL2 <- list(fot = c(10,Inf), age = seq(0,150,5))
  
  BL3 <- list(fot = c(5, Inf), per = c(1900, 2100), age = c(25,100))
  
  BL4 <- list(fot = 0:10)
  
  BL5 <- list(fot = 5:10)
  
  BL6 <- list(per = 1990:2000, age = c(50,70))
  
  
  BL <- list(BL1, BL2, BL3, BL4, BL5, BL6)
  
  x <- Lexis(data=sire2[dg_date < ex_date], entry=list(fot=0, per=dg_yrs, age=dg_age),
             exit=list(per=ex_yrs), merge=TRUE, exit.status=1L, entry.status = 0L)
  setDT(x)
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  
  #   x2 <- splitMulti(x, breaks = BL[[6]], drop = F)
  #   x3 <- splitMultiEpi(x, breaks = BL[[6]], drop = F)
  #   
  #   x2d <- splitMulti(x, breaks = BL[[6]], drop = T)
  #   x3d <- splitMultiEpi(x, breaks = BL[[6]], drop = T)
  #   x2d <- intelliDrop(x2, breaks = BL[[6]])
  #   x3d <- intelliDrop(x3, breaks = BL[[6]])
  #   compareSMWithEpi(x, BL[[6]], drop = F)
  
  
  # one row per id ---------------------------------------------------------------
  
  test_that("splitMulti and splitLexis congruent with one row per id", {
    for (sc in seq_along(BL)) {
      compareSMWithEpi(x, BL[[sc]])
    }
  })
  
  
  
  # multiple rows per id ---------------------------------------------------------
  
  sire2 <- sire2[rep(1:.N, each = 2)]
  
  x <- Lexis(data=sire2[dg_date < ex_date], entry=list(fot=0, per=dg_yrs, age=dg_age),
             exit=list(per=ex_yrs), merge=TRUE, exit.status=1L, entry.status = 0L, id = id)
  setDT(x)
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  
  for (sc in seq_along(BL)) {
    test_that(paste0("splitLexisDT and splitLexis congruent with multiple rows per id with breaks no. ", sc), {
      compareSMWithEpi(x, BL[[sc]])
    })
  }
  
  # multistate using Lexis -----------------------------------------------------
  
    
  sire2[, EX := suppressWarnings(factor(status, levels = 0:2, labels = c("ok", "dead", "dead"), ordered = T))]    
  sire2[, EN := suppressWarnings(factor(0L, levels = 0:2, labels = c("ok", "dead", "dead"), ordered = T))]
  
  x <- Lexis(data=sire2[dg_date < ex_date & !duplicated(id)], entry=list(fot=0, per=bi_yrs, age=0),
             exit=list(per=ex_yrs), merge=TRUE, exit.status=EX, entry.status = EN, id = id)
  
  x <- cutLexis(x, cut = x$dg_yrs, timescale = "per", new.state = "sick", precursor.state = "ok")
  setDT(x)   
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  
  BL[[1L]] <- NULL ## this would drop all rows in split data
  
  for (sc in seq_along(BL)) {
    test_that(paste0("splitLexisDT and splitLexis congruent with multiple Lexis states per id using breaks list no. ", sc), {
      compareSMWithEpi(x, BL[[sc]])
    })
  }
  
  # multistate using mstate ----------------------------------------------------
  
#   test_that("splitMulti works with data produced by mstate", {
#     
#   }
  
})



