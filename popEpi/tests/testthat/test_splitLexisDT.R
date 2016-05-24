context("Compare splitLexisDT results with splitLexis results")

test_that("splitLexisDT and splitLexis are congruent", {
  skip_on_cran()
  skip_on_travis()
  library(Epi)
  
  sire2 <- copy(sire)
  sire2[, dg_yrs := get.yrs(dg_date, "actual")]
  sire2[, ex_yrs := get.yrs(ex_date, "actual")]
  sire2[, bi_yrs := get.yrs(bi_date, "actual")]
  sire2[, id := 1:.N]
  
  BL <- list(fot = 0:5, fot = 2:6, fot = c(10,Inf), fot = c(0, Inf),
             per = 1990:1995, per = seq(1950,2010,10), per = c(1900, 2100), 
             age = seq(0,150,5), age = c(25,100), age = c(0, 60))
  
  x <- Lexis(data=sire2[dg_date < ex_date], entry=list(fot=0, per=dg_yrs, age=dg_age),
             exit=list(per=ex_yrs), merge=TRUE, exit.status=1L, entry.status = 0L)
  setDT(x)
  forceLexisDT(x, breaks = list(fot=NULL,per=NULL,age=NULL),
               allScales = c("fot", "per", "age"))
  
  
  # x2 <- splitLexis(x, breaks = BL[[3]], time.scale = "fot")
  # x3 <- splitLexisDT(x, breaks = BL[[3]], timeScale = "fot", drop = FALSE)
  # x2 <- intelliDrop(setDT(x2),  breaks = list(fot = BL[[3]]))
  # x3 <- intelliDrop(x3,  breaks = list(fot = BL[[3]]))
  # x4 <- splitLexisDT(x, breaks = BL[[3]], timeScale = "fot", drop = TRUE)
  
  
  # one row per id ---------------------------------------------------------------
  
  test_that("splitLexisDT and splitLexis congruent with one row per id", {
    for (sc in seq_along(BL)) {
      test_that(paste0("results congruent using breaks ", sc), {
        popEpi:::compareSLDTWithEpi(data = x, breaks = BL[[sc]], timeScale = names(BL)[sc])
      })
    }
  })
  
  
  
  # multiple rows per id ---------------------------------------------------------
  
  sire2 <- sire2[rep(1:.N, each = 2)]
  
  
  x <- Lexis(data=sire2[dg_date < ex_date], entry=list(fot=0, per=dg_yrs, age=dg_age),
             exit=list(per=ex_yrs), merge=TRUE, exit.status=1L, entry.status = 0L, id = id)
  
  
  test_that("splitLexisDT and splitLexis congruent with multiple rows per id", {
    for (sc in seq_along(BL)) {
      test_that(paste0("results congruent using breaks ", sc), {
        popEpi:::compareSLDTWithEpi(data = x, breaks = BL[[sc]], timeScale = names(BL)[sc])
      })
    }
  })
  
  
  # multistate using Lexis -----------------------------------------------------
  
  sire2[, EX := suppressWarnings(factor(status, levels = 0:2, labels = c("ok", "dead", "dead"), ordered = T))]    
  sire2[, EN := suppressWarnings(factor(0L, levels = 0:2, labels = c("ok", "dead", "dead"), ordered = T))]
  
  x <- Lexis(data=sire2[dg_date < ex_date & !duplicated(id)], entry=list(fot=0, per=bi_yrs, age=0),
             exit=list(per=ex_yrs), merge=TRUE, exit.status=EX, entry.status = EN, id = id)
  
  x <- cutLexis(x, cut = x$dg_yrs, timescale = "per", new.state = "sick", precursor.state = "ok")
  setDT(x)
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  
  test_that("splitLexisDT and splitLexis congruent with multiple Lexis states per id", {
    for (sc in seq_along(BL)) {
      test_that(paste0("results congruent using breaks ", sc), {
        popEpi:::compareSLDTWithEpi(data = x, breaks = BL[[sc]], timeScale = names(BL)[sc])
      })
    }
  })
  
  # using mstate package -------------------------------------------------------
  test_that("splitLexisDT and splitLexis congruent using mstate data", {
    
    ## taken directly from ?msprep examples
    library(mstate)
    tmat <- trans.illdeath()
    # some data in wide format
    tg <- data.frame(stt=rep(0,6),sts=rep(0,6),
                     illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
                     dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
                     x1=c(1,1,1,2,2,2),x2=c(6:1))
    tg$x1 <- factor(tg$x1,labels=c("male","female"))
    tg$patid <- factor(2:7,levels=1:8,labels=as.character(1:8))
    # define time, status and covariates also as matrices
    tt <- matrix(c(rep(NA,6),tg$illt,tg$dt),6,3)
    st <- matrix(c(rep(NA,6),tg$ills,tg$ds),6,3)
    keepmat <- data.frame(gender=tg$x1,age=tg$x2)
    # data in long format using msprep
    df <- msprep(time=tt,status=st,trans=tmat,keep=as.matrix(keepmat))
    
    x <- Lexis(data=df, entry=list(TT = Tstart), exit=list(TT = Tstop), 
               merge=TRUE, exit.status=to, entry.status = from, id = id)
    setDT(x)
    setattr(x, "class", c("Lexis", "data.table", "data.frame"))
    BL2 <- list(TT = 0:4, TT = 1:2, TT = 0:1, TT = c(1,Inf))
    
    for (sc in seq_along(BL2)) {
      test_that(paste0("results congruent using breaks ", sc), {
        popEpi:::compareSLDTWithEpi(data = x, breaks = BL2[[sc]], timeScale = "TT")
      })
    }
  })
  
})




