library(DNAprofiles)
context("RMP is KI with UN vs ID")

data(freqsNLsgmplus,freqsNLngm)
set.seed(123)
x1 <- sample.profiles(N = 1e3,freqs = freqsNLngm)

test_that(desc = "RMP equals KI with hyp.1=UN versus hyp.2=ID",{
  tol <- 1e-6
  
  # rmp is ki with ID vs UN
  r1 <- ki(x1 = x1,x2 = x1,hyp.1 = "UN",hyp.2 = "ID")
  r2 <- rmp(x1)
  expect_true(all.equal(r1,r2, tolerance = tol,scale = r1),TRUE)
  
  # with subset of markers
  r1 <- ki(x1 = x1,x2 = x1,hyp.1 = "UN",hyp.2 = "ID",markers = names(freqsNLsgmplus))
  r2 <- rmp(x1,markers = names(freqsNLsgmplus))
  expect_true(all.equal(r1,r2, tolerance = tol,scale = r1),TRUE)  
  
  # with theta>0
  r1 <- ki(x1 = x1,x2 = x1,hyp.1 = "UN",hyp.2 = "ID",markers = names(freqsNLsgmplus),theta = 1e-1)
  r2 <- rmp(x1,markers = names(freqsNLsgmplus),cmp = TRUE,theta = 1e-1)
  expect_true(all.equal(r1,r2, tolerance = tol,scale = r1),TRUE)  
  
  # with missing values
  i <- sample(seq(nrow(x1)),size = 1e2)
  x1[i,"VWA.1"] <- NA
  x1[i,"VWA.2"] <- NA
  x1[i,"TH01.1"] <- NA
  x1[i,"TH01.2"] <- NA
  
  r1 <- ki(x1 = x1,x2 = x1,hyp.1 = "UN",hyp.2 = "ID",markers = names(freqsNLsgmplus),theta = 1e-1)
  r2 <- rmp(x1,markers = names(freqsNLsgmplus),cmp = TRUE,theta = 1e-1)
  expect_true(all.equal(r1,r2, tolerance = tol,scale = r1),TRUE)
  
  # with NAs and per marker
  r1 <- ki(x1 = x1,x2 = x1,hyp.1 = "UN",hyp.2 = "ID",markers = names(freqsNLsgmplus),theta = 1e-1,ret.per.marker = TRUE)
  r2 <- rmp(x1,markers = names(freqsNLsgmplus),cmp = TRUE,theta = 1e-1,ret.per.marker = TRUE)
  expect_true(all.equal(r1,r2, tolerance = tol,scale = r1),TRUE)
})