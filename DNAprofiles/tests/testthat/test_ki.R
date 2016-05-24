library(DNAprofiles)
context("KI computation")

test_that(desc = "KIs can be computed",{
  data(freqsNLngm,freqsNLsgmplus)
  
  set.seed(123)
  
  for (theta in c(0,0.1,1)){
    # one-to-one (single pair)
    x <- sample.profiles(N = 1,freqsNLngm)
    y <- sample.profiles(N = 1,freqsNLngm)  
    r <- ki(x,y,hyp.1 = "FS",theta = theta)  
    expect_equal(1L,length(r))
    expect_equal(class(r),"matrix")
    
    # one-to-one (many pairs)
    x <- sample.profiles(N = 1e2,freqsNLngm)
    y <- sample.profiles(N = 1e2,freqsNLngm)  
    r <- ki(x,y,hyp.1 = "FS",theta = theta)  
    expect_equal(100L,length(r))
    expect_equal(class(r),"matrix")
    
    # one-to-many (single profile vs db)
    x <- sample.profiles(N = 1,freqsNLngm)
    y <- sample.profiles(N = 1e2,freqsNLngm)    
    expect_error(r <- ki(x,y,hyp.1 = "FS",theta = theta)  )
    r <- ki.db(x,y,"FS",theta = theta)
    expect_equal(100L,length(r))
    expect_equal(class(r),"matrix")
    
    # one-to-many (many profiles vs db)
    x <- sample.profiles(N = 50,freqsNLngm)
    y <- sample.profiles(N = 1e2,freqsNLngm)    
    R <- ki.db(x,y,"FS",theta = theta)
    expect_equal(dim(R),c(100L,50L))
    expect_equal(class(R),"matrix")    
  }
  
  # with return per marker and a subset of markers
  
  for (theta in c(0,0.1,1)){
    markers <- names(freqsNLsgmplus)
    # one-to-one (single pair)
    x <- sample.profiles(N = 1,freqsNLngm)
    y <- sample.profiles(N = 1,freqsNLngm)  
    r <- ki(x,y,hyp.1 = c(0,1,0),markers = markers,theta = theta,ret.per.marker = TRUE)  
    
    expect_equal(dim(r),c(1,length(markers)))
    expect_equal(class(r),"matrix")
    expect_false(any(is.na(r)))
    
    # one-to-one (many pairs)
    x <- sample.profiles(N = 1e2,freqsNLngm)
    y <- sample.profiles(N = 1e2,freqsNLngm)  
    r <- ki(x,y,hyp.1 = c(0,1,0),markers = markers,theta = theta,ret.per.marker = TRUE)  
    
    expect_equal(dim(r),c(100L,length(markers)))
    expect_equal(class(r),"matrix")
    expect_false(any(is.na(r)))
        
    # one-to-many (single profile vs db)
    x <- sample.profiles(N = 1,freqsNLngm)
    y <- sample.profiles(N = 1e2,freqsNLngm)    
    expect_error( ki(x,y,hyp.1 = c(0,1,0),markers = markers,theta = theta,ret.per.marker = TRUE))
    r <- ki.db(x,y,hyp.1 = c(0,1,0),markers = markers,theta = theta,ret.per.marker = TRUE)
    expect_equal(dim(r),c(100L,10L))
    expect_equal(class(r),"matrix")
    expect_false(any(is.na(r)))
    
    # one-to-many (many profiles vs db)
    x <- sample.profiles(N = 50,freqsNLngm)
    y <- sample.profiles(N = 1e2,freqsNLngm)    
    expect_error(r <- ki.db(x,y,hyp.1 = c(0,1,0),markers = markers,theta = theta,ret.per.marker = TRUE))
  }
  
  
})