library(EcoSimR)
context("Niche overlap null model tests")

nicheTmat <- matrix(rnorm(100),ncol=10,nrow=10)

nicheTmat[c(1,19,13,23,42,90)] <- 0
uval <- length(unique(nicheTmat))



test_that("ra1 algorithm works:",{
  expect_true(is.matrix(ra1()))
  expect_true(is.matrix(ra1(dataMacWarb[,2:5])))
  
}
)


test_that("ra2 algorithm works:",{
  expect_true(is.matrix(ra2()))
  expect_true(is.matrix(ra2(as.matrix(dataMacWarb[,2:5]))))
  expect_equal(length(which(ra2(nicheTmat)==0)),6)
}
)

test_that("ra3 algorithm works:",{
  expect_true(is.matrix(ra3()))
  expect_true(is.matrix(ra3(dataMacWarb[,2:5])))
  expect_equal(length(unique(ra3(nicheTmat))),uval)  
}
)


test_that("ra4 algorithm works:",{
  expect_true(is.matrix(ra4()))
  expect_true(is.matrix(ra4(dataMacWarb[,2:5])))
  expect_equal(sum(apply(ra4(nicheTmat),2,function(x){sum(x==0)}) - apply(nicheTmat,2,function(x){sum(x==0)})), 0)
  expect_equal(sum(apply(ra4(nicheTmat),1,function(x){sum(x==0)}) - apply(nicheTmat,1,function(x){sum(x==0)})), 0)
}
)

test_that("pianka metric works:",{
expect_true(is.numeric(pianka(ra1())))
expect_true(is.numeric(pianka(ra4(dataMacWarb[,2:5]))))
expect_true(pianka(ra4(dataMacWarb[,2:5])) > .5 )

})

test_that("czekanowski metric works:",{
  expect_true(is.numeric(czekanowski(ra1())))
  expect_true(is.numeric(czekanowski(ra4(dataMacWarb[,2:5]))))
  expect_true(czekanowski(ra4(dataMacWarb[,2:5])) > .5 )
  
})

test_that("pianka_var metric works:",{
  expect_true(is.numeric(pianka_var(ra1())))
  expect_true(is.numeric(pianka_var(ra4(dataMacWarb[,2:5]))))
  expect_true(pianka_var(ra4(dataMacWarb[,2:5])) < .5 )
  
})


test_that("czekanowski_var metric works:",{
  expect_true(is.numeric(czekanowski_var(ra1())))
  expect_true(is.numeric(czekanowski_var(ra4(dataMacWarb[,2:5]))))
  expect_true(czekanowski_var(ra4(dataMacWarb[,2:5])) < .5 )
  
})

test_that("pianka_skew metric works:",{
  expect_true(is.numeric(pianka_skew(ra1())))
  expect_true(is.numeric(pianka_skew(ra4(dataMacWarb[,2:5]))))
  
})


test_that("czekanowski_skew metric works:",{
  expect_true(is.numeric(czekanowski_skew(ra1())))
  expect_true(is.numeric(czekanowski_skew(ra4(dataMacWarb[,2:5]))))  
})


test_that("niche_null_model works with all combinations of metrics and algorithms",{
  ### Test that proper object is returned
  expect_is(niche_null_model(dataMacWarb,metric ="pianka" ,algo = "ra1",nRep=10),"nichenullmod")
  expect_is(niche_null_model(dataMacWarb,metric ="pianka_var" ,algo = "ra1",nRep=10),"nichenullmod")
  expect_is(niche_null_model(dataMacWarb,metric ="pianka_skew" ,algo = "ra1",nRep=10),"nichenullmod")
  expect_is(niche_null_model(dataMacWarb,metric ="czekanowski" ,algo = "ra1",nRep=10),"nichenullmod")
  expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_var" ,algo = "ra1",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_skew" ,algo = "ra1",nRep=10),"nichenullmod")
 
expect_is(niche_null_model(dataMacWarb,metric ="pianka" ,algo = "ra2",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="pianka_var" ,algo = "ra2",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="pianka_skew" ,algo = "ra2",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski" ,algo = "ra2",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_var" ,algo = "ra2",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_skew" ,algo = "ra2",nRep=10),"nichenullmod")


expect_is(niche_null_model(dataMacWarb,metric ="pianka" ,algo = "ra3",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="pianka_var" ,algo = "ra3",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="pianka_skew" ,algo = "ra3",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski" ,algo = "ra3",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_var" ,algo = "ra3",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_skew" ,algo = "ra3",nRep=10),"nichenullmod")


expect_is(niche_null_model(dataMacWarb,metric ="pianka" ,algo = "ra4",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="pianka_var" ,algo = "ra4",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="pianka_skew" ,algo = "ra4",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski" ,algo = "ra4",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_var" ,algo = "ra4",nRep=10),"nichenullmod")
expect_is(niche_null_model(dataMacWarb,metric ="czekanowski_skew" ,algo = "ra4",nRep=10),"nichenullmod")


nmod <- niche_null_model(dataMacWarb,metric ="czekanowski" ,algo = "ra1",nRep=10)

expect_output(summary(nmod),"Metric:  czekanowski")
expect_true(is.null(plot(nmod,type="hist")))
expect_true(is.list(plot(nmod,type="niche")))
})


  
test_that("all text data frames are handled proprely",{
  textMat <- cbind(letters[1:dim(nicheTmat)[1]],nicheTmat)
  expect_is(niche_null_model(textMat,metric ="czekanowski_skew" ,algo = "ra4",nRep=10),"nichenullmod")
  
})





