require(testthat)
context("Error Localization: numerical data")

test_that("errorLocalizer for numerical data",{
    bt <- errorLocalizer(
            E = editmatrix("x + y == z"),
            x = c(x=1,y=1,z=2))
    expect_true( all(bt$searchNext()$adapt==c(FALSE,FALSE,FALSE)) )
    expect_true( is.null(bt$searchNext()))
    expect_true(1 == errorLocalizer(editmatrix("x+y==z"),c(x=1,y=1,z=3))$searchNext()$w)
    expect_true(1 == errorLocalizer(editmatrix("x+y==z"),c(x=1,y=1,z=3,u=5))$searchNext()$w)
    expect_true(is.null(errorLocalizer(editmatrix("x+y==z"),c(x=1,y=1,z=3))$searchNext(maxduration=-1)) )
    expect_true(is.null(errorLocalizer(editmatrix("x+y==z"),c(x=1,y=1,z=3),maxadapt=0)$searchNext()) )
    expect_true(is.null(errorLocalizer(editmatrix("x+y==z"),c(x=1,y=1,z=3),maxweight=0)$searchNext()) )
    expect_that(errorLocalizer(editmatrix("x+y==z"),c(x=1,y=NA,z=3),weight=c(1,NA,1))$searchNext(),throws_error())
})

test_that("weight calculation when checkDatamodel is activated",{
   expect_equal(
      localizeErrors(
         editmatrix(expression(x+y+z==w,x>0)),
         data.frame(x=-1,y=1,z=0,w=0)
      )$status$weight, 2
   )
   expect_equal(
      localizeErrors(
         editmatrix(expression(x+y+z==w,x>0)),
         data.frame(x=-1,y=1,z=0,w=0),
         method='mip'
      )$status$weight, 2
   )

})



#d <- "../../../pkg/R"
#for ( b in file.path(d,dir(d)) ) dmp <- source(b,echo=FALSE)

context("Error localization: categorical data")
test_that("errorLocalizer for categorical data",{
    E <- editarray(c(
        "positionInHouseHold %in% c('marriage partner','child','other')",
        "age %in% c('under aged','adult')",
        "maritalStatus %in% c('unmarried','married','widowed')",
        "if (age == 'under aged') maritalStatus == 'unmarried'",
        "if (maritalStatus != 'unmarried' ) !positionInHouseHold %in% c('child','marriage partner')"
        ))
    r <- c(age='under aged',maritalStatus='married',positionInHouseHold='child')
    e <- errorLocalizer(E,r)
    expect_equivalent(e$searchBest()$adapt,c(FALSE,TRUE,FALSE))
    expect_true(e$degeneracy==1)
    e$reset()
    expect_true(e$searchBest()$w==1)
    expect_true(is.null(e$searchNext()))
    expect_true(is.null(errorLocalizer(E,r,maxweight=0)$searchBest()) )
    expect_true(is.null(errorLocalizer(E,r,maxadapt=0)$searchBest()) )
    expect_true(is.null(errorLocalizer(E,r)$searchBest(maxduration=-1)) )
    expect_that(errorLocalizer(E,r,weight=c(1,NA,1)),throws_error())
})

test_that("errorlocalizer.editarray runs when blocks reduces the datamodel",{
# Thanks to Bikram Maharjan for reporting a bug in version 2.0-1
    ageGrps=rep(1,5)
    genderIndex=rep(1,5)
    icdIndex=rep(103,5)
    diagnosisIndex=c('TOT', 'A00', 'A02', 'C53','C54')
     
     
    e.diag=editarray(c(
      "ageGrps %in% c(1,2,3,4)",
      "genderIndex %in% c(0,1,2)",
      "icdIndex %in% c(103)",
      "diagnosisIndex %in% c('TOT', 'A00', 'A02', 'C53','C54')", 
      "if(genderIndex=='1' & icdIndex=='103') diagnosisIndex %in% c('C53','C54')"
    ))
    dat.diag=data.frame(ageGrps,genderIndex,icdIndex,diagnosisIndex)
    err=localizeErrors(e.diag,dat.diag)
})

test_that("if-else edits are parsed to MIP correctly",{
  e <- editset("if (x > 0 ) y> 0")
  d <- data.frame(x = 1,y=0)
  localizeErrors(e,d,method='mip')
})




