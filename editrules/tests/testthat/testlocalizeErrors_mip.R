require(testthat)

context("Localize errors using MIP")

test_that("localizeError.mip",{
  Et <- editmatrix(expression( p + c == t
                             , c - 0.6*t >= 0
                             , c>=0
                             , p >=0
                             )
                  )
  
  x <- c(p=755,c=125,t=200)
  
  sol <- errorLocalizer_mip(Et, x)
  expect_equal(sol$w, 1)
  expect_equivalent(sol$adapt, c(TRUE, FALSE, FALSE))
  expect_equal(unname(sol$x_feasible[1:3]), c(75, 125, 200), tolerance=1e-10)
})

test_that('localizeErrors works without specified weight',{
  E <- editmatrix(c('x+y==z','x < 1'))
  dat <- data.frame(
    x = c(1,0,2),
    y = c(1,1,1),
    z = c(1,1,1)
  )
  
  loc <- localizeErrors( E = E
                       , dat = dat
                       , method="mip"
                       )
  
#   print(loc)
  expect_equivalent( loc$adapt
                   , matrix(c(
                      TRUE , FALSE, FALSE,
                      FALSE , FALSE, FALSE,
                      TRUE , FALSE, FALSE),
                           nrow=3,
                           byrow=TRUE
                           )
                  )
  
  
})



test_that('localizeErrors works with single specified weight',{
  
  expect_equivalent(localizeErrors(
    E       = editmatrix(c('x+y==z','x<1')),
    dat     = data.frame(
      x = c(1,1,1),
      y = c(1,1,1),
      z = c(1,1,1)
      ),
    weight  = c(1,2,2),
    method="mip"
    )$adapt,
                    matrix(c(
                      TRUE , FALSE, FALSE,
                      TRUE , FALSE, FALSE,
                      TRUE , FALSE, FALSE),
                           nrow=3,
                           byrow=TRUE
                           )
                    )
  
  
})


test_that('localizeErrors works with different weights per record',{
  le <- localizeErrors(
    E       = editmatrix('x+y==z'),
    dat     = data.frame(
      x = c(1,1,1),
      y = c(1,1,1),
      z = c(1,1,1)
      ),
    weight  = matrix(c(
      1,2,2,
      2,1,2,
      2,2,1),
                     nrow=3,
                     byrow=TRUE
                     ),
    method="mip"
    )
    expect_equivalent(le$adapt,
                    matrix(c(
                      TRUE , FALSE, FALSE,
                      FALSE, TRUE , FALSE ,
                      FALSE, FALSE, TRUE),
                           nrow=3,
                           byrow=TRUE
                           )
                    )
  expect_that(localizeErrors(
    E = editmatrix('x +y==z'),
    dat     = data.frame(
      x = c(1,1,1),
      y = c(1,1,1),
      z = c(1,1,1)
      ),
    weight  = matrix(c(
      1,2,2,
      2,2,1),
                     nrow=3,
                     byrow=TRUE
                     ),
    method="mip"
    ),
              throws_error()
              )
})

test_that('localizeErrors handles single edits with mip method',{
    expect_true(
        localizeErrors(editmatrix("x>0"),data.frame(x=-1),method='mip')$adapt[1,1]
    )
    loc <- localizeErrors(editmatrix("x>0"),data.frame(x=NA),method='mip')
    expect_true(loc$adapt[1,1])
})

test_that('localizeErrors handles trivial single edits with mip method',{
  expect_false( 
    localizeErrors( editmatrix("z == 100")
                              , data.frame(z=100)
                              , method='mip'
                              )$adapt[1,1])
})

test_that('localizeErrors handles a ">" edits correctly.',{
  loc <- localizeErrors( editmatrix("x > 1")
                       , data.frame(x=1)
                       , method='mip'
                       )
  #print(loc)
  expect_true( localizeErrors( editmatrix("x < 1")
                             , data.frame(x=1)
                             , method='mip'
                             )$adapt[1,1]
             )
})


## editset
test_that("localizeError_mip editset",{
  Et <- editset(expression( p + c == t
                               , c - 0.6*t >= 0
                               , c>=0
                               , p >=0
                               )
                   )
  
  x <- c(p=755,c=125,t=200)
  
  sol <- errorLocalizer_mip(Et, x)
  expect_equal(sol$w, 1)
  expect_equivalent(sol$adapt, c(TRUE, FALSE, FALSE))
  expect_equivalent(sol$x_feasible, c(75, 125, 200))
})

test_that('localizeErrors works without specified weight, editset',{
  E <- editset(c('x+y==z','x < 1'))
  dat <- data.frame(
    x = c(1,0,2),
    y = c(1,1,1),
    z = c(1,1,1)
    )
  
  loc <- localizeErrors( E = E
                         , dat = dat
                         , method="mip"
                         )
  
  #print(loc)
  expect_equivalent( loc$adapt
                     , matrix(c(
                       TRUE , FALSE, FALSE,
                       FALSE , FALSE, FALSE,
                       TRUE , FALSE, FALSE),
                              nrow=3,
                              byrow=TRUE
                              )
                     )
  
  
})



test_that('localizeErrors works with single specified weight, editset',{
  
  expect_equivalent(localizeErrors(
    E       = editset(c('x+y==z','x<1')),
    dat     = data.frame(
      x = c(1,1,1),
      y = c(1,1,1),
      z = c(1,1,1)
      ),
    weight  = c(1,2,2),
    method="mip"
    )$adapt,
                    matrix(c(
                      TRUE , FALSE, FALSE,
                      TRUE , FALSE, FALSE,
                      TRUE , FALSE, FALSE),
                           nrow=3,
                           byrow=TRUE
                           )
                    )
  
  
})


test_that('localizeErrors works with different weights per record, editset',{
  expect_equivalent(localizeErrors(
    E       = editset('x+y==z'),
    dat     = data.frame(
      x = c(1,1,1),
      y = c(1,1,1),
      z = c(1,1,1)
      ),
    weight  = matrix(c(
      1,2,2,
      2,1,2,
      2,2,1),
                     nrow=3,
                     byrow=TRUE
                     ),
    method="mip"
    )$adapt,
                    matrix(c(
                      TRUE , FALSE, FALSE,
                      FALSE, TRUE , FALSE ,
                      FALSE, FALSE, TRUE),
                           nrow=3,
                           byrow=TRUE
                           )
                    )
  expect_that(localizeErrors(
    E = editset('x +y==z'),
    dat     = data.frame(
      x = c(1,1,1),
      y = c(1,1,1),
      z = c(1,1,1)
      ),
    weight  = matrix(c(
      1,2,2,
      2,2,1),
                     nrow=3,
                     byrow=TRUE
                     ),
    method="mip"
    ),
              throws_error()
              )
})

test_that('localizeErrors handles single edits containing NA with mip method, editset',{
  expect_true(
    localizeErrors(editset("x>0"),data.frame(x=-1),method='mip')$adapt[1,1]
    )
  loc <- localizeErrors(editset("x>0"),data.frame(x=NA),method='mip')
  expect_true(loc$adapt[1,1])
})

test_that('localizeErrors handles trivial single edits with mip method, editset',{
  expect_false( 
    localizeErrors( editset("z == 100")
                    , data.frame(z=100)
                    , method='mip'
                    )$adapt[1,1])
})

test_that('localizeErrors handles a "<" edits correctly., editset',{
  loc <- localizeErrors( editset("x > 1")
                         , data.frame(x=1)
                         , method='mip'
                         )
  #print(loc)
  expect_true( localizeErrors( editset("x < 1")
                               , data.frame(x=1)
                               , method='mip'
                               )$adapt[1,1]
               )
})

test_that('localizeErrors',{
  E <- editset(c(
    "age %in% c('under aged','adult')",
    "maritalStatus %in% c('unmarried','married','widowed','divorced')",
    "positionInHousehold %in% c('marriage partner', 'child', 'other')",
    "if( age == 'under aged' ) maritalStatus == 'unmarried'",
    "if( maritalStatus %in% c('married','widowed','divorced')) !positionInHousehold %in% c('marriage partner','child')"
    ))
  record <- data.frame(age='under aged', maritalStatus='married', positionInHousehold='child')
  expect_equivalent(
    localizeErrors(E,record, method="mip")$adapt
    ,
    array(c(FALSE,TRUE,FALSE),dim=c(1,3))
    )
})

test_that('localizeErrors with outofrange error',{
  E <- editset(c(
    "age %in% c('under aged','adult')",
    "maritalStatus %in% c('unmarried','married','widowed','divorced')",
    "positionInHousehold %in% c('marriage partner', 'child', 'other')",
    "if( age == 'under aged' ) maritalStatus == 'unmarried'",
    "if( maritalStatus %in% c('married','widowed','divorced')) !positionInHousehold %in% c('marriage partner','child')"
    ))
  record <- data.frame(age='under aged', maritalStatus='unmarried', positionInHousehold='out-of-range')
  expect_equivalent(
    localizeErrors(E,record, method="mip")$adapt
    ,
    array(c(FALSE,FALSE,TRUE),dim=c(1,3))
    )
})

test_that("localizeErrors works with TRUE/FALSE",{
  E <- editset(expression(
    A %in% c(TRUE,FALSE),
    B %in% letters[1:4],
    if ( !A ) B %in% letters[1:2]
    ))
  
  # should run without errors...
  localizeErrors(E,data.frame(A=c(TRUE,FALSE),B=c('c',"d")), method="mip")
})


context("MIP-based error localization numerical stability tests")

test_that("Records of range 1-1e9",{
   p <- 1:9
   x <- 10^(sqrt(p))
   
   names(x) <- paste("x",p,sep="")
   e <- editmatrix(expression(
      x1 + x2 == x3,
      x4 + x5 + x6 == x7,
      x7 + x3 + x8 == x9
   ))
   expect_equal(
      errorLocalizer_mip(e,x)$w,
      errorLocalizer(e,x)$searchBest()$w
   )
   
   x[] <- 10^(p)
   # MIP is sensitive to (very) large differences in values (mvdl/ejne 11.08.2012)
   expect_true(
     errorLocalizer_mip(e,x)$w != errorLocalizer(e,x)$searchBest()$w
   )
})


test_that("Consistency with B&B algorithm",{
# generated from smoketest. Used to generate bug in weight count.
   r <- c(
      x1 = -0.556503362667066,
      x2 = -0.0839342133749159,
      x3 = 0.775427129507229,
      x4 = -1.94883105542653,
      x5 = 1.11931659004076,
      x6 = -1.453437500466,
      x7 = 0.47278745357947
   )
   e <- editmatrix(expression(
      x1 + x2 == x3
      ,x3 + x4 == x5
      ,x3 + x5 == x7
      ,0 <= x1
      ,0 <= x2
      ,0 <= x3
      ,0 <= x4
      ,0 <= x5
      ,0 <= x6
      ,0 <= x7
   ))

   expect_equal(
      errorLocalizer(e,r)$searchBest()$w, 
      errorLocalizer_mip(e,r)$w
   )
   X <- as.data.frame(t(r))
   expect_equal(
      localizeErrors(e,X)$status$weight[1],
      localizeErrors(e,X,method="mip")$status$weight[1]
   )
})


test_that("strict inequalities are treated correctly", {
  E <- editmatrix("x > 0")
  expect_that(errorLocalizer_mip(E, c(x=0))$adapt, is_equivalent_to(TRUE))
  expect_that(errorLocalizer_mip(E, c(x=1e-4))$adapt, is_equivalent_to(TRUE))

  expect_that(errorLocalizer_mip(E, c(x=1e-3))$adapt, is_equivalent_to(FALSE))
})











