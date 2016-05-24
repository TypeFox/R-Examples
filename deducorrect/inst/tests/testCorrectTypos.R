library(testthat)
library(editrules)


context("Correction of typos")

E <- editmatrix( c("x1 + x2 == x3"
                  ,"x2 == x4"
                  ,"x5 + x6 + x7 == x8"
                  ,"x3 + x8 == x9"
                  ,"x9 - x10 == x11"
                  )
               )

dat <- read.csv(textConnection(
"    , x1, x2 , x3  , x4 , x5 , x6, x7, x8 , x9   , x10 , x11
4  , 1452, 116, 1568, 116, 323, 76, 12, 411,  1979, 1842, 137
4.1, 1452, 116, 1568, 161, 323, 76, 12, 411,  1979, 1842, 137
4.2, 1452, 116, 1568, 161, 323, 76, 12, 411, 19979, 1842, 137
4.3, 1452, 116, 1568, 161,   0,  0,  0, 411, 19979, 1842, 137
4.4, 1452, 116, 1568, 161, 323, 76, 12,   0, 19979, 1842, 137"
))

test_that("correctTypos works",{   
   cor <- correctTypos(E,dat)
   corrected <- cor$corrected
   expect_equal(corrected[1,], dat[1,])
   
   expect_equal(as.integer(corrected[2,]), as.integer(dat[1,]))
   
   expect_equal(as.integer(corrected[3,]), as.integer(dat[1,]))
})


E2 <- editmatrix( c("x1 + x2 == x3"
                  )
               )

dat2 <- read.csv(textConnection(
"    , x2, x1 , x3  
     , 1452, 116, 1568
     , 1452, 161, 1568"
))


test_that("correctTypos reorder works",{   
   cor <- correctTypos(E2,dat2)
   corrected <- cor$corrected
   expect_equivalent(corrected[2,], dat2[1,])
})

test_that("correctTypos fixation works",{   
   # fixate x2, no real fixation
   cor <- correctTypos(E2,dat2, fixate="x2")
   corrected <- cor$corrected
   expect_equivalent(corrected[2,], dat2[1,])

   # fixate x1, no solution
   cor <- correctTypos(E2,dat2, fixate="x1")
   expect_equal(nrow(cor$corrections), 0)
})

test_that("correctTypos with inequality contraint record fails",{   
      # overconstraint editmatrix (only works for x1=0, and x2=0)
      E <- editmatrix( c( "x1 == x2"
                        , "x2 < 10"
                        )
                     )
      data <- data.frame(
       x1 = c(10,29),
       x2 = c(1,92))

      cor <- correctTypos(E,data)
      expect_equal(as.character(cor$status$status), c("invalid", "invalid"))
})


test_that("correctTypos with noncorrectable record works",{   
      # overconstraint editmatrix (only works for x1=0, and x2=0)
      E <- editmatrix( c( "x1 == x2"
                        , "9*x1 == x2"
                        )
                     )
      data <- data.frame(
       x1 = 10,
       x2 = 99)

      cor <- correctTypos(E,data)

      expect_equal(as.character(cor$status$status[1]), "invalid")
})

test_that("correctTypos with noncorrectable record works",{   
      # overconstraint editmatrix (only works for x1=0, and x2=0)
      E <- editmatrix( c( "x1 <= x2"
                        , "9*x1 == x2"
                        )
                     )
      data <- data.frame(
       x1 = 10,
       x2 = 99)

      cor <- correctTypos(E,data)

      expect_equal(as.character(cor$status$status[1]), "invalid")
})


test_that("correctTypos with missing variable works",{   
   # valid edit matrix (but missing x4)
   E <- editmatrix(
    "x1 == x2 + x3 + x5 + x6")

   #print(E)
   data <- data.frame(
    x1 = 42280000,
    x2 = 11289000,
    x3 = 4328000,
    x4 = 361300,
    x5 = 11201000,
    x6 = 11849000)

   # fail:
   cor <- correctTypos(E,data)
   #print(cor)
   expect_equal(as.character(cor$status$status[1]), "invalid")
})






test_that("correctTypos.editset works with pure numerical",{

    v <- correctTypos(
        editset(expression(
        x + y == z))$num,
        data.frame(
        x=1,
        y=1,
        z=1
        )
    )
    w <- correctTypos(
        editset(expression(
        x + y == z)),
        data.frame(
        x=1,
        y=1,
        z=1
        )
    )
    expect_equal(v$corrected,w$corrected)
    expect_equal(v$corrections, w$corrections) 
})


test_that("correctTypos.editset works with pure categorical",{

    v <- correctTypos(
        editset(expression(
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b'
        )),
        data.frame(
            A = 'a',
            B = NA)
    )
    expect_equal( nrow(v$corrections), 0 )

})

test_that("correctTypos.editset works with unconnected numeric/categorical",{

    E <- editset(expression(
        x + y == z,
        x >= 0,
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a' ) B == 'b'
    ))

    x <- data.frame(
        x = 1,
        y = 12,
        z = 1,
        A = 'a',
        B = NA
    )
    v <- correctTypos(E,x)
    w <- correctTypos(E$num,x)
    expect_equal(v$corrected,w$corrected)
    expect_equal(v$corrections, w$corrections)
})


test_that("correctTypos.editset works with unconnected numeric/categorical",{
    
    E <- editset(expression(
        x + y == z,
        x >= 0,
        y > 0,
        y < 2,
        z > 1,
        z < 3,
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a' ) B == 'b',
        if ( B == 'b' ) x > 0
    ))
    
    x <- data.frame(
        x = 10,
        y = 1,
        z = 2,
        A = 'a',
        B = NA
    )
    
    
    v <- correctTypos(E,x)
    w <- correctTypos(E$num,x)
    expect_equal(v$corrected,w$corrected)
    expect_equal(v$corrections, w$corrections)
     
    ## with revert
    E <- editset(expression(
        x + y == z,
        x >= 0,
        y > 0,
        y < 2,
        z > 1,
        z < 3,
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a' ) B == 'b',
        if ( B == 'b' ) x > 3
    ))
    
    x <- data.frame(
        x = 10,
        y = 1,
        z = 2,
        A = 'a',
        B = 'b'
    )
    
    
    v <- correctTypos(E,x)
    expect_equal(nrow(v$corrections),0)

})


