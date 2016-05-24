

context("Correction of typos")

v <- validate::validator( 
   x1 + x2 == x3
  ,x2 == x4
  ,x5 + x6 + x7 == x8
  ,x3 + x8 == x9
  ,x9 - x10 == x11
)

dat <- read.csv(textConnection(
"    , x1, x2 , x3  , x4 , x5 , x6, x7, x8 , x9   , x10 , x11
4  , 1452, 116, 1568, 116, 323, 76, 12, 411,  1979, 1842, 137
4.1, 1452, 116, 1568, 161, 323, 76, 12, 411,  1979, 1842, 137
4.2, 1452, 116, 1568, 161, 323, 76, 12, 411, 19979, 1842, 137
4.3, 1452, 116, 1568, 161,   0,  0,  0, 411, 19979, 1842, 137
4.4, 1452, 116, 1568, 161, 323, 76, 12,   0, 19979, 1842, 137"
))

test_that("correct_typos works",{   
   cor <- correct_typos(dat,v)
   expect_equal(cor[1,], dat[1,])
   expect_equal(as.integer(cor[2,]), as.integer(dat[1,]))
   expect_equal(as.integer(cor[3,]), as.integer(dat[1,]))
})


v2 <- validate::validator( x1 + x2 == x3 )


dat2 <- read.csv(textConnection(
" x2, x1 , x3  
  1452, 116, 1568
  1452, 161, 1568"
))


test_that("correct_typos reorder works",{   
   cor <- correct_typos(dat2, v2)
   expect_equivalent(cor[2,], dat2[1,])
})

test_that("correct_typos fixation works",{   
   # fixate x2, no real fixation
   cor <- correct_typos(dat2, v2, fixate="x2")
   expect_equivalent(cor[2,], dat2[1,])

   # fixate x1, no solution
   cor <- correct_typos(dat2, v2, fixate="x1")
   expect_equal(sum(abs(cor-dat2)),0)
})

test_that("correct_typos with inequality contraint record fails",{   
      # overconstraint editmatrix (only works for x1=0, and x2=0)
      v <- validate::validator( x1 == x2, x2 < 10 )

      dat <- data.frame(
       x1 = c(10,29),
       x2 = c(1,92))

      cor <- correct_typos(dat,v)
      expect_equivalent(rowSums(values(confront(cor,v)))==2,c(FALSE,FALSE))
})


test_that("correct_typos with noncorrectable record works",{   
      # overconstraint editmatrix (only works for x1=0, and x2=0)
      v <- validate::validator(x1 == x2, 9*x1 == x2)
                    
      dat <- data.frame(
       x1 = 10,
       x2 = 99)

      cor <- correct_typos(dat,v)
      expect_true(rowSums(values(confront(cor,v)))[1] < 2)
})

test_that("correct_typos with noncorrectable record works",{   
      # overconstraint editmatrix (only works for x1=0, and x2=0)
      v <- validate::validator(x1 <= x2, 9*x1 == x2)
                     
      dat <- data.frame(
       x1 = 10,
       x2 = 99)

      cor <- correct_typos(dat,v)

      expect_true(rowSums(values(confront(cor,v)))[1]<2)
})


test_that("correct_typos with missing variable works",{   
   # valid edit matrix (but missing x4)
   v <- validate::validator(x1 == x2 + x3 + x5 + x6)

   dat <- data.frame(
    x1 = 42280000,
    x2 = 11289000,
    x3 = 4328000,
    x4 = 361300,
    x5 = 11201000,
    x6 = 11849000)

   # fail:
   cor <- correct_typos(dat,v)
   #print(cor)
   expect_true(rowSums(values(confront(cor,v)))[1]<2)

})

