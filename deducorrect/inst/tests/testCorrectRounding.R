library(testthat)
library(editrules)

context("Correction of rounding errors")

test_that("correctRounding works",{
   R <- editmatrix(c("a == 1"))
   dat <- data.frame(a=2)
   sol <- correctRounding(R,dat)
   expect_equivalent(sol$corrected[1,], data.frame(a=1))
})



test_that("correctRounding with fixate works",{
   R <- editmatrix(c("a + b == 2"))
   dat <- data.frame(a=2,b=1)
   sol <- correctRounding(R,dat, fixate="b")
   expect_equivalent(sol$corrected[1,], data.frame(a=1,b=1))
})

test_that("correctRounding with Q works",{
   set.seed(1)
   R <- editmatrix(c("a + b == 2", "b>0"))
   dat <- data.frame(a=2,b=1)
   sol <- correctRounding(R,dat)
   expect_equivalent(sol$corrected[1,], data.frame(a=1,b=1))
})

test_that("correctRounding works with Scholtus 2008 example",{
   E <- editmatrix(c( "x1 + x2 == x3"
                    , "x2 == x4"
                    , "x5 + x6  + x7 == x8"
                    , "x3 + x8 == x9"
                    , "x9 - x10 == x11"
                    )
                  )

   dat <- data.frame( x1=12
                    , x2=4
                    , x3=15
                    , x4=4
                    , x5=3
                    , x6=1
                    , x7=8
                    , x8=11
                    , x9=27
                    , x10=41
                    , x11=-13
                    )
   sol <- correctRounding(E,dat)
   expect_equal(as.character(sol$status$status), "corrected")
})


# smoke test
context("Smoke test of correctrounding")
s <- sample(.Machine$integer.max,1)
cat("Randseed is ",s,":")
set.seed(s)
E <- editmatrix(c("x+y==z","x>=0","y>=0"))
for ( i in 1:10 ){
    cat(i); flush.console()
    x <- sample(0:2,10,replace=TRUE)
    y <- sample(0:2,10,replace=TRUE)
    z <- x + y + sample(c(-1,1),10,replace=TRUE) 
    dat <- data.frame(x=x,y=y,z=z)
    v <- correctRounding(E,dat)
    test_that("No extra inequalities are generated",{
        f1 <- violatedEdits(E,v$corrected)
        f2 <- violatedEdits(E,dat)
        expect_true(all(which(f1[,2]) %in% which(f1[,2])))
        expect_true(all(which(f1[,3]) %in% which(f2[,3])))
    })
}

test_that("correctRounding.editset works with pure numerical",{
    v <- correctRounding(
        editset("x + y == z"),
        data.frame(x=1,y=1,z=1)
    )
    expect_equal(nrow(v$corrections),1)
})


test_that("correctRounding.editset works with pure categorical",{
    
    v <- correctRounding(
        editset(expression(
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b')
        ),
        data.frame(
            A = 'a',
            B = NA)
    )
    expect_equal(nrow(v$corrections),0)
    
})


test_that("correctRounding.editset works with unconnected numeric/categorical",{

    v <- correctRounding(
        editset(expression(
            x + y == z,
            x >= 0,
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b'
        )),
        data.frame(
            x = 1,
            y = 1,
            z = 1,
            A = 'a',
            B = NA
        )
    )
    expect_equal(nrow(v$corrections),1)
})

test_that("correctRounding.editset works with connected numeric/categorical",{
    # with NA (uncheckable)
    v <- correctRounding(editset(expression(
            x + y == z,
            x >= 0,
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b',
            if ( B == 'b' ) x > 0
        )),
        data.frame(
            x = NA,
            y = 1,
            z = 1,
            A = 'a',
            B = NA
        )
    )
    # without revert
    v <- correctRounding(editset(expression(
            x + y == z,
            x >= 0,
            y > 0,
            y < 1,
            z > 1,
            z < 3,
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b',
            if ( B == 'b' ) x > 0
        )),
        data.frame(
            x = 0,
            y = 1,
            z = 2,
            A = 'a',
            B = 'b'
        )
    )
    expect_equal(nrow(v$corrections),1)
    # with revert
    v <- correctRounding(
        editset(expression(
            x + y == z,
            x >= 0,
            y > 0,
            y < 2,
            z > 1,
            z < 3,
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b',
            if ( B == 'b' ) x < 1
        )), 
        data.frame(
            x = 0,
            y = 1,
            z = 2,
            A = 'a',
            B = 'b'
        )
    )
    expect_equal(nrow(v$corrections),0)
})














