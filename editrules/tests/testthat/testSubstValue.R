
context("Value substitution")

test_that("value substitution for numerical data",{
    E <- editmatrix(c("x+y==z","x-u==1"))
    expect_true(
        sum(abs(
            getAb(substValue(E,c('x','y'),c(1,2)))-
            matrix(c(
                0, 0,-1, 0,-3,
                0, 0, 0,-1, 0),nrow=2,byrow=TRUE)
        )) == 0
    )
    expect_true(
        sum(abs(
            getAb(substValue(E,c('x','y'),c(1,2),reduce=TRUE))-
            matrix(c(
                -1, 0,-3,
                 0,-1, 0),nrow=2,byrow=TRUE)
        )) == 0
    )
})


test_that("value substitution for categorical data",{

    E <- editarray(c(
        "v %in% letters[1:3]",
        "w %in% letters[4:6]",
        "x %in% letters[7:9]",
        "if ( v=='a' ) w !='d'",
        "if ( w %in% c('d','e')) x != 'g'"
    ))
    expect_true(
        all(
            getArr(substValue(E,c('v','w'),c('a','d')))
            ==
            matrix(c(
                TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,  TRUE,  TRUE,
                TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE
            ),nrow=c(2),byrow=TRUE)
        )
    )
    expect_true(
        all(
            getArr(substValue(E,c('v','w'),c('a','d'),reduce=TRUE))
            ==
            matrix(c(
                TRUE, TRUE, TRUE,  TRUE,  TRUE,
                TRUE, TRUE, TRUE, FALSE, FALSE
            ),nrow=c(2),byrow=TRUE)
        )
    )
    

})

test_that("value substitution works for boolean data",{
    E <- editarray(c(
    "g %in% c('m','f')",
    "p %in% c(TRUE, FALSE)"))
   expect_true(nrow(getArr(substValue(E,'p',TRUE)))==0)
})

test_that("Bug reported by Sander Scholtus is fixed",{
  # this used to cause a crash after upgrading to R>=3.x.x
  E1 <- editset(expression(
    A %in% letters[1:2],
    B %in% letters[2:3],
    if ( (A == 'a') ) x > y
  ))
  E1 <- substValue(E1, var = "x", val = 1)
})


