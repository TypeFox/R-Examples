
context("Echelon")


test_that("Matrix reduces to echelon form",{
    expect_equal(
        round(echelon(
            matrix(c(
                1,3,1,4,
                2,7,3,-9,
                1,5,3,1,
                1,2,0,8), byrow=TRUE, nrow=4
            )
        )),
        matrix(c(
            1,0,-2,0,
            0,1,1,0,
            0,0,0,1,
            0,0,0,0),byrow=TRUE,nrow=4
        )
    )
    expect_equal(
        round(echelon(
            matrix(c(
                2,1,-1,8,
               -3,-1,2,-11,
               -2,1,2,-3), byrow=TRUE, nrow=3
            )
        )),
        matrix(c(
            1,0,0,2,
            0,1,0,3,
            0,0,1,-1),
            byrow=TRUE,nrow=3
        )
    )
})


