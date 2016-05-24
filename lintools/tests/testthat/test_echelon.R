
context("Echelon")


test_that("Matrix reduces to echelon form",{
    expect_equivalent(
        echelon(
            A = matrix(c(
                1,3,1,
                2,7,3,
                1,5,3,
                1,2,0), byrow=TRUE, nrow=4)
            , b = c(4,-9,1,8)
            , neq=4
            )
        ,
        list(
          A = matrix(c(
                1,0,-2,
                0,1, 1,
                0,0, 0),byrow=TRUE,nrow=3)
          , b = c(0,0,1)
          , neq = 3
          , nleq = 0
        )
        
    )
    expect_equivalent(
        echelon(
            A = matrix(c(
                2,1,-1,
               -3,-1,2,
               -2,1,2), byrow=TRUE, nrow=3)
            , b = c(8,-11,-3)
            , neq=3
        )
        , list(
          A = diag(rep(1,3))
          , b = c(2,3,-1)
          , neq=3
          , nleq=0
        )
    )
    # with an inequality present
    expect_equivalent(
        echelon(
            A = matrix(c(
                2,1,-1,
               -3,-1,2,
               -2,1,2,
                1,2,3), byrow=TRUE, nrow=4)
            , b = c(8,-11,-3,0)
            , neq=3
        )
        , list(
          A = rbind(diag(rep(1,3)),c(1,2,3))
          , b = c(2,3,-1,0)
          , neq=3
          , nleq=0
        )
    )
})


