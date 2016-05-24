
context('Check datamodel')

test_that("checkDatamodel.editmatrix works",{
    expect_true(
        checkDatamodel(
            editmatrix("x > 0"),
            data.frame(x=-1)
        )$adapt[1,1]
    )
    # test with NA
    expect_true(
        checkDatamodel(
            editmatrix("x > 0"),
            data.frame(x=NA)
        )$adapt[1,1]
    )
    # test with valid entry
    expect_false(
        checkDatamodel(
            editmatrix("x < 0"),
            data.frame(x=-1)
        )$adapt[1,1]
    )
    # test with no single-variable edits
    expect_false(
        checkDatamodel(
            editmatrix("x + y ==  1"),
            data.frame(x=-1,y=2)
        )$adapt[1,1]
    )
    
})


test_that('checkDatamodel.editarray works',{
    # dat has column also in E
    expect_equivalent(    
        checkDatamodel(
            E = editarray('x %in% 1:2'),
            dat = data.frame(x=1:3)
        )$adapt[,1],
        c(FALSE,FALSE,TRUE)
    )
    # dat has column not specified by E
    expect_equivalent(    
        checkDatamodel(
            E = editarray('x %in% 1:2'),
            dat = data.frame(y=1:3,x=1:3)
        )$adapt[,1],
        c(FALSE,FALSE,FALSE)
    )
    # dat misses a variable, specified in E
    expect_error(    
        checkDatamodel(
            E = editarray('x %in% 1:2'),
            dat = data.frame(y=1:3)
        )
    )
    # dat computes correct weights
    expect_equivalent(
        checkDatamodel(
            E = editarray(c('x %in% 1:2','y %in% c("a","b")')),
            dat = data.frame(x=1:4,y=c('a','c','b','c'))
        )$status$weight,
        c(0,1,1,2)
    )
})


test_that("checkDatamodel.editset works with pure numerical edits",{
    v <- checkDatamodel(
        editset(expression(
            x > 0,
            x + y == 3
        )),
        data.frame(x=c(-1,2),y=c(1,1))
    )
    expect_equivalent(v$adapt,array(c(TRUE,FALSE,FALSE,FALSE),dim=c(2,2)))
})


test_that("checkDatamodel.editset works with pure categorical edits",{
    v <- checkDatamodel(
        editset(expression(
            A %in% letters[1:3],
            B %in% 1:3
        )),
        data.frame(A=c('q','c'),B=c(1,10))
    )
    expect_equivalent(v$adapt,array(c(TRUE,FALSE,FALSE,TRUE),dim=c(2,2)))
})



test_that("checkDatamodel.editset works with conditional numerical edits",{
    v <- checkDatamodel(
        editset(expression(
            x + y == 3,
            x > 0,
            if ( x > 2 ) y < 1,
            v %in% letters[1:3])),
        data.frame(
            x = c(-1,3),
            y = c(0, 1),
            v = c("a","out-of-range")
    ))
    expect_equivalent(
        v$adapt,
        array(c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE),dim=c(2,3))
    )
})


test_that("checkDatamodel.editset works with conditional categorical/numerical edits",{
    v <- checkDatamodel(
        editset(expression(
            x + y == 3,
            x > 0,
            v %in% letters[1:3],
            if ( v == 'a' ) y > 0
            )),
        data.frame(
            x = c(-1,3),
            y = c(0, 1),
            v = c("a","out-of-range")
    ))

})




