
require(testthat)

context("Detect violated edits")

test_that("Numerical edit violations are detected",{
    expect_equivalent(
        violatedEdits(
            editmatrix(c( "x+3*y==2*z", "x==z")),
                data.frame( 
                   x = c(0,2,1),
                   y = c(0,0,1),
                   z = c(0,1,1))
        )[,,drop=FALSE],
        matrix(c(FALSE,FALSE,TRUE,FALSE,TRUE,FALSE),nrow=3)
    )
    # with a tolerance
    expect_equivalent(
        violatedEdits(
            editmatrix(c( "x+3*y==2*z", "x==z")),
            data.frame( 
                x = c(0,2,1),
                y = c(0,0,1),
                z = c(0,1,1)),
            tol=100
        )[,,drop=FALSE],
        matrix(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),nrow=3)
    )
})
				 
test_that("An empty editmatrix is always valid",{
    expect_equivalent(
        violatedEdits(
            editmatrix("x==1")[0],
                data.frame( 
                   x = c(0,2,1),
                   y = c(0,0,1),
                   z = c(0,1,1))
        )[,,drop=FALSE],
        matrix(nrow=3, ncol=0)
    )
    expect_equivalent(
        violatedEdits(
            editmatrix("x==1")[0],
                c(x = 0,y = 0,z = 0)
        )[,,drop=FALSE],
        matrix(nrow=1, ncol=0)
    )
})

test_that("violation of inequalities are detected",{
    expect_true(violatedEdits(editmatrix("x<0"),c(x=1))[1])
    expect_true(violatedEdits(editmatrix("x<0"),c(x=0))[1])
    expect_true(violatedEdits(editmatrix("x<0"),c(x=0,tol=1e-8))[1])
    expect_false(violatedEdits(editmatrix("x<0"),c(x=-1))[1])
    expect_true(all(
        violatedEdits(
        editmatrix("x < y"),
        data.frame(x=c(1,2),y=c(-1,3)))==c(TRUE,FALSE))
    )
})

test_that("summary.violatedEdits works when no edits are violated",{
    
   expect_equal(summary(violatedEdits(
        editmatrix("x+y==z"),
        c(x=1,y=2,z=3)
    )),NULL)

})

test_that("summary.violatedEdits works with NA's",{
    
   expect_equal(summary(violatedEdits(
        editmatrix("x+y==z"),
        c(x=1,y=2,z=NA)
    )),NULL)

})

test_that("NA's in data are handled correctly by violatededits.editmatrix",{
    expect_identical(
        violatedEdits(editmatrix(c('x==0','y==0')),c(x=NA,y=1))[,],
        c(num1=NA,num2=TRUE)
    )
    expect_identical(
        violatedEdits(editmatrix(c('x>0','y==0')),c(x=NA,y=1))[,],
        c(num1=NA,num2=TRUE)
    )
    expect_identical(
        violatedEdits(editmatrix(c('x>0','y==0','x+y>=1')),c(x=NA,y=1))[,],
        c(num1=NA,num2=TRUE,num3=NA)
    )
    expect_identical(
        violatedEdits(editmatrix(c('x==0','y==0')),data.frame(x=c(NA,1),y=c(1,NA)))[,],
        array(c(NA,TRUE,TRUE,NA),dim=c(2,2),dimnames=list(record=1:2,edit=c('num1','num2')))
    )

    expect_identical(
        violatedEdits(editmatrix(c('x==0','y==0')),c(x=NA,y=1),tol=0.1)[,],
        c(num1=NA,num2=TRUE)
    )
    expect_identical(
        violatedEdits(editmatrix(c('x>0','y==0')),c(x=NA,y=1),tol=0.1)[,],
        c(num1=NA,num2=TRUE)
    )
    expect_identical(
        violatedEdits(editmatrix(c('x>0','y==0','x+y>=1')),c(x=NA,y=1),tol=0.1)[,],
        c(num1=NA,num2=TRUE,num3=NA)
    )
    expect_identical(
        violatedEdits(editmatrix(c('x==0','y==0')),data.frame(x=c(NA,1),y=c(1,NA)),tol=0.1)[,],
        array(c(NA,TRUE,TRUE,NA),dim=c(2,2),dimnames=list(record=1:2,edit=c('num1','num2')))
    )

})

test_that("categorical edit violations are detected",{
    E <-  editarray(c(
        "gender %in% c('male','female')",
        "pregnant %in% c(TRUE, FALSE)",
        "if( gender == 'male' ) !pregnant"))    
    dat <- data.frame(
        gender=c('male','male','female','cylon'), 
        pregnant=c(TRUE,FALSE,TRUE,TRUE)
    )
    expect_equivalent(
        violatedEdits(E,dat)[,,drop=FALSE], 
        matrix(c(
            FALSE, FALSE,  TRUE,
            FALSE, FALSE, FALSE,
            FALSE, FALSE, FALSE,
            TRUE,  FALSE, FALSE),byrow=TRUE,nrow=4)
    )
    expect_equivalent(
        violatedEdits(E,dat,datamodel=FALSE)[,,drop=FALSE], 
        matrix(c(
             TRUE,
            FALSE,
            FALSE,
            FALSE),byrow=TRUE,nrow=4)
    )
})



test_that("Empty editarray is handled correctly",{
    expect_equivalent(
        violatedEdits(editarray(expression()),data.frame(x="a",y="b"))[,,drop=FALSE],
        array(logical(0),dim=c(1,0))
    )
    expect_equivalent(
        violatedEdits(editarray(expression()),data.frame(x=c("a","a"),y=c("b","b")))[,,drop=FALSE],
        array(logical(0),dim=c(2,0))
    )

})

#for ( d in dir("../../../pkg/R",full.names=TRUE) ) dmp <- source(d)
test_that("Pure numerical editset",{
    expect_true(
        violatedEdits(
            E = editset("x + y == z"),
            dat = data.frame(x=1,y=1,z=1)
        )[1,1]
    )
})


