require(editrules)

context("Correction of sign errors")

test_that("correctSigns",{
    expect_identical(
        correctSigns(
            editmatrix(c("x>0","y>0","z>0","x+y==z")),
            data.frame(x=-1,y=1,z=2))$corrected,
        data.frame(x=1,y=1,z=2)
    )
    expect_identical(
        correctSigns( # nothing can be done without specifying swap
            editmatrix(c("x>0","y>0","z>0","x+y==z")),
            data.frame(x=1,y=2,z=1))$corrected,
       data.frame(x=1,y=2,z=1)
    )
    expect_identical(
        correctSigns( # specify swap
            editmatrix(c("x>0","y>0","z>0","x+y==z")),
            data.frame(x=1,y=2,z=1),
            swap=list("z","y"))$corrected,
       data.frame(x=1,y=1,z=2)
    )
   expect_identical(
      correctSigns(    
         editmatrix(expression(
            y >= 0,
            z >= 0,
            x + y == z
         )),
         data.frame(x=1,y=1,z=0),
         swap=list(c('z','y')),
         flip=c()
      )$corrected,
      data.frame(x=1,y=0,z=1)
   )
    expect_identical( # flip has variable not in E
        correctSigns(
            editmatrix(c("y>0","z>0","x+y==z")),
            data.frame(x=1,y=2,z=1,u=3),
            flip=c("x","u"))$corrected,
       data.frame(x=-1,y=2,z=1,u=3)
    )
    expect_identical( # swap has variable not in E
        correctSigns(
            editmatrix(c("x>0","y>0","z>0","x-y==z")),
            data.frame(x=1,y=2,z=1,u=3),
            flip=c(),swap=list(c("x","y"),c("u","z")))$corrected,
       data.frame(x=2,y=1,z=1,u=3)
    )
    
}) 


test_that("correctSigns.editset works with pure numerical",{
    
    v <- correctSigns(
        editset(expression(
        x + y == z)),
        data.frame(
            x=-1,
            y=1,
            z=2
        )
    )

    w <- correctSigns(
        editset(expression(
        x + y == z))$num,
        data.frame(
            x=-1,
            y=1,
            z=2
        )
    )
    expect_equal(v$corrected,w$corrected)
    expect_equal(v$corrections,w$corrections)

})

test_that("correctSigns.editset works with pure numerical",{

    v <- correctTypos(editset(expression(
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a' ) B == 'b')
        ),
        data.frame(
            A = 'a',
            B = NA
    ))

    expect_equal(nrow(v$corrections),0)
})


test_that("correctSigns.editset works with unconnected numerical/categorical",{

    v <- correctSigns(
        editset(expression(
            x + y == z,
            x >= 0,
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b'
        )),
        data.frame(
            x = -1,
            y = 1,
            z = 2,
            A = 'a',
            B = NA
        )
    )
    w <- correctSigns(
        editset(expression(
            x + y == z,
            x >= 0,
            A %in% c('a','b'),
            B %in% c('c','d'),
            if ( A == 'a' ) B == 'b'
        ))$num,
        data.frame(
            x = -1,
            y = 1,
            z = 2,
            A = 'a',
            B = NA
        )
    )
    expect_equal(v$corrected, w$corrected)
    expect_equal(v$corrections,w$corrections)
})



test_that("correctSings.editset works with connected numerical/categorical",{


## without revert
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
        x = -1,
        y = 1,
        z = 2,
        A = 'a',
        B = NA
    )
    v <- correctSigns(E,x)
    w <- correctSigns(E$num,x)
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
        if ( B == 'b' ) x < 1
    ))

    x <- data.frame(
        x = -1,
        y = 1,
        z = 2,
        A = 'a',
        B = 'b'
    )

    v <- correctSigns(E,x)
    expect_equal(nrow(v$corrections),0)

})

# fixes issue nr XXXX
test_that("constant vector is accounted for",{
   expect_equal(correctSigns(editmatrix("x + y == 2"), data.frame(x=-1,y=1))$corrected$x[1],1)
})





