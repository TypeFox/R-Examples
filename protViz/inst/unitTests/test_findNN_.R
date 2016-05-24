#R


test_findNN_ <-
function(){
# TEST 1
    # the drawback of the findNN implementation. Use findNN_!
    checkEqualsNumeric(findNN(3.5, 1:5), findNN(3.5, 1:6), tolerance=1.0)

    checkEqualsNumeric(findNN_(3.5, 1:5), findNN_(3.5, 1:6), tolerance=0.0)

# TEST 2
    DB<-sort(rnorm(100, mean=100, sd=10))
    checkEqualsNumeric(unique(DB[findNN(DB,DB)] - DB), 0, tolerance=0.0)
    checkEqualsNumeric(unique(DB[findNN_(DB,DB)] - DB), 0, tolerance=0.0)

# TEST 3 -- testing lower and upper index
    DB <- seq(-1,1,length=101)
    query <- c(-1000,0,0.001,10,10000)
    result <- c(1, 51, 51, 101, 101)

    checkEqualsNumeric(findNN(q=query,  vec=DB), result)
    checkEqualsNumeric(findNN_(q=query,  vec=DB), result)
    checkEqualsNumeric(findNN(q=query,  vec=DB), findNN_(q=query,  vec=DB))


}

test_findNN_()
