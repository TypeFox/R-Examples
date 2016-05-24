context("Checking list2df")

test_that("list2df gived intended output",{
	
    lst1 <- list(x=c("foo", "bar"), y=1:5)
    expected1 <- structure(list(X1 = c("foo", "bar", "1", "2", "3", "4", "5"), 
            X2 = c("x", "x", "y", "y", "y", "y", "y")), .Names = c("X1", 
        "X2"), row.names = c(NA, -7L), class = "data.frame")
    
    expect_equivalent(list2df(lst1), expected1)
    
    lst2 <- list(a=c("hello", "everybody"), b = mtcars[1:6, 1])
    expected2 <- structure(list(`col 1` = c("hello", "everybody", "21", "21",     
         "22.8", "21.4", "18.7", "18.1"), `col 2` = c("a", "a", "b", "b",      
         "b", "b", "b", "b")), .Names = c("col 1", "col 2"), row.names = c(NA, 
         -8L), class = "data.frame") 
    
    expect_equivalent(list2df(lst2, "col 1", "col 2"), expected2)
  
})


test_that("matrix2df gived intended output",{
	
    expect_true(is.data.frame(matrix2df(mtcars)))
    expect_true(all(matrix2df(mtcars)[, 1] == rownames(mtcars)))
    expect_true(all(dim(matrix2df(mtcars)) == c(32, 12)))

    expected4 <- structure(list(var1 = c("1", "2", "3"), `1` = 1:3, `2` = 4:6, 
            `3` = 7:9), .Names = c("var1", "1", "2", "3"), row.names = c(NA, 
        -3L), class = "data.frame")

    expect_equivalent(matrix2df(matrix(1:9, ncol=3)), expected4)
})


test_that("vect2df gived intended output",{
	
    expected5 <- structure(list(X1 = structure(1:10, .Label = c("x01", "x02", 
    "x03", "x04", "x05", "x06", "x07", "x08", "x09", "x10"), class = "factor"), 
        X2 = 1:10), .Names = c("X1", "X2"), row.names = c(NA, -10L
    ), class = "data.frame")

    expect_equivalent(vect2df(1:10), expected5)
    expect_true(is.data.frame(vect2df(1:10)))
    expect_equivalent(vect2df(1:10)[, 1], factor(paste0("x", pad(1:10))))
})

test_that("list_df2df gived intended output",{
	
    expect_equivalent(
        list_df2df(list(mtcars, mtcars)),
        data.frame(X1=rep(paste0("L", 1:2), each=nrow(mtcars)), 
            rbind(mtcars, mtcars), row.names=NULL, stringsAsFactors = FALSE)
    )
})

test_that("list_vect2df gived intended output",{
	
    L1 <- list(a=1:10, b=1:6, c=5:8)
    expect_true(is.data.frame(list_vect2df(L1)))

    expected6 <- structure(list(X1 = c("a", "a", "a", "a", "a", "a", "a", "a", 
    "a", "a", "b", "b", "b", "b", "b", "b", "c", "c", "c", "c"), 
        X2 = structure(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 
        11L, 12L, 13L, 14L, 15L, 16L, 15L, 16L, 17L, 18L), .Label = c("x01", 
        "x02", "x03", "x04", "x05", "x06", "x07", "x08", "x09", "x10", 
        "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"), class = "factor"), 
        X3 = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 1L, 2L, 3L, 
        4L, 5L, 6L, 5L, 6L, 7L, 8L)), .Names = c("X1", "X2", "X3"
    ), row.names = c(NA, -20L), class = "data.frame")

    expect_equivalent(list_vect2df(L1), expected6)
	
    L2 <- list(
        months=setNames(1:12, month.abb),
        numbers=1:6,
        states=setNames(factor(state.name[1:4]), state.abb[1:4])
    )

    outcome <- list_vect2df(L2, order = FALSE)
    expect_true(is.data.frame(outcome))
    expect_true(all(dim(outcome)==c(22, 3)))	
})

test_that("counts2list gived intended output",{
	
    expect_true(is.list(counts2list(mtcars[1:10, 8:10])))
    expect_true(length(counts2list(mtcars[1:10, 8:10])) == 10)
    expect_true(all(names(counts2list(mtcars[1:10, 8:10])) == rownames(mtcars)[1:10]))
    expect_true(all(unique(unlist(counts2list(mtcars[1:10, 8:10]))) %in% colnames(mtcars)[8:10]))

})

test_that("vect2list gived intended output",{
	
    i <- vect2list(LETTERS[1:10])
    j <- vect2list(LETTERS[1:10], numbered.names = TRUE)
    x <- setNames(LETTERS[1:4], paste0("Element_", 1:4))
    k <- vect2list(x)
    l <- vect2list(x, FALSE)
    m <- vect2list(x, FALSE, TRUE)
    
    expect_true(all(sapply(list(i, j, k, l, m), is.list)))
    expect_true(all(names(i) == i))
    expect_false(all(names(j) == j))
    expect_true(all(names(j) == pad(1:10)))
    expect_true(all(names(k) == names(x)))
    expect_false(all(names(l) == names(x)))
    expect_true(all(names(l) == LETTERS[1:4]))
    expect_false(all(names(m) == names(x)))


    expected7 <- structure(list(A = "A", B = "B", C = "C", 
        D = "D"), .Names = c("A", "B", "C", "D"))


    expect_equivalent(l, expected7)


    expected8 <- structure(list(`1` = "A", `2` = "B", `3` = "C", 
        `4` = "D"), .Names = c("1", "2", "3", "4"))

    expect_equivalent(m, expected8)

})

test_that("df2matrix gived intended output",{
	
    cnts <- structure(list(month = c("January", "February", "March", "April", 
        "May", "June"), X1 = c(1L, 0L, 1L, 1L, 2L, 2L), X2 = c(2L, 1L, 
        2L, 1L, 2L, 2L), X3 = c(0L, 1L, 2L, 2L, 0L, 2L)), .Names = c("month", 
        "X1", "X2", "X3"), row.names = c(NA, -6L), class = "data.frame")
    
    df2matrix(cnts)
    m <- df2matrix(cnts)
    expect_equivalent(rownames(m), cnts[, 1])
    expect_true(is.matrix(m))
    expect_true(mode(m) == "numeric")

    m2 <- df2matrix(cnts, 2)
    expect_equivalent(as.integer(rownames(m2)), cnts[, 2])
    expect_true(is.matrix(m2))
    expect_true(mode(m2) == "character")

    m3 <- df2matrix(cnts, "X2")
    expect_equivalent(as.integer(rownames(m3)), cnts[, 3])
    expect_true(is.matrix(m3))
    expect_true(mode(m3) == "character")
	
})

test_that("matrix2long gived intended output",{
    
    exmat <- structure(list(cols = c("1", "1", "1", "2", "2", "2", "3", "3", 
        "3"), rows = c("1", "2", "3", "1", "2", "3", "1", "2", "3"), 
            vals = 1:9), .Names = c("cols", "rows", "vals"), row.names = c(NA, 
        -9L), class = "data.frame")
    
    mat <- matrix(1:9, ncol=3)
    expect_equivalent(matrix2long(mat), exmat)
    
})
