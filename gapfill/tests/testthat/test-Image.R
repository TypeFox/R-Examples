#require(testthat); library("gapfill", lib.loc = "../../../lib/")
context("test-Image")

test_that("Image",{
   expect_equal(class(Image(ndvi)), c("gg", "ggplot"))
   expect_equal(class(Image(ndvi[,,,1])), c("gg", "ggplot"))
   expect_equal(class(Image(ndvi[,,1,1])), c("gg", "ggplot"))
   expect_error(Image(ndvi[,1,1,1]),
                "2 <= length\\(dim\\(x\\)\\) && length\\(dim\\(x\\)\\) <= 4 is not TRUE")
   expect_equal(class(Image(unname(ndvi))), c("gg", "ggplot"))
   expect_equal(class(Image(ndvi, zlim = c(1,2))), c("gg", "ggplot"))
   expect_error(Image(ndvi, zlim = TRUE),
                "is.numeric\\(zlim\\) is not TRUE")
   expect_error(Image(ndvi, zlim = c(1)),
                "identical\\(length\\(zlim\\), 2L\\) is not TRUE")
})

