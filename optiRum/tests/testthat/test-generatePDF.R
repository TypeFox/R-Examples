context("generatePDF")

# initial values
wd <- getwd()
basepath <- wd

testpath <- file.path(basepath, "temp")
dir.create(testpath)

test_that("generatePDF - correct behaviour, DATED=FALSE,CLEANUP=TRUE", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = FALSE)
    
    expect_true(file.exists(file.path(testpath, "basic.pdf")))
    expect_true(file.exists(file.path(testpath, "basic.tex")))
    expect_false(file.exists(file.path(testpath, "basic.log")))
    expect_false(file.exists(file.path(testpath, "basic.aux")))
    expect_false(file.exists(file.path(testpath, "basic.out")))
    expect_false(file.exists(file.path(testpath, "basic.toc")))
})

test_that("generatePDF - correct behaviour, DATED=FALSE,CLEANUP=TRUE, compiler=xelatex", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = FALSE)
    
    expect_true(file.exists(file.path(testpath, "basic.pdf")))
    expect_true(file.exists(file.path(testpath, "basic.tex")))
    expect_false(file.exists(file.path(testpath, "basic.log")))
    expect_false(file.exists(file.path(testpath, "basic.aux")))
    expect_false(file.exists(file.path(testpath, "basic.out")))
    expect_false(file.exists(file.path(testpath, "basic.toc")))
})

test_that("generatePDF - correct behaviour, quick generate", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "datatabletest")
    
    expect_true(file.exists(file.path(basepath, "datatabletest.pdf")))
    expect_true(file.exists(file.path(basepath, "datatabletest.tex")))
    expect_false(file.exists(file.path(basepath, "datatabletest.log")))
    expect_false(file.exists(file.path(basepath, "datatabletest.aux")))
    expect_false(file.exists(file.path(basepath, "datatabletest.out")))
    expect_false(file.exists(file.path(basepath, "datatabletest.toc")))
    
    dtfiles <- dir(path = basepath, pattern = "datatabletest*", full.names = TRUE)
    dtfiles <- dtfiles[!dtfiles %in% file.path(basepath, c("datatabletest.Rnw", "datatabletest.pdf"))]
    file.remove(dtfiles)
    
})

test_that("generatePDF - correct behaviour, quick & quiet generate", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "datatabletest", QUIET = TRUE)
    
    expect_true(file.exists(file.path(basepath, "datatabletest.pdf")))
    expect_true(file.exists(file.path(basepath, "datatabletest.tex")))
    expect_false(file.exists(file.path(basepath, "datatabletest.log")))
    expect_false(file.exists(file.path(basepath, "datatabletest.aux")))
    expect_false(file.exists(file.path(basepath, "datatabletest.out")))
    expect_false(file.exists(file.path(basepath, "datatabletest.toc")))
    
    dtfiles <- dir(path = basepath, pattern = "datatabletest*", full.names = TRUE)
    dtfiles <- dtfiles[!dtfiles %in% file.path(basepath, c("datatabletest.Rnw", "datatabletest.pdf"))]
    file.remove(dtfiles)
    
})

test_that("generatePDF - correct behaviour, DATED=FALSE,CLEANUP=TRUE, compiler=xelatex", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "datatabletest", destpath = testpath, destname = "datatabletest", DATED = FALSE, compiler = "xelatex")
    
    expect_true(file.exists(file.path(testpath, "datatabletest.pdf")))
    expect_true(file.exists(file.path(testpath, "datatabletest.tex")))
    expect_false(file.exists(file.path(testpath, "datatabletest.log")))
    expect_false(file.exists(file.path(testpath, "datatabletest.aux")))
    expect_false(file.exists(file.path(testpath, "datatabletest.out")))
    expect_false(file.exists(file.path(testpath, "datatabletest.toc")))
})

test_that("generatePDF - correct behaviour, DATED=TRUE,CLEANUP=TRUE", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = TRUE)
    
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".pdf"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".tex"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".log"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".aux"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".out"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".toc"))))
})


test_that("generatePDF - correct behaviour, DATED=FALSE,CLEANUP=FALSE", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = FALSE, CLEANUP = FALSE)
    
    expect_true(file.exists(file.path(testpath, "basic.pdf")))
    expect_true(file.exists(file.path(testpath, "basic.tex")))
    expect_true(file.exists(file.path(testpath, "basic.log")))
    expect_true(file.exists(file.path(testpath, "basic.aux")))
    expect_true(file.exists(file.path(testpath, "basic.toc")))
})


test_that("generatePDF - correct behaviour, DATED=TRUE,CLEANUP=FALSE", {
    skip_on_cran() 
    file.remove(dir(testpath, full.names = TRUE))
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = TRUE, CLEANUP = FALSE)
    
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".pdf"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".tex"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".log"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".aux"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".toc"))))
})


test_that("generatePDF - multiple calls still performs as expected correct behaviour,CLEANUP=TRUE", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = FALSE)
    
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = TRUE)
    
    expect_true(file.exists(file.path(testpath, "basic.pdf")))
    expect_true(file.exists(file.path(testpath, "basic.tex")))
    expect_false(file.exists(file.path(testpath, "basic.log")))
    expect_false(file.exists(file.path(testpath, "basic.aux")))
    expect_false(file.exists(file.path(testpath, "basic.toc")))
    
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".pdf"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".tex"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".log"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".aux"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".toc"))))
})

test_that("generatePDF - multiple calls still performs as expected correct behaviour,CLEANUP=FALSE", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = FALSE, CLEANUP = FALSE)
    
    generatePDF(srcpath = basepath, srcname = "basic", destpath = testpath, destname = "basic", DATED = TRUE, CLEANUP = FALSE)
    
    expect_true(file.exists(file.path(testpath, "basic.pdf")))
    expect_true(file.exists(file.path(testpath, "basic.tex")))
    expect_true(file.exists(file.path(testpath, "basic.log")))
    expect_true(file.exists(file.path(testpath, "basic.aux")))
    expect_true(file.exists(file.path(testpath, "basic.toc")))
    
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".pdf"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".tex"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".log"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".aux"))))
    expect_true(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".toc"))))
})


test_that("generatePDF - errors - DATED=FALSE, source file does not exist", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    expect_error(generatePDF(srcpath = basepath, srcname = "basic1", destpath = testpath, destname = "basic1", DATED = FALSE))
    
    expect_false(file.exists(file.path(testpath, "basic1.pdf")))
    expect_false(file.exists(file.path(testpath, "basic1.tex")))
    expect_false(file.exists(file.path(testpath, "basic1.log")))
    expect_false(file.exists(file.path(testpath, "basic1.aux")))
    expect_false(file.exists(file.path(testpath, "basic.toc")))
})


test_that("generatePDF - errors - DATED=TRUE, source file does not exist", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    expect_error(generatePDF(srcpath = basepath, srcname = "basic1", destpath = testpath, destname = "basic1", DATED = TRUE))
    
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".pdf"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".tex"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".log"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".aux"))))
    expect_false(file.exists(file.path(testpath, paste0("basic", format(Sys.Date(), "%Y%m%d"), ".toc"))))
})


test_that("generatePDF - errors - missing inputs", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    expect_error(generatePDF(srcpath = basepath, destpath = testpath, destname = "basic"))
    
    expect_false(file.exists(file.path(testpath, "basic1.pdf")))
    expect_false(file.exists(file.path(testpath, "basic1.tex")))
    expect_false(file.exists(file.path(testpath, "basic1.log")))
    expect_false(file.exists(file.path(testpath, "basic1.aux")))
    expect_false(file.exists(file.path(testpath, "basic.toc")))
})

test_that("generatePDF - errors - missing inputs", {
    skip_on_cran()
    file.remove(dir(testpath, full.names = TRUE))
    expect_error(generatePDF())
    
    expect_false(file.exists(file.path(testpath, "basic1.pdf")))
    expect_false(file.exists(file.path(testpath, "basic1.tex")))
    expect_false(file.exists(file.path(testpath, "basic1.log")))
    expect_false(file.exists(file.path(testpath, "basic1.aux")))
    expect_false(file.exists(file.path(testpath, "basic.toc")))
})

# file.remove(testpath) 
