context("file system utility functions")

test_that('abs2rel works',{
  testans=abs2rel("/Volumes/JData/JPeople/Sebastian/images","/Volumes/JData/JPeople/")
  realans="Sebastian/images"
  expect_equal(testans,realans)
  
  testans=abs2rel(file.path("/Volumes/JData/JPeople/Sebastian", LETTERS), 
                  "/Volumes/JData/JPeople/")
  realans=file.path("Sebastian", LETTERS)
  expect_equal(testans, realans)
  
  testans=abs2rel("/Volumes/JData/JPeople/Sebastian/images","/Volumes/JData/JPeople")
  realans="Sebastian/images"
  expect_equal(testans,realans)
  
  expect_error(abs2rel("/some/other/path","/Volumes/JData/JPeople/",
                       StopIfNoCommonPath=TRUE))
  
  expect_error(abs2rel(c("/some/path", "/someother/path"), "/some",
                       StopIfNoCommonPath=TRUE), 
               "stempath.*/some.*is not present in .* /someother/path")
    
  expect_warning(testans<-abs2rel("/some/other/path","/Volumes/JData/JPeople/"))
  realans="/some/other/path"
  expect_equal(testans,realans)
})

test_that('touch works',{
  
  if(.Platform$OS.type!="unix")
    skip("touch only supported on unix platforms")
    
  tf=replicate(2,tempfile())
  on.exit(unlink(tf))
  
  expect_error(touch(tf[1],Create=FALSE),
      info="Throws exception if Create=FALSE and file does not exist")
  
  expect_true(touch(tf[1]))
  expect_true(file.exists(tf[1]),"touching a file without other argument creates it")
  expect_true(touch(tf[2]))
  
  t1=ISOdatetime(2001, 1, 1, 12, 12, 12)
  t2=ISOdatetime(2011, 1, 1, 12, 12, 12)
  expect_true(touch(tf[1],t1))
  expect_true(touch(tf[2],t2,Create=FALSE),
      "Check no error when Create=FALSE and target file exists")
  fis=file.info(tf)
  fis$mtime
  expect_equivalent(fis$mtime[1],t1,"Change to a specific time")
  expect_equivalent(fis$mtime[2],t2,"Change to a specific time")
  
  # Change modification time to that of a reference file, leaving access intact
  expect_true(touch(tf[2],reference = tf[1],timestoupdate = "modification"))
  fis2=file.info(tf[2])
  expect_equal(fis2$mtime,fis$mtime[1],
               info="Change mtime to that of a reference file")
  expect_equal(fis2$atime,fis$atime[2],info="Leave atime intact")
})

test_that("common_path works",{
  pp=c("/a/b/c/d", "/a/b/c")
  # no normalisation  
  expect_equal(common_path(pp), "/a/b/c")
  expect_equal(common_path(c("a","b")), "")
  expect_equal(common_path(c("","")), "")
  expect_equal(common_path(c("","/a")), "")
  expect_equal(common_path(c("/a","/b")), "/")
  expect_equal(common_path(c("/a/b/d","/b/c/d")), "/")
  expect_equal(common_path(c("/a/b/","/a/b")), "/a/b")
  expect_equal(common_path(c("/a/b/d","/a/b/c/d")), "/a/b/")
  
  # with normalisation
  np<-function(x) normalizePath(x, winslash = .Platform$file.sep, mustWork = F)
  expect_equal(common_path(pp, normalise = T), np("/a/b/c"))
  expect_equal(common_path(c("",""), normalise = T), "")
  expect_equal(common_path(c("","/a"), normalise = T), np(""))
  expect_equal(common_path(c("/a","/b"), normalise = T), np("/"))
  expect_equal(common_path(c("/a/b/","/a/b"), normalise = T), np("/a/b"))
  expect_equal(common_path(c("/a/b/d","/a/b/c/d"), normalise = T), np("/a/b/"))
  
  # expansion required
  expect_equal(common_path(c("~","~/"), normalise = F), "~")
  expect_equal(common_path(c("~","~/"), normalise = T), np("~"))
  expect_equal(common_path(c("~/a/b/d","~/a/b/c/d"), normalise = F), "~/a/b/")
  expect_equal(common_path(c("~/a/b/d","~/a/b/c/d"), normalise = T), np("~/a/b/"))
})

test_that("split_path works",{
  p1="/a/b/c"
  p1c=c("a", "b", "c")
  expect_equal(split_path(p1), p1c)
  expect_equal(split_path("a/b/c"), p1c)
  expect_equal(parts<-split_path(p1, include.fseps=TRUE),
               c("/", "a", "/", "b", "/","c"))
  # join parts back up again
  expect_equal(paste(parts, collapse = ""), p1)
  expect_equal(split_path("//a/b//c", include.fseps=TRUE, omit.duplicate.fseps=TRUE),
               c("/", "a", "/", "b", "/","c"))
  # Windows style
  expect_equal(split_path("C:\\a\\b\\c", fsep="\\"), c("C:", p1c))
  
  # errors
  expect_error(split_path(c(p1, p1)), 'one path')
})

test_that("file.swap works", {
  tf=paste0(tempfile(), 1:2)
  cat("1", file = tf[1])
  cat("2", file = tf[2])
  contents<-lapply(tf, scan, quiet=TRUE)
  file.swap(tf[1],tf[2])
  expect_equal(lapply(tf, scan, quiet=TRUE), rev(contents))
})
