# tests of file importing in various formats
require(networkDynamic)
require(testthat)
pkgpath<-path.package('networkDynamic') # for accessing test files
# below is a work around for local-specific sort order bug in R CMD check 
# that causes order of attributes in the list to be sorted differently
Sys.setlocale("LC_COLLATE", "C")


# ----- read.son import tests ---
# try reading basic file containing alpha ids and with start time but no end time
# also includes a fake cluster, just to mess with things
test_that("file with alphaIds",{
  expect_warning(alphaIdNet<-read.son(paste(pkgpath,'/extdata/alphaIdTest.son',sep='')),regexp="Unable to locate an arc column for 'EndTime'")
  expect_equal(as.data.frame(alphaIdNet)$onset,c(0,0))
  expect_equal(as.data.frame(alphaIdNet)$terminus,c(0,0))
  expect_equal(as.data.frame(alphaIdNet)$tail,c(1,3))
  expect_equal(as.data.frame(alphaIdNet)$head,c(2,2))
  # did it load pid
  expect_equal(alphaIdNet%n%'vertex.pid','vertex.names')
  expect_equal(alphaIdNet%v%'vertex.names',c('1','2','three'))
})

test_that('reading attributes as tea',{
  # try reading the extra var test
  varsNet<-read.son(paste(pkgpath,'/extdata/extraVarTest.son',sep=''))
  # "unchanging" should be static, others active
  
  expect_equal(list.vertex.attributes(varsNet),c("Happiness.active","Label.active","Unchanging","active","na","vertex.names"))
  # ArcWeight should be loaded in as a static attriube             
  expect_equal(list.edge.attributes(varsNet),c("ArcWeight","active","na"))
})

# check behavior of guess.TEA
test_that('guess.TEA flag works',{
  varsNet2<-read.son(paste(pkgpath,'/extdata/extraVarTest.son',sep=''),guess.TEA=FALSE)
  # unchanging should be active
  expect_equal(list.vertex.attributes(varsNet2),c("Happiness.active","Label.active","Unchanging.active","active","na","vertex.names"))
  # ArcWeight should be loaded in as a TEA attriube             
  expect_equal(list.edge.attributes(varsNet2),c("ArcWeight.active","active","na"))
})

# try reading mcfarland classroom with attributes
test_that('mcfarland classroom file',{
cls33<-read.son(paste(pkgpath,'/extdata/cls33_10_16_96.son',sep=''))
expect_equal(network.size(cls33),20)
expect_equal(range(get.change.times(cls33)),c(0,49))
expect_equal(list.vertex.attributes(cls33),c("BorderColor", "BorderWidth", "ColorName",   "Label", "NodeShape", "NodeSize","active","na","vertex.names"))

# check that vertex attributes were parsed correctly
expect_equal(get.vertex.attribute(cls33,'BorderColor'),rep('black',20))
expect_equal(get.vertex.attribute(cls33,'BorderWidth'),rep(1.5,20))
expect_equal(get.vertex.attribute(cls33,'ColorName'),c("gray","gray","gray","darkGray","darkGray","gray","orange" ,  "gray","darkGray","gray","darkGray","gray","darkGray","orange" ,  "gray"  ,   "gray"   ,  "gray"   ,  "gray" ,    "gray"     ,"gray" ))
expect_equal(get.vertex.attribute(cls33,'Label'),c(122658 ,129047, 129340, 119263, 122631, 144843,   1003, 113433, 131642, 139722, 139195, 133105, 116749,      3, 146757, 121402, 127265, 121829, 113140, 128065))
expect_equal(get.vertex.attribute(cls33,'NodeShape'),c("ellipse", "rect",    "rect",    "rect",    "ellipse", "rect",    "rect",    "rect",    "rect",    "ellipse", "ellipse", "rect", "rect",    "rect",    "ellipse", "rect",    "ellipse", "ellipse", "rect",    "ellipse"))
expect_equal(get.vertex.attribute(cls33,'NodeSize'),rep(5,20))

# check that edge attributes parsed correctly
expect_equal(list.edge.attributes(cls33),c("ArcWeight.active","ArcWidth.active","ColorName.active","active","na"))

# check that a few specific edges have the correct values
# and check edge attrs read correct values
eid<-get.edgeIDs(cls33,v=14, alter=12)
expect_equal(get.edge.activity(cls33,e=eid)[[1]][,1], c(0.125,1.167,4.667,9.964,21.0,37.737,41.553))
expect_equal(unlist(get.edge.attribute(cls33,'ColorName.active',unlist=FALSE)[[eid]][[1]]),c("blue",  "black", "red",   "black", "red",   "black", "black"))

eid<-get.edgeIDs(cls33,v=4, alter=5)
expect_equal(get.edge.activity(cls33,e=eid)[[1]][,1], c(2.5,7.286,9.321,15.964,16.393,17.679,18.964,19.393,27.977,30.628,35.368,36.421))

expect_equal(unlist(get.edge.attribute(cls33,'ColorName.active',unlist=FALSE)[[eid]][[1]]),c("blue",  "blue",  "black", "blue",  "blue",  "black", "blue",  "blue",  "black", "black", "black", "blue" ))
eid<-get.edgeIDs(cls33,v=17, alter=6)
expect_equal(get.edge.activity(cls33,e=eid)[[1]][,1], 44.0)
expect_equal(unlist(get.edge.attribute(cls33,'ColorName.active',unlist=FALSE)[[eid]][[1]]),c("blue"))
})


# test file with no edges
# no arcs header
testContent<-"NodeId\tLabel\tStartTime\tEndTime\n1\tA\t0\t1\n2\tB\t0\t1\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
expect_warning(testNet<-read.son(tempFile),regexp='Unable to locate an arc header line')
expect_true(network.size(testNet)==2 & network.edgecount(testNet)==0)

# has arcs header, but no arc rows
testContent<-"NodeId\tLabel\tStartTime\tEndTime\n1\tA\t0\t1\n2\tB\t0\t1\nFromId\tToId\tStartTime\tEndTime\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
testNet<-read.son(tempFile)
expect_true(network.size(testNet)==2 & network.edgecount(testNet)==0)

# test conversion of Label attribute to vertex.names
testContent<-"NodeId\tLabel\tStartTime\tEndTime\n1\tA\t0\t1\n2\tB\t0\t1\nFromId\tToId\tStartTime\tEndTime\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
testNet<-read.son(tempFile)
expect_equal(network.vertex.names(testNet),c('A','B'))

testContent<-"NodeId\tLabel\tStartTime\tEndTime\n1\tA\t0\t1\n1\tB\t1\t2\nFromId\tToId\tStartTime\tEndTime\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
testNet<-read.son(tempFile)
expect_equal(network.vertex.names(testNet),1) # Label is TEA so not copied to vertex names


# test file with multiple attribute rows doesn't create multiple edges. 
testContent<-"NodeId\tLabel\tStartTime\tEndTime\n1\tA\t0\t10\n1\tB\t0\t10\nFromId\tToId\tStartTime\tEndTime\tValue\n1\t2\t0\t1\tA\n1\t2\t1\t5\tB\n1\t2\t5\t10\tC\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
testNet<-read.son(tempFile)
expect_equal(network.edgecount(testNet),1) # only one edge created despite 3 rows
expect_equal(get.edge.attribute.active(testNet,'Value',at=5),"C")

# test file with missing header
testContent<-"1\tA\t0\t10\n1\tB\t0\t10\nFromId\tToId\tStartTime\tEndTime\tValue\n1\t2\t0\t1\tA\n1\t2\t1\t5\tB\n1\t2\t5\t10\tC\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
expect_error(read.son(tempFile),regexp='Unable to locate the header line')

# test missing start time
testContent<-"NodeId\tLabel\tEndTime\n1\tA\t10\n1\tB\t10\nFromId\tToId\tStartTime\tEndTime\tValue\n1\t2\t0\t1\tA\n1\t2\t1\t5\tB\n1\t2\t5\t10\tC\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
expect_error(read.son(tempFile),regexp="Unable to locate a node column for 'StartTime'")

# test missing end time (start time should be duplicated)
testContent<-"NodeId\tLabel\tStartTime\n1\tA\t10\n2\tB\t10\nFromId\tToId\tStartTime\tEndTime\tValue\n1\t2\t0\t1\tA\n1\t2\t1\t5\tB\n1\t2\t5\t10\tC\n"
tempFile<-tempfile()
cat(testContent,file=tempFile)
expect_message(testNet<-read.son(tempFile),regexp="Unable to locate a node column for 'EndTime'")
expect_equal(unlist(get.vertex.activity(testNet)),c(10,10,10,10))
