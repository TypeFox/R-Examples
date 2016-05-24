library(repijson)

context("Constructors")

#data for simple start tests
#later would be good maybe to read this data in from the package
dF <- data.frame(id=c("A","B","3D"),
                 name=c("tom","andy","ellie"),
                 dob=c("1984-03-14","1985-11-13","1987-06-16"),
                 gender=c("male","male","female"),
                 rec1contact=c(2,1,5),
                 rec1date=c("2014-12-28","2014-12-29","2015-01-03"),
                 rec1risk=c("high","high","low"),  
                 rec1temp=c(39,41,41),
                 rec2contact=c(4,1,1),
                 rec2date=c("2015-01-02","2015-01-12","2015-01-09"),
                 rec2risk=c("high","low","high"),stringsAsFactors=FALSE)

#create objects here for various tests after

#create attribute objects
attributeTst1 <- create_ejAttribute(name="name", type="string", value=dF$name[1])
attributeTst2 <- create_ejAttribute(name="name", type="string", value=dF$name[2])
#create metadata object
metadataTst <- create_ejMetadata( attributes=list(attributeTst1, attributeTst2) )


test_that("constructing an Attribute object doesn't modify value", {
  
  #expect the value isn't changed
  expect_equal( attributeTst1$value,
                dF$name[1] )  
  
})


test_that("constructing a Metadata object doesn't modify value", {
  
  #expect the value isn't change in 2 attribute objects
  expect_equal( metadataTst[[1]]$value,
                dF$name[1] ) 
  expect_equal( metadataTst[[2]]$value,
                dF$name[2] )  
  
})





