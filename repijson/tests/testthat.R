library(testthat)
library(repijson)


#TODO: Tests that we need to do
#1) construct an attribute of each type
#2) construct an event a) +date -location b) -date + location c) +date +location
#3) construct a record (with the events from above)
#4) construct an ejObject
#5) create an ejObject from a) dataframe b) spatialdataframe c) obkdata d) EpiJSON
#6) from and ejObject create a) dataframe b) spataildataframe c) obkdata e) EpiJSON

## and when we have finished the code
#7) add two datasets together ejObject + ejObject
#8) add attributes to a) metadata (ejObject) b) a record c) an event
#9) subsetting of ejObjects (should return a new ejObject with only the subset records)

test_check("repijson")
