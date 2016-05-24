library(RNetLogo)

# path to NetLogo installation folder
# PLEASE FILL IN THE PATH TO YOUR NetLogo INSTALLATION FOLDER !!!!
path.to.NetLogo <- "C:/Program Files/NetLogo 5.3/app"

# start NetLogo, if it is not started yet
if (!exists("nl.test1", -1)) 
{
  # an reference name
  nl.test1 <- "nl.test1"
  NLStart(path.to.NetLogo, gui=TRUE, nl.obj=nl.test1)
}

# load a sample model
# PLEASE FILL IN THE PATH TO THE SAMPLE NetLogo FILE !!!!
path.to.local.file <- 'C:/Users/jthiele/Documents/R/win-library/3.2/RNetLogo/examples/code_samples/9-NLDfToList/dftest.nlogo'
NLLoadModel(path.to.local.file, nl.obj=nl.test1)

# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# create a data.frame with four sample vectors
list1 <- c(1,2,3,4,5)
list2 <- c(6,7,8,9,10)
list3 <- c('test1','test2','test3','test4','test5')
list4 <- c(TRUE,FALSE,TRUE,FALSE,TRUE)

sample.df <- data.frame(list1,list2,list3,list4)

# fill automaticlly the four NetLogo list with the names "list1", "list2", "list3" and "list4" (defined in "globals")
# with the values of the corresponding columns of the data.frame 
NLDfToList(sample.df, nl.obj=nl.test1)

# show content of the NetLogo list in NetLogo Command Center
NLCommand("show-lists", nl.obj=nl.test1)


# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
