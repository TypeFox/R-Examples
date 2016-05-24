library(kinship2)
source("ped2trio.r")

# Here is a function that recursively prints the descendants of a subject
# using eval to get the trio of each child of the current subject
find.descendants = function(trio,trio.list)
# arguments:
# trio: the subject for which we want to print the descendants
# trio.list: a list of trios representing the pedigree
{
for (i in 1:length(trio$offspring))
  {
  print(trio$offspring[[i]])
  # If the child is a parent in a trio
  if( length(grep("trio",trio$offspring[[i]])>0) )
    {
    # Here is where I get the next trio
    child.trio = eval(parse(text=trio$offspring[[i]]),envir=trio.list)
    recurtrio (child.trio,trio.list)
    }
}
}

# Creating a list of trios representing a pedigree
test = ped2trio(ped.list[[47]])
# Testing the function
find.descendants(test$trio4,test)
