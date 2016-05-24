remnodes<-function(left,right,list){
#Removes branches from a tree
#
#list is a vector, branches whose root is mentioned in the list
#  will be removed
#
num<-length(list)
for (i in 1:num){
  left[list[i]]<-0
  right[list[i]]<-0
}
return(list(left=left,right=right))
}


