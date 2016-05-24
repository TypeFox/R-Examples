f.haptable.list <- function(object){
##
## APPLIES haptable TO EACH ELEMENT OF object
## object IS A list OF haplin OBJECTS
## 
#
#
## CREATE STANDARD haptable FOR EACH ELEMENT
.tab <- lapply(object, haptable)
#
## STACK TABLES (EFFICIENTLY)
.tab <- toDataFrame(.tab)
#
return(.tab)
}
