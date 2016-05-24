apply.tbin <-
function(idv,train.output,test){
 
         xyz.pbin = function(tab,y,z) {
 trname = row.names(t(tab))
 tname = row.names(tab)
 testname = levels(y)
 new_levels(testname,tname,z)
 for (j in 1:length(tname)) {
 for (k in 1:length(trname)) {
if (tab[j,k]>0){
     levels(y)[levels(y)==tname[j]] <- trname[k]
 }
}
}
return(y) 
}
# error handling function for new levels in test data
new_levels <- function(a,b,z){
if( !(all(a %in% b)) ) stop(sprintf("%s in the test data contains new levels",z))
 }
  
  
# temproary data creation 
temp = data.frame(matrix(rep(0,length(idv)*nrow(test)),nrow(test),length(idv)))
name = names(test)
sub_name = idv
for (j in 1:length(sub_name)){
sub_name[j] = paste(sub_name[j],"cat",sep = "_")
}
colnames(temp) <- c(sub_name)

 k = 1
for (i in 1:length(idv)) {
     tab = table(train.output[,idv[i]],train.output[,sub_name[i]])
 temp[,i] = xyz.pbin(tab,test[,idv[i]],sub_name[i])
}
     test = cbind(test,temp)
return(test)
}
