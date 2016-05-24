sumry_categorical <-
function(categ_var_test){
names(categ_var_test)
nr=0
for (i in 1:ncol(categ_var_test)){
  x=length(unique(categ_var_test[,i]))
  nr=nr+x
}
nr=nr+ncol(categ_var_test)
freq_mat=matrix(0,nrow=nr ,ncol=4)

#Assigning Column name
colnames(freq_mat)<-c("Var_Desc","Levels","Frequency","Percentage")

#Assigning VARIABLE name and values
var_name=names(categ_var_test)
freq_mat[1,1]=var_name[1]

ln=ncol(categ_var_test)-1
t=0
for (i in 1:ln){
  x=length(unique(categ_var_test[,i])) +1
  t=t+x
  freq_mat[t+1,1]=var_name[i+1]
}

##
stpt=1
enpt=0
for(i in 1:ncol(categ_var_test)){
  x=length(unique(categ_var_test[,i]))
  enpt=enpt+x+1
  
  freq_mat[stpt:enpt,2]=rownames(freq(categ_var_test[,i]))
  freq_mat[stpt:enpt,3:4]=(freq(categ_var_test[,i]))[,1:2]
  stpt=enpt+1
}
final_cat_data=as.data.frame(freq_mat)
final_cat_data[,1]=as.character(final_cat_data[,1])
for(i in 1:nrow(final_cat_data)){
  
if(final_cat_data[i,1]==0)  {
  final_cat_data[i,1]=" "
}
}
return(final_cat_data)
}
