onko<-function(rivi,j){
#Checks whether j is in rivi
#
#rivi is vector where beginning is positive integers, rest NA
#j is positive inetger
#
#Returns TRUE is j is in rivi, FALSE otherwise
#
len<-length(rivi)
res<-FALSE
i<-1
while ((!is.na(rivi[i])) && (i<=len) && (rivi[i]<=j)){
  if (rivi[i]==j) res<-TRUE
  i<-i+1
}
return(res)
}

