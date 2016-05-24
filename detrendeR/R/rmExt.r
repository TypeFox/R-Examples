rmExt= function(string){
for (i in 1:nchar(string)){
strFind <- substr(string, nchar(string)+1-i, nchar(string)+1-i)
if (strFind=="."){
break
}
nChr=i+1
}
string<-substr(string, 1, nchar(string)-nChr)
return(string)
}
