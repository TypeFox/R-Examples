getFileName = function(fullpath){
fullpath=rmExt(fullpath)
fullpath=basename(fullpath)
return(fullpath)
}
#getFileName("C:\teste\teste\teste.rwl")