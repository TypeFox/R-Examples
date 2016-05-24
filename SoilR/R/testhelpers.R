#
# vim:set ff=unix expandtab ts=2 sw=2:
namedPlot=function(lexp,env){
  #get the name of the caller
  print(sys.calls())
  fileName=paste(as.character(sys.calls()[[sys.nframe()-1]]),"pdf",sep=".")
  plotAndCheck(fileName,lexp,env)
}
#####################################################################################################
plotAndCheck=function(fileName,lexp,env){
pdf(file=fileName)
eval(lexp,env)
dev.off()
res=system(command=paste("qpdf --check ",fileName,sep=""))
checkEquals(attr(res,"status"),NULL)
}
#####################################################################################################
pp=function(# print out an 
### This function is used to print out a value of a variable together with its name and is helpful for debugging
string,env){
print(paste(string,"="))
print(get(string,env))
}
#####################################################################################################
pe=function(string,env){
print(string)
print(eval(string,env))
}
