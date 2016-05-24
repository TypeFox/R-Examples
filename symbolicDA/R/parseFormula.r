.parseFormula<-function(form,sdt){
t1<-t<-as.character(sdt$variables[,"label"])
dim(t)<-c(1,length(t))
t<-as.data.frame(t)
names(t)<-t1
z1<-terms(formula(form),data=t)
variableSelection<-attr(z1,"term.label")
variables<-as.matrix(sdt$variables)
indivN<-as.matrix(sdt$indivN)
vsNrs<-NULL
for(v in variableSelection){
vsNrs<-c(vsNrs,as.numeric(variables[variables[,"label"]==v,"num"]))
}
#eval(attr(z1,"variables")
#attach(t)
#vc<-as.character(eval(attr(z1,"variables"))[[1]])
#detach(t)
vc<-with (t,as.character(eval(attr(z1,"variables"))[[1]]))
categorialVariable<-as.numeric(variables[variables[,"label"]==vc,"num"])
classes<-as.numeric(indivN[indivN[,"variable"]==categorialVariable,"value"])
resul<-list(classesVariable=categorialVariable,classes=classes,variableSelection=vsNrs)
}
