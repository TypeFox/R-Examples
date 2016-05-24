# global script for CA2
if(inherits(x, "data.frame")){
  nomData=nomData
  newdata=x
  colonnesup=NULL
  lignesup=NULL
  catsup=NULL
  axe1=1
  axe2=2
  Invisible=NULL
  selec1="no"
  selec2="no"
  valueselec1=NULL
  valueselec2=NULL
  size=1
  title1="CA factor map"
}

if(inherits(x, "CAshiny")){
nomData=x$nomData
newdata=x$data
colonnesup=x$a
lignesup=x$b
catsup=x$c
axe1=x$d
axe2=x$e
Invisible=x$f
selec1=x$type1
selec2=x$type2
valueselec1=x$selec1
valueselec2=x$selec2
size=x$taille
title1=x$title1
}

if(inherits(x, "CA")){
  nomData=as.character(x$call$call[2])
#  nomData=unlist(strsplit(nomData, split='[', fixed=TRUE))[1]
  newdata=x$call$Xtot
  colonnesup=rownames(x$col.sup$coord)
  lignesup=rownames(x$row.sup$coord)
  catsup=NULL
  axe1=1
  axe2=2
  Invisible=NULL
  selec1="no"
  selec2="no"
  valueselec1=NULL
  valueselec2=NULL
  size=1
  title1="CA factor map"
}

withna=c()
rowna=c()
nomrow=c()
for (i in 1:dim(newdata)[2]){
  if(any(is.na(newdata[,i])==TRUE)){
    if(is.numeric(newdata[,i])==TRUE){
      withna=c(withna,colnames(newdata)[i])
    }
  }
}

for (i in 1:dim(newdata)[1]){
  if(any(is.na(newdata[i,])==TRUE)){
      rowna=c(rowna,i)
      nomrow=c(nomrow,rownames(newdata)[i])
  }
}

quanti=names(which(sapply(newdata,is.numeric)))
quali=names(which(!(sapply(newdata,is.numeric))))
VariableChoice=quanti
noms=rownames(newdata)
nums=c(1:length(noms))
QualiChoice=quali
IdChoice=c(1:length(VariableChoice))
Idqualisup=c(1:length(QualiChoice))
sup=c()
for(i in IdChoice){
  if(VariableChoice[i]%in%withna){
    sup=c(sup,i)
  }
}
sup2=c()
if(!(is.null(sup))){
  VariableChoices=VariableChoice[-sup]
}
if(is.null(sup)){
  VariableChoices=VariableChoice
}
IdChoices=c(1:length(VariableChoices))
for(i in nums){
  if(noms[i]%in%nomrow){
    sup2=c(sup2,i)
  }
}
if(!(is.null(sup2))){
  nom=noms[-sup2]
}
if(is.null(sup2)){
  nom=noms
}
num=c(1:length(nom))