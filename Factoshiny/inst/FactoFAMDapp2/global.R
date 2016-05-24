#global script for PCA2
if(inherits(x, "data.frame")){
  nomData=nomData
  newdata=x
  quantisup=NULL
  qualisup=NULL
  indsupl=NULL
  axe1=1
  axe2=2
  labind=TRUE
  labvar=TRUE
  habillageind=NULL
  selection="NONE"
  selection2=NULL
  selection3="NONE"
  selection4=NULL
  selection5="NONE"
  selection6=NULL
  size=1
  size2=1
  size3=1
  title1="Individual factor map"
  title2="Variables representation"
  title3="Correlation circle"
}

if(inherits(x, "FAMDshiny")){
  nomData=x$nomData
  newdata=x$data
  quantisup=x$b
  qualisup=x$c
  indsupl=x$d
  axe1=x$e
  axe2=x$f
  habillageind=x$g
  selection=x$h
  selection2=x$i
  selection3=x$j
  selection4=x$k
  selection5=x$o
  selection6=x$p
  size=x$l
  size2=x$m
  size3=x$n
  title1=x$title1
  title2=x$title2
  title3=x$title3
  labind=x$labind
  labvar=x$labvar
}
if(inherits(x, "FAMD")){
  nomData=as.character(x$call$call[2])
  newdata=x$call$X
  quantisup=colnames(x$call$quanti.sup)
  qualisup=colnames(x$call$quali.sup$quali.sup)
  indsupl=rownames(x$ind.sup$coord)
  axe1=1
  axe2=2
  habillageind=NULL
  selection="NONE"
  selection2=NULL
  selection3="NONE"
  selection4=NULL
  selection5="NONE"
  selection6=NULL
  size=1
  size2=1
  size3=1
  title1="Individual factor map"
  title2="Variables representation"
  title3="Correlation circle"
  labind=TRUE
  labvar=TRUE
}  
all=colnames(newdata)
quanti=names(which(sapply(newdata,is.numeric)))
quali=names(which(!(sapply(newdata,is.numeric))))
VariableChoices=quanti
nom=rownames(newdata)
num=c(1:length(nom))
QualiChoice=quali
IdChoices=c(1:length(VariableChoices))
Idqualisup=c(1:length(QualiChoice))
Idall=c(1:length(all))