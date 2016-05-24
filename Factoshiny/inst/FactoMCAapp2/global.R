#global script for MCA
if(is.data.frame(x)==TRUE){
  newdata=x
  nomData=nomData
  axe1=1
  axe2=2
  varsup=c("suplquali","suplquanti")
  indvar=c("Ind","Mod","Modsup","Indsup")
  labvar=c()
  habillageind=NULL
  selection="NONE"
  selection2=NULL
  selection3="NONE"
  selection4=NULL
  
  quanti=names(which(sapply(newdata,is.numeric)))
  quali=names(which(!(sapply(newdata,is.numeric))))
  supquali=NULL
  quantiS=NULL
  varsup=c("act")
  indvar=c("Ind","Mod","Indsup","Modsup")
  labvar=c("Ind","Mod","Indsup","Modsup")
  indsupl=NULL
  title1="MCA factor map"
  title2=""
  title3="Supplementary variables on the MCA factor map"
  
}
if(inherits(x, "MCA")||inherits(x, "MCAshiny")){

if(inherits(x, "MCAshiny")){
  ###
  newdata=x$data
  nomData=x$nomData
  quantiS=x$b
  supquali=x$c
  varsup=x$z
  indvar=x$y
  indsupl=x$d
  axe1=x$e
  axe2=x$f
  habillageind=x$g
  selection=x$h
  selection2=x$i
  selection3=x$j
  selection4=x$k
  title1=x$title1
  title2=x$title2
  title3=x$title3
}
if(inherits(x, "MCA")){
  nomData=as.character(x$call$call[2])
#  nomData=unlist(strsplit(nomData, split='[', fixed=TRUE))[1]
  newdata=x$call$X
  quantiS=rownames(x$quanti.sup$coord)
  supquali=rownames(x$quali.sup$coord)
  indsupl=rownames(x$ind.sup$coord)
  axe1=1
  axe2=2
  varsup=c("suplquali","suplquanti")
  indvar=c("Ind","Mod","Modsup","Indsup")
  labvar=c("Ind","Mod","Indsup","Modsup")
  habillageind=NULL
  selection="NONE"
  selection2=NULL
  selection3="NONE"
  selection4=NULL
  title1="MCA factor map"
  title2=""
  title3="Supplementary variables on the MCA factor map"
}  }


###
quanti=names(which(sapply(newdata,is.numeric)))
quali=names(which(!(sapply(newdata,is.numeric))))
VariableChoices=quali

nom=rownames(newdata)
num=c(1:length(nom))

QuantiChoice=quanti

IdChoices=c(1:length(VariableChoices))
Idquantisup=c(1:length(QuantiChoice))