triangle.test <- function (design,answer,preference=NULL){

answer = gsub("(\\w)", "\\U\\1", as.character(answer), perl=TRUE)
labprod = levels(as.factor(c(as.character(design[,1]),as.character(design[,2]),as.character(design[,3]))))
nbprod = length(labprod)
nb.answer = nb.good = pref = matrix(0,nbprod,nbprod)
for (i in 1:nrow(design)){
  for (j in 1:nbprod){
     if (labprod[j] == design[i,1]) i1 = j
     if (labprod[j] == design[i,2]) i2 = j
     if (labprod[j] == design[i,3]) i3 = j
  }
  if (i1!=i2) nb.answer [i1,i2] = nb.answer[i1,i2]+1
  if (i1==i2) nb.answer [i1,i3] = nb.answer[i1,i3]+1
  if ((i1==i2)&(answer[i]=="Z")){
    nb.good[i1,i3]=nb.good[i1,i3]+1
    if (length(preference)>0){
      if (preference[i]!="Z") pref[i3,i1] = pref[i3,i1] +1
      if (preference[i]=="Z") pref[i1,i3] = pref[i1,i3] +1
    }
  }
  if ((i1==i3)&(answer[i]=="Y")){
    nb.good[i1,i2]=nb.good[i1,i2]+1
    if (length(preference)>0){
      if (preference[i]!="Y") pref[i2,i1] = pref[i2,i1] +1
      if (preference[i]=="Y") pref[i1,i2] = pref[i1,i2] +1
    }
  }
  if ((i2==i3)&(answer[i]=="X")){
    nb.good[i1,i2]=nb.good[i1,i2]+1
    if (length(preference)>0){
      if (preference[i]!="X") pref[i1,i2] = pref[i1,i2] +1
      if (preference[i]=="X") pref[i2,i1] = pref[i2,i1] +1
    }
  }
}
nb.good = nb.good + t(nb.good)
nb.answer = nb.answer + t(nb.answer)

diag(nb.answer)=1
prob = pbinom(nb.good-1,nb.answer,1/3,lower.tail=FALSE)
maxML = recognize = minimum = matrix(NA,nbprod,nbprod)
for (i in 1: (nbprod-1)){
  for (j in (i+1):nbprod){
    aux = matrix(0,nb.good[i,j]+1,1)
    for (k in 0:nb.good[i,j]) aux[k] = dbinom(nb.good[i,j]-k,nb.answer[i,j]-k,1/3)
    maxML[i,j] = maxML[j,i] = max(aux)
    recognize[i,j] = recognize[j,i] = rev(order(aux))[1]-1
    mini = 0
    for (k in round(nb.answer[i,j]/3):nb.answer[i,j]) if ((mini==0)&(dbinom(k,nb.answer[i,j],1/3)<0.05)) mini=k
    minimum[i,j]=minimum[j,i]=mini
  }
}

confusion = (nb.answer-recognize) / nb.answer
diag(nb.answer)=diag(recognize)=0
diag(maxML)=diag(confusion)=1

rownames(nb.answer) = colnames(nb.answer) = rownames(nb.good) = colnames(nb.good) = labprod
rownames(prob) = colnames(prob)= rownames(confusion) = colnames(confusion)= labprod
rownames(maxML) = colnames(maxML) = rownames(minimum) = colnames(minimum) = rownames(recognize) = colnames(recognize) = labprod
if (length(preference)>0) rownames(pref) = colnames(pref) = labprod

res = list()
res$nb.comp = nb.answer
res$nb.ident = nb.good
res$p.value = prob
res$nb.recognition = recognize
res$maxML = maxML
res$confusion = confusion
res$minimum = minimum
if (length(preference)>0) res$pref = pref
##res$complete = result
return(res)
}
