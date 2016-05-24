

epigen<- function(response, marker, kin, numkeep=floor(length(response)*.5), selectvar, corbnd=0.5, mafb=0.04, method="complete") {
### This function construct the genetic effect matrix including both main and epistatic effects.
### y is the vector of the phenotye values, marker is the matrix of main effects.
### kin is relationship matrix of all lines.
### numkeep and selectvar are passed to amltest() in selecting main effects used to construct epistatic effects
### First the function fit adaptive mixed lasso using amltest() with numkeep and selectvar as parameters
### The selected markers (the number is given by selectvar) are used to construct two-way epistatic effects.
### The resulting epistatic matrix is then cleaned with the Hclust package to remove highly correlated columns
### If several columns (genetic effects) are highly correlated, only one will be retained
### hcb and minor are parameters for SNPclust(), see its documentation.
### The value of this function is a list of three items.
###    sel.effects is the matrix including main and epistatic effects after removing highly correlated entries
###    marker1 is the one of the party to each of the genetic effects given as marker name (column names in marker)
###    marker2 is the other party.  if marker1[i]=marker2[i], it means the ith effect is a main effect.
### For example, suppose marker1[41]="wPt.3695",  marker2[41]="wPt.0959"
###    it means the 41th colum in sel.effects represents epistatic effect of markers "wPt.3695" and "wPt.0959"
### Since marker names are taken from input marker matrix, it helps to give meaningful column names



  tem <- amltest(response, marker, kin, numkeep=numkeep, selectvar=selectvar)
  epiX<- marker[, tem$est[,1]]
  mkname<- colnames(epiX)
  prek<- dim(epiX)[2]
  co1<-1:prek
  co2<-1:prek

  for(i in 1:(prek-1)){

        co1<-c(co1,rep(i, prek-i))
        co2<-c(co2,(i+1):prek)
        interact<- epiX[,i]*epiX[,(i+1):prek, drop=FALSE]
        colnames(interact)<- paste(mkname[i],"..",mkname[(i+1):prek])
        epiX<-cbind(epiX, interact)

##        cat(i,j,"\n")

   }


sel.effects<- cleanclust(epiX,  mafb=mafb, corbnd=corbnd, method=method)


return(list(effects=sel.effects$newmarker, marker1=mkname[co1[sel.effects$tagged]], marker2=mkname[co2[sel.effects$tagged]])) 


}



