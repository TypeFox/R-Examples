ConvertPed=function(ped,persons=NULL){
#Columns of ped are: ID FID MID SEX AFF(not used) locus1.1 locus1.2 locus2.1 locus2.2 ...
if(is.null(persons)) persons=ped[,1]
if(class(ped)!="data.frame") stop("First argument, ped, should be a data.frame.")
if(dim(ped)[1]!=length(unique(persons))) stop("Number of rows in ped should equal number of persons.")
ped[ped==0]=NA
gender=c("male","female")
ped2=pedigree(id=persons[ped[,1]], dadid=persons[ped[,2]], momid=persons[ped[,3]],sex=gender[ped[,4]]) 
if(dim(ped)[2]>5) {
datamatrix=ped[,-(1:5)]
rownames(datamatrix)=persons
}
else
datamatrix=NULL
list(ped=ped2,datamatrix=datamatrix)
}

