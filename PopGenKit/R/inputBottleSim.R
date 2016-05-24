inputBottleSim <-
function(datafile, ndigit=3)
{

input=read.delim(file=datafile, sep='\r')
noext= gsub("[.].*$","",datafile)

popcount=rbind(which(input == 'pop', arr.ind=T),which(input=='Pop',arr.ind=T),which(input=='POP',arr.ind=T))
popcount=popcount[,1]
npops=NROW(popcount)

popcount=sort(popcount)
nloci=popcount[1]-1

popsizes=vector(length=npops)

for (h in 1:npops)
{
if (h<npops) popsizes[h]=popcount[h+1]-popcount[h]-1
if (h==npops) popsizes[h]=nrow(input)-popcount[h]
}

n.individuals=sum(popsizes)

bottlesim=matrix(NA, n.individuals+3, nloci*2+1)

bottlesim[1,1]=n.individuals
bottlesim[2,1]=nloci
for (a in 1:nloci){bottlesim[3,a]=as.character(input[a,])}

input2=input[-popcount,]
input3=input2[-(1:nloci)]

#cette boucle pour chaque individu 
for (b in 1:n.individuals)
{
inddata=input3[b]
inddata=as.vector(inddata)
indID= gsub("[ ]","",inddata)
indID= gsub("[\t]|[,].*$","",indID)

bottlesim[(3+b),1]=substr(indID,1,9) 

inddata=gsub( "[^[:alnum:]]", "", inddata)
#ceci pour savoir où commencent les allèles, au cas où nom de lind inclus alphanum
total=nchar(inddata)
diff=total-(ndigit*2*nloci)

#cette boucle pour chaque locus par individu
for (x in 1:nloci)
{
bottlesim[(3+b),(1+((2*x)-1))]=substr(inddata,diff+(2*ndigit*(x-1))+(ndigit-(ndigit-1)), diff+(2*ndigit*(x-1))+ndigit)
bottlesim[(3+b),(1+(2*x))]=substr(inddata,diff+(2*ndigit*(x-1))+ndigit+1, diff+(2*ndigit*(x-1))+(2*ndigit))
}

}

if (ndigit==2) missing='00'
if (ndigit==3) missing='000'


#mettre les cellules NA vides pour lesthétisme !
for (r in 1:nrow(bottlesim))
{
for (c in 1:ncol(bottlesim))
{
if (is.na(bottlesim[r,c])==T) bottlesim[r,c]=''
}
} 

#remplacer les données manquantes
for (x in 1:nrow(bottlesim)) 
{
for (y in 1:ncol(bottlesim))
{
if  (bottlesim[x,y]==missing) bottlesim[x,y]='?'
}}

bottlesim2=bottlesim

alleles.bsim=c(0,1,2,3,4,5,6,7,8,9, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z')

#déterminer le nb dallèles par locus, modifier les données, etc
final.vec=NULL
kvec=rep(NA,times=nloci)
for (i in 1:nloci)
{
alleles=bottlesim[4:(nrow(bottlesim)),(i*2):((i*2)+1)]
alleles=as.vector(alleles)
vec=unique(alleles)
if (length(which(vec[]=='?')>0)) vec=vec[-which(vec[]=='?')]
kvec[i]=length(vec)

sort.vec=sort(vec)

for (j in 1:kvec[i])
{
for (x in 4:(nrow(bottlesim)) )
{
for (y in (i*2):((i*2)+1))
{
if (bottlesim[x,y]==sort.vec[j]) bottlesim2[x,y] = alleles.bsim[j]
}}}
final.vec=c(final.vec,sort.vec)
}

filename=paste(noext,'_Bsim.txt',sep='')
write.table(bottlesim2, file=filename, quote=F, row.names=F, col.names=F, sep='\t')
save(final.vec, kvec, file='alleles')

return(list('final.vec'=final.vec, 'kvec'=kvec))

}

