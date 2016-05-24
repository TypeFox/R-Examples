######################################
##MULTIPLE-SITE BIODIVERSITY INDICES##
##MBI                               ##
##YOUHUA CHEN					    ##
##10/13/2012						##
##CPL>=2.0						    ##
######################################


#
#Baselga's full index (nestdness+turnover)
#rows of the data are sites
cfull<-function(data)
{
st<-dim(data)[1]
deno<-0
sp<-dim(data)[2]
cn<-seq(1,sp,1)
#
index<-t(combn(c(1:st),2))
for(k in 1:dim(index)[1])
{
i<-index[k,1]
j<-index[k,2]
cn1<-cn[which(data[i,]>0)]
cn2<-cn[which(data[j,]>0)]
bij<-length(cn1[!cn1%in%cn2])
bji<-length(cn2[!cn2%in%cn1])
#
deno<-deno+min(bij,bji)+max(bij,bji)
}#k
#sites
si<-0
for(i in 1:sp)
{
si<-si+length(which(data[,i]>0))
}
sim<-deno/(deno+2*(si-sp))
return(sim)
}#end


#Baselga's turnover index
#rows of the data are sites
ct<-function(data)
{
st<-dim(data)[1]
deno<-0
sp<-dim(data)[2]
cn<-seq(1,sp,1)
#
index<-t(combn(c(1:st),2))
for(k in 1:dim(index)[1])
{
i<-index[k,1]
j<-index[k,2]
cn1<-cn[which(data[i,]>0)]
cn2<-cn[which(data[j,]>0)]
bij<-length(cn1[!cn1%in%cn2])
bji<-length(cn2[!cn2%in%cn1])
#
deno<-deno+min(bij,bji)
}#k
#sites
si<-0
for(i in 1:sp)
{
si<-si+length(which(data[,i]>0))
}
sim<-deno/(deno+si-sp)
return(sim)
}#end

#Baselga's nestedness
#rows of the data are sites
cn<-function(data)
{
st<-dim(data)[1]
dmin<-0
dmax<-0
sp<-dim(data)[2]
cnn<-seq(1,sp,1)
#
index<-t(combn(c(1:st),2))
for(k in 1:dim(index)[1])
{
i<-index[k,1]
j<-index[k,2]
cn1<-cnn[which(data[i,]>0)]
cn2<-cnn[which(data[j,]>0)]
bij<-length(cn1[!cn1%in%cn2])
bji<-length(cn2[!cn2%in%cn1])
#
dmin<-dmin+min(bij,bji)
dmax<-dmax+max(bij,bji)
}#k
#sites
si<-0
for(i in 1:sp)
{
si<-si+length(which(data[,i]>0))
}
sim<-(dmax-dmin)/(2*(si-sp)+dmax+dmin)*(si-sp)/(si-sp+dmin)
return(sim)
}

#lennon's richness
cl<-function(data)
{
st<-dim(data)[1]
dmin<-0
dmax<-0
sp<-dim(data)[2]
cn<-seq(1,sp,1)
#
index<-t(combn(c(1:st),2))
for(k in 1:dim(index)[1])
{
i<-index[k,1]
j<-index[k,2]
cn1<-cn[which(data[i,]>0)]
cn2<-cn[which(data[j,]>0)]
bij<-length(cn1[!cn1%in%cn2])
bji<-length(cn2[!cn2%in%cn1])
#
dmin<-dmin+min(bij,bji)
dmax<-dmax+max(bij,bji)
}#k
#sites
si<-0
for(i in 1:sp)
{
si<-si+length(which(data[,i]>0))
}
sim<-2*(dmax-dmin)/(2*(si-sp)+dmax+dmin)
return(sim)
}



#Ulrich's nesetdness index for both abundance and occurrence analysis
#rows of the data are sites
wnodf<-function(data)
{
n<-dim(data)[2]
m<-dim(data)[1]
if(length(which(data>1))!=0)
{
type=1
}else
{
type=2
}
cindex<-t(combn(c(1:n),2))
#
fc<-0
for(k in 1:dim(cindex)[1])
{
i<-cindex[k,1]
j<-cindex[k,2]
if(type==1)
{
num<-length(which(data[,j]-data[,i]<0 & data[,j]>0)) #must be non-zero cells
nj<-length(which(data[,j]>0))
fc<-fc+num/nj
}else if(type==2)
{
if(sum(data[,j])>=sum(data[,i]))
{
po=0
}else
{
that<-which(data[,i]==1)
ind<-which(data[that,j]==1)
po=length(ind)/sum(data[,j])
}
fc<-fc+po
}
}#k
#
rindex<-t(combn(c(1:m),2))
fr<-0
for(k in 1:dim(rindex)[1])
{
i<-rindex[k,1]
j<-rindex[k,2]
if(type==1)
{
num<-length(which(data[j,]-data[i,]<0 & data[j,]>0)) #must be non-zero cells
nj<-length(which(data[j,]>0))
fr<-fr+num/nj
}else if(type==2)
{
if(sum(data[j,])>=sum(data[i,]))
{
po=0
}else
{
that<-which(data[i,]==1)
ind<-which(data[j,that]==1)
po=length(ind)/sum(data[j,])
}
fr<-fr+po
}
}#k
#
#total
ff<-2*(fc+fr)*100/(m*(m-1)+n*(n-1))
#
return(list(total=ff,column=fc*200/(n*(n-1)),row=fr*200/(m*(m-1))))
}#end

#replacement of community
#rows of the data are sites
crep<-function(data)
{
st<-dim(data)[1]
deno<-0
dem<-0
sp<-dim(data)[2]
cn<-seq(1,sp,1)
#
index<-t(combn(c(1:st),2))
for(k in 1:dim(index)[1])
{
i<-index[k,1]
j<-index[k,2]
cn1<-cn[which(data[i,]>0)]
cn2<-cn[which(data[j,]>0)]
bij<-length(cn1[!cn1%in%cn2])
bji<-length(cn2[!cn2%in%cn1])
#
deno<-deno+min(bij,bji)
dem<-dem+bij+bji
}#k
#sites
si<-0
for(i in 1:sp)
{
si<-si+length(which(data[,i]>0))
}
sim<-2*deno/(dem+si-sp)
return(sim)
}#end


#richness difference
#rows of the data are sites
crich<-function(data)
{
st<-dim(data)[1]
dmin<-0
dmax<-0
sp<-dim(data)[2]
cn<-seq(1,sp,1)
index<-t(combn(c(1:st),2))
for(k in 1:dim(index)[1])
{
i<-index[k,1]
j<-index[k,2]
cn1<-cn[which(data[i,]>0)]
cn2<-cn[which(data[j,]>0)]
bij<-length(cn1[!cn1%in%cn2])
bji<-length(cn2[!cn2%in%cn1])
#
dmin<-dmin+min(bij,bji)
dmax<-dmax+max(bij,bji)
}#k
#sites
si<-0
for(i in 1:sp)
{
si<-si+length(which(data[,i]>0))
}
sim<-(dmax-dmin)/(2*(si-sp)+dmax+dmin)
return(sim)
}


#whittaker's beta
wbeta<-function(data)
{
dat<-data
st<-dim(dat)[1]
sp<-dim(dat)[2]
dat[which(dat>1)]=1
si<-rowSums(dat)
return(sp/mean(si))
}

#Harrison's dissimilarity
harrison<-function(data)
{
term<-wbeta(data)
st<-dim(data)[1]
return((term-1)/(st-1))
}

#Diserud-Odegaard similarity index
do<-function(data)
{
return(1-harrison(data))
}

#Harrison's turnover
ht<-function(data)
{
sp<-dim(data)[2]
data[which(data>1)]<-1
si<-rowSums(data)
tt<-dim(data)[1]
#
return((sp/max(si)-1)/(tt-1))
}

#William's turnover
wt<-function(data)
{
sp<-dim(data)[2]
data[which(data>1)]<-1
si<-rowSums(data)
return(1-max(si)/sp)
}

#calculate mean Jaccard distance index
#rows of the data are sites
mjaccard<-function(data)
{
st<-dim(data)[1]
index<-t(combn(c(1:st),2))
#
mm<-vector()
#
for(i in 1:dim(index)[1])
{
aa<-which(data[index[i,1],]>0)
bb<-which(data[index[i,2],]>0)
#
onlya<-length(aa[!aa%in%bb])
onlyb<-length(bb[!bb%in%aa])
ab<-length(aa[aa%in%bb])
mm[i]<-(onlya+onlyb)/(ab+onlya+onlyb)
}
#
#calculate mean of these pairwise sites
res<-mean(mm,na.rm=TRUE)
return(res)
}#end


#calculate mean Sorensen distance index
#rows of the data are sites
msorensen<-function(data)
{
st<-dim(data)[1]
index<-t(combn(c(1:st),2))
#
mm<-vector()
#
for(i in 1:dim(index)[1])
{
aa<-which(data[index[i,1],]>0)
bb<-which(data[index[i,2],]>0)
#
onlya<-length(aa[!aa%in%bb])
onlyb<-length(bb[!bb%in%aa])
ab<-length(aa[aa%in%bb])
mm[i]<-(onlya+onlyb)/(2*ab+onlya+onlyb)
}
#
#calculate mean of these pairwise sites
res<-mean(mm,na.rm=TRUE)
return(res)
}#end


#calculate mean Baselga's turnover index
#rows of the data are sites
mt<-function(data)
{
st<-dim(data)[1]
index<-t(combn(c(1:st),2))
#
mm<-vector()
#
for(i in 1:dim(index)[1])
{
aa<-which(data[index[i,1],]>0)
bb<-which(data[index[i,2],]>0)
#
onlya<-length(aa[!aa%in%bb])
onlyb<-length(bb[!bb%in%aa])
ab<-length(aa[aa%in%bb])
mm[i]<-min(onlya,onlyb)/(ab+min(onlya,onlyb))
}
#
#calculate mean of these pairwise sites
res<-mean(mm,na.rm=TRUE)
return(res)
}#end

#calculate mean Baselga's nestedness index
#rows of the data are sites
mn<-function(data)
{
st<-dim(data)[1]
index<-t(combn(c(1:st),2))
#
mm<-vector()
#
for(i in 1:dim(index)[1])
{
aa<-which(data[index[i,1],]>0)
bb<-which(data[index[i,2],]>0)
#
onlya<-length(aa[!aa%in%bb])
onlyb<-length(bb[!bb%in%aa])
ab<-length(aa[aa%in%bb])
mm[i]<-ab/(ab+min(onlya,onlyb))*abs(onlya-onlyb)/(2*ab+onlya+onlyb)
}
#
#calculate mean of these pairwise sites
res<-mean(mm)
return(res)
}#end



#calculate mean compositional richness difference index
#rows of the data are sites
mrich<-function(data)
{
st<-dim(data)[1]
index<-t(combn(c(1:st),2))
#
mm<-vector()
#
for(i in 1:dim(index)[1])
{
aa<-which(data[index[i,1],]>0)
bb<-which(data[index[i,2],]>0)
#
onlya<-length(aa[!aa%in%bb])
onlyb<-length(bb[!bb%in%aa])
ab<-length(aa[aa%in%bb])
mm[i]<-abs(onlya-onlyb)/(ab+abs(onlya-onlyb))
}
#
#calculate mean of these pairwise sites
res<-mean(mm,na.rm=TRUE)
return(res)
}#end



#calculate mean species replacement index
#rows of the data are sites
mrep<-function(data)
{
st<-dim(data)[1]
index<-t(combn(c(1:st),2))
#
mm<-vector()
#
for(i in 1:dim(index)[1])
{
aa<-which(data[index[i,1],]>0)
bb<-which(data[index[i,2],]>0)
#
onlya<-length(aa[!aa%in%bb])
onlyb<-length(bb[!bb%in%aa])
ab<-length(aa[aa%in%bb])
mm[i]<-2*min(onlya,onlyb)/(ab+onlya+onlyb)
}
#
#calculate mean of these pairwise sites
res<-mean(mm,na.rm=TRUE)
return(res)
}#end

#calculate mean Lennon richness index
#rows of the data are sites
ml<-function(data)
{
st<-dim(data)[1]
index<-t(combn(c(1:st),2))
#
mm<-vector()
#
for(i in 1:dim(index)[1])
{
aa<-which(data[index[i,1],]>0)
bb<-which(data[index[i,2],]>0)
#
onlya<-length(aa[!aa%in%bb])
onlyb<-length(bb[!bb%in%aa])
ab<-length(aa[aa%in%bb])
mm[i]<-2*(max(onlya,onlyb)-min(onlya,onlyb))/(2*ab+max(onlya,onlyb)+min(onlya,onlyb))
}
#
#calculate mean of these pairwise sites
res<-mean(mm,na.rm=TRUE)
return(res)
}#end


 #calculate diversity indices, batch handling
 #mat is a list of matrices with varying sizes
 batch.calculation<-function(mat)
 {
  if(is.list(mat)==TRUE)
 {
 res<-matrix(0,nrow=length(mat),ncol=20)
 colnames(res)<-c("cn","ct","crep","crich","wnodfT","wnodfC","wnodfR","wbeta","harrison","ht","wt","do","cl"
                  ,"msorensen","mjaccard","mn","mt","mrep","mrich","ml","cfull")
 for(i in 1:length(mat))
 {
 res[i,1]<-cn(mat[[i]])
 res[i,2]<-ct(mat[[i]])
 res[i,3]<-crep(mat[[i]])
 res[i,4]<-crich(mat[[i]])
 wf<-wnodf(mat[[i]])
 res[i,5:7]<-c(wf$total, wf$column,wf$row)
 res[i,8]<-wbeta(mat[[i]])
 res[i,9]<-harrison(mat[[i]])
 res[i,10]<-ht(mat[[i]])
 res[i,11]<-wt(mat[[i]])
 res[i,12]<-do(mat[[i]])
 res[i,13]<-cl(mat[[i]])
 res[i,14]<-msorensen(mat[[i]])
 res[i,15]<-mjaccard(mat[[i]])
 res[i,16]<-mn(mat[[i]])
 res[i,17]<-mt(mat[[i]])
 res[i,18]<-mrep(mat[[i]])
 res[i,19]<-mrich(mat[[i]])
 res[i,20]<-ml(mat[[i]])
 res[i,21]<-cfull(mat[[i]])
 #
 #
 }#
 }else
 {
 res<-vector()
 res[1]<-cn(mat)
 res[2]<-ct(mat)
 res[3]<-crep(mat)
 res[4]<-crich(mat)
 wf<-wnodf(mat)
 res[5:7]<-c(wf$total, wf$column,wf$row)
 res[8]<-wbeta(mat)
 res[9]<-harrison(mat)
 res[10]<-ht(mat)
 res[11]<-wt(mat)
 res[12]<-do(mat)
 res[13]<-cl(mat)
 res[14]<-msorensen(mat)
 res[15]<-mjaccard(mat)
 res[16]<-mn(mat)
 res[17]<-mt(mat)
 res[18]<-mrep(mat)
 res[19]<-mrich(mat)
 res[20]<-ml(mat)
 res[21]<-cfull(mat)
 names(res)<-c("cn","ct","crep","crich","wnodfT","wnodfC","wnodfR","wbeta","harrison","ht","wt","do","cl"
 ,"msorensen","mjaccard","mn","mt","mrep","mrich","ml","cfull")
 }
 return(res)
 }#end batch


#rarity
#rows of the data are sites
rarity<-function(data,percent=.3)
{
sp<-dim(data)[2]
st<-dim(data)[1]
maxa<-max(colSums(data))
rc<-0
dc<-0
for(i in 1:sp)
{
if(length(which(data[,i]>0))<=percent*st)
{
rc=rc+1
}
if(sum(data[,i])<=percent*maxa)
{
dc=dc+1
}
}
return(list(rangerarity=rc/st,densityrarity=dc/maxa))
}#end
