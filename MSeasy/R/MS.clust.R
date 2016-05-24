MS.clust <-
function(data_tot, quant=FALSE, clV, ncmin, ncmax, Nbc, varRT = 0.1, disMeth="euclidean", linkMeth="ward", clustMeth="hierarchical")
{
#Rprof()
if(missing(data_tot)){
data_tot<-tk_choose.files(default=file.path(getwd(),"*.txt"),caption="Please, select a file for MS.clust (initial_DATA.txt)", multi=FALSE, filters=matrix(c("your file","*.txt"),ncol=2))
data_tot<-read.table(data_tot, header=TRUE, check.names=FALSE)
}	
### MSeasy v 1.3 March 2011 with improved errorhandling
st<-strsplit(date(), " ")[[1]]
stBis<-strsplit(st[4], ":")[[1]]
Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
Date<-paste(st[1], st[2], st[3], Hour, sep="_")
Mypath<-paste("output_MSclust", "_", "result", Date, sep="")
dir.create(Mypath)
library(cluster)
data_tot<-as.data.frame(data_tot)
if (quant==TRUE) {Beg<-5}
if (quant==FALSE) {Beg<-3}

if (clV==FALSE) #check for error on Nbc or ncmax value
{
	if (as.numeric(dim(data_tot)[1])<as.numeric(Nbc)) 
		{
			
			stop(paste(cat("\n")," Error Nbc too high. Max value is: "),dim(data_tot)[1],cat("\n"))
		}
}else{
	if (as.numeric(dim(data_tot)[1])<as.numeric(ncmax)) 
		{
			
			stop(paste(cat("\n")," Error ncmax too high. Max value is: "),dim(data_tot)[1],cat("\n"))
		}
}

## if clV ==FALSE, the number(s) of clusters in the dataset is already determined; Nclust is a vector with the chosen number(s) of clusters

   if (clV==FALSE)
   {
   Nclust<-length(Nbc)
   answer<-Nbc
   }
   
## if clV==TRUE, determine the optimal number(s) of clusters
               
   else
   {
   	if (clustMeth=="hierarchical")
       {
	   require(amap)
		Dis<-amap::Dist(data_tot[,Beg:ncol(data_tot)],method = disMeth, diag = FALSE, upper = FALSE)
		gc()
		Hc<- hcluster(data_tot[,Beg:ncol(data_tot)],method = disMeth,link = linkMeth, nbproc= 2, doubleprecision = TRUE)
		gc()
		si<-c()
		l=1
		for (i in ncmin:ncmax){
			si[l]<-mean(silhouette(cutree(Hc, k=i), Dis)[,3])
			l=l+1
			gc()
			}
		plot(x=ncmin:ncmax, y=si, type="l", main = "silhouette width" , xlab = "Cluster", ylab = "silhouette width")
		abline(v=as.numeric(as.numeric(ncmin)+as.numeric(match(max(si),si)))-1, col="red")
		text(x= as.numeric(as.numeric(ncmin)+as.numeric(match(max(si),si)))+2, y=as.numeric(max(si)), col="red", label=as.character(as.numeric(ncmin)+as.numeric(match(max(si),si))-1))
		
		res_clValid<-cbind(ncmin:ncmax,si)
		print(res_clValid)
		write.table(res_clValid, file=paste(Mypath, "/", "res_clValid.txt", sep=""))
		#rm(dist1, hc1)
		gc()
	   }
	   else{
		library(clValid)
		clV<-clValid(data_tot[,Beg:ncol(data_tot)], nClust=ncmin:ncmax, validation="internal", clMethods = clustMeth, metric = disMeth, method = linkMeth, maxitems=nrow(data_tot))

		plot(x=ncmin:ncmax, y=measures(clV)[3,,1], type="l", main="silhouette width", xlab="number of clusters", ylab="silhouette width")
		abline(v=as.numeric(as.vector(optimalScores(clV)[3,3])), col="red")
		text(x= as.numeric(as.vector(optimalScores(clV)[3,3]))+2, y=as.numeric(as.vector(optimalScores(clV)[3,1])), col="red", label=as.character(as.vector(optimalScores(clV)[3,3])))

		res_clValid<-measures(clV)[3,,1]
		print(res_clValid)
		write.table(res_clValid, file=paste(Mypath, "/", "res_clValid.txt", sep=""))
		gc()
		}
		
   cat("How many clustering separations? \n If you want just one clustering separation, type 1, \n otherwise type the number of times you want to separate the dataset", "\n")
   Nclust<-readLines(n=1)

   cat("clustering method : disMeth= ",disMeth, "linkMeth= ", linkMeth, "clustMeth= ",clustMeth ,"\n")
   cat("How many clusters? \n If more than one clustering separation, type each number of clusters followed by ENTER successively", "\n")
   answer<-readLines(n=Nclust)
   }


   for (NC in 1:Nclust)
   {
              if (clustMeth=="hierarchical")
       {
		#in order to avoid multiple calculation of Dist
		if (clV==TRUE){
			
			cl<-cutree(Hc, k=as.numeric(answer[NC]))
		}
		else{
			require(amap)
			Dis<-amap::Dist(data_tot[,3:ncol(data_tot)],method=disMeth, diag = FALSE, upper = FALSE)
			Hc<-amap::hcluster(data_tot[,3:ncol(data_tot)], method = disMeth,link = linkMeth, nbproc= 2, doubleprecision = TRUE)
			cl<-cutree(Hc, k=as.numeric(answer[NC]))
		}
       }
	   
	   #if (clustMeth=="agnes")
       #{
       #Dis<-dist(data_tot[,Beg:ncol(data_tot)],method=disMeth)
       #Hc<-agnes(Dis, diss=TRUE, method=linkMeth)
       #cl<-cutree(Hc, k=as.numeric(answer[NC]))
       #}
	   
	   if (clustMeth=="diana")
       {
       Dis<-as.dist((1-cor(t(data_tot[,Beg:ncol(data_tot)]), method="spearman")))
       Hc<-diana(Dis, diss=TRUE)
       cl<-cutree(Hc, k=as.numeric(answer[NC]))
       }
	   
       if (clustMeth=="pam")
       {
       Dis<-dist(data_tot[,Beg:ncol(data_tot)],method=disMeth)
       cl<-pam(x=data_tot[,Beg:ncol(data_tot)], k=as.numeric(answer[NC]), metric=disMeth, cluster.only=TRUE)
       }
		
	if (clustMeth=="kmeans")
       {
       Dis<-dist(data_tot[,Beg:ncol(data_tot)],method=disMeth)
       cl<-kmeans(x=data_tot[,Beg:ncol(data_tot)], as.numeric(answer[NC]))$cluster
       }
	   
       

data_tot_temp<-cbind(cl, data_tot)
Sil.fin<-silhouette(cl, dist=Dis)
data_tot_temp_new<-cbind(data_tot_temp[1:3],Sil.fin[,3],Sil.fin[,2], data_tot_temp[4:ncol(data_tot_temp)])
colnames(data_tot_temp_new)<-c("cluster", "analyses","RT", "Silhouette_indiv","neighbor_cluster", as.numeric(colnames(data_tot_temp)[4:ncol(data_tot_temp)]))

resPclus<-vector()

        for (i in 1:as.numeric(answer[NC]))
        {
        resPclus[i]<-mean(Sil.fin[Sil.fin[,1]==i,3])
        }

clusters<-levels(as.factor(cl))
res<-matrix(nrow=length(clusters), ncol= 6 + ncol(data_tot_temp[,(Beg+1):ncol(data_tot_temp)]))
rownames(res)<-as.vector(clusters)
colnames(res)<-c("number_of_distinct_analyses_in_the_cluster", "number_of_analyses_in_the_cluster", "mean_retention_time", "range_retention_time", "mean_silhouette_width", "8_first_mz", colnames(data_tot_temp)[(Beg+1):ncol(data_tot_temp)])


        for (i in 1:length(clusters))
        {
        temp<-subset(data_tot_temp, data_tot_temp[,1]==clusters[i])
        res[i,1]<-length(levels(as.factor(as.vector(temp[,2])))) #one line per cluster from 1 to the optimal number of clusters entered by the user
        res[i,2]<-nrow(temp)									 #number of lines (analyses) for the tested number of clusters [i]
        res[i,3]<-mean(as.numeric(as.vector(temp[,3])))			 #calculation of the mean RT for the number of clusters [i]
        res[i,4]<-max(as.numeric(as.vector(temp[,3])))-min(as.numeric(as.vector(temp[,3])))

# determine the 8 most abundant m/z in the mass spectrum of the cluster [i]

       temp2<-as.data.frame(temp[,(Beg+1):ncol(temp)])

            for (j in 1:ncol(temp2))
            {
            temp2[,j]<-as.numeric(as.vector(temp2[,j]))
            }
            
            ms_moyen<-apply(temp2, MARGIN=2, mean)
            res[i,7:ncol(res)]<-ms_moyen

            ms_moyen<-sort(ms_moyen, decreasing=TRUE)
            ms_moyen_8<-substr(names(ms_moyen)[1:8], 1,4)

            MOY<-ms_moyen_8[1]
                for (mo in 2:8)
                {
                MOY<-paste(MOY, "/",ms_moyen_8[mo])
                }
                
            res[i,6]<-MOY
        }
        
        
res[,5]<-resPclus

## check which clusters are homogeneous

clus_ok<-vector()
l<-1

    for (k in 1:nrow(res))
        if (as.numeric(res[k,4])<varRT)
        {
        clus_ok[l]<-rownames(res)[k]
        l<-l+1
        }
print(paste("cluster ok length=",length(clus_ok),sep=" "))
if (length(clus_ok)==0) stop(cat("No cluster OK found, please increase Nbc value or varRT \n"), call.=FALSE)

m<-match(clus_ok,clusters)
mn<-1:length(clusters)
clus_pb<-mn[-m]

dat<-subset(data_tot_temp, data_tot_temp[,1]==clus_ok[1])

        for (iii in 2:length(clus_ok))
        {
        dat<-rbind(dat, subset(data_tot_temp, data_tot_temp[,1]==clus_ok[iii]))
		
        }
        
colnames(dat)[1]<-"num_cluster"


## construct the fingerprinting/profiling matrix

mol<-levels(as.factor(dat[,1]))
In<-levels(as.factor(dat[,2]))

mat_dat<-matrix(nrow=length(In), ncol=length(mol))
rownames(mat_dat)<-In
colnames(mat_dat)<-mol

mat_dat2<-matrix(nrow=length(In), ncol=length(mol))
rownames(mat_dat2)<-In
colnames(mat_dat2)<-mol


## modified Patch for quantification Florence & Cyrille 280611:
       for (i in 1:length(In))
       {
        tempo<-subset(dat, dat[,2]==In[i])
		tempMol<-as.factor(as.vector(tempo[,1]))

        vect<-is.na(match(mol, tempMol))# les pics présents dans l'échantillon
                if (quant ==FALSE)
                for (j in 1:length(vect))
                    if (vect[j]==TRUE)
                    mat_dat[i,j]<-0
                    else (mat_dat[i,j]<-1)
                    
                if (quant ==TRUE)
                {
                for (j in 1:length(vect))
                    if (vect[j]==TRUE)
                    mat_dat[i,j]<-0
                    else ( if (length(tempo[which(tempo$num_cluster==mol[j]),4])>1)
				mat_dat[i,j]<-NA
			   else(mat_dat[i,j]<-tempo[which(tempo$num_cluster==mol[j]),4])
			  )
                    
                for (j in 1:length(vect))
                    if (vect[j]==TRUE)
                    mat_dat2[i,j]<-0
                    else ( if (length(tempo[which(tempo$num_cluster==mol[j]),5])>1)
				mat_dat2[i,j]<-NA
			   else(mat_dat2[i,j]<-tempo[which(tempo$num_cluster==mol[j]),5])
			  )    
                }
}               



## construct a pdf with the histograms of retention times for each inhomogeneous cluster
if (length(clus_pb)>1)
{
pdf(file.path(Mypath,paste("Hist_cluster_problem_RT", answer[NC], ".pdf", sep="")))
for (i in 1:length(clus_pb))
{
hist(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_pb[i]])), xlab="retention time", ylab="number of analyses", main=paste("Distribution of retention times for the cluster", clus_pb[i], sep=" "))
mtext(paste("mean retention time =", round(mean(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_pb[i]]))), digits=3), " ; ", "range of retention times =", round((max(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_pb[i]])))-min(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_pb[i]])))), digits=3)))
}
dev.off()
}

## construct a pdf with the histograms of silhouette widths for each inhomogeneous cluster
if (length(clus_pb)>1)
{
pdf(file.path(Mypath,paste("Hist_cluster_problem_silhouette", answer[NC], ".pdf", sep="")))
for (i in 1:length(clus_pb))
{
hist(as.numeric(as.vector(data_tot_temp_new$Silhouette_indiv[data_tot_temp$cl==clus_pb[i]])), xlab="silhouette width", ylab="number of analyses", main=paste("Distribution of silhouette widths for the cluster", clus_pb[i], sep=" "))
mtext(paste("mean silhouette width =", round(mean(as.numeric(as.vector(data_tot_temp_new$Silhouette_indiv[data_tot_temp$cl==clus_pb[i]]))), digits=3), " ; ", "minimum silhouette =", round(min(as.numeric(as.vector(data_tot_temp_new$Silhouette_indiv[data_tot_temp$cl==clus_pb[i]]))), digits=3)))
}
dev.off()
}

## construct a pdf with the histograms of retention times for each homogeneous cluster
if (length(clus_ok)>1)
{
pdf(file.path(Mypath,paste("Hist_cluster_ok_RT", answer[NC], ".pdf", sep="")))
for (i in 1:length(clus_ok))
{
hist(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_ok[i]])), xlab="retention time", ylab="number of analyses", main=paste("Distribution of retention time for the cluster", clus_ok[i], sep=" "))
mtext(paste("mean retention time =", round(mean(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_ok[i]]))), digits=3), " ; ", "range of retention times =", round((max(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_ok[i]])))-min(as.numeric(as.vector(data_tot_temp_new$RT[data_tot_temp$cl==clus_ok[i]])))), digits=3)))
}
dev.off()
}

## construct a pdf with the histograms of silhouette widths for each homogeneous cluster

if (length(clus_ok)>1)
{
pdf(file.path(Mypath,paste("Hist_cluster_ok_silhouette", answer[NC], ".pdf", sep="")))
for (i in 1:length(clus_ok))
{
hist(as.numeric(as.vector(data_tot_temp_new$Silhouette_indiv[data_tot_temp$cl==clus_ok[i]])), xlab="silhouette width", ylab="number of analyses", main=paste("Distribution of silhouette widths for the cluster", clus_ok[i], sep=" "))
mtext(paste("mean silhouette width =", round(mean(as.numeric(as.vector(data_tot_temp_new$Silhouette_indiv[data_tot_temp$cl==clus_ok[i]]))), digits=3), " ; ", "minimum silhouette =", round(min(as.numeric(as.vector(data_tot_temp_new$Silhouette_indiv[data_tot_temp$cl==clus_ok[i]]))), digits=3)))
}
dev.off()
}

## construct the cluster matrix with 1 for homogeneous clusters and 0 for inhomogeneous clusters

clus_ok<-cbind(clus_ok, rep(1, length(clus_ok)))
clus_pb<-cbind(clus_pb, rep(0, length(clus_pb)))
clus_check<-as.matrix(rbind(clus_ok, clus_pb))
clus_check<-clus_check[order(as.numeric(clus_check[,1])),]
res<-cbind(clus_check[,1:2], res)
colnames(res)[1:2]<-c("cluster_number", "check_quality_cluster")



## add the 8 first mz for each cluster
##Cyr
assign("data_tot_temp_new_bug",data_tot_temp_new,envir=parent.frame())


ms_sort<-sort(data_tot_temp_new[1,8:ncol(data_tot_temp_new)], decreasing=TRUE)#8 replaced 6
ms_8<-substr(names(ms_sort)[1:8], 1,4)

for (li in 2:nrow(data_tot_temp_new))
{
ms_sort_temp<-sort(data_tot_temp_new[li,8:ncol(data_tot_temp_new)], decreasing=TRUE)#8 replaced 6
ms_8<-rbind(ms_8, substr(names(ms_sort_temp)[1:8], 1,4))
}

##Cyr
assign("ms8_bug",ms_8,envir=parent.frame())

## write the outfiles in Mypath

data_final_all<-cbind(data_tot_temp_new[,1:7], ms_8, data_tot_temp_new[8:ncol(data_tot_temp_new)])
if (quant ==TRUE)
{
	colnames(data_final_all)[6]<-"quantification1"
	colnames(data_final_all)[7]<-"quantification2"
}
write.table(res, file=file.path(Mypath,paste("output_cluster", answer[NC], ".txt", sep="")), row.names=FALSE)
write.table(data_final_all, file=file.path(Mypath,paste("output_peak", answer[NC], ".txt", sep="")), row.names=FALSE)

if (quant ==FALSE)
{
mat_dat<-cbind(rownames(mat_dat), mat_dat)
colnames(mat_dat)[1]<-"analysis"
write.table(mat_dat, file=file.path(Mypath,paste("output_fingerprintingmatrix", answer[NC], ".txt", sep="")), row.names=FALSE)
}
if (quant ==TRUE)
{
mat_dat<-cbind(rownames(mat_dat), mat_dat)
colnames(mat_dat)[1]<-"analysis"
write.table(mat_dat, file=file.path(Mypath,paste("output_profilingmatrix_quantification1_", answer[NC], ".txt", sep="")), row.names=FALSE)

mat_dat2<-cbind(rownames(mat_dat2), mat_dat2)
colnames(mat_dat2)[1]<-"analysis"
write.table(mat_dat2, file=file.path(Mypath,paste("output_profilingmatrix_quantification2_", answer[NC], ".txt", sep="")), row.names=FALSE)

}


#Rprof(NULL)
#summaryRprof(filename = "Rprof.out")$sampling.time
print(paste("Your working directory is :", eval(getwd())))
print(paste("Text and pdf files have been generated in the folder:", Mypath))
print("Done")

}
}

