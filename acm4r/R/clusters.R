#' Clustering
#'
#' Clustering for exact matching and BABLBS matching
#'
#' @param input is the result of the call_acm function in the format of ULI ULI numbands numbands num_matches
#' @param bablbs is the result of the call_bablbs, call_gd1 or  call_gd2 function in the format of ULI ULI
#' @param work_dir is where the datasets should be stored
#' @param type value ("exact", "bablbs", "gd1", "gd2") indicating if you want exact matching or BABLBS/Genetic Distance. Default is "exact" i.e. exact matching. 
#' @return A list with 7 components.  The first (SINGLE) is a 3 column matrix of all fingerprint IDs that do not belong to a cluster. The columns correspond to the Cluster_Number, the Cluster_size, and the Fingerprint ID. The second (CLUSTERED) is a matrix of all fingerprint IDs that do belong to a cluster. The columns correspond to the Cluster_Number, the Cluster_size, and the Fingerprint IDs that belong to that cluster. The third (BOTH) combines the others into one matrix. The fourth and fifth calculate RTIN and RTIn-1.  The last two are used for the histograms that are produced by a call to this function. 
#' @references 
#' Salamon et. al (1998) Accommodating Error Analysis in Comparison and Clustering of Molecular Fingerprints. Emerging Infectious Diseases Vol. 4, No. 2, April-June 1998
#' @references Abasci LLC. JAMES v1.0 User Documentation. 2002. 
#' @author Andrea Benedetti \email{andrea.benedetti@@mcgill.ca}
#' @author Sahir Rai Bhatnagar
#' @author XiaoFei Zhao
#' @examples 
#' \dontshow{
#' (WD <- getwd())
#' if (!is.null(WD)) setwd(WD)
#' data(replicates.in)
#' write.table(replicates.in,  "replicates.in", quote=FALSE, row.names=FALSE)
#' data(experiments.in)
#' write.table(experiments.in, "experiments.in", quote=FALSE, row.names=FALSE)
#' call_erra("replicates.in",  dnum=1, sd=1, delete=TRUE)
#' res1<-call_acm("experiments.in")
#' res_bab<-call_bablbs(res1)
#' res_gd1<-call_gd1(res1)
#' }
#' #synthesize the results
#' # Exact matching clusters
#' exact<-clusters(input=res1,type="exact")
#' names(exact)
#' exact$RTIN
#' exact$RTIN1
#' # Clustering based on BABLBS
#' bablbs<-clusters(input=res1, bablbs=res_bab,type="bablbs")
#' names(bablbs)
#' bablbs$RTIN
#' bablbs$RTIN1
#' # Clustering based on GD1
#' gd1<-clusters(input=res1, bablbs=res_gd1,type="gd1")
#' names(gd1)
#' gd1$RTIN
#' gd1$RTIN1
#' @export 
clusters<-function(input, bablbs, type="exact", work_dir="") {
#returns the subjects in their clusters.
#type means it is for type matching and takes a result file with uli1 uli2 numbands1 
#numbands2 nummatching bands as input
#else takes a file with pairs of ulis
	
	if(type=="exact"){
		
		#first get fpts with identical matches
		idents<-input[input[,3]==input[,4] & input[,3]==input[,5],];idents
		
		#list of all patient ids
		patient_ids<-unique(matrix(cbind(as.character(input[,1]),as.character(input[,2]))),ncol=1)
	
		#get the ids of fpts with an identical match
		has.ident<-unique(matrix(cbind(as.character(idents[,1]),as.character(idents[,2])),ncol=1))
		max.size<-nrow(has.ident)
		resu<-matrix(NA,ncol=(max.size+1),nrow=(max.size+nrow(patient_ids)))
		
		#keep a count as well
		count<-matrix(0, nrow=(max.size+nrow(patient_ids)),ncol=1)
		count[1]<-2
			}
	else {
		
		#list of all patient ids
		patient_ids<-unique(matrix(cbind(as.character(input[,1]),as.character(input[,2]))),ncol=1)
		
		#input has differnt format ULI1 ULI2 -- result from call_bablbs
		idents<-bablbs
		has.ident<-unique(matrix(cbind(as.character(idents[,1]),as.character(idents[,2])),ncol=1));has.ident
		max.size<-nrow(has.ident)
		resu<-matrix(NA,ncol=(max.size+1),nrow=(max.size+nrow(has.ident)))
		
		#keep a count as well
		count<-matrix(0, nrow=(max.size+nrow(has.ident)),ncol=1)
		count[1]<-2
		}
	
	j<-1
	#start it off
	
	#cluster number
	resu[1,1]<-1
	
	#### Incomplete Fix: if every patient is in a singleton then do not continue ####
	if (dim(idents)[1] < 1 | dim(idents)[2] < 1) {
		print ("No matching is found, RTI(n) = 0, RTI(n-1) = 0")
		print ("No recent transmission occurred.")
		return (0)
	}
	
	#first two members
	resu[1,2]<-as.character(idents[1,1])
	resu[1,3]<-as.character(idents[1,2])
	
	##### Fix: RTIn-1 is based on number of clusters of size of at least 2 #####
	zxf.num.of.clusters <- 0

	for(i in 2:nrow(idents)){
		#if they are both already in there than just carry on
		if( !( idents[i,1] %in% resu[1:j,2:(max.size+1)] &  idents[i,2] %in% resu[1:j,2:(max.size+1)])) {
			
			if(  idents[i,1] %in% resu[1:j,2:(max.size+1)]) {
				#get the row where it it
				pos.row<-which(resu==idents[i,1],arr.ind=TRUE)[,1]
				resu[pos.row,count[pos.row]+2]<-as.character(idents[i,2])
				count[pos.row]<-count[pos.row]+1
			}
			else if( idents[i,2] %in% resu[1:j,2:(max.size+1)]) {
				pos.row<-which(resu==idents[i,2],arr.ind=TRUE)[,1]
				resu[pos.row,count[pos.row]+2]<-as.character(idents[i,1])
				count[pos.row]<-count[pos.row]+1
			}
			else if(!  idents[i,1] %in% resu[1:j,2:(max.size+1)] & !  idents[i,2] %in% resu[1:j,2:(max.size+1)]) {
				j<-j+1
				count[j]<-2
				resu[j,1]<-j
				resu[j,2]<-as.character(idents[i,1])
				resu[j,3]<-as.character(idents[i,2])
				##### Fix: RTIn-1 is based on number of clusters of size of at least 2 #####
				zxf.num.of.clusters <- zxf.num.of.clusters + 1
			}}
		
	}
		
	resu2<-cbind(count,resu)
	max.clus<-max(count)
	

	#get rid of empty rows
	resu3<-resu2[resu2[,1]>0,1:(max.clus+2)];resu3
	
	#re-arrange columns
	a<-resu3[,2]; b<-resu3[,1]; c<-resu3[,3:(max.clus+2)]
	resu4<-cbind(resu3[,2],resu3[,1],resu3[,3:(max.clus+2)])
	
	#------------data frame of everyone who is in a cluster------------------------------------------------
	clustered<-cbind(a,b,c)
	emptyy<-rep(0,ncol(clustered)-2)
	uli_names<-for(i in 1:(ncol(clustered)-2)){emptyy[i]<-paste("ULI_",i,sep="")}
	colnames(clustered)<-c("Cluster_Number","Cluster_Size",emptyy)
	
	datpath1 = paste(work_dir, ifelse(type=="exact","exact_matching_",ifelse(substr(type,1,2)=="gd",
							paste("Genetic_Distance_",noquote(substr(type,3,3)),sep=""),"BABLBS_")), "clustered",sep="")
	write.matrix(clustered, file = datpath1, sep = "\t")
	file.show(datpath1, title = ifelse(type=="exact","Exact Matching Clusters - ULI's that are in a Cluster",
					ifelse(substr(type,1,2)=="gd",
							paste("Genetic Distance ",noquote(substr(type,3,3))," Clusters - ULI's that are in a Cluster",sep=""),"BABLBS Clusters - ULI's that are in a Cluster")))
	
	#--------------to create matrix of Singletons------------------------------------------------------
	
	uli_unique<-matrix(setdiff(patient_ids,has.ident),ncol=1)	
			
	max_cluster<-max(as.integer(clustered[,1]));max_cluster
	row_singletons<-nrow(uli_unique);row_singletons
	mat_na<-matrix(NA, nrow = row_singletons, ncol = (ncol(clustered)-3));mat_na

	num_rows<-c((max_cluster+1):(max_cluster+row_singletons))
	singletons<-(cbind(c(num_rows),c(rep(1,row_singletons)),uli_unique))
	colnames(singletons)<-c("Cluster_Number","Cluster_Size","ULI");singletons
	
	datpath3 = paste(work_dir, ifelse(type=="exact","exact_matching_",ifelse(substr(type,1,2)=="gd",paste("Genetic_Distance_",noquote(substr(type,3,3)),sep=""),"BABLBS_")), 
			"Singletons",sep="")
	write.matrix(singletons, file = datpath3, sep = "\t")
	file.show(datpath3, title = ifelse(type=="exact","Exact Matching Clusters - ULI's that are NOT in a Cluster",
					ifelse(substr(type,1,2)=="gd",paste("Genetic Distance ",noquote(substr(type,3,3))," Clusters - ULI's that are NOT in a Cluster",sep=""),
							"BABLBS Clusters - ULI's that NOT are in a Cluster")))
			
	
	#--------------to create matrix of Singletons with Clustered ULI's----------------------------------
	jj<-cbind(singletons,mat_na);jj
	safe<-rbind(resu4,jj);safe
	clus_single<-safe
	
	emptyy2<-rep(0,ncol(clustered)-2)
	uli_names2<-for(i in 1:(ncol(clustered)-2)){emptyy2[i]<-paste("ULI_",i,sep="")}
	colnames(clus_single)<-c("Cluster_Number","Cluster_Size",emptyy2)
	
	datpath2 = paste(work_dir, ifelse(type=="exact","exact_matching_",ifelse(substr(type,1,2)=="gd",paste("Genetic_Distance_",noquote(substr(type,3,3)),sep=""),"BABLBS_")), 
					"clustered_and_single",sep="")
	write.matrix(clus_single, file = datpath2, sep = "\t")
	file.show(datpath2, title = ifelse(type=="exact","Exact Matching Clusters - All ULI's",
					ifelse(substr(type,1,2)=="gd",paste("Genetic Distance ",noquote(substr(type,3,3))," Clusters - All ULI's",sep=""),"BABLBS Clusters - All ULI's")))
	
	#-----------Matrix of ID, Cluster ID, Belongs to a cluter? (0=FALSE, 1=TRUE)------------------------
	
	max.size1<-nrow(patient_ids);max.size1
	max.size2<-nrow(clus_single);max.size2
	
	mat_logical<-matrix(NA,ncol=3,nrow=max.size1);colnames(mat_logical)<-c("ID","Cluster_ID","Clustered");mat_logical
	
	for (i in 1:max.size1){mat_logical[i,1]<-patient_ids[i]};mat_logical
	
	
	for (i in 1:max.size1)
	{
		ID_temp<-mat_logical[i,1]
		for (j in 1:max.size2)
		{
			#if its in the jth row
			if (ID_temp %in% clus_single[j,]) 
			{
				mat_logical[i,2]<-clus_single[j,1]
				ifelse(clus_single[j,2]>1,mat_logical[i,3]<-1,mat_logical[i,3]<-0)
			}
		}
	}
	
	datpath4 = paste(work_dir, ifelse(type=="exact","exact_matching_",ifelse(substr(type,1,2)=="gd",paste("Genetic_Distance_",noquote(substr(type,3,3)),sep=""),"BABLBS_")), 
			"ID distribution",sep="")
	write.matrix(mat_logical, file = datpath4, sep = "\t")
	file.show(datpath4, title = ifelse(type=="exact","Exact Matching Clusters - ID Distribution",
					ifelse(substr(type,1,2)=="gd",paste("Genetic Distance ",noquote(substr(type,3,3))," Clusters - ID Distribution",sep=""),"BABLBS Clusters - ID Distribution")))
	
	#--------------RTI(n) and RTI (n-1)------------------------------------
	
	rtin<-matrix(nrow(has.ident)/nrow(patient_ids), ncol=1,dimnames=list(paste("RTIn = #clustered/#patients = ",
							nrow(has.ident),"/",nrow(patient_ids)," = ",sep=""),NULL))
	
	rtin1<-matrix((nrow(has.ident)-max.clus)/nrow(patient_ids), ncol=1,
			dimnames=list(paste("RTIn-1 = (#clustered - clust_size>1)/#patients = (",nrow(has.ident),"-",max.clus,")/",nrow(patient_ids)," = ",NULL)))
	
	##### Fix: RTIn-1 is based on number of clusters of size of at least 2 #####
	rtin1<-matrix((nrow(has.ident)-zxf.num.of.clusters)/nrow(patient_ids), ncol=1,
			dimnames=list(paste("RTIn-1 = (#clustered - clust_size>1)/#patients = (",nrow(has.ident),"-",max.clus,")/",nrow(patient_ids)," = ",NULL)))
	
	#------Histograms------------------------------------------------------------------------
	#-------------ULI's that belong to a Cluster---------------------------------------
	dev.new()
	par(mfrow=c(1,2))
	freq_clus<-as.numeric(clustered[,2])
		
	hist(freq_clus, col="lightblue", breaks=c(-1:max.clus+1), 
		main = ifelse(type=="exact","Exact Matching Clustering:\n ULI's that belong to a Cluster",
				ifelse(substr(type,1,2)=="gd",paste("Genetic Distance ",noquote(substr(type,3,3))," Clusters - \n ULI's that belong to a Cluster",sep=""),
						"BABLBS Clustering\n ULI's that belong to a Cluster")),
		xlab="Cluster Size", font.lab=2)


	#-------------ULI's that belong to a Cluster with Singletons---------------------------------------
	clus_vector<-as.numeric(clustered[,2])
	freq_both<-append(clus_vector,as.numeric(singletons[,2]))
	
	hist(freq_both, col="lightblue", breaks=c(-1:max.clus+1),
			main = ifelse(type=="exact","Exact Matching Clustering:\n Clustered and Singletons",
					ifelse(substr(type,1,2)=="gd",paste("Genetic Distance ",noquote(substr(type,3,3))," Clustering:\n Clustered and Singletons",sep=""),"BABLBS Clustering:\n Clustered and Singletons")),
			xlab="Cluster Size", font.lab=2, plot=TRUE)
	legend("top",legend=c(paste("RTI(n) = ",round(rtin, digits=4),sep=""),paste("RTI(n-1) = ",round(rtin1,digits=4),sep="")), 
			bg="white", bty = "n", text.width = 2, xjust = 1, yjust = 1)

	
	#---------Print Histogram to PDF---------------------------------------------- 
	pdf(paste(work_dir,"hist_both.pdf",sep=""),7,7,pointsize=10)
	par(lend=1, ljoin=1, lmitre=4)
	hist(freq_both, col="lightblue", breaks=c(-1:max.clus+1),
			main = ifelse(type=="exact","Exact Matching Clustering:\n Clustered and Singletons",
					ifelse(substr(type,1,2)=="gd",paste("Genetic Distance ",noquote(substr(type,3,3))," Clustering:\n Clustered and Singletons",sep=""),"BABLBS Clustering:\n Clustered and Singletons")),
			xlab="Cluster Size", font.lab=2)
	legend("top",legend=c(paste("RTI(n) = ",round(rtin, digits=4),sep=""),paste("RTI(n-1) = ",round(rtin1,digits=4),sep="")), 
			bg="white", bty = "n", text.width = 2, xjust = 1, yjust = 1)
	dev.off()
	
	#---------Print Histogram to PDF---------------------------------------------- 
	pdf(paste(work_dir,"hist_clust.pdf",sep=""),7,7,pointsize=10)
	par(lend=1, ljoin=1, lmitre=4)
	hist(freq_clus, col="lightblue", breaks=c(-1:max.clus+1), 
			main = ifelse(type=="exact","Exact Matching Clustering:\n ULI's that belong to a Cluster",
					ifelse(substr(type,1,2)=="gd",paste("Genetic Distance ",noquote(substr(type,3,3))," Clustering\n ULI's that belong to a Cluster",sep=""),
							"BABLBS Clustering\n ULI's that belong to a Cluster")),
			xlab="Cluster Size", font.lab=2)
	legend("top",legend=c(paste("RTI(n) = ",round(rtin, digits=4),sep=""),paste("RTI(n-1) = ",round(rtin1,digits=4),sep="")), 
			bg="white", bty = "n", text.width = 2, xjust = 1, yjust = 1)
	dev.off()
	
		
	#--------------output-------------------------------------------------------	
	list1<-list(SINGLE=singletons,CLUSTERED=clustered,BOTH=clus_single, RTIN=rtin, RTIN1=rtin1, DISTRIBUTION=mat_logical, HIST_BOTH=freq_both)

			
}







