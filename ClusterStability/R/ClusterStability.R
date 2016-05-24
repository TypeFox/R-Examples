#
# Functions for ClusterStability
# Authors: Etienne Lord, Matthieu Willems, Vladimir Makarenkov
# Since: December 2015-July 2015
#

# Function to return the Stirling numbers of the second kind
Stirling2nd<-function(n,k) {
	total=0;
	somme=0;
	for (j in 0:k) {
		
		tt=(-1)^(k-j)*choose(k,j)*j^n;
		somme=somme+tt;
	}
	total=(1/factorial(k))*somme;
	return(total)
}

# Internal function to return the p_n_k probability
# Note: Now use the recurence version from the R copula package
# to calculate the Stirling numbers of the second kind
p_n_k<-function(n,k) {
	no=copula::Stirling2(n-1,k);
	de=copula::Stirling2(n,k);
	node=no/de;
	if (is.na(node)||is.infinite(node)) return (1/k);
	return (node);
}

# Internal function to return the p_t_n_k probability
# Note: Now use the recurence version from the R copula package
# to calculate the Stirling numbers of the second kind
p_tilde_n_k<-function(n,k){
	no=copula::Stirling2(n-1,k-1);
	de=copula::Stirling2(n,k);
	node=no/de;
	if (is.na(node)||is.infinite(node)) return (0);
	return (node);
}

#Function to return if the member in group (index) are in the same partition
is_partition_group<-function(partition, group) {
	par<-partition[group[1]];
	for (i in group) {
		if (partition[i]!=par) return (FALSE);
	}
	return (TRUE);
}

#Function to calculate the singleton indice
calculate_singleton<-function(indices, partition, indice, total_indice) {
	total_singleton<-array(0, c(length(indices)));	
	#cat(total_indice);
	#cat(indices);
	#if (total_indice==0) total_indice=1;
	for (k in 1:length(partition)) {
		part<-as.vector(partition[[k]]);						
		a<-table(part);		
		#cat(part);
		#cat ("\nTable:", a[1],a[2],a[3],(a[1]+a[2]+a[3]),"\n");
		for (i in 1:length(a)) {
			if (a[i]==1) {
				#find the corresponding element in alone group
				for (j in 1:length(indices)) {
					if (part[j]==i) {
						if (!is.finite(indice[j])||is.na(indice[j])) indice[j]=0.0;						
						total_singleton[j]<-total_singleton[j]+indice[j];
					}
				}
			}
		}		
	}
	for (j in 1:length(indices)) {
		total_singleton[j]=total_singleton[j]/total_indice;		
	}	
	for (j in 1:length(indices)) {
		total_singleton[j]=max(total_singleton[j], 1-total_singleton[j]);
	}	
	return (total_singleton);
}


#Main ClusterStability function (approximative)
ClusterStability<-function(dat, k=3, replicate=1000, type='kmeans') {
	
	mylist<-list();
	dat=as.matrix(dat);
	len<-nrow(dat);
	
	partitions<-list();
	indices<-list();
	indice_list_ch<-array(0, c(replicate));
	indice_list_dunn<-array(0, c(replicate));
	indice_list_db<-array(0, c(replicate));
	indice_list_sil<-array(0, c(replicate));
	
	starts<-sample(1:10000000, replicate, replace=FALSE);
	total_calinski_harabasz=0;	
	total_silhouette=0;
	total_dunn=0;	
	total_db=0;
	
	for (i in 1:replicate) {
		if (type=='kmeans') {
			cluster<-Reorder(kmeans(dat,centers=k, nstart=1, iter.max=100, algorithm="MacQueen")$cluster);
		} else {		
			cluster<-Reorder(wcKMedoids(dist(dat),k,npass=0,cluster.only=TRUE));		
		}		
		indice_kmeans<-intCriteria(dat, cluster,c("Calinski_Harabasz","Dunn","Davies_Bouldin"))	
		total_calinski_harabasz<-total_calinski_harabasz+indice_kmeans[1]$calinski_harabasz;	
		total_dunn<-total_dunn+indice_kmeans[2]$dunn;	
		total_db<-total_db+indice_kmeans[3]$davies_bouldin;			
		ind<-summary(silhouette(cluster,dist(dat)))$avg.width;			
		if (is.nan(ind)) { 
			ind=0.0;			
		} else if (ind<0) {
			ind=(ind+1)/2;
		}
		set.seed(starts[i])
		total_silhouette<-total_silhouette+ind;				
		partitions[[i]]<-as.vector(cluster);
		indice_list_ch[i]<-indice_kmeans[1]$calinski_harabasz;
		indice_list_sil[i]<-ind;		
		indice_list_db[i]<-indice_kmeans[3]$davies_bouldin;	
		indice_list_dunn[i]<-indice_kmeans[2]$dunn;
		
	}
	
	r<-list("partition"=partitions, "calinski_harabasz"=indice_list_ch, "silhouette"=indice_list_sil,"total_calinski_harabasz"=total_calinski_harabasz, "total_silhouette"=total_silhouette, "dunn"=indice_list_dunn, "db"=indice_list_db,"total_dunn"=total_dunn, "total_db"=total_db);
	indices=1:nrow(dat);
	combinations=Kcombination(indices, 2);
	total_combination_ch<-calculate_indices(combinations, r$partition, r$calinski_harabasz, r$total_calinski_harabasz);	
	total_combination_sil<-calculate_indices(combinations, r$partition, r$silhouette, r$total_silhouette);
	total_combination_dunn<-calculate_indices(combinations, r$partition, r$dunn, r$total_dunn);	
	total_combination_db<-calculate_indices(combinations, r$partition, r$db, r$total_db);
	total_singletons_ch<-calculate_singleton(indices, r$partition,r$calinski_harabasz,r$total_calinski_harabasz);
	total_singletons_sil<-calculate_singleton(indices, r$partition,r$silhouette,r$total_silhouette);
	total_singletons_dunn<-calculate_singleton(indices, r$partition,r$dunn,r$total_dunn);
	total_singletons_db<-calculate_singleton(indices, r$partition,r$db,r$total_db);
	
	global_PSG_ch<-0;
	global_PSG_sil<-0;
	global_PSG_dunn<-0;
	global_PSG_db<-0
	total_PSG_dunn<-calculate_individual_PSG_approximative(k,combinations, total_singletons_dunn,total_combination_dunn, indices);	
	total_PSG_db<-calculate_individual_PSG_approximative(k,combinations, total_singletons_db,total_combination_db, indices);
	total_PSG_ch<-calculate_individual_PSG_approximative(k,combinations, total_singletons_ch, total_combination_ch, indices);	
	total_PSG_sil<-calculate_individual_PSG_approximative(k,combinations, total_singletons_sil,total_combination_sil, indices);
	global_PSG_dunn<-mean(total_PSG_dunn);
	global_PSG_db<-mean(total_PSG_db);
	global_PSG_ch<-mean(total_PSG_ch);
	global_PSG_sil<-mean(total_PSG_sil);
	
	return(list("ST_ch"=total_PSG_ch,
				"ST_sil"=total_PSG_sil,
				"ST_dunn"=total_PSG_dunn,
				"ST_db"=total_PSG_db,
				"ST_global_ch"=global_PSG_ch,
				"ST_global_sil"=global_PSG_sil, 
				"ST_global_dunn"=global_PSG_dunn,
				"ST_global_db"=global_PSG_db
				)
			); 
}

#Main ClusterStability function (exact)
ClusterStability_exact<-function(dat, k=3, replicate=1000, type='kmeans') {
	
	mylist<-list();
	dat=as.matrix(dat);
	len<-nrow(dat);
	
	partitions<-list();
	indices<-list();
	indice_list_ch<-array(0, c(replicate));
	indice_list_dunn<-array(0, c(replicate));
	indice_list_db<-array(0, c(replicate));
	indice_list_sil<-array(0, c(replicate));
	
	starts<-sample(1:10000000, replicate, replace=FALSE);
	total_calinski_harabasz=0;	
	total_silhouette=0;
	total_dunn=0;	
	total_db=0;
	for (i in 1:replicate) {
		if (type=='kmeans') {
			cluster<-Reorder(kmeans(dat,centers=k, nstart=1, iter.max=100, algorithm="MacQueen")$cluster);
		} else {		
			cluster<-Reorder(wcKMedoids(dist(dat),k,npass=0,cluster.only=TRUE));		
		}		
		indice_kmeans<-intCriteria(dat, cluster,c("Calinski_Harabasz","Dunn","Davies_Bouldin"))	
		total_calinski_harabasz<-total_calinski_harabasz+indice_kmeans[1]$calinski_harabasz;	
		total_dunn<-total_dunn+indice_kmeans[2]$dunn;	
		total_db<-total_db+indice_kmeans[3]$davies_bouldin;			
		ind<-summary(silhouette(cluster,dist(dat)))$avg.width;			
		if (is.nan(ind)) { 
			ind=0.0;			
		} else if (ind<0) {
			ind=(ind+1)/2;
		}
		set.seed(starts[i])
		total_silhouette<-total_silhouette+ind;				
		partitions[[i]]<-as.vector(cluster);
		indice_list_ch[i]<-indice_kmeans[1]$calinski_harabasz;
		indice_list_sil[i]<-ind;		
		indice_list_db[i]<-indice_kmeans[3]$davies_bouldin;	
		indice_list_dunn[i]<-indice_kmeans[2]$dunn;
		
	}
	r<-list("partition"=partitions, "calinski_harabasz"=indice_list_ch, "silhouette"=indice_list_sil,"total_calinski_harabasz"=total_calinski_harabasz, "total_silhouette"=total_silhouette, "dunn"=indice_list_dunn, "db"=indice_list_db,"total_dunn"=total_dunn, "total_db"=total_db);
	indices=1:nrow(dat);
	combinations=Kcombination(indices, 2);
	total_combination_ch<-calculate_indices(combinations, r$partition, r$calinski_harabasz, r$total_calinski_harabasz);	
	total_singletons_ch<-calculate_singleton(indices, r$partition,r$calinski_harabasz,r$total_calinski_harabasz);
	
	total_combination_sil<-calculate_indices(combinations, r$partition, r$silhouette, r$total_silhouette);
	total_singletons_sil<-calculate_singleton(indices, r$partition,r$silhouette,r$total_silhouette);
	
	total_combination_dunn<-calculate_indices(combinations, r$partition, r$dunn, r$total_dunn);	
	total_singletons_dunn<-calculate_singleton(indices, r$partition,r$dunn,r$total_dunn);
	
	total_combination_db<-calculate_indices(combinations, r$partition, r$db, r$total_db);
	total_singletons_db<-calculate_singleton(indices, r$partition,r$db,r$total_db);
	
	global_PSG_ch<-0;
	global_PSG_sil<-0;
	global_PSG_dunn<-0;
	global_PSG_db<-0;
	
	#k, r_combinations, total_indices, indices
	#Warning message 
	if (is.na(Stirling2(nrow(dat), k))||is.infinite(Stirling2(nrow(dat), k))) {
		msg<-paste("Warning, the current values of (n=",nrow(dat),") and (k=",k,") are greater than the supported values for this function. The use of its approximate version 'ClusterStability' is recommended in this case.");
		warning(msg)
	}
	pnk=p_n_k(nrow(dat), k);
	pnktilde=p_tilde_n_k(nrow(dat), k);
	total_PSG_dunn<-calculate_individual_PSG_exact(k,combinations, total_singletons_dunn,total_combination_dunn, indices,pnk,pnktilde);	
	total_PSG_db<-calculate_individual_PSG_exact(k,combinations, total_singletons_db,total_combination_db, indices,pnk,pnktilde);
	total_PSG_ch<-calculate_individual_PSG_exact(k,combinations,  total_singletons_ch,total_combination_ch, indices,pnk,pnktilde);	
	total_PSG_sil<-calculate_individual_PSG_exact(k,combinations, total_singletons_sil,total_combination_sil, indices,pnk,pnktilde);
	global_PSG_dunn<-mean(total_PSG_dunn);
	global_PSG_db<-mean(total_PSG_db);
	global_PSG_ch<-mean(total_PSG_ch);
	global_PSG_sil<-mean(total_PSG_sil);
	
	return(list("ST_ch"=total_PSG_ch,
				"ST_sil"=total_PSG_sil,
				"ST_dunn"=total_PSG_dunn,
				"ST_db"=total_PSG_db,
				"ST_global_ch"=global_PSG_ch,
				"ST_global_sil"=global_PSG_sil, 
				"ST_global_dunn"=global_PSG_dunn,
				"ST_global_db"=global_PSG_db
				)
			); 
}
