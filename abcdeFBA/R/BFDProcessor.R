### processes the output files from batch simulations. removes single gene deletion fatalities
BFD_Processor<-function(fba_object,EXSDR=list())
{
if(length(EXSDR)==0)
	{
	EXSDR<-Exhaustive_single_deletion(fba_object)
	flux_vec=EXSDR$lethal_dels
	}

flux_vec=EXSDR$lethal_dels

print("Put the batch files named in the format results(num) in a folder called BKO and put in the system directory before continuing")
sim4=matrix(c(0,0,0),nrow=1)
num_files<-as.numeric(readline("Enter the number of input files"))
for(i in 1:num_files)
	{
	sim4<-rbind(sim4,read.table(paste("./BKO/results",i,sep="")))
	}
	Sim4<-sim4[which(sim4[,3]!=1),]
for(i in 1:length(flux_vec))
	{
	Sim4<-Sim4[which(Sim4[,1]!=flux_vec[i]),]
	Sim4<-Sim4[which(Sim4[,2]!=flux_vec[i]),]
	}

Sim5<-as.matrix(Sim4)
Sim5<-rbind(Sim5,cbind(Sim5[,2],Sim5[,1],Sim5[,3]))
pairwise_fatality_redundant<-cbind(fba_object$reaction_list[Sim5[,1]],fba_object$reaction_list[Sim5[,2]])
pairwise_fatality_unique<-cbind(fba_object$reaction_list[Sim4[,1]],fba_object$reaction_list[Sim4[,2]])
write.table(pairwise_fatality_redundant,file="Fatal_Knockout_pairs_redundant.xls",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(pairwise_fatality_unique,file="Fatal_Double_knockouts_unique.xls",sep="\t",row.names=FALSE,col.names=FALSE)

set1_uniq<-unique(pairwise_fatality_unique[,1])
set2_uniq<-unique(pairwise_fatality_unique[,2])
freq_reac1=0
freq_reac2=0

for(i in 1:length(set1_uniq))
	{freq_reac1[i]=length(which(set1_uniq[i]==pairwise_fatality_unique))
	+length(which(set2_uniq[i]==pairwise_fatality_unique))}

for(i in 1:length(set2_uniq))
	{freq_reac2[i]=length(which(set2_uniq[i]==pairwise_fatality_unique))
	+length(which(set1_uniq[i]==pairwise_fatality_unique))}

png("Freq_curve_bko_log.png")
plot(rev(sort(freq_reac1)),log="xy",type="p",col="blue",ylab="Synthetic Double Lethality Pairing Frequency",xlab="Reaction ID",main="Log-log plot of Lethal Interactions")
dev.off()

#source("../system/dko_freq_anlyz.R")
################################################# ANALYTICAL SECTION #######################################################
subsystem_wise_fatality<-cbind(fba_object$sub_system[Sim5[,1]],fba_object$sub_system[Sim5[,2]])
write.table(subsystem_wise_fatality,file="subsys_pairs_redundant",sep="\t",row.names=FALSE,col.names=FALSE,quote=F)
print("Creating adjacency matrix using Subsystem Classifications")
################ Check the file for duplicate entries caused due to Case sensitivity and remove #################
subsystem_wise_fatality<-read.delim("subsys_pairs_redundant",header=F)

freq_subsys_reac1=0
freq_subsys_reac2=0

set1_subsys_uniq<-sort(unique(subsystem_wise_fatality[,1]))
set2_subsys_uniq<-sort(unique(subsystem_wise_fatality[,2]))

subsys_adjacency=matrix(0,length(set1_subsys_uniq),length(set2_subsys_uniq))

for(i in 1:length(set1_subsys_uniq))
	{
	for(j in 1:length(set2_subsys_uniq))
		{
		subsys_adjacency[i,j]=length(which(subsystem_wise_fatality
		[which(subsystem_wise_fatality[,1]==set1_subsys_uniq[i]),2]
		==set2_subsys_uniq[j]))
		}
	}

write.table(set1_subsys_uniq,file="Subsys Legend",quote=F,col.names=F,sep="\t")

################################################ HEATMAPPER AND CORRELATION PLOTTERS #######################################
message("Normalizing Data for Heatmaps and Correlation Plots")
unique_subsystems<-sort(unique(fba_object$sub_system))

cross_reac<-matrix(0,length(set1_subsys_uniq),length(set2_subsys_uniq))

for(i in 1:length(set1_subsys_uniq))
	{
	for(j in 1:length(set2_subsys_uniq))
		{
		cross_reac[i,j]<-length(which(fba_object$sub_system==
		set1_subsys_uniq[i]))*length(which(fba_object$sub_system
		==set2_subsys_uniq[j]))
		}
	}

normalized_subsys_adjacency<-subsys_adjacency/cross_reac
pdf("Heatmaps.pdf")
levelplot(normalized_subsys_adjacency)
heatmap(normalized_subsys_adjacency)
dev.off()

norm_subsys_cor<-cor(normalized_subsys_adjacency)
rownames(norm_subsys_cor)=set1_subsys_uniq
colnames(norm_subsys_cor)=set1_subsys_uniq
pdf("Correlation_Plot.pdf")
corrplot(norm_subsys_cor,addgrid.col=NULL,order="hclust")
dev.off()
}
############################################################################################################################
