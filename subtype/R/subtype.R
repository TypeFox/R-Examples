subtype <-
function(
   GEset,
   outcomeLabels, 
   treatment = NULL,
   Npermutes = 10,  #how many times to permute the whole probe list
   Nchunks = 25,    #how many parts of the probe list
   minClusterSizeB = 20, #min no. of patients per subtype 
   NclustersASet=100, #cut tree branches of variable size:
   FDRpermutation=TRUE,
   nFDRperm=50,
   seed=NULL,
   testMode = "quick", #only measure p-value distribution from t-test, no cross-validation
   survivaltimes = NULL,
   method = "penalized", # or "dlda" or "lasso"
   top_best_probes=100, #by t-test, to be input for penalized / dlda
   Niter = 20, #iterations of { (TrainingSet, TestSet)->training->test->recordResults }
   showMovie = 0, #display RUC/Surv curves and heatmaps
   redefineSubtypeMembers = 0, #detect subtype members after every hold-out
   holdOut = 10 #out of the subtype, i.e. Nsubtype - holdOut = Ntraining_set
)
{

#require(rsmooth);
library(penalized);
library(ROCR);
if (!is.null(seed)) {set.seed(seed)}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Record ALL/SUCCESSFUL sets of probes and patients:

GenesDefiningSubtypes_allMembers<-c(1); #probes
GenesDefiningSubtypes<-list(empty=c("")); #probes
SubtypeSpecificPredictors<-list(empty=c("")); #probes
SubtypeSpecificSignatures<-list(empty=c(""));
geneSetSignatures<-c(1); #merged probe lists
SubtypePatients<-list(empty=c("")); #patients
SubtypeFDR01<-list(empty=c("")); #permutated FDRs
SubtypeFDR02<-list(empty=c("")); #permutated FDRs
SubtypeFDR03<-list(empty=c("")); #permutated FDRs
SubtypeBestByPvalue<-list(empty=c("")); #probes, only in 'quick' mode
FakeID<-NULL;


# and results overall:
resultsAll<-NULL;  #t-test only
subtResults<-NULL; #t-test AND cross-validation results
trtResults<-list(empty=c(""));
subtypeSpecificProbes<-list(empty=c(""));
BestSubtypeGenes<-list(empty=c(""));

chunkSize = floor(nrow(GEset) / Nchunks);
Npatients = ncol(GEset);


if (is.null(rownames(GEset)) | is.null(colnames(GEset)) )  {stop("The row and column of the dataset should have their names.")}
if (!is.null(rownames(GEset))) {probe_list<-rownames(GEset)}
if (!is.null(colnames(GEset))) {Patients<-colnames(GEset)}



relapse_labels <- outcomeLabels; names(relapse_labels)<-Patients;


if (testMode=="normal"){
relapse_times <- survivaltimes; names(relapse_times)<-Patients
trt<-treatment; names(trt)<-Patients
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for (permute in c(1:Npermutes)) {


orderunif<-order(runif(c(1:nrow(GEset))))
probe_list<-rownames(GEset)[orderunif]


for (chunk in c(1:Nchunks)) {

# FIND GENE (=probe) SETS THAT CAN DELINEATE 1 or, better, 2 subtypes in the patient population:
# start from the entire original set of probes:

classifying_set<-probe_list[c((chunkSize * (chunk - 1) + 1):(chunk  * chunkSize))]
if (chunk==Nchunks) {classifying_set<-probe_list[c((chunkSize * (chunk - 1) + 1):nrow(GEset))]}


a1<-hclust(dist(GEset[classifying_set,],method="euclidean"))


#cut tree branches of variable size:


for (NclustersA in NclustersASet) {

ctreeA<-cutree(a1, k=NclustersA);
tcut<-2*(chunkSize/NclustersA)


if (length(which(table(ctreeA)> tcut ))>0){
uuctreeA<-as.numeric(names(table(ctreeA)[which(table(ctreeA)> tcut )]))
}
if (length(which(table(ctreeA)> tcut ))==0){
uuctreeA<-as.numeric(names(table(ctreeA)[which.max(table(ctreeA))]))
}

for (clA in uuctreeA) {
subtype_genes_set<-names(ctreeA)[which(ctreeA == clA)];

if (length(subtype_genes_set) > tcut) {   ## 10: require (at least) 10 genes, but arbitrary
fingerprint = paste(sort(subtype_genes_set), collapse="");

if (is.na(geneSetSignatures[fingerprint])) {
geneSetSignatures[fingerprint] = 1;
b1<-hclust(dist(t(GEset[subtype_genes_set,]),method="euclidean"));



#NclustersB = 1;
#clusSizes<-c(minClusterSizeB+11,minClusterSizeB+10);
#while (clusSizes[2] > minClusterSizeB & (NclustersB < (Npatients - holdOut))) {
#NclustersB = NclustersB + 1;
#ctreeB<-cutree(b1, k=NclustersB);
#clusSizes<-sort(table(ctreeB), decreasing = T);
#if (clusSizes[2]<minClusterSizeB) {ctreeB<-cutree(b1, k=NclustersB-1);}
#}

Bcut<-ncol(GEset)/minClusterSizeB;
ctreeB<-cutree(b1, k= (Bcut+1));
clusSizes<-sort(table(ctreeB), decreasing = T);


nsubt<-sum(table(ctreeB)>=minClusterSizeB);
ctreeB12<-names(sort(table(ctreeB), decreasing = T))[1:nsubt];


#subtype1<-names(ctreeB)[which(ctreeB == ctreeB12[1])]; ## We will compare the patients within each subtype (1 or 2).
#subtype2<-names(ctreeB)[which(ctreeB == ctreeB12[2])];


###############################3
for (subt in nsubt) {

#if (subt == 1) {subtypeN = subtype1} else {subtypeN = subtype2}
subtypeN<-names(ctreeB)[which(ctreeB == ctreeB12[subt])]


if (length(unique(relapse_labels[as.character(subtypeN)]))>1 ) {

#tOC<-array(rep(NA, times=length(rownames(GEset))));
tOC<-tstatistics2(GEset[,subtypeN], relapse_labels[as.character(subtypeN)])$tstat;
names(tOC)<-rownames(GEset);
dfree<-length(c(subtypeN)) - 2;
pt1<-2*pt(abs(tOC), dfree, lower.tail = FALSE);
FDR.Res<-ELF(xdat=GEset[,subtypeN], grp=relapse_labels[as.character(subtypeN)],nperm=0);  ##### crude FDR for filtering out !
pNFDR01 = sum(FDR.Res$ELF<0.1)  


pValue_Area = low_p_value_integral(pt1);



#ginteraction(xdat=GEset[,subtypeN], grp=relapse_labels[as.character(subtypeN)], trt=trt_labels[as.character(subtypeN)])


promising = 0;

#if (pNFDR01>2) {promising=1}
promising= 1;
#if (testMode == "quick" ) {promising = 1;
#} else {if (pValue_Area[1] > 1.3 ) {promising = 1}}


if (promising > 0) {

subtypeID = paste("Pm", permute, chunk, NclustersA, clA, subt, sep="_");
GenesDefiningSubtypes[[subtypeID]] = c(subtype_genes_set);
SubtypePatients[[subtypeID]] = c(subtypeN);


if (abs(FDR.Res$b)>0.2){
FDR.Res<-ELF(xdat=GEset[,subtypeN], grp=relapse_labels[as.character(subtypeN)] , nperm=50, normalize=TRUE);
}
NFDR01 = sum(FDR.Res$ELF<0.1) 
NFDR02 = sum(FDR.Res$ELF<0.2) 
NFDR03 = sum(FDR.Res$ELF<0.3)
#print(NFDR02)


NFDR01_p<-NULL
NFDR02_p<-NULL
NFDR03_p<-NULL

if (FDRpermutation==TRUE){

NFDR01_p<-rep(0,nFDRperm)
NFDR02_p<-rep(0,nFDRperm)
NFDR03_p<-rep(0,nFDRperm)

for (j in 1:nFDRperm){ 

numcc<-1
while (numcc==1){
pseudo_relapse_labels<-sample(relapse_labels); 
numcc<-length(unique(pseudo_relapse_labels[as.character(subtypeN)]))
}

#pseudo_relapse_labels<-sample(relapse_labels[as.character(subtypeN)]); 
#names(pseudo_relapse_labels)<-subtypeN


names(pseudo_relapse_labels)<-Patients


#sum(pseudo_relapse_labels[as.character(subtypeN)]==1)
#sum(pseudo_relapse_labels[as.character(subtypeN)]==2)


FDR.Res_p<-ELF(xdat=GEset[,subtypeN], grp=pseudo_relapse_labels[as.character(subtypeN)] , nperm=50, normalize=FALSE);

NFDR01_p[j] = sum(FDR.Res_p$ELF<0.1) 
NFDR02_p[j] = sum(FDR.Res_p$ELF<0.2)
NFDR03_p[j] = sum(FDR.Res_p$ELF<0.3)

#print(NFDR02_p[j] )
}

}

SubtypeFDR01[[subtypeID]] = c(NFDR01_p);
SubtypeFDR02[[subtypeID]] = c(NFDR02_p);
SubtypeFDR03[[subtypeID]] = c(NFDR03_p);

if (showMovie > 0) {
heatmap(t(GEset[classifying_set, ]), labCol=ifelse(classifying_set %in% subtype_genes_set, "x", ""), xlab="GENES", ylab="PATIENTS", scale="none");
}


#####################################
##AUC_and_SDiff<-NULL;
#if (redefineSubtypeMembers > 0) {
#centroid<-apply(GEset[subtype_genes_set, subtypeN], 1, mean);
#AUC_and_SDiff = test_subt_massively(subtype_genes_set, centroid, top_best_probes, method);
#} else {
#subtype_individuals = subtypeN;
#if (testMode == "quick") {
#} else {
#AUC_and_SDiff = test_subtype_only(subtype_individuals, top_best_probes, method, subtypeID, holdOut, Niter, relapse_labels, relapse_times, GEset, trt, showMovie);
#}}
#######################################

subtype_individuals = subtypeN;
if (testMode == "quick") {
} else {
AUC_and_SDiff = test_subtype_only(subtype_individuals, top_best_probes, method, subtypeID, holdOut, Niter, relapse_labels, relapse_times, GEset, trt, showMovie);
}


promising2 = 0;
if (testMode == "quick") { promising2 = 1
} else { 
if (AUC_and_SDiff$AUC>0.6 | AUC_and_SDiff$SDiff>4) {promising2 = 1;}
}


if (promising*promising2==0) {FakeID<-c(FakeID,subtypeID)}

if (promising2 > 0) {
for (pr in subtype_genes_set) {
GenesDefiningSubtypes_allMembers[pr] = ifelse(is.na(GenesDefiningSubtypes_allMembers[pr]),1, GenesDefiningSubtypes_allMembers[pr] + 1);
}

if (testMode == "quick" ) {
SubtypeBestByPvalue[[subtypeID]]<-sort(pt1)[1:top_best_probes];
resultsAll <- rbind(resultsAll, c(permute, chunk, subtypeID, NclustersA, clA, length(subtype_genes_set), subt, length(subtypeN), pValue_Area, NFDR01,NFDR02,NFDR03,pNFDR01));
} else {
#features = c("Relapse", "Adj.chemo", "pT", "pN", "pM", "Type", "Prev.cancer_any");
#assocs<-NULL;
#for (fe in features) {
#assocs = c(assocs, association(fe, subtypeN)); #associations with clinical features
#}
#subtResults<- rbind(subtResults, c(subtypeID, NclustersA, clA, length(subtype_genes_set), subt, length(subtypeN), pValue_Area, AUC_and_SDiff$AUC, AUC_and_SDiff$SDiff, assocs));

trtResults[[subtypeID]] = AUC_and_SDiff$results3D;
BestSubtypeGenes[[subtypeID]] = AUC_and_SDiff$best_in_subtype;


subtypeSpecificProbes[[subtypeID]] = AUC_and_SDiff$SubtypeSpecificSignatures;

subtResults<- rbind(subtResults, c(subtypeID, NclustersA, clA, length(subtype_genes_set), subt, length(subtypeN), pValue_Area, NFDR01,NFDR02,NFDR03, AUC_and_SDiff$AUC, AUC_and_SDiff$SDiff, pNFDR01));
#print(c(subtypeID, NclustersA, clA, length(subtype_genes_set), subt, length(subtypeN), pValue_Area, NFDR02, AUC_and_SDiff));
}


} ### promising2
} ### promising
} ### the uniqueness of relapse labels
} ### subt (subtype 1 or 2)
} ### fingerprint
} ## condition: length (unique ( class of subtype ) ) >1
} ### clA
} ### NclustersA
} ### chunk
} ### permute

if (testMode == "quick" && !is.null(resultsAll)) {
colnames(resultsAll)[1:14] = c("Permute", "Chunk", "SubtypeID", "NclustersInitial", "Cluster", "Genes_in_cluster", "Subtype", "Patients_in_subtype", "low_pValue_Area", "AUCCDF","NFDR01","NFDR02","NFDR03","pNFDR01");
} 
 
if (testMode == "normal" && !is.null(subtResults)) {
colnames(subtResults)[1:14] = c("SubtypeID", "NclustersInitial", "Cluster", "Genes_in_cluster", "Subtype", "Patients_in_subtype", "low_pValue_Area", "AUCCDF","NFDR01","NFDR02","NFDR03", "AUC", "SurvDiffChiSq","pNFDR01");
NoFakeID<-setdiff(names(GenesDefiningSubtypes[-1]),FakeID); NoFakeIDN<-which(subtResults[,"SubtypeID"] %in% NoFakeID);
}

if (testMode == "normal" && is.null(subtResults)) {NoFakeID<-NULL; NoFakeIDN<-NULL;}

if (testMode == "quick") {
return(list(resultsAll=resultsAll,GenesDefiningSubtypes=GenesDefiningSubtypes[-1],SubtypePatients=SubtypePatients[-1],SubtypeFDR01=SubtypeFDR01[-1],SubtypeFDR02=SubtypeFDR02[-1],SubtypeFDR03=SubtypeFDR03[-1]));
} else {
return(list(subtResults=subtResults[NoFakeIDN,],GenesDefiningSubtypes=GenesDefiningSubtypes[NoFakeID],SubtypePatients=SubtypePatients[NoFakeID],SubtypeFDR01=SubtypeFDR01[NoFakeID],SubtypeFDR02=SubtypeFDR02[NoFakeID],SubtypeFDR03=SubtypeFDR03[NoFakeID],trtResults=trtResults[NoFakeID],BestSubtypeGenes=BestSubtypeGenes[NoFakeID],SubtypeSpecificSignatures=SubtypeSpecificSignatures[NoFakeID]));
}

}
