test_subtype_only <-
function(subtype_individuals, top_best_probes, method, subtypeID, holdOut, Niter, relapse_labels, relapse_times, GEset, trt, showMovie,SubtypeSpecificSignatures,SubtypeSpecificPredictors,NclustersA,clA) {

minTrainClassSize = 10
if ((length(subtype_individuals) - holdOut)< minTrainClassSize) {print(paste(subtypeID,":","Training data has too few samples."))}


#par(mfrow=c(2,1));

results2D<-NULL;results3D<-NULL; 

probeNames<-NULL;aucs<-NULL; predictor_genes<-NULL;

#if (showMovie > 0) {
#heatmap_plus(actual_GE[subtype_genes_set,], clindat[,c(25,33,34)], covar=3, na.rm=TRUE);
#}

for (i in 1:Niter) {

# create training and test sets completely separated from each other:
train<-NULL; test<-NULL;
train<-sample(subtype_individuals, (length(subtype_individuals) - holdOut));
test<-subtype_individuals[which(!(subtype_individuals %in% train))];

dfree<-length(train) - 2; # df for t-test
#table(as.character(train) %in% as.character(names(ctreeB)))
subtype_train<-train;
subtype_test<-test;
toTrain = sort(table(relapse_labels[subtype_train]), decreasing = T)[2];
toTest = sort(table(relapse_labels[subtype_train]), decreasing = T)[2];


if ((!is.na(toTrain) & toTrain > minTrainClassSize) & (!is.na(toTest) & toTest > -1)) { #the training and test sets must have enough cases of both categories

tOC<-array(rep(NA, times=length(rownames(GEset))));
tOC<-tstatistics2(GEset[,subtype_train], relapse_labels[subtype_train])$tstat;
names(tOC)<-rownames(GEset);
pt1<-pt(abs(tOC), dfree, lower.tail = F)
#pt2<-p.adjust(pt1, method="fdr");
best_in_subtype<-names(sort(pt1))[1:top_best_probes];
probeNames = c(probeNames, best_in_subtype);


border = 1.5;

#################################################################
##produce the predictor:
#if (method == "lasso") {
#l1<-lars(t(GEset[best_in_subtype, subtype_train]),relapse_labels[subtype_train],type = "lasso")
#cont_predictions<-predict(l1, t(GEset[best_in_subtype, subtype_test]))$fit
#cont_predictions <- cont_predictions[,ncol(cont_predictions)];
#dicr_predictions = ifelse((cont_predictions>border), 2, 1)
#}


#if (method == "dlda") {
#dicr_predictions<-stat.diag.da(
#t(GEset[best_in_subtype, subtype_train]),
#relapse_labels[subtype_train],
#t(GEset[best_in_subtype, subtype_test]),
#pool=1)$pred;
#cont_predictions = local_dlda(
#t(GEset[best_in_subtype, subtype_train]),
#relapse_labels[subtype_train],
#t(GEset[best_in_subtype, subtype_test]),
#1);
#}
#################################################################


if (method == "penalized") {
#best_in_subtype<-best10_stN_0; #subtype_test=clindat$Patid;
p1<-penalized(factor(relapse_labels[subtype_train]), t(GEset[best_in_subtype, subtype_train]),standardize=TRUE, lambda1=1);
#cont_predictions<-predict(p1, t(GEset[best_in_subtype, subtype_test]))
cont_predictions<-predict(p1, t(GEset[best_in_subtype, subtype_test])) #[,1]
dicr_predictions = ifelse((cont_predictions>0.5), 2, 1)
predictor_genes<-c(predictor_genes,  names(p1@penalized)[which(p1@penalized != 0)]);
SubtypeSpecificSignatures[[subtypeID]][[i]]<<-p1@penalized[which(p1@penalized != 0)];
}



results3D<-rbind(
results3D, cbind(
dicr_predictions,
relapse_labels[subtype_test] ,
relapse_times[subtype_test], trt[subtype_test]));
colnames(results3D)<-c("subtype","event","times","treatment")

results2D<-rbind(results2D, cbind(
cont_predictions,
relapse_labels[subtype_test]));
#if (showMovie > 0 & i == 10) {
#heatmap_plus(actual_GE[names(p1@penalized)[which(p1@penalized != 0)],], clindat[,c(25,33,34)], covar=3, na.rm=TRUE);
#}
} ### if condition
} ### iter: 1 ~ Niter


#heatmap(t(actual_GE[names(p1@penalized)[which(p1@penalized != 0)], train]), labRow=ifelse(names(ctreeB) %in% subtype_train, "x", ""), xlab="GENES", ylab="PATIENTS", scale="none");


if (!is.null(results2D)) {
if (nrow(results2D) >= 25) {
SubtypeSpecificPredictors[[subtypeID]]<<-table(predictor_genes);
pred1<-prediction(results2D[,1], results2D[,2]);
AUC<-as.numeric(performance(pred1, measure="auc")@y.values)
SDiff<-survdiff(Surv(results3D[,3],results3D[,2])~ factor(results3D[,1]))$chisq;
#SDiff_trt<-coxph(Surv(results3D[,3],results3D[,2])~ factor(results3D[,1])+factor(results3D[,4])+factor(results3D[,1]*results3D[,4]))


if (showMovie > 0) {
plot(performance(pred1, measure="tpr", x.measure="fpr"));
fit =survfit(Surv(results3D[,3],results3D[,2])~ factor(results3D[,1]));
plot(fit, mark.time=F, lty=1:3, xlab='Years since surgery', ylab='Relapse-free survival');
title(paste(NclustersA, "in.clust; cl.#", clA), line=2, cex.main=0.9);
title(paste(top_best_probes, " probes;", "hold-out=", holdOut, ""), line=1,  cex.main=0.75);
}

#mean(results2D[,1]); SDiff; AUC;
return(list(AUC=round(AUC, digits=4), SDiff=round(SDiff, digits=2),results3D=results3D,best_in_subtype=best_in_subtype,SubtypeSpecificSignatures=SubtypeSpecificSignatures));
}
if (nrow(results2D) < 25) {return(list(AUC=0,SDiff=0,results3D=0,best_in_subtype=0,SubtypeSpecificSignatures=0))}
}

if (is.null(results2D)) {return(list(AUC=0,SDiff=0,results3D=0,best_in_subtype=0,SubtypeSpecificSignatures=0))}

}
