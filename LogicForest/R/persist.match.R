persist.match <-
function(fit, pred.nms, preds, match.list)
{
 if (missing(pred.nms)) {
   prednames<-colnames(fit$Xs)
   predcomps<-paste("!", prednames, sep="")
   pred.names<-c(prednames,predcomps)
   }
 else {
   predcomps<-paste("!", pred.nms, sep="")
   pred.names<-c(pred.nms, predcomps)
   }
 if (class(fit)=="logforest") {allPI.list<-names(unlist(fit$PI.importance))}
 if (class(fit)=="LBoost") {allPI.list<-names(fit$PI.frequency)}
 PI.freq<-fit$PI.frequency
 for (i in 1:length(allPI.list))
   {
   if (allPI.list[i]=="") {allPI.list[i]<-"Null"}
   else {allPI.list[i]<-allPI.list[i]}
   }
 if (!missing(match.list))
   {
   PI.list<-match.list
   for (j in 1:length(PI.list))
     {
     if (PI.list[j]=="") {PI.list[j]<-"Null"}
     else {PI.list[j]<-PI.list[j]}
     }
   PI.mtch<-list.PIs(match.list=PI.list, allPIs=allPI.list, preds=preds)
   all.PIs<-PI.mtch$PIsetlist
   allPIs<-PI.mtch$allPI.list
   minPI<-min(PI.mtch$lengths)
   maxPI<-max(PI.mtch$lengths)
   lgth<-PI.mtch$lengths
   all.lgth<-PI.mtch$all.lengths
   }
 else{ 
   PI.list<-unique(c(pred.names, allPI.list))
   PI.info<-list.PIs(match.list=PI.list, allPIs=allPI.list, preds=preds)
   all.PIs<-PI.info$PIsetlist
   allPIs<-PI.info$allPI.list
   minPI<-min(PI.info$lengths)
   maxPI<-max(PI.info$lengths)
   lgth<-PI.info$lengths
   all.lgth<-PI.info$all.lengths
   }
 PIsets<-matrix(0, nrow=length(lgth), ncol=length(all.PIs))
 whlPI.set<-matrix(0, nrow=length(all.lgth), ncol=length(allPIs))
 colnames(PIsets)<-PI.list
 rownames(PIsets)<-paste("size", lgth, sep=" ")
 colnames(whlPI.set)<-allPI.list
 rownames(whlPI.set)<-paste("size", all.lgth, sep=" ")
 for (i in 1:length(all.PIs))
   {
   PIset<-all.PIs[[i]]
   for (j in 1:length(lgth))
     {
     sz<-lgth[j]
     PIsets[j,i]<-ifelse(length(PIset)==sz, 1, 0)
     }
   }
 for (i in 1:length(allPIs))
   {
   PIs<-allPIs[[i]]
   for (j in 1:length(all.lgth))
     {
     sz2<-all.lgth[j]
     whlPI.set[j,i]<-ifelse(length(PIs)==sz2, 1, 0)
     }
   }
 numsizes<-rowSums(whlPI.set)
 names(lgth)<-paste("size", lgth, sep=" ")
 names(all.lgth)<-paste("size", all.lgth, sep=" ")
 names(numsizes)<-paste("size", all.lgth, sep=" ")
 PI.sets<-vector("list", length(lgth))
 allPI.sets<-vector("list", length(all.lgth))
 for (i in 1:length(lgth))
   {
   PIids<-which(PIsets[i,]==1)
   PInms<-vector("list", length(PIids))
   if(length(PIids)==0) {PInms<-"NA"}
   else{
     for(j in 1:length(PIids))
       {
       PInms[[j]]<-all.PIs[[PIids[j]]]
       }
     }
   PI.sets[[i]]<-PInms
   if (length(PIids)==1) {names(PI.sets[[i]])<-colnames(PIsets)[PIids]}
   else {names(PI.sets[[i]])<-colnames(PIsets[,PIids])}
   }
 for (i in 1:length(all.lgth))
   {
   allPIds<-which(whlPI.set[i,]==1)
   allPInms<-vector("list", length(allPIds))
   if (length(allPIds)==0) {allPInms<-"NA"}
   else {
     for(g in 1:length(allPIds))
       {
       allPInms[[g]]<-allPIs[[allPIds[g]]]
       } 
     }
   allPI.sets[[i]]<-allPInms
   if (length(allPIds)==1) {names(allPI.sets[[i]])<-colnames(whlPI.set)[allPIds]}
   else {names(allPI.sets[[i]])<-colnames(whlPI.set[,allPIds])}
   }
 matches<-vector("list", (length(lgth)-1))
 PIlgth<-c()
 num.matches<-vector("list", (length(lgth)-1))
 for(i in 1:length(lgth))
   {
   sz<-lgth[i]
   locat<-which(all.lgth==sz)
   all<-length(all.lgth)
   if (length(locat)==0) {comp.szs<-all.lgth}
   else {comp.szs<-all.lgth[locat:all]}
   sz.notmtch<-length(setdiff(all.lgth, comp.szs))
   extras<-length(comp.szs)
   PIs<-PI.sets[[i]]
   PIlgth<-append(PIlgth, length(PIs))
   match.mat<-matrix(0, nrow=length(PIs), ncol=extras)
   rownames(match.mat)<-names(PI.sets[[i]])
   colnames(match.mat)<-paste("size",comp.szs, sep=" ")
   PI.extras<-vector("list", extras)
   for(k in 1:extras)
     {
     PI.extras[[k]]<-allPI.sets[[k+sz.notmtch]]
     }
   numPIs<-length(PIs)
   PI.match<-vector("list", numPIs)
   names(PI.match)<-names(PIs)
   for (j in 1:numPIs)
     {
     match.names<-vector("list", extras)
     for (m in 1:extras)
       {
       PI.comp<-PI.extras[[m]]
       tr.match<-c()
       for (d in 1:length(PI.comp))
         {
         mch<-length(which(PI.comp[[d]]%in%PIs[[j]]))
         mtch.nm<-ifelse(mch==sz, names(PI.extras[[m]][d]), "Null")
         true<-ifelse(mch==sz, 1, 0)
         tr.match<-append(tr.match, true)
         if (mtch.nm!="Null") {match.names[[m]]<-append(match.names[[m]], mtch.nm)}
         }
       match.mat[j,m]<-sum(tr.match)
       }
     PI.match[[j]]<-match.names
     } 
   matches[[i]]<-PI.match
   num.matches[[i]]<-match.mat
   }
 names(num.matches)<-paste("size", lgth, sep=" ")
 ans<-list(numsizes=numsizes, uniq.sz=lgth, all.sz=all.lgth, matches=matches, num.matches=num.matches)
}
