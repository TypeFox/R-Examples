gbm2sas <-
function(
gbmobject,
data=NULL,
sasfile=NULL,
ntrees=NULL,
mysasdata="mysasdata",
treeval="treeval",
prefix="do_"
) {
if(is.null(ntrees)) ntrees<-gbmobject$n.trees
maxhmmt<-0
hasprefix<-prefix!="do_"
hasmysasdata<-mysasdata!="mysasdata"
hastreeval<-treeval!="treeval"
prepwords<-"data mysasdata; set mysasdata;"
if(hasmysasdata) prepwords<-gsub("mysasdata", mysasdata, prepwords)
write.table(prepwords, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE)
numtrees<-ntrees
for(treeloop in 1:numtrees) {
pgt<-pretty.gbm.tree(gbmobject,i.tree = treeloop)[1:7]
hmmt<-dim(pgt)[1]
maxhmmt<-max(maxhmmt, hmmt)
wordsa<-"do_x=0;"
for(loop in 0:(hmmt-1)) {
if(loop>0) {
wordsb<-gsub("x", loop, wordsa)
} else {
wordsb<-"do_0=1;"
}
if(hasprefix) wordsb<-gsub("do_", prefix, wordsb)
write.table(wordsb, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
} 
words0<-"if missing(V A R1) then do_V A R5=1; else do;" 
words1<-"if V A R1 lt V A R2 then do_V A R3=1; else do_V A R4=1; end;"
words2<-"if V A R1 in (V A R2) then do_V A R3=1; else do_V A R4=1; end;"
words2b<-"do_V A R4=1; end;"
words3<-"end;"
if(hasprefix) {
words0<-gsub("do_", prefix, words0)
words1<-gsub("do_", prefix, words1)
words2<-gsub("do_", prefix, words2)
words2b<-gsub("do_", prefix, words2b)
words3<-gsub("do_", prefix, words3)
}
thevarnames<-gbmobject$var.names
types<-lapply (lapply(data[,gbmobject$var.names],class), function(i) ifelse (strsplit(i[1]," ")[1]=="ordered","ordered",i))
levels<-lapply(data[,gbmobject$var.names],levels)
for(loop in 1:hmmt) {
prepwords<-paste("if do_", (loop-1), ">0 then do;", sep="")
if(hasprefix) prepwords<-gsub("do_", prefix, prepwords)
write.table(prepwords, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
splitvar<-1+as.numeric(as.vector(pgt[loop,]$SplitVar))
splitcodepred<-as.numeric(as.vector(pgt[loop,]$SplitCodePred))
leftnode<-as.numeric(as.vector(pgt[loop,]$LeftNode))
rightnode<-as.numeric(as.vector(pgt[loop,]$RightNode))
missingnode<-as.numeric(as.vector(pgt[loop,]$MissingNode))
if(splitvar>0) {
words0a<-gsub("V A R1", thevarnames[splitvar], words0)
words1a<-gsub("V A R1", thevarnames[splitvar], words1)
words2a<-gsub("V A R1", thevarnames[splitvar], words2)
words0a<-gsub("V A R5", missingnode, words0a)
words1a<-gsub("V A R3", leftnode, words1a)
words2a<-gsub("V A R3", leftnode, words2a)
words1a<-gsub("V A R4", rightnode, words1a)
words2a<-gsub("V A R4", rightnode, words2a)
words2ab<-gsub("V A R4", rightnode, words2b)
thistype<-types[[splitvar]]
leftstring<-" "
rightstring<-" "
write.table(words0a, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
if(thistype=="numeric") {
words1a<-gsub("V A R2", splitcodepred, words1a)
write.table(words1a, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
} else {
if(thistype=="ordered") {
splitcodepred<-ceiling(splitcodepred)
if(splitcodepred>=1) {
theleft<-c(levels[[splitvar]][1:splitcodepred], NA)
} else {
theleft<-rep(NA, 2)                   
}
} else {
describer<-unlist(gbmobject$c.splits[1+splitcodepred])
theleft<-c(levels[[splitvar]][describer==-1], NA)
}
logic<-!is.na(theleft)
if(sum(as.numeric(logic))>0) {
theleft<-theleft[logic]
hmmt2<-length(theleft) 
leftstring<-NULL
for(loop2 in 1:hmmt2) {
leftstring<-paste(leftstring, "'", theleft[loop2], "'", sep="")
if(loop2<hmmt2) leftstring<-paste(leftstring, ", ", sep="")
}
} else {
leftstring<-"blah"
}
if(leftstring!="blah") {
words2a<-gsub("V A R2", leftstring, words2a)
write.table(words2a, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
} else {
write.table(words2ab, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
}
}
} else {
prepwords<-paste("treeval", treeloop, "=", splitcodepred, ";", sep="")
if(hastreeval) prepwords<-gsub("treeval", treeval, prepwords)
write.table(prepwords, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
}
write.table(words3, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
}
}
wordsa<-"drop do_x;"
for(loop in 0:(maxhmmt-1)) {
if(loop>0) {
wordsb<-gsub("x", loop, wordsa)
} else {
wordsb<-"drop do_0;"
}
if(hasprefix) wordsb<-gsub("do_", prefix, wordsb)
write.table(wordsb, sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
}
write.table("run;", sasfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
}
