orddom_p <- function (x,y,alpha=.05,paired=FALSE,sections="1234a4b5a5b",header="Y",sorted="XY",outfile="orddom_csv_tab.txt",appendfile=FALSE,show=1,description="") { 
#prints tab-delimited output of character matrix from dominance matrix with signs
#General settings
library(psych)
if(sorted==""){sorted="@@"}
if(paired==TRUE){pd<-TRUE}else{pd<-FALSE}
err<-getOption("show.error.messages")
options("show.error.messages"=FALSE)
sink(file=outfile, append=appendfile, type="output", split=FALSE)
if(paired==FALSE){
 x_head<-"INDEPENDENT GROUPS"
 x_name<-paste("X/Group1/Baseline (",colnames(x),")",sep="")
 y_name<-paste("Y/Group2/Treatment (",colnames(y),")",sep="")
} else {#paired=TRUE 
 x_head<-"DEPENDENT GROUPS / PAIRED/REPEATED DATA"
 x_name<-paste("X/Pre (",colnames(x),")",sep="")
 y_name<-paste("Y/Post (",colnames(y),")",sep="")
} 
 x_name<-gsub("(^ +)|( +$)", "", x_name)
 y_name<-gsub("(^ +)|( +$)", "", y_name)

if (header=="Y"){
if (description!="") { cat("***** ",description," *****\n",sep="",append=TRUE) }
cat("===== ORDINAL DOMINANCE ANALYSIS FOR ",x_head," =====",sep="",append=TRUE)
}
if(any(grep("[1]",sections))){
#=========Output Section 1 Raw Data=====================================================================================================
 if (header=="Y"){
  cat("\n\n1. Raw Data",sep="",append=TRUE)
  if (any(grep("[XY]",sorted,TRUE))==TRUE) {cat (" (sorted:",sorted,")")}}
 cat ("\n(i)\t",x_name,sep="",append=TRUE)
 if(paired==FALSE) {
 #indep---------------------------------------------------------------------------------
  cat ("\n",sep="",append=TRUE)
  #Ausgabe Rohdaten Gruppe 1
  if (any(grep("[X]",sorted,TRUE))) {#sorted
  write.table(cbind(matrix(seq(1:length(x)),ncol=1),x)[order(x),],append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  } else {#unsorted
  write.table(data.frame(x),append=TRUE,quote=FALSE,sep="\t",col.names=FALSE)}
  #Ausgabe Rohdaten Gruppe 2  
  cat ("\n(j)\t",y_name,"\n",sep="",append=TRUE)
  if (any(grep("[Y]",sorted,TRUE))) {#sorted
  write.table(cbind(matrix(seq(1:length(y)),ncol=1),y)[order(y),],append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  } else {#unsorted
  write.table(data.frame(y),append=TRUE,quote=FALSE,sep="\t",col.names=FALSE)}
  } 
else 
{
 #paired---------------------------------------------------------------------------------
  if (any(grep("[XY]",sorted,TRUE))) {#sorted
  if (any(grep("[X]",sorted,TRUE))) {#sorted X->Y
  table<-cbind(matrix(seq(1:length(x)),ncol=1),x,y)[order(x,y),]} else {
  table<-cbind(matrix(seq(1:length(x)),ncol=1),x,y)}
  if (any(grep("[Y]",sorted,TRUE))) {#sorted Y->X
  table<-cbind(table,cbind(matrix(seq(1:length(x)),ncol=1),x,y)[order(y,x),])}else{
  table<-cbind(matrix(seq(1:length(x)),ncol=1),x,y,table)}
  cat("\t",y_name,"\t(j)\t",x_name,"\t",y_name,"\n",sep="",append=TRUE)
  write.table(table,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  } else {#unsorted
  cat("\t",y_name,"\ni,j=1..n\t",sep="",append=TRUE)
  write.table(cbind(x,y),quote=FALSE,sep="\t")}#Ausgabe Rohdaten Pre/post
}}
if(any(grep("[2]",sections))){
#=========Output Section 2 Descriptives=====================================================================================================
 if (header=="Y"){ cat("\n\n2. Descriptives\n",sep="",append=TRUE)}
  cat("\t",x_name,"\t",y_name,"\n",sep="",append=TRUE)
  write.table(t(rbind(describe(x),describe(y))),append=TRUE,quote=FALSE,sep="\t",col.names=FALSE)} #Ausgabe Descriptives print
if(any(grep("[3]",sections))){
#=========Output Section 3 Metric Differences t ============================================================================================
 if (header=="Y"){ cat("\n\n3. Metric t-test for ",tolower(x_head),"\n",sep="",append=TRUE) }
 write.table(metric_t(x,y,alpha,paired),append=TRUE,quote=FALSE,sep="\t",col.names=FALSE) #Ausgabe t-test results
}
if(any(grep("[4][a]",sections))){
#=========Output Section 4a Difference Matrix =============================================================================================
 if (header=="Y"){ cat("\n\n4a. Difference Matrix (Rows:",x_name,"  /  Cols:",y_name,")\n",sep="",append=TRUE) }
 if(length(y)>50) {cat ("Number of cases in ",y_name," exceeds 50. Printout skipped.\n",sep="",append=TRUE)} else {
 cat("i\tX\\Y\t",sep="",append=TRUE)
 table<-cbind(t(colnames(dm(x,y))),"di.")
 write.table(table,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
 if (any(grep("[X]",sorted,TRUE))) {#sorted
 table<-cbind(matrix(seq(1:length(x)),ncol=1),x,-dm(x,y,diff=TRUE))[order(x),]} else {
 table<-cbind(matrix(seq(1:length(x)),ncol=1),x,-dm(x,y,diff=TRUE))} #unsorted
 table<-cbind(table,t(t(rowMeans(table[,3:(length(y)+2)])))) #di.
 table<-rbind(table,c("","d.j",t(colMeans(table[,3:(length(y)+3)])))) #d.j, dij
 write.table(table,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)}}
if(any(grep("[4][b]",sections))){
#=========Output Section 4b Difference Matrix =============================================================================================
 if (header=="Y"){ cat("\n\n4b. Difference Matrix (Rows:",y_name,"  /  Cols:",x_name,")\n",sep="",append=TRUE) }
 if(length(x)>50) {cat ("Number of cases in ",x_name," exceeds 50. Printout skipped.\n",sep="",append=TRUE)} else {
 cat("j\tY\\X\t",sep="",append=TRUE)
 table<-cbind(t(colnames(dm(y,x))),"d.j")
 write.table(table,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
 if (any(grep("[Y]",sorted,TRUE))) {#sorted
 table<-cbind(matrix(seq(1:length(y)),ncol=1),y,t(-dm(x,y,diff=TRUE)))[order(y),]} else {
 table<-cbind(matrix(seq(1:length(y)),ncol=1),y,t(-dm(x,y,diff=TRUE)))} #unsorted
 table<-cbind(table,t(t(rowMeans(table[,3:(length(x)+2)])))) #d.j
 table<-rbind(table,c("","di.",t(colMeans(table[,3:(length(x)+3)])))) #di., dij
 write.table(table,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)}}
if(any(grep("[5][a]",sections))){
#=========Output Section 5a Dominance Matrix }=============================================================================================
 if (header=="Y"){ cat("\n\n5a. Dominance Matrix (Rows:",x_name,"  /  Cols:",y_name,")\n",sep="",append=TRUE)
 cat("\t[+",sep="",append=TRUE)
 if (paired==TRUE) {cat("/>",sep="",append=TRUE)}
 cat(" indicate higher values for ",y_name,"]\n",sep="",append=TRUE) }
 if(length(y)>50) {cat ("Number of cases in ",y_name," exceeds 50. Printout skipped.\n",sep="",append=TRUE)} else {
 cat("i\tX\\Y\t",sep="",append=TRUE)
 table<-cbind(t(colnames(dm(x,y))),"di.")
 write.table(table,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
 if(paired==FALSE){
  dcol<-t(t(rowMeans(-dm(x,y))))
  drow<-t(t(colMeans(-dm(x,y))))
  last<-mean(-dm(x,y))} else {
  dcol<-t(t((rowSums(-dm(x,y))-diag(-dm(x,y)))/(length(y)-1)))
  drow<-t(t((colSums(-dm(x,y))-diag(-dm(x,y)))/(length(x)-1)))
  last<-mean(dcol)}
 if (any(grep("[X]",sorted,TRUE))) {#sorted
 table<-cbind(matrix(seq(1:length(x)),ncol=1),x,dms(-dm(x,y),pd),dcol)[order(x),]} else {
 table<-cbind(matrix(seq(1:length(x)),ncol=1),x,dms(-dm(x,y),pd),dcol)} #unsorted
 table<-rbind(table,c("","d.j",drow,last))
 write.table(table,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)} }
if(any(grep("[5][b]",sections))){
#=========Output Section 5b Dominance Matrix }=============================================================================================
 if (header=="Y"){ cat("\n\n5b. Dominance Matrix (Rows:",y_name,"  /  Cols:",x_name,")\n",sep="",append=TRUE) 
 cat("\t[+",sep="",append=TRUE)
 if (paired==TRUE) {cat("/>",sep="",append=TRUE)}
 cat(" indicate higher values for ",y_name,"]\n",sep="",append=TRUE) }
 if(length(x)>50) {cat ("Number of cases in ",x_name," exceeds 50. Printout skipped.\n",sep="",append=TRUE)} else {
 cat("j\tY\\X\t",sep="",append=TRUE)
 table<-cbind(t(colnames(dm(y,x))),"d.j")
 write.table(table,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
 if(paired==FALSE){
  dcol<-t(t(colMeans(-dm(x,y))))
  drow<-t(t(rowMeans(-dm(x,y))))
  last<-mean(-dm(x,y))} else {
  dcol<-t(t((colSums(-dm(x,y))-diag(-dm(x,y)))/(length(y)-1)))
  drow<-t(t((rowSums(-dm(x,y))-diag(-dm(x,y)))/(length(x)-1)))
  last<-mean(dcol)}
 if (any(grep("[Y]",sorted,TRUE))) {#sorted
 table<-cbind(matrix(seq(1:length(x)),ncol=1),x,t(dms(-dm(x,y),pd)),dcol)[order(x),]} else {
 table<-cbind(matrix(seq(1:length(x)),ncol=1),x,t(dms(-dm(x,y),pd)),dcol)} #unsorted
 table<-rbind(table,c("","d.j",drow,last))
 write.table(table,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)}}
sink()
closeAllConnections()
if (show==1) {file.show(outfile)}
options("show.error.messages"=err) 
}
