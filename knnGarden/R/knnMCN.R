knnMCN <-
function
(TrnX,OrigTrnG,TstX=NULL,K=1,ShowObs=F)
{

OrigTrnG=as.factor(OrigTrnG)
TrnG=as.numeric(OrigTrnG)
CodeMeaning=data.frame(TrnG,OrigTrnG)

TK=sort(as.matrix(table(TrnG)),decreasing=F)
if(K>TK[1])
{
stop(c("
NOTES: 
sorry, the value of K ","(K=",K,") ",
"you have selected is bigger than the capacity of one class in your training data set",
"(","the capacity is ",TK[1],")",",","please choose a less value for K"))
}

if(is.null(TstX)==T) 
{
IsTst=1
TstX<-as.matrix(TrnX)
}else
{
IsTst=0
}

if(is.matrix(TstX)==F) 
{
TstX<-as.matrix(TstX)
}

TrnX<-as.matrix(TrnX)
ElmTrnG=union(TrnG,TrnG)
LevTrnG=length(ElmTrnG)
TrnTotal=cbind(TrnG,TrnX)

if(abs(det(cov(TrnX[which(TrnTotal[,1]==ElmTrnG[1]),])))<1e-07)
{
stop("
Warnings:
sample variance-covariance matrix is singular,
and larger class sample capacity is required ,
or you can try other methods in knnWXM of this package")
}else
{
MaDisList=list(solve(cov(TrnX[which(TrnTotal[,1]==ElmTrnG[1]),]),LINPACK=T))
}

if(LevTrnG>1)
{
for(i in (1+1):LevTrnG)
{
if(abs(det(cov(TrnX[which(TrnTotal[,1]==ElmTrnG[i]),])))<1e-07)
{
stop("
Warnings:
sample variance-covariance matrix is singular,
and larger class sample capacity is required ,
or you can try other methods in knnWXM of this package")
}else
{
MaDisNode=list(solve(cov(TrnX[which(TrnTotal[,1]==ElmTrnG[i]),]),LINPACK=T))
MaDisList=c(MaDisList,MaDisNode)
}
}
}


NTstX=nrow(TstX)
NTrnTotal=nrow(TrnTotal)

VoteResult=NULL
VoteResultList=NULL
for(i in 1:nrow(TstX)){
RankBoardI<-NULL
RankBoardIJ<-NULL
for(j in 1:LevTrnG){
TempTrnXI=TrnX[which(TrnTotal[,1]==ElmTrnG[j]),]
TempCovJ=as.matrix(MaDisList[[j]])
TempTstXI=NULL
for(k in 1:nrow(TempTrnXI)){
TempTstXI=rbind(TempTstXI,TstX[i,])
}
TempMadisBoardI<-sqrt(diag((TempTstXI-TempTrnXI)%*%TempCovJ%*%t(TempTstXI-TempTrnXI)))
MadisBoardI<-as.matrix(TempMadisBoardI)
GBoardI<-as.matrix(rep(ElmTrnG[j],nrow(TempTrnXI)))
RankBoardI<-cbind(GBoardI,MadisBoardI)
RankBoardIJ<-rbind(RankBoardIJ,RankBoardI)
}

VoteAndWeight=RankBoardIJ[sort(RankBoardIJ[,2],index.return=T)$ix[1:k],1:2]
TempVote4TstXI=RankBoardIJ[sort(RankBoardIJ[,2],index.return=T)$ix[1:k],1]
ElmVote=union(TempVote4TstXI,TempVote4TstXI)

CountVote=as.matrix(sort(table(TempVote4TstXI),decreasing=T))
TempWinner=as.numeric(rownames(CountVote))

if(length(CountVote)==1|K==1){
Winner=TempWinner[1]
TstXIBelong=union(CodeMeaning$OrigTrnG[which(CodeMeaning$TrnG==Winner)],
CodeMeaning$OrigTrnG[which(CodeMeaning$TrnG==Winner)])
VoteResultNode=data.frame(TstXIBelong)
VoteResultList=rbind(VoteResultList,VoteResultNode)
}else{
NumOfTie=CountVote[1]
FinalList=NULL
j=1
TempWeight=sum(VoteAndWeight[which(VoteAndWeight[,1]==TempWinner[j]),2])
FinalList=data.frame(TempWinner[j],TempWeight)
while(CountVote[j]==CountVote[j+1]&j<length(CountVote)){
TempWeight=sum(VoteAndWeight[which(VoteAndWeight[,1]==TempWinner[j+1]),2])
FinalListNode=c(TempWinner[j+1],TempWeight)
FinalList=rbind(FinalList,FinalListNode)
j=j+1  
}

FinalList=FinalList[sort(FinalList$TempWeight,index.return=T)$ix[1],]
TstXIBelong=union(CodeMeaning$OrigTrnG[which(CodeMeaning$TrnG==FinalList[1,1])],
CodeMeaning$OrigTrnG[which(CodeMeaning$TrnG==FinalList[1,1])])
VoteResultNode=data.frame(TstXIBelong)
VoteResultList=rbind(VoteResultList,VoteResultNode)
}
}

if(IsTst==1){
CheckT=as.matrix(table(data.frame(VoteResultList,OrigTrnG)))
AccuStat=1-sum(CheckT-diag(diag(CheckT)))/length(TrnG)
cat("test results","\n")
print(CheckT)
cat("the classification accuracy of this algorithm on this training dataset is: ",
AccuStat*100,"%","\n\n\n")
}

if(IsTst==1&ShowObs==F){
result=data.frame(VoteResultList,OrigTrnG)
}else{
if(IsTst==1&ShowObs==T){
result=data.frame(TstX,VoteResultList,OrigTrnG)
}else{
if(ShowObs==F){
result=data.frame(VoteResultList)
}else{
result=data.frame(TstX,VoteResultList)
}
}
}
return(result)
}
