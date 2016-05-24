knnVCN <-
function
(TrnX,OrigTrnG,TstX=NULL,K=1,ShowObs=F, method = "euclidean",p = 2)
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

NTstX=nrow(TstX)
NTrnTotal=nrow(TrnTotal)

VoteResult=NULL
VoteResultList=NULL

for(i in 1:nrow(TstX))
{

RankBoardI<-NULL
RankBoardIJ<-NULL 

Total=rbind(TstX[i,],TrnX)
RankBoardI=as.matrix(dist(Total, method = method ,p=p)[1:nrow(TrnX)])
RankBoardIJ=cbind(TrnG,RankBoardI)

VoteAndWeight=RankBoardIJ[sort(RankBoardIJ[,2],index.return=T)$ix[1:K],1:2]
TempVote4TstXI=RankBoardIJ[sort(RankBoardIJ[,2],index.return=T)$ix[1:K],1]
ElmVote=union(TempVote4TstXI,TempVote4TstXI)

CountVote=as.matrix(sort(table(TempVote4TstXI),decreasing=T))
TempWinner=as.numeric(rownames(CountVote))

if(length(CountVote)==1|K==1)
{
Winner=TempWinner[1]
TstXIBelong=union(CodeMeaning$OrigTrnG[which(CodeMeaning$TrnG==Winner)],
CodeMeaning$OrigTrnG[which(CodeMeaning$TrnG==Winner)])
VoteResultNode=data.frame(TstXIBelong)
VoteResultList=rbind(VoteResultList,VoteResultNode)

}else   
{
NumOfTie=CountVote[1]
FinalList=NULL

j=1
TempWeight=sum(VoteAndWeight[which(VoteAndWeight[,1]==TempWinner[j]),2])
FinalList=data.frame(TempWinner[j],TempWeight)
while(CountVote[j]==CountVote[j+1]&j<length(CountVote))
{
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

if(IsTst==1)
{
CheckT=as.matrix(table(data.frame(VoteResultList,OrigTrnG)))
AccuStat=1-sum(CheckT-diag(diag(CheckT)))/length(TrnG)
print(CheckT)
cat("the classification accuracy of this algorithm on this training dataset is: ",
AccuStat*100,"%","\n\n\n")

}

if(IsTst==1&ShowObs==F){
result=data.frame(VoteResultList,OrigTrnG)
}else
{
if(IsTst==1&ShowObs==T){
result=data.frame(TstX,VoteResultList,OrigTrnG)
}else
{
if(ShowObs==F){
result=data.frame(VoteResultList)
}else{
result=data.frame(TstX,VoteResultList)
}
}
}
return(result)
}
