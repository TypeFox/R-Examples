Margin.vote <-
function(vote,type="unsupervised",observed)
{
         myDim<-dim(vote)
         n_classes<-myDim[1]
         n_instances<-myDim[2]
         myMargin<-rep(-99,n_instances)
    if(type=="unsupervised")
     {

         
         for( i in 1 : n_instances)
         {
              myMaxVote<-max(vote[,i])
              whichMax<-which.max(vote[,i])
              mySecVote<-max(vote[-whichMax,i])
              myMargin[i]<-myMaxVote-mySecVote
              
             
         }
         
     }
    else
      {
           for(i in 1 : n_instances)
            {
                  myLabelVote<-vote[observed[i],i]
                  myMaxvote<-max(vote[-observed[i],i])
                  myMargin[i]<-myLabelVote-myMaxvote
            }
      }

      return (myMargin)
}
