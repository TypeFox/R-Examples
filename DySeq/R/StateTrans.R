#'StateTrans
#'
#'Produces a state transition table for dyadic binary sequences.
#'
#'That is, the behavior of interest in interval t, 
#'is mapped against the combination of the observed 
#'behaviors in the preceding interval (t - 1). 
#'Hence, the total absolute frequency equals the number
#'of time intervals minus 1. And the number of obtained 
#'tables is equal the number of sequence-pairs.   
#'
#'
#'printing the output will display mean frequencies. For inspecting individual cases use [[original-rownumber]].
#'
#'For an extensive overview see Kenny, Kashy and Cook (2006).
#'The original idea stems from (to our knowledge) Bakeman and Gottman (1997).
#'
#'
#'
#'@param x Dataframe or matix containing combined sequences, see help(StateExpand)
#'@param first logical value indicating if the first sequence should used as dependend variable (TRUE) or the second (FALSE)
#'@param dep.lab two-element string vector with labels for dependend variable (first entry corresponds to the value zero, the second to one)
#'@param indep.lab four-element string vector with labels for the combined variable (order corresponds to the order of the StateExpand function)
#'
#'@references 
#'\itemize{
#'  \item Bakeman, R., & Gottman, J. M. (1997) <DOI: 10.1017/cbo9780511527685 >
#'  \item Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006) <DOI: 10.1177/1098214007300894>
#'}
#' 
#'
#'@examples
#'
#'# Example 1: Sequences from couples cope
#'
#'data(CouplesCope)
#'my.s<-StateExpand(CouplesCope, 2:49, 50:97)
#'
#'# First sequence is dependend variable 
#'# - what behavior preceeds stress signals?
#'StateTrans(my.s) 
#'
#'# Second sequence is dependend variable 
#'# - what behavior preceeds dyadic coping signals?
#'StateTrans(my.s, FALSE) 
#'
#'# investigating a single case
#'StateTrans(my.s, FALSE)[[41]] 
#'
#'@export


StateTrans<-function(x,first=TRUE, dep.lab=c("1","0"), indep.lab=c("1-1","1-0","0-1","0-0")) {

  output<-list()
  comb<-x

  for(case in 1:length(x[,1])){

    y<-matrix(rep(0,8),4,2)
    colnames(y)<-dep.lab
    rownames(y)<-indep.lab

    if(first){
      for(i in 2:length(comb[case,])){
        if(comb[case,(i-1)]==0 & (comb[case,i]==1 |comb[case,i]==3)) y[4,1]<-y[4,1]+1
        if(comb[case,(i-1)]==1 & (comb[case,i]==1 |comb[case,i]==3)) y[2,1]<-y[2,1]+1
        if(comb[case,(i-1)]==2 & (comb[case,i]==1 |comb[case,i]==3)) y[3,1]<-y[3,1]+1
        if(comb[case,(i-1)]==3 & (comb[case,i]==1 |comb[case,i]==3)) y[1,1]<-y[1,1]+1

        if(comb[case,(i-1)]==0 & (comb[case,i]==0 |comb[case,i]==2)) y[4,2]<-y[4,2]+1
        if(comb[case,(i-1)]==1 & (comb[case,i]==0 |comb[case,i]==2)) y[2,2]<-y[2,2]+1
        if(comb[case,(i-1)]==2 & (comb[case,i]==0 |comb[case,i]==2)) y[3,2]<-y[3,2]+1
        if(comb[case,(i-1)]==3 & (comb[case,i]==0 |comb[case,i]==2)) y[1,2]<-y[1,2]+1
      }
    }


    if(!first){
      for(i in 2:length(comb[case,])){
        if(comb[case,(i-1)]==0 & (comb[case,i]==2 |comb[case,i]==3)) y[4,1]<-y[4,1]+1
        if(comb[case,(i-1)]==1 & (comb[case,i]==2 |comb[case,i]==3)) y[2,1]<-y[2,1]+1
        if(comb[case,(i-1)]==2 & (comb[case,i]==2 |comb[case,i]==3)) y[3,1]<-y[3,1]+1
        if(comb[case,(i-1)]==3 & (comb[case,i]==2 |comb[case,i]==3)) y[1,1]<-y[1,1]+1

        if(comb[case,(i-1)]==0 & (comb[case,i]==0 |comb[case,i]==1)) y[4,2]<-y[4,2]+1
        if(comb[case,(i-1)]==1 & (comb[case,i]==0 |comb[case,i]==1)) y[2,2]<-y[2,2]+1
        if(comb[case,(i-1)]==2 & (comb[case,i]==0 |comb[case,i]==1)) y[3,2]<-y[3,2]+1
        if(comb[case,(i-1)]==3 & (comb[case,i]==0 |comb[case,i]==1)) y[1,2]<-y[1,2]+1
      }
    }

    output[[case]]<-y
    class(output)[2]<-"state.trans"

  }
  return(output)
}



