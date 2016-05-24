#' Weighted binary set for multiclass data
#' 
#' A function to get the weighted supersets and subsets for multiclass data
#' 
#' @param input a vector of multiclass data
#' @param data a matrix of multiclass data as training data
#' @return a list of the following components:
#' \item{superset}{the weighted supersets of the input data from the training data}
#' \item{subset}{the weighted subsets of the input data from the training data}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' data(tyner_interaction_multiclass)
#' sets<-getBinary1(tyner_interaction_multiclass[1,], tyner_interaction_multiclass)
getBinary1<-function(input, data){
  #example:
  #input<-c(0,1,2,3)
  #data<-sample(0:5,40, T)
  #data<-array(data,dim=c(10,4))
  input1<-matrix(rep(input, nrow(data)), ncol=ncol(data), byrow=TRUE)
  res<-input1-data
  drugs<-c(1:nrow(data))
  neg_val<-which(res<0, arr.ind=TRUE)
  neg_row<-unique(neg_val[,1])
  #zeros<-apply(res,1, function(x) if(all(x==0)) return(1) else return(0))
  #identical<-which(zeros==1)
  zeros<-rowSums(abs(res))
  identical<-which(zeros==0)
  # subset no negative values
  sub<-setdiff(drugs, neg_row)
  sub<-setdiff(sub, intersect(sub, identical))
  # weight for every subset to the input data
  if(length(sub)!=0){
    sub_w<-rowSums(matrix(res[sub,], ncol=ncol(data)))
  }else{sub_w<-0}
  
  
  pos_val<-which(res>0, arr.ind=TRUE)
  pos_row<-unique(pos_val[,1])
  # superset no positive values
  sup<-setdiff(drugs, pos_row)
  # weight for every superset to the input data
  if(length(sup)!=0){
    sup_w<-abs(rowSums(matrix(res[sup,], ncol=ncol(data))))
  }else{sup_W<-0}
  
  #res1<-data-input1
  #minones<-apply(res,1, function(x) if(length(which(x<0))==0) return(1) else return(0))
  # get the index in the data for subset
  #sub<-which(minones==1)
  #sub<-data[which(minones==1),]    
  #positiveval<-apply(res,1, function(x) if(length(which(x>0))==0) return(1) else return(0))
  # get the index in the data for subset
  #sup<-which(positiveval==1)
  #sup<-data[which(positiveval==1),]
  return(list(subset=sub, superset=sup, subw=sub_w, supw=sup_w))
}
