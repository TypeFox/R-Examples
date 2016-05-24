#' fit crr model
#' 
#' Wrapper function to fit fine and gray competing risk model using function crr
#' from package cmprsk
#' 
#' @param f formula for the model. Currently the formula only works by using the name
#' of the column in a dataframe. It does not work by using $ or [] notation.
#' @param data dataframe containing data
#' @keywords model
#' @export
crrRx<-function(f,data){
  k<-as.character(f)[3]
  covs<-removedollar(k)
  ff<-modelmatrix(f,data)
  m1<-crr(ff[[1]][,1],ff[[1]][,2],ff[[2]])
  m1$call<-paste("~",covs)
  return(m1)
}

#' fit box cox transformed linear model
#' 
#' Wrapper function to fit fine and gray competing risk model using function crr
#' from package cmprsk
#' 
#' @param f formula for the model. Currently the formula only works by using the name
#' of the column in a dataframe. It does not work by using $ or [] notation.
#' @param data dataframe containing data
#' @param lambda boolean indicating if you want to output the lamda used in the boxcox transformation. If so the function will return a list of length 2 with the model as the first element and a vector of length 2 as the second.
#' @keywords model
#' @export
boxcoxfitRx<-function(f,data,lambda=F){
  x<-as.character(f)[3]
  y<-as.character(f)[2]
  time<- gsub("\\s","",unlist(strsplit(y,"+",fixed=T))[1])
  covs<-removedollar(x)
  tempindexboxcoxfitRx<-seq_len(nrow(data))
  df1<-data.frame(tempindexboxcoxfitRx,data)  
  f2<-as.formula(paste("tempindexboxcoxfitRx+",y,"~",x))    
  temp<-modelmatrix(f2,df1)
  ff<-list(temp[[1]][,-1,drop=F],temp[[2,drop=F]])
  temp<-temp[[1]][,1,drop=F]
  lambda1<-unlist(unlist(boxcoxfit(ff[[1]],ff[[2]],lambda2=T))[1:2])                     
  ff[[1]]<-((ff[[1]]+lambda1[2])^lambda1[1]-1)/lambda1[1]
  df<-merge(df1,temp,by="tempindexboxcoxfitRx")[,-1,drop=F]
  df[,time]<-ff[[1]]
  out<-lm(f,data=df)
  out$call<-paste("~",covs)
  if(lambda)  return(list(out,lambda1))
  return(out)
}


# bwselect<-function(data,response,covs,strata=1,type,force=NULL,p=0.05,test=F,boxcox=F){
#   while(T){
#     
#     if(type=="logistic"){
#       pvalues<-sapply(covs,function(cov){
#         m0<-glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                 data=subset(data,!is.na(data[,cov])),family="binomial")
#         m1<-update(m0,as.formula(paste(".~.-",cov,sep="")))
#         1-pchisq(summary(m1)$deviance-summary(m0)$deviance,
#                  summary(m1)$df.residual-summary(m0)$df.residual)
#       })          
#       if(test) print(pvalues)
#       if(length(covs)==1){
#         if(max(pvalues)>p){
#           return(glm(as.formula(paste(response,"~1",sep="")),
#                      data=data,family="binomial"))
#         }else{
#           return(glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                      data=data,family="binomial"))
#         }}else{
#           if(max(pvalues)>p){
#             covs<-covs[-which.max(pvalues)]
#           }else{
#             return(glm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                        data=data,family="binomial"))
#           }}
#     }else if(type=="linear"){
#       pvalues<-sapply(covs,function(cov){
#         if(!boxcox){
#           m0<-lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                  data=subset(data,!is.na(data[,cov])))
#           m1<-update(m0,as.formula(paste(".~.-",cov,sep="")))
#         }else{
#           m0<-boxcoxlm(data[,response],data[,covs])[[1]]
#           m1<-boxcoxlm(data[!is.na(data[,cov]),response],data[!is.na(data[,cov]),setdiff(covs,cov)])[[1]]
#         }
#         1-pchisq(2*(logLik(m0)-logLik(m1)),
#                  summary(m0)$df[1]-summary(m1)$df[1])
#       })          
#       if(test) print(pvalues)
#       if(length(covs)==1){
#         if(max(pvalues)>p){
#           return(lm(as.formula(paste(response,"~1",sep="")),
#                     data=data))
#         }else{
#           return(lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                     data=data))
#         }}else if(length(covs)==2 & boxcox & max(pvalues)>p){
#           return(boxcoxlm(data[,response],data[,covs[which.min(pvalues)]]))
#           
#         }
#       else{
#         if(max(pvalues)>p){
#           covs<-covs[-which.max(pvalues)]
#         }else{
#           if(!boxcox){
#             return(lm(as.formula(paste(response,"~",paste(covs,collapse="+"),sep="")),
#                       data=data))
#           }else{
#             return(boxcoxlm(data[,response],data[,covs]))
#           }
#         }}}
#     
#   }}