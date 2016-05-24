rhg<-function(X,selection)
{
if(missing(selection)) stop("'selection' is missing")
if(!is.data.frame(X)) X=as.data.frame(X)
if(ncol(X)<2) stop("the input is incorrect")
if(any(is.na(X))) stop("the input should not contain missing values")
if(is.null(names(X))) stop("specify the names of columns in X")
if(!all(c("ID_unit","status") %in% names(X))) stop("ID_unit and/or status are not columns of X") 
m=match(selection,names(X),nomatch = 0)
if(sum(m)==0) stop("'selection' should be the name of one of the X's columns")
X1=cbind.data.frame(X$status,X[,m])
names(X1)=c("status",names(X)[m])
x=unique(X1[,2])
if(ncol(X1)>=3)
for(i in 3:ncol(X1))
  x=list(x,unique(X1[,i]))
x=expand.grid(x)
ng=1
prob=rhgroup=numeric(nrow(X1))
for (i in 1:nrow(x)) {
            expr=rep(FALSE, nrow(X1))
            for(j in 1:nrow(X1)) {
                                  expr[j] = all(X1[j,2:ncol(X1)] == x[i, ])
                                  if(expr[j]) rhgroup[j]=ng
                                 }
            if(any(expr)) ng=ng+1 
                     }
gr=unique(rhgroup)
if(is.data.frame(X1))
X1=cbind.data.frame(X1,rhgroup)
else X1=cbind(X1,rhgroup)
for(i in 1:length(gr))
  {l=nrow(X1[X1[,ncol(X1)]==gr[i],]) 
   lr=nrow(X1[X1[,ncol(X1)]==gr[i] & X1[,1]==1,])
   for(j in 1:length(prob)) 
             if(rhgroup[j]==gr[i] & X1[j,1]==1) prob[j]=lr/l
   }
result=cbind.data.frame(X$ID_unit,X1,prob)
names(result)=c("ID_unit",names(X1),"prob_resp")
res = NULL
mm = match(names(X), names(result), nomatch = 0)
if(0 %in% mm)
{index = (1:ncol(X))[mm == 0]
res = cbind.data.frame(X[X$ID_unit==result$ID_unit, index], result)
names(res)[1:length(index)] = names(X)[index]
}
else res=result
res
}
