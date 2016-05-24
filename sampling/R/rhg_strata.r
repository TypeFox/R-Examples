rhg_strata<-function(X,selection)
{
if(is.matrix(X)) X=as.data.frame(X)
m=match(selection,names(X),nomatch=0)
if(sum(m)==0) stop("the 'selection' should be the name of one the X columns")
if(!("Stratum" %in% names(X))) stop("the column 'Stratum' is missing")
result=NULL
u=unique(X$Stratum)
for(i in 1:length(u))
  {si=X[X$Stratum==u[i],]
   x=cbind.data.frame(si$ID_unit,si$status,si[,m])
   names(x)=c("ID_unit","status",names(X)[m])
   result=rbind.data.frame(result,rhg(x,selection))
   }
res = NULL
mm = match(names(X), names(result), nomatch = 0)
index = (1:ncol(X))[mm == 0]
if (length(index) > 0) {
res = cbind.data.frame(X[X$ID_unit==result$ID_unit, index], result)
names(res)[1:length(index)] = names(X)[index]
                        }
res
}

