RAI <-
function(x,conf.levels = 0.95){
 d = dim(x) 
if(is.null(d)) stop("argument is not a matrix")
if(conf.levels>1 || conf.levels<0){stop("magnitude of the returned confidence interval must be a single number between 0 and 1")}
nsm = 1000
catg = sort(unique(c(x))) #categories
C = length(catg) # dim categories
psj = matrix(0,nsm,C) #proportion of agreement specific to category 
po = matrix(0,1,nsm) #proportion of observed agreement
K = d[1] #rated cases

ii = 1:K
for (sm in 1:nsm){
njk = matrix(0,K,C) #number of actual agreements
for(j in 1:C){
for (k in 1:K){
njk[k,j] = length(which(x[ii[k],] == catg[j]))
}
}

nk =  apply(njk,1,sum) # total number of ratings made on case k

Sj = matrix(0,1,C) #total number of agreements specifically on rating level
Spossj = matrix(0,1,C) # number of possible agreements on category
for(j in 1:C){
for (k in 1:K){
Sj[,j] = Sj[,j] + njk[k,j]*(njk[k,j]-1)
Spossj[,j] = Spossj[,j] + njk[k,j]*(nk[k]-1)
}
}

psj[sm,] = Sj/Spossj

O = sum(Sj) # total number of actual agreements
Oposs = sum(nk*(nk-1)) #total number of possible agreements

po[sm] = O/Oposs
 ii = round(runif(d[1],1,d[1]))
}
 q1 = (1-conf.levels)/2
 q2 = 1-q1
 qs = matrix(0,C,2)
 for (i in 1:C){
 qs[i,] = quantile(psj[,i],c(q1,q2),na.rm=TRUE)
}
 qo = quantile(po,c(q1,q2))
 dfo = data.frame(po[1],qo[1],qo[2],row.names=NULL)
 colnames(dfo) = c('value',paste('(',q1*100,"%"),paste('-',q2*100,"%",')'))

 dfs = data.frame(psj[1,],qs[,1],qs[,2],row.names=catg)
 colnames(dfs) = c('value',paste('(',q1*100,"%"),paste('-',q2*100,"%",')'))


return(list(Subjects=K,Observers=max(nk),Overall_agreement= dfo, Categories = catg,Specific_agreement = dfs ))
}
