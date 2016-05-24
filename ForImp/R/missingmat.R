missingmat <-
function(mat, nummissing, pattern="r", nk=1, p=0.1, w=3) 
{
n<-nrow(mat)

m<-ncol(mat)

matmiss<-mat

###
# missing totally at random
###

if (pattern=="r")
{

posmissings<-sample(n*m, nummissing, replace=FALSE)

rowvalues<-matrix(data = matmiss, ncol=1)

rowvalues[posmissings]=NA

matmiss<-matrix(data = rowvalues, nrow=n, ncol=m, byrow=FALSE)
}

# missing on the lowest value (rate p) on variables nk

else if (pattern=="l")
{
minimum<-min(mat[,nk])

lowest<-which(mat[,nk]==minimum)

l<-sample(lowest, round(length(lowest)*p), replace=FALSE)

matmiss[l,nk]<-NA
}

# "block" missing

else if(pattern=="b")
{
posmissings<-sample(n, nummissing, replace=FALSE)

matmiss[posmissings,nk]<-NA
}

# missing MNAR

else if (pattern=="n")
{
weights<-rep(1,n)

minimum<-min(mat[,nk])

weights[matmiss[,nk]==minimum]<-w

pik<-weights/sum(weights)*nummissing

# UPpivotal: see package sampling

posmissingsindex<-UPpivotal(pik)

posmissings<-(1:length(weights))[posmissingsindex>0.5]

matmiss[posmissings, nk]<-NA
}
matmiss
}

