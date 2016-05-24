#HINoV.metric<-function(x, type="metric", s = 2, u, #distance=NULL, method = "kmeans", Index = "cRAND")
#{
#	HINoV.MOD(x, type="metric", s, u, distance, #method, Index)
#}
#HINoV.nonmetric<-function(x, type="nonmetric", s = 2, u, #distance=NULL, method = "kmeans", Index = "cRAND")
#{
#	HINoV.MOD(x, type="nonmetric", s, u, distance, #method, Index)
#}

HINoV.Mod<-function(x, type="metric", s = 2, u=NULL, distance=NULL, method = "kmeans", Index = "cRAND")
{
short2LongName <-function(value,fullName=FALSE)
{
longMethods<- c("single","complete","average","mcquitty","pam","ward","centroid","median","k-means","diana")
longDistances<-c("manhattan","minkowski","maximum","euclidean","gdm1","canberra","bc","gdm2","sm")
fullDistances <-c("Manhattan","Euclidean","Chebyschev","Squared Euclidean","GDM1","Canberra","Bray-Curtis","GDM2","Sokal-Michener")

longBinaryDistances <-c("1","2","3","4","5","6","7","8","9","10")
fullBinaryDistances <-c("Binary 1","Binary 2","Binary 3","Binary 4","Binary 5","Binary 6","Binary 7","Binary 8","Binary 9","Binary 10")

if (value=="") l2SN<-value
else
{
type=substr(value,1,1)
index=as.integer(substr(paste(value," ",sep=""),2,3))

if (type=="n")
	l2SN<-value
if (type=="m")
{
	l2SN<-longMethods[as.integer(index)]
}
if (type=="d")
{
	if (fullName==TRUE)
		l2SN<-fullDistances[index]
	else
	{	
		l2SN<-longDistances[index]
	}
}
if (type=="b")
{
	if (fullName==TRUE)
		l2SN<-fullBinaryDistances[index]
	else
		l2SN<-longBinaryDistances[index]
}
}
l2SN
}
	if(is.null(dim(x))){
    stop("Number of variables have to be greater than one")
	}
	z<-x
	if (is.null(u) && (type!="nonmetric")) stop ("for metric and mixed data number of classes must be set")
	#if(!require("cluster")) stop ("Please install cluster package")
	#if(!require("e1071")) stop ("Please install e1071 package")
	#if(!require("ade4")) stop ("Please install ade4 package")
	if (is.null(distance) && method != "kmeans" && method!="pam" && method != "diana") stop("For hierarchical methods parameter distance cannot be NULL")
	if (Index != "RAND" && Index!="cRAND") stop("Wrong index type, only RAND or cRAND are allowed")
	if((length(type)>1) && (length(type)< ncol(z))) stop ("Wrong length of type parameter, must be equal to number of variables")
	if (!is.null(distance))
	{
		if (sum(c("d1","d2","d3","d4","d5","d6","d7")==distance)==0)
			stop("wrong distance")
		if ((s!=1) && ((distance=="d6")||(distance=="d7")))
			stop("wrong distance for non-ratio data")
	}
	liczba_klas=u
	klasyfikacje<-NULL
	cl<-NULL
	for(i in 1:ncol(z))
	{
		x<-z[,i]
		y<-as.matrix(x)
		if (length(type)>1)
		{
			if (type[i]=="m") ttype<-"metric"
			else ttype<-"nonmetric"
		}
		else
		{
			ttype<-type
		}
		if (!is.null(distance))
		{
			if (distance=="d5")
			{
				d<-dist.GDM(y)
			}
			else
			if (distance=="d7")
			{
				d<-dist.BC(y)
			}
			else
			{
				
				d <- dist(y, method=short2LongName(distance))				
			}
		}
		if (ttype=="metric")
		{
			if (method=="kmeans")
			{
				cl<-kmeans(y,y[initial.Centers(y,u),])$cluster
			}
			else if(method=="pam")
			{
				if (is.null(distance))
				{
					cl<-pam(y,u,diss=FALSE)$clustering
				}
				else
				{
					cl<-pam(d,u,diss=TRUE)$clustering
				}
			}
      else if(method=="diana"){
          cl<-cutree(as.hclust(diana(d)),k=u)

      }
			else
			{
				cl<-cutree(hclust(d,method=method),u)
			}
		}
		else
			cl<-as.integer(y)
		klasyfikacje<-cbind(klasyfikacje,cl)
	}
	#print("DEBUG: A6")
	liczba_zmiennych=ncol(z)
	wynik<-array(0,c(liczba_zmiennych,liczba_zmiennych))
	for (i in 1:liczba_zmiennych)
	for (j in 1:liczba_zmiennych)
	{
		#print(paste(1,j))
		if(i==j)
		{
			wynik[i,j]=1
		}
		else
		{
	#print("DEBUG: A8")
			#print(klasyfikacje[,i])
			#print(klasyfikacje[,j])
			tk=table(klasyfikacje[,i],klasyfikacje[,j])
	#print("DEBUG: A9")
			w<-classAgreement(tk)
			#print(w)
			if (Index=="cRAND")
			wynik[i,j]<-w$crand
			else
			wynik[i,j]<-w$rand
		}
	}
	#print("DEBUG: A10")
	posortowane<-array(0,c(2,liczba_zmiennych))
	for (i in 1:liczba_zmiennych)
	{
		posortowane[1,i]<-i
		posortowane[2,i]<-sum(wynik[i,])-1
	}
	topri<-posortowane
	if(liczba_zmiennych>1){
    for (i in 1:liczba_zmiennych)
    for (j in 1:(liczba_zmiennych-1))
    {
      if(posortowane[2,j]<posortowane[2,j+1])
      {
        p1<-posortowane[1,j+1]
        p2<-posortowane[2,j+1]
        posortowane[1,j+1]<-posortowane[1,j]
        posortowane[2,j+1]<-posortowane[2,j]
        posortowane[1,j]<-p1
        posortowane[2,j]<-p2
      }
    }
  }
	resul<-list(parim=wynik,topri=t(topri),stopri=t(posortowane))
	resul
}





