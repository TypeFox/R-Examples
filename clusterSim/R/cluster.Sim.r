initial.Centers<-function(x,k)
{
	ic<-1:k
	d<-as.matrix(dist(x,"minkowski",p=2))
	for(i in (k+1):nrow(x))
	{
		d2<-as.vector(d[ic,ic])
		dmn<-min(d2[d2!=0])
		for(j in ic)
		{
			if (sum(d[j,ic]==dmn)!=0)
			{
				n<-j
				for (l in ic)
				{
					if(d[n,l]==dmn)
						m<-l
				}
			}	
		}
		if(min(d[i,ic])>dmn)
		{
			if(d[i,m]<d[i,n]) 
			{
				#print(paste("zmiana",n,i))
				ic[ic==n]<-i
			}
			else
			{
				#print(paste("zmiana",m,i))
				ic[ic==m]<-i
			}
		}
		else
		{
			dq=min(d[i,ic])
			for( l in ic)
			{
				if (d[i,l]==dq)
					q<-l
			}
			d2<-as.vector(d[q,ic])
			miniq<-min(d2[d2!=0])
			#print(paste("i:",i," ||| miniq" ,miniq, " ||| q: ", q, " ||| dq:",dq))
			if(dq > miniq)
			{
				#print(paste("zmiana",q,i))
				ic[ic==q]<-i
			}
		}
	}
	ic	
}



cluster.Sim<-function(x,p=1,minClusterNo,maxClusterNo,icq="S",outputHtml="",outputCsv="",outputCsv2="",normalizations=NULL,distances=NULL,methods=NULL)
{
#require("R2HTML")
#require("e1071")
#require("ade4")
short2LongName <-function(value,fullName=FALSE)
{
longMethods<- c("single","complete","average","mcquitty","pam","ward","centroid","median","k-means")
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
  dim(x)<-c(length(x),1)
}
maxint=100000
#if(!require("cluster")) stop ("Please install cluster")
#if(!require("ade4")) stop ("Please install ade4")
#if(!require("R2HTML")) stop ("Please install R2HTML")
if ((p<1) || (p>9))  stop ("Path must be one of the values:\n 1-ratio data, 2-interval or mixed (ratio & interval), 3-ordinal data,\n 4-nominal data, 5-binary data, 6-ratio data without normalization \n7 - interval or mixed (ratio & interval) data without normalization,\n8 - ratio data with k-means,\n9 - interval or mixed (ratio & interval) data with k-means")
if(icq=="G1" && p>=3 && p<=5) stop ("Calinski Pseudo F statistic may be used only with metric data")
if (minClusterNo<2 || minClusterNo>nrow(x)-1) stop (paste("Number of classes must be in <2 ,",nrow(x)-1,">" ))
if (maxClusterNo<2 || maxClusterNo>nrow(x)-1) stop (paste("Number of classes must be in <2 ,",nrow(x)-1,">" ))
if ((icq == "G3") && (maxClusterNo==nrow(x)-1) ) stop (paste("Number of classes for G3 index can not be equal",nrow(x)-1 ))
if ((icq == "KL") && (maxClusterNo==nrow(x)-1) ) stop (paste("Number of classes for KL index can not be equal",nrow(x)-1 ))
if (minClusterNo>maxClusterNo) stop ("minClusterNo cannot be greater than maxClusterNo")
if (p==5) for (t in x) if(t!=0 && t!=1) stop("Path number five is for binary data only")
if ((p==8 || p==9) && (icq!="G1") && (icq!="KL")) stop("For clustering methods not based on distance matrix like k-means only G1 and K-L indexes can be used")
if ((p==3 || p==4 || p==5) && (icq=="G1" || icq=="KL")) stop ("G1 and K-L indexes may be used only with metric (ordinal or ratio) data")
if (p==1 || p==6 || p==8)
	if (sum(as.matrix(x)<0)!=0)
		stop(paste("for path number",p,"all variables in data matrix have to be measured on ratio scale"))
sciezka=c("Ratio data","Interval or mixed (ratio & interval)","Ordinal data","Nominal data","Binary data","Ratio data without normalization","Interval or mixed (ratio & interval) data without normalization","Ratio data with K-means","Interval or mixed (ratio & interval) data with K-means")
v_norm=c("n1","n2","n3","n4","n5","n6","n7","n8","n9","n10","n11")    # n1-n5 - piec odleglosci
v_norm=cbind(v_norm,c("n1","n2","n3","n4","n5","","","","","",""))
v_norm=cbind(v_norm,c("n0","","","","","","","","","",""))
v_norm=cbind(v_norm,c("n0","","","","","","","","","",""))
v_norm=cbind(v_norm,c("n0","","","","","","","","","",""))
v_norm=cbind(v_norm,c("n0","","","","","","","","","",""))
v_norm=cbind(v_norm,c("n0","","","","","","","","","",""))
v_norm=cbind(v_norm,c("n1","n2","n3","n4","n5","n6","n7","n8","n9","n10","n11"))
v_norm=cbind(v_norm,c("n1","n2","n3","n4","n5","","","","","",""))


s_dist =c("d1","d2","d3","d4","d5","d6","d7","","","")
s_dist =cbind(s_dist,c("d1","d2","d3","d4","d5","","","","",""))
s_dist =cbind(s_dist,c("d8","","","","","","","","",""))
s_dist =cbind(s_dist,c("d9","","","","","","","","",""))
s_dist =cbind(s_dist,c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10"))
s_dist =cbind(s_dist,c("d1","d2","d3","d4","d5","d6","d7","","",""))
s_dist =cbind(s_dist,c("d1","d2","d3","d4","d5","","","","",""))
s_dist =cbind(s_dist,c("","","","","","","","","",""))
s_dist =cbind(s_dist,c("","","","","","","","","",""))

v_dist =c("manhattan","minkowski","maximum","euclidean","gdm1","canberra","bc","","","")
v_dist =cbind(v_dist,c("manhattan","minkowski","maximum","euclidean","gdm1","","","","",""))
v_dist =cbind(v_dist,c("gdm2","","","","","","","","",""))
v_dist =cbind(v_dist,c("sm","","","","","","","","",""))
v_dist =cbind(v_dist,c("1","2","3","4","5","6","7","8","9","10"))
v_dist =cbind(v_dist,c("manhattan","minkowski","maximum","euclidean","gdm1","canberra","bc","","",""))
v_dist =cbind(v_dist,c("manhattan","minkowski","maximum","euclidean","gdm1","","","","",""))
v_dist =cbind(v_dist,c("N.A.","","","","","","","","",""))
v_dist =cbind(v_dist,c("N.A.","","","","","","","","",""))

v_distDesc =c("Manhattan","Euclidean","Chebyschev","Squared Euclidean","GDM1","Canberra","Bray-Curtis","","","")
v_distDesc =cbind(v_distDesc,c("Manhattan","Euclidean","Chebyschev","Squared Euclidean","GDM1","","","","",""))
v_distDesc =cbind(v_distDesc,c("GDM2","","","","","","","","",""))
v_distDesc =cbind(v_distDesc,c("Sokal-Michener","","","","","","","","",""))
v_distDesc =cbind(v_distDesc,c("Binary 1","Binary 2","Binary 3","Binary 4","Binary 5","Binary 6","Binary 7","Binary 8","Binary 9","Binary 10"))
v_distDesc =cbind(v_distDesc,c("Manhattan","Euclidean","Chebyschev","Squared Euclidean","GDM1","Canberra","Bray-Curtis","","",""))
v_distDesc =cbind(v_distDesc,c("Manhattan","Euclidean","Chebyschev","Squared Euclidean","GDM1","","","","",""))
v_distDesc =cbind(v_distDesc,c("N.A.","","","","","","","","",""))
v_distDesc =cbind(v_distDesc,c("N.A.","","","","","","","","",""))

v_method=c("single","complete","average","mcquitty","pam","ward","centroid","median")
v_method=cbind(v_method,c("single","complete","average","mcquitty","pam","ward","centroid","median"))
v_method=cbind(v_method,c("single","complete","average","mcquitty","pam","ward","centroid","median"))
v_method=cbind(v_method,c("single","complete","average","mcquitty","pam","ward","centroid","median"))
v_method=cbind(v_method,c("single","complete","average","mcquitty","pam","ward","centroid","median"))
v_method=cbind(v_method,c("single","complete","average","mcquitty","pam","ward","centroid","median"))
v_method=cbind(v_method,c("single","complete","average","mcquitty","pam","ward","centroid","median"))
v_method=cbind(v_method,c("k-means","","","","","","",""))
v_method=cbind(v_method,c("k-means","","","","","","",""))

l_norm=nrow(v_norm)
l_dist=nrow(v_dist)
l_method=nrow(v_method)
max_norm=0
max_dist=0
max_method=0
max_classes=0
max_cl=NULL
result<-array(0,c(400*(maxClusterNo-minClusterNo+1),7))
i_result<-0
norms_str<-""
for (i in 1:l_norm)
if (v_norm[i,p]!="") {
	if (norms_str=="") norms_str<- paste("\"",v_norm[i,p],"\"",sep="")
	else norms_str<- paste(norms_str,",\"", v_norm[i,p],"\"",sep="")}

dist_str<-""
for (i in 1:l_dist)
if (s_dist[i,p]!="") {
	if (dist_str=="") dist_str<- paste("\"",s_dist[i,p],"\"",sep="")
	else dist_str<- paste(dist_str,",\"", s_dist[i,p],"\"",sep="")}


method_str<-""
if (p==8 || p==9)
	method_str="\"m9\""
else
for (i in 1:l_method)
if (v_method[i,p]!="") {
	if (method_str=="") method_str<- paste("\"m",i,"\"",sep="")
	else method_str<- paste(method_str,",\"m",i,"\"",sep="")}
if(!missing(normalizations))
{
	if (!is.vector(normalizations)) stop (paste("Optional parameter \"normalizations\" for path",p,"(",sciezka[p],")","must be a vector containing only values from the following list: \n",norms_str))
	for (i in 1:length(normalizations))
	{
	if (sum(v_norm[,p]==short2LongName(normalizations[i]))!=1)
		stop (paste("Optional parameter \"normalizations\" for path",p,"(",sciezka[p],")","must be a vector containing only values from the following list: \n",norms_str))
	if (sum(normalizations==normalizations[i])>1)
		stop (paste("Optional parameter \"normalizations\" for path",p,"(",sciezka[p],")","contains doubled values (",normalizations[i],")"))
	}
	for (i in 1:l_norm)
	{
		if (i<=length(normalizations)) v_norm[i,p]=short2LongName(normalizations[i])
		else v_norm[i,p]=""
	}
	norms<-normalizations
}
if(!missing(distances))
{
	if (p==8 || p==9) stop(paste("For path",p,"(",sciezka[p],")","optional parameter \"distances\" is not applicable"))
	if (!is.vector(distances)) stop (paste("Optional parameter \"distances\" for path",p,"(",sciezka[p],")","must be a vector containing only values from the following list: \n",dist_str))
	for (i in 1:length(distances))
	{
	if (sum(v_dist[,p]==short2LongName(distances[i]))!=1)
		stop (paste("Optional parameter \"distances\" for path",p,"(",sciezka[p],")","must be a vector containing only values from the following list: \n",dist_str))
	if (sum(distances==distances[i])>1)
		stop (paste("Optional parameter \"distances\" for path",p,"(",sciezka[p],")","contains doubled values (",distances[i],")"))
	}
	t_dist<-v_dist[,p]
	for(i in 1:l_dist) 
		v_dist[i,p]<-""
	for (i in 1:l_dist)
	{
		for(j in 1:length(distances))
		if(t_dist[i]==short2LongName(distances[j]))
			v_dist[i,p]<-short2LongName(distances[j])
		
	}
}
if(!missing(methods))
{
	if (!is.vector(methods)) stop (paste("Optional parameter \"methods\" for path",p,"(",sciezka[p],")","must be a vector containing only values from the following list: \n",method_str))
	for (i in 1:length(methods))
	{
	if (sum(v_method[,p]==short2LongName(methods[i]))!=1)
		stop (paste("Optional parameter \"methods\" for path",p,"(",sciezka[p],")","must be a vector containing only values from the following list: \n",method_str))
	if (sum(methods==methods[i])>1)
		stop (paste("Optional parameter \"methods\" for path",p,"(",sciezka[p],")","contains doubled values (",methods[i],")"))

	}
	for (i in 1:l_method)
	{
		if (i<=length(methods)) v_method[i,p]<-short2LongName(methods[i])
		else v_method[i,p]<-""
	}

}


if (icq=="G3")
{
	wynik_global=maxint
}
else
{
	wynik_global=-maxint
}
stepG2=0
if (icq=="G2")
	print("Calculations of G2 internal cluster quality index may be sligthly slow. Iteration step informations will be diplayed")
d<-array(0,c(1))
czas1<-Sys.time()
for (i_norm in 1 : l_norm)
{
	if(v_norm[i_norm,p]!="")
	{
	z<-data.Normalization(x, v_norm[i_norm,p])
	if(getRversion() >= '2.14'){
	zz<-as.matrix(z)
	}
	else{
	zz<-z
	}
	if(sum(is.nan(zz))!=0 || sum(is.na(zz))!=0 || sum(is.infinite(zz))!=0){
    stop(paste("Na/NaN/Infinite after \"", v_norm[i_norm,p] , "\" normalization \nPlease exclude this normalization from simulation \nexplicitly use normalizations = c(...)",sep=""))
	}
	z<-as.data.frame(z)
	for (i_dist in 1:l_dist)
	{
		if(v_dist[i_dist,p]!="")
		if(!(((v_dist[i_dist,p]=="bc" || v_dist[i_dist,p]=="canberra")&&(v_norm[i_norm,p]=="n1" 
		|| v_norm[i_norm,p]=="n2"  || v_norm[i_norm,p]=="n3" || v_norm[i_norm,p]=="n4" || v_norm[i_norm,p]=="n5"))))  
		{
		if (v_dist[i_dist,p]=="gdm1")
		{
			y<-as.matrix(z)
			d<-dist.GDM(y)
		}
		else
		if (v_dist[i_dist,p]=="gdm2")
		{
			y<-as.matrix(z)
			d<-GDM2(y)
		}
		else
		if (v_dist[i_dist,p]=="bc")
		{
			y<-as.matrix(z)
			d<-dist.BC(y)
		}
		else
		if (v_dist[i_dist,p]=="sm")
		{
			y<-as.matrix(z)
			d<-dist.SM(y)
		}
		else
		if (v_dist[i_dist,p]=="1" || v_dist[i_dist,p]=="2"
|| v_dist[i_dist,p]=="3" || v_dist[i_dist,p]=="4" || v_dist[i_dist,p]=="5"
|| v_dist[i_dist,p]=="6" || v_dist[i_dist,p]=="7" || v_dist[i_dist,p]=="8"
|| v_dist[i_dist,p]=="9" || v_dist[i_dist,p]=="10" )
		{
			y<-as.matrix(z)
			d<-dist.binary(z,v_dist[i_dist,p])
		}
		else
		if (v_dist[i_dist,p]=="N.A.")
		{
			y<-as.matrix(z)
		}
		else
		{
			if (v_dist[i_dist,p]=="minkowski")
			{
				d <- dist(z, method=v_dist[i_dist,p],p=2)			
			}
			else
			{
				d <- dist(z, method=v_dist[i_dist,p])			
				if (v_dist[i_dist,p]=="euclidean"){
          d<-d^2
				}
			}
		}
		if(v_dist[i_dist,p]=="euclidean")
			lmax_method=l_method
		else
			lmax_method=l_method-3
		if (is.nan(d) || is.nan(d) || is.infinite(d))
		{
			(paste("Nieskonczonosc",i_norm,i_dist))
			print(d)
			stop()
		}	
		#d[is.infinite(d)]<-maxint
		#d[is.nan(d)]<-maxint
		if(sum(is.na(d))!=0 || sum(is.nan(d))!=0 || sum(is.infinite(d))!=0 ){
      stop(paste("Na/NaN/Infinite in distance \"", s_dist[i_dist,p] , "\" calculation\nPlease exclude this distance from simulation \nexplicitly use distances = c(...)",sep="")) 
    }
		for (i_method in 1:lmax_method)
		{		
			if(v_method[i_method,p]!="")
			{
			#print(paste(p,Sys.time(),v_norm[i_norm,p],v_dist[i_dist,p],v_method[i_method,p]))
			if(v_method[i_method,p]=="pam")
			{
				if (icq=="S")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<- pam(d, liczba_klas,diss=TRUE)
					t<-index.S(d,c$clustering)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c$clustering
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}
				

				if (icq=="G2")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<- pam(d, liczba_klas,diss=TRUE)
					t<-index.G2(d,c$clustering)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c$clustering
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
					stepG2=stepG2+1
					#print(paste("I've accomplished",stepG2,"step. Current time is:",Sys.time()))
				}}



				if (icq=="G3")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<- pam(d, liczba_klas,diss=TRUE)
					t<-index.G3(d,c$clustering)
					if(t<wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_classes=liczba_klas
						max_method=i_method
						max_cl<-c$clustering
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t

				}}

				if (icq=="G1")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<- pam(d, liczba_klas,diss=TRUE)
					t<-index.G1(z,c$clustering)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c$clustering
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}

				if (icq=="KL")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c1<- pam(d, liczba_klas-1,diss=TRUE)
					c2<- pam(d, liczba_klas,diss=TRUE)
					c3<- pam(d, liczba_klas+1,diss=TRUE)
					clall<-cbind(c1$clustering,c2$clustering,c3$clustering)
					t<-index.KL(z,clall)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c2$clustering
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}





			}		
			else # KMEANS
			if (v_method[i_method,p]=="k-means")
			{
				if (icq=="S")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<-kmeans(y,y[initial.Centers(y,liczba_klas),],100)	# pomyslec o seeds
					t<-index.S(d,c$cluster)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c$cluster
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t

				}}

				if (icq=="G2")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<-kmeans(y,y[initial.Centers(y,liczba_klas),],100)	# pomyslec o seeds
					t<-index.G2(d,c$cluster)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c$cluster
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
					stepG2=stepG2+1
					#print(paste("I've accomplished",stepG2,"step. Curretnt time is:",Sys.time()))
				}}


				if (icq=="G3")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<-kmeans(y,y[initial.Centers(y,liczba_klas),],100)	# pomyslec o seeds
					t<-index.G3(d,c$cluster)
					if(t<wynik_global)
					{
						wynik_global<-t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c$cluster
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}

 				if (icq=="G1")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					#print(y)
					c<-kmeans(y,y[initial.Centers(y,liczba_klas),],100)	# pomyslec o seeds
					t<-index.G1(z,c$cluster)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c$cluster
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}

 				if (icq=="KL")
				{
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					if (liczba_klas==2)
						c1cl<-array(1,c(nrow(x)))
					else
					{
						c1<-kmeans(y,y[initial.Centers(y,liczba_klas-1),],100)	# pomyslec o seeds
						c1cl=c1$cluster
					}
					c2<-kmeans(y,y[initial.Centers(y,liczba_klas),],100)	# pomyslec o seeds
					c3<-kmeans(y,y[initial.Centers(y,liczba_klas+1),],100)	# pomyslec o seeds
					clall<-cbind(c1cl,c2$cluster,c3$cluster)
					t<-index.KL(z,clall)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c2$cluster
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}

			}
			else
			{
				if (icq=="S")
				{
				hc <- hclust(d, method = v_method[i_method,p])
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<-cutree(hc,k=liczba_klas)
					t<-index.S(d,c)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c
					}
					if(v_dist[i_dist,p]=="bc")
					{
						
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t

				}}

				if (icq=="G2")
				{
				hc <- hclust(d, method = v_method[i_method,p])
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<-cutree(hc,k=liczba_klas)
					t<-index.G2(d,c)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
					stepG2=stepG2+1
				}}


				if (icq=="G3")
				{
				hc <- hclust(d, method = v_method[i_method,p])
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<-cutree(hc,k=liczba_klas)
					t<-index.G3(d,c)
					if(t<wynik_global)
					{
						wynik_global<-t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}

 				if (icq=="G1")
				{
				hc <- hclust(d, method = v_method[i_method,p])
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c<-cutree(hc,k=liczba_klas)
					t<-index.G1(z,c)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}


 				if (icq=="KL")
				{
				hc <- hclust(d, method = v_method[i_method,p])
				for (liczba_klas in minClusterNo:maxClusterNo)
				{
					c1<-cutree(hc,k=liczba_klas-1)
					c2<-cutree(hc,k=liczba_klas)
					c3<-cutree(hc,k=liczba_klas+1)
					clall<-cbind(c1,c2,c3)
					t<-index.KL(z,clall)
					if(t>wynik_global)
					{
						wynik_global=t
						max_norm=i_norm
						max_dist=i_dist
						max_method=i_method
						max_classes=liczba_klas
						max_cl<-c2
					}
					i_result<-i_result+1
					result[i_result,1]<-i_result
					result[i_result,2]<-liczba_klas
					result[i_result,3]<-i_norm
					result[i_result,4]<-i_dist
					result[i_result,5]<-i_method
					result[i_result,6]<-t
				}}


			}}
		}}
	}}
}
czas2<-Sys.time()
result<-result[1:i_result,]
if (icq=="G3")
	t<-order(result[,6],decreasing=FALSE)
else
	t<-order(result[,6],decreasing=TRUE)
for (i in 1:i_result)
{
	result[t[i],7]<-as.double(i)
}
sorted<-result
sorted<-sorted[t,]

#Uciecie posortowanych danych
if (i_result>10)
	ile<-10
else
	ile<-i_result
while(ile<i_result && sorted[ile,6]==sorted[ile+1,6]) ile<-ile+1
#sorted<-sorted[1:ile,]


if(p>=3 && p<=5)
	v_norm[,]<-"N.A."
if(p==6 || p==7)
	v_norm[,]<-"Without normalization"

tt<-sorted
for(i in 1:nrow(sorted))
{
	sorted[i,3]=v_norm[as.integer(tt[i,3]),p]
	sorted[i,4]=v_distDesc[as.integer(tt[i,4]),p]
	sorted[i,5]=v_method[as.integer(tt[i,5]),p]
}
tt<-result
for(i in 1:nrow(result))
{
	result[i,3]=v_norm[as.integer(tt[i,3]),p]
	result[i,4]=v_distDesc[as.integer(tt[i,4]),p]
	result[i,5]=v_method[as.integer(tt[i,5]),p]
}

sorted<-cbind(sorted[,7],sorted[,1:6])
sorted<-as.data.frame(sorted)
result<-as.data.frame(result)

indexName<-switch(icq,"S"="Silhouette","G2"="G2 index","G3"="G3 index","G1"="Calinski-Harabasz index ","KL"="Krzanowski-Lai index")
names(result)<-c(" No. "," No. of clusters "," Normalization formula "," Distance measure "," Clustering method ",indexName,"Rank")
names(sorted)<-c(" Rank "," No. "," No. of clusters "," Normalization formula "," Distance measure "," Clustering method ",indexName)

if (outputHtml!="")
{
	target <- HTMLInitFile(".",filename=outputHtml)
	options("R2HTML.format.decimal.mark"=",")
	HTMLChangeCSS("Pastel")
	HTML.title("RESULTS OF CLASSIFICATIONS ",file=target,align="center",HR=3,color="00FF00",align="center")
	HTML.title(paste("PATH = ",p," (",sciezka[p],")"),HR=4)
	HTML.title(paste("INDEX = ",indexName),HR=4)
	HTML.title(paste("NO. OF CLUSTERS = <",minClusterNo,";",maxClusterNo,">"),HR=4)
	HTML.title("Unsorted results",HR=3)
	HTML(result,file=target,innerborder=1,Border=1,align="center")
	HTML.title("Sorted results",HR=3)
	HTML(sorted,file=target,innerborder=1,Border=1,align="center")
	HTML.title("Estimated calculation time :",HR=5)
	HTML.title(czas2-czas1,HR=5)
	#browseURL(target)
	HTMLEndFile(target)
}
if (outputCsv!="")
{	write.table(result,paste(outputCsv,".csv",sep=""),dec=".",sep=",",row.names=FALSE)
}
if (outputCsv2!="")
{	write.table(result,paste(outputCsv2,".csv",sep=""),dec=",",sep=";",row.names=FALSE)
}
resul<-list(path=p,result=wynik_global,normalization=v_norm[max_norm,p],distance=v_distDesc[max_dist,p],method=v_method[max_method,p],classes=as.vector(sorted[1,3]),optClustering=max_cl,optClusteringDescription=cluster.Description(as.matrix(x),max_cl),time<-(czas2-czas1))
resul
}
