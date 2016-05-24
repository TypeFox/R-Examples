
HINoV.SDA<-function(table.Symbolic, u=NULL, distance="H", Index="cRAND", method = "pam",...)
{
	z<-table.Symbolic
	if (is.null(u)) stop ("for symbolic data number of classes must be set")
	if (is.null(distance)) stop("For hierarchical methods parameter distance cannot be NULL")
	if (Index != "RAND" && Index!="cRAND") stop("Wrong index type, only RAND or cRAND are allowed")
	if (!is.null(distance))
	{
		if (sum(!any(c(paste("U",2:4,sep="_"),"C_1",paste("SO",1:5,sep="_"),"H","L_1","L_2")==distance)))
			stop("wrong distance")
	}
	liczba_klas=u
	klasyfikacje<-NULL
	cl<-NULL
	for(i in 1:nrow(z$variables))
	{
		#print(i)
		#print(x)
		d <- dist.SDA(z, type=distance,variableSelection=i,...)	
		#print("po dist.symbolic")
		if(method=="pam")
		{
			cl<-pam(d,u,diss=TRUE)$clustering
		}
		else if(casefold(method)=="sclust"){
      print("SCLUST")
      cl<-SClust(z,u,variableSelection=i)
		}
		else if(casefold(method)=="dclust"){
      print("DCLUST")
      cl<-DClust(d,u)
    }
		else
		{
			cl<-cutree(hclust(d,method=method),u)
		}
		klasyfikacje<-cbind(klasyfikacje,cl)
#		#kl<<-klasyfikacje
	}
	liczba_zmiennych<-nrow(z$variables)
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
	if(liczba_zmiennych!=1)
	{
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

#data<-parse.SO("car")
#print(HINoV.SDA(data,u=3,dist="U_2",method="sclust",probType="J"))