#***********************************************************************************************************************************************
#*  
#*  (C) 2010     Andrzej Dudek    Uniwersytet Ekonomiczny we Wrocławiu
#*  
#*  Metody klasyfikacji dynamicznej dla danych symbolicznych
#*  Skrypt do książki:
#*  "Podejście symboliczne w analizie danych ekonomicznych" powstałej w ramach projektu badawczego habilitacyjnego N N111 105 234
#*  
#*  Kod poniższy może być modyfikowany, kopiowany i rozprowadzany na warunkach li-cencji GPL 2 (http://gnu.org.pl/text/licencja-gnu.html), 
#*  a w szczególności pod warunkiem umieszczenia w zmodyfikowanym pliku widocznej informacji o dokonanych zmianach, wraz z datą ich dokonania. 
#*  
#***********************************************************************************************************************************************

.medoid<-function(dist,cl,klasa)
{
	min_i<-NA
	min_dist<-NA
	for (i in 1:nrow(dist))
	{

		if (klasa!=0)
		{
			i_sum<-sum(dist[i,cl==klasa])
			if (is.na(min_dist) || ((i_sum<min_dist) && (cl[i]==klasa)))
			{
				min_dist<-i_sum
				min_i<-i			
			}
		}
		else
		{
			i_sum<-sum(dist[i,])
			if (is.na(min_dist) || i_sum<min_dist)
			{
				min_dist<-i_sum
				min_i<-i			
			}
		}
	}
	min_i
}


DClust<-function(dist,cl,iter=100)
{
  medoids<-1:cl
	dist<-as.matrix(dist)
	#print(ncol(dist))
	#print(cl)
	if(length(cl)!=ncol(dist)){
    cl<-cl_result<-rep(1:cl,ncol(dist)%/%cl+1)[1:ncol(dist)]
	}
	#print(cl_result)
	test<-TRUE
	k<-max(cl_result)
	test_licznik<-0
	while(test)
	{
		test<-FALSE
		test_licznik<-test_licznik+1
		for(i in 1:length(cl_result))
		{
      distMedoids<-dist[i,medoids]
      min_j<-which.min(distMedoids)
      #print(distMedoids)
      #print(paste("which.min",i,min_j))
			if (min_j!=cl_result[i])
			{
				test<-TRUE
				cl_result[i]<-min_j
			}
		}
    for(i in 1:max(cl)){
      medoids[i]<-.medoid(dist,cl_result,i)
    }
		if (test_licznik>iter) test<-FALSE
	}
	if (test_licznik<iter)
	{
		if (sum(cl==cl_result)==length(cl))
		{
			
			#print (paste("DCLUST SKONCZONY PRAWIDLOWO ALE NIC NIE ZOSTALO ZMIENIONE"))
			#print(cl)
			print(cl_result)
		}
		else
		{		
			#print (paste("DCLUST SKONCZONY PRAWIDLOWO PO ",test_licznik," ITERACJACH"))
		}
	}
	else
	{ 
	#print ( paste("DCLUST PRZEKROCZYL "),iter,"iteracji")
	}
	cl_result
}

.SDist<-function(table.Symbolic,variableSelection,objectSelection){
  resul<-array(0,c(length(objectSelection),length(objectSelection)))
  for(k in variableSelection){
    if(table.Symbolic$variables[k,"type"]=="IC"){
          t<-table.Symbolic$indivIC[objectSelection,k,]
          dim(t)<-c(length(objectSelection),1,2)
          resul<-resul+as.matrix(dist.Symbolic(t,type="H")^2)
          #print(resul)
          #stop("test")
    }
    if(table.Symbolic$variables[k,"type"]=="NM"){
          resul<-resul+as.matrix(dist.SDA(table.Symbolic,probType="CHI",power=1,variableSelection=k))[objectSelection,objectSelection]
          #print(resul)
    }
    if(table.Symbolic$variables[k,"type"]=="N"){
          resul<-resul+as.matrix(dist.SDA(table.Symbolic,type="L_2",power=1,variableSelection=k))[objectSelection,objectSelection]
          #print(resul)
    }
    if(table.Symbolic$variables[k,"type"]=="MN"){
          #print(as.matrix(dist.SDA(table.Symbolic,type="L_2",power=1,variableSelection=k)))
          resul<-resul+as.matrix(dist.SDA(table.Symbolic,type="L_2",power=1,variableSelection=k))[objectSelection,objectSelection]
          #print(resul)
    }
  }
  resul
}

SClust<-function(table.Symbolic,cl,iter=100,variableSelection=NULL,objectSelection=NULL){
  if(is.null(variableSelection)){
    variableSelection<-1:nrow(table.Symbolic$variables)
  }
  if(is.null(objectSelection)){
    objectSelection<-1:nrow(table.Symbolic$individuals)
  }
    dist<-.SDist(table.Symbolic,variableSelection,objectSelection)
    DClust(dist,cl,iter)
}


