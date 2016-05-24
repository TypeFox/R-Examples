#***********************************************************************************************************************************************
#*  
#*  (C) 2010     Andrzej Dudek    Uniwersytet Ekonomiczny we Wrocławiu
#*  
#*  Metoda Ichino wyboru zmiennych w procesie klasyfikacji
#*  Skrypt do książki:
#*  "Podejście symboliczne w analizie danych ekonomicznych" powstałej w ramach projektu badawczego habilitacyjnego N N111 105 234
#*  
#*  Kod poniższy może być modyfikowany, kopiowany i rozprowadzany na warunkach licencji GPL 2 (http://gnu.org.pl/text/licencja-gnu.html), 
#*  a w szczególności pod warunkiem umieszczenia w zmodyfikowanym pliku widocznej informacji o dokonanych zmianach, wraz z datą ich dokonania. 
#*  
#***********************************************************************************************************************************************


IchinoFS.SDA<-function(table.Symbolic){
liczbaZmiennych<-nrow(table.Symbolic$variables)
liczbaObiektow<-nrow(table.Symbolic$individuals)
liczbaKombinacji<-2^liczbaZmiennych-1
wynik<-array(0,c(liczbaKombinacji))
kombinacje<-array(0,c(liczbaKombinacji))
dlugosci<-array(0,c(liczbaKombinacji))
for(i in 1:liczbaKombinacji)
{
	t<-i
	kombinacja=c()
  #print("debug 1")

	for(j in 1:liczbaZmiennych)
	{
		if((t %% 2)==1) kombinacja<-cbind(kombinacja,c(j))
		t<-t%/%2
	}
  #print("debug 2")
#  if(i%%100==1){
    #print(paste("krok",i,"kombinacja",toString(kombinacja)))
#	}
	kombinacje[i]<-toString(kombinacja)
	dlugosci[i]<-length(kombinacja)
	skrzyzowania<-array(FALSE,c(liczbaObiektow,liczbaObiektow))
	for(k in 1:(liczbaObiektow-1))
	for(l in (k+1):liczbaObiektow){
		for(m in kombinacja){
		  #print(paste("przed",k,l,m))
      if(.gprod(table.Symbolic,k,l,m)!=0){
				skrzyzowania[l,k]<-TRUE
				skrzyzowania[k,l]<-TRUE
        }
		  #print(paste("po",k,l,m))
		}
	}
	wektor<-rep(0,liczbaObiektow)
	istniejaNiepowiazane<-TRUE
	skupienie<-0
	while(istniejaNiepowiazane){
		istniejaNiepowiazane<-FALSE
		for(j in 1:liczbaObiektow){
			if(wektor[j]==0){
				istniejaNiepowiazane<-TRUE
				break
			}
		}
		if(istniejaNiepowiazane){
			skupienie<-skupienie+1
			#print(paste("skupienie:",skupienie))
			wektor[j]<-skupienie
			#print(wektor)
			istniejawSkupieniu<-TRUE
			while (istniejawSkupieniu && j!=liczbaObiektow){
				istniejawSkupieniu<-FALSE
				for(k in j:(liczbaObiektow-1)){
				for(l in (k+1):liczbaObiektow){
					if(wektor[k]==skupienie && wektor[l]==0 && skrzyzowania[l,k]){
						wektor[l]<-skupienie
						istniejawSkupieniu<-TRUE
					}
				}
        }
      }
    }
    #print("debug 9")

  }
  #print(paste("liczba skupien",skupienie))
  #print("debug 10")
  wynik[i]<-0
  for(j in 1:max(wektor)){
	wynik[i]<-wynik[i]+sum(wektor==j)*(sum(wektor==j)-1)/2
  }
  #print("debug 11")
}
maksymalneWyniki<-rep(0,liczbaZmiennych)
maksymalneWynikiKombinacje<-rep("",liczbaZmiennych)
maksymalneWynikiRoznice<-rep(0,liczbaZmiennych)
for(i in 2:liczbaZmiennych){
	maksymalneWyniki[i]<-min(wynik[dlugosci==i])
	maksymalneWynikiKombinacje[i]<-kombinacje[((1:liczbaKombinacji)[dlugosci==i])[which.min(wynik[dlugosci==i])]]
	if(i>2)
		maksymalneWynikiRoznice[i]<-maksymalneWyniki[i]-maksymalneWyniki[i-1]
}
print(paste("Kombinacja  ",maksymalneWynikiKombinacje[which.max(maksymalneWynikiRoznice)]))
plot(NULL,xlim=c(2,liczbaZmiennych), ylim=c(min(maksymalneWyniki)-5,max(maksymalneWyniki)+5),ylab="Maksymalne Wyniki")
for(i in 2:(liczbaZmiennych-1)){
  segments(i,maksymalneWyniki[i],i+1,maksymalneWyniki[i])
  segments(i+1,maksymalneWyniki[i],i+1,maksymalneWyniki[i+1])
  text(i,maksymalneWyniki[i],paste("{",maksymalneWynikiKombinacje[i],"}",sep=""))
}
resul<-list(kombinacja=maksymalneWynikiKombinacje[which.max(maksymalneWynikiRoznice)],maksymalneWynikiKombinacje=maksymalneWynikiKombinacje,maksymalneWynikiRoznice=maksymalneWynikiKombinacje,wyniki=wynik)
resul
}

