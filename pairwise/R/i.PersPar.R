### ------------------------- PersPar -----------------------------
PersPar = function(np,sb,m, imax=30,limit=0.0001, k=length(m),mscs=sum(m-1) ){ # gew. gesch. Personenpar.
  # und Standardfehler se(fw). Ausdruck, wenn iter <= pri
  # argumente:
  # np (number person) vector der länge mscs+1 (anzahl der rohwertgruppen) Anzahl der Personen mit Scoresumme r
  # sb matrix mit  bedingten Schwierigkeitsparametern (sb) aus tau: t(apply(U$tau,1,function(x){cumsum(x)}))
  # imax Anzahl der Iterationen
  # limit genauigkeit
  # pri weggelassen 
  # ergaenzt um : 
  # mscs integer maximale rohwertsumme 
  # m vector der länge k mit anzahl der kategorien (ehemals k)
  # k integer anzahl der items (ehemals m)
  #################### start der funktion ##############
  fw = np*0 # startwerte bei 0
  # fw <- (scale(0:mscs)[,1]) # startwerte zentrierter rohwert vector
  iterproc<-list()
  for(iter in 1:imax){ U = GewLL(np=np,sb=sb,fw=fw, m=m,k=k,mscs=mscs); wLL = U$wLL
                       verb = U$verb; vmax = max(abs(verb)); itmax = iter
                       # if(iter <= pri) cat(iter,vmax,wLL,round(fw,3),sep=" ",fill=T)
                       iterproc[[iter]]<-list(iter=iter,vmax=vmax,wLL=wLL,estimates=fw )
                       if(vmax>10)return(list(itmax=iter,vmax=vmax,wLL=wLL,fw=fw,se.fw=U$f2))# war fehler: se.fw=f2  d.h. das U$ fehlte
                       if(vmax>1)verb = verb/vmax; fw = fw-verb # NR-Schritt
                       if(vmax<limit)break() }                    # Konvergenz
  erg<-list(iterproc=iterproc,itmax=itmax,vmax=vmax,wLL=wLL,fw=fw,se.fw=1/sqrt(-U$f2)) 
  return(erg)
}