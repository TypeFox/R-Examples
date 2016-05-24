vottrans <-
function(Ro,Rn,v=1,nw=FALSE){
X2<-Ro; Y2<-Rn
# Ro beinhaltet die Ergebnisse der Vergleichswahl in Absolutzahlen; Rn die der aktuellen Wahl; die erste Spalte muss die Anzahl der Wahlberechtigten enthalten
N1<-nrow(X2); I<-ncol(X2)-1; J<- ncol(Y2)-1 
# N= Anzahl der Gemeinden/Sprengel, I= Anzahl der Parteien bei der Vergleichswahl; I=Anzahl der Parteien bei er aktuellen Wahl

if(nw==TRUE){Y2[,(I+1)]<-Y2[,(I+1)]+Y2[,1]-X2[,1]}
#Option, welche die Differenz der Wahlberechtigten zur Anzahl der NW bei der aktuellen Wahl addiert

loeschen<-function(A){k<-1; b<-ncol(A); while(k<(N1+1)){if(A[k,b]<0){A<-A[-k,]; N1<-N1-1} else{k<-k+1}};return(A)} 

# Loescht Eintraege mit negativen Werten fuer Nichtwaehler

X<-loeschen(X2)
Y<-loeschen(Y2)
N<-nrow(X)

Xabs<-X[,-1];Yabs<-Y[,-1];

#Erstellt Matrix, welche die Ergebnisse ohne die erste Spalte enthaelt

Diag<-function(X){E<-X[,1]; DD<- matrix(0,N,N) ; for(i in 1:N){DD[i,i]<-E[i]};return(DD)}


# erstellt eine Diagonalmatrix deren Eintraege die Anzahl der Wahlberechtigten sind
DD<-Diag(Y)
Xpro<-solve(Diag(X))%*% Xabs
Ypro<-solve(Diag(Y))%*% Yabs

#berechnet die Ergebnismatrix in Prozent


teilcode<-function(A,i){r<-nrow(A); teilcode<-matrix(7,r,1); k<-1; while(k<i){teilcode<-cbind(teilcode, matrix(0,r,I)); k<-k+1}; teilcode<-cbind(teilcode,A); k<-0 ; while(k<(J-i)){teilcode<-cbind(teilcode,matrix(0,r,I));k<-k+1}; teilcode<-teilcode[,-1];return(teilcode)}

#erstellt zeilenweise Teile der geblockten Matrix DD 

geblockt<-function(A){geblockt<-matrix(7,1,ncol(A)*J); k<-1; while(k<(J+1)){geblockt<-rbind(geblockt,teilcode(A,k)); k<-k+1};geblockt<-geblockt[-1,]; return(geblockt)}

#erstellt eine grosze geblockte Matrix mit den Elementen von 'teilcode'

Ygeb<-function(C){Ygeb<-c(1); k<-1; while(k<(J+1)){Ygeb<-matrix(rbind(Ygeb,matrix(c(C[,k]))),ncol=1); k<-k+1}; Ygeb<-Ygeb[-1,]; return(Ygeb)}
Y1<-Ygeb(Ypro)

# blockt die aktuellen Ergebnisse der einzelnen Parteien zu einem Spaltenvektor



summe<-function(F){summe<-c(1); ges<-sum(F[,1]); k<-2; while(k< (ncol(F)+1)){summe<-cbind(summe,sum(F[,k])/ges); k<-k+1}; summe<-summe[,-1]; return(summe)}
summealt<-summe(X)
summeneu<-summe(Y)

#berechnet die Gesamtergebnisse der Vorwahl der einzelnen Parteien


teila<-function(){teila<-matrix(7,I,1); k<-0; while(k<J){teila<-cbind(teila,diag(I)); k<-k+1}; teila<-teila[,-1]; return(teila)}
teilcode1<-function(i){teilcode1<-c(7); k<-1 ; while(k<i){teilcode1<-cbind(teilcode1, matrix(0,1,I)); k<-k+1}; teilcode1<-cbind(teilcode1, matrix(summealt,1)); k<-0; while(k<J-i){teilcode1<-cbind(teilcode1,matrix(0,1,I)); k<-k+1}; teilcode1<-teilcode1[,-1]; return(teilcode1)}

#erstellt ersten Teil der Nebenbedingungsmatrix Amat


teilb<-function(){teilb<-matrix(7,1,I*J); k<-1; while(k<(J+1)){teilb<-rbind(teilb,teilcode1(k)); k<-k+1}; teilb<-teilb[-1,]; return(teilb)}
#erstellt zweiten Teil der Nebenbedingungsmatrix Amat

Amat<-rbind(teila(),teilb(),diag(I*J))

bvec<-rbind(matrix(1,I,1),matrix(summeneu,ncol=1),matrix(0,I*J,1))

if(v ==1){vottrans<-solve.QP(geblockt(t(Xpro)%*%Xabs)*1e-9,t(Y1%*%geblockt(Xabs))*1e-9,t(Amat),bvec,meq=2*I)}
if(v ==2){vottrans<-solve.QP(geblockt(t(Xpro)%*%Xpro)*1e-9,t(Y1%*%geblockt(Xpro))*1e-9,t(Amat),bvec,meq=2*I)}
if(v ==3){vottrans<-solve.QP(geblockt(t(Xabs)%*%Xabs)*1e-9,t(Y1%*%geblockt(Xabs))*1e-9,t(Amat),bvec,meq=2*I)}
solution1<-vottrans[1]

teil<-function(A,n){teil<-c(1); for(k in (I*(n-1)+1):(I*n)){teil<-rbind(teil,A[[1]][k]); k<-k+1}; teil<-teil[-1,]}
gesamtf<-function(A){gesamtf<-matrix(7,I,1); for(k in 1:J){gesamtf<-cbind(gesamtf, teil(A,k)); k<-k+1}; gesamtf<-gesamtf[,-1]; return(gesamtf)}
gesamt<-gesamtf(solution1)

#waehlt die gewuenschte Loesung aus der Ausgabe von solve.QP aus und formt diese in eine Matrix mit den Uebergangswahrscheinlichkeiten um.

return(gesamt)
}
