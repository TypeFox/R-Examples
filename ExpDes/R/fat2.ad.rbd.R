fat2.ad.rbd <-
function(factor1, factor2, block, resp, respAd, quali=c(TRUE,TRUE), mcomp='tukey', fac.names=c('F1','F2'), sigT=0.05, sigF=0.05) {

cat('------------------------------------------------------------------------\nLegend:\n')
cat('FACTOR 1: ',fac.names[1],'\n')
cat('FACTOR 2: ',fac.names[2],'\n------------------------------------------------------------------------\n\n')


fatores<-cbind(factor1,factor2)
Fator1<-factor(factor1)
Fator2<-factor(factor2)
Block<-factor(block)
nv1<-length(summary(Fator1))   #Diz quantos niveis tem o fator 1.
nv2<-length(summary(Fator2))   #Diz quantos niveis tem o fator 2.
lf1<-levels(Fator1)
lf2<-levels(Fator2)
J=length(respAd)
n.trat2<-nv1*nv2

#ANAVA do fatorial 2
anavaF2<-summary(aov(resp~ Block + Fator1*Fator2))

SQa<-anavaF2[[1]][2,2]
SQb<-anavaF2[[1]][3,2]
SQab<-anavaF2[[1]][4,2]

#Anava de todos os tratamentos do experimento (fatorial 2 + adicional)
col1<-numeric(0)
for(i in 1:n.trat2) {
col1<-c(col1, rep(i,J))     
                    }
col1<-c(col1,rep('ad',J))
col2<-c(block,rep(1:J))
col3<-c(resp,respAd)
tabF2ad<-data.frame("TRAT2"=col1, "BLOCK"=col2, "RESP2"=col3)
TRAT2<-factor(tabF2ad[,1])
BLOCK<-factor(tabF2ad[,2])
anava<-aov(tabF2ad[,3]~ BLOCK + TRAT2)
anavaTr<-summary(anava)

SQB<-anavaTr[[1]][1,2]
SQad<-anavaTr[[1]][2,2] - (SQa+SQb+SQab)
SQE<-anavaTr[[1]][3,2]
SQT<-anavaTr[[1]][1,2]+anavaTr[[1]][2,2]+anavaTr[[1]][3,2]

glB=J-1
gla=nv1-1
glb=nv2-1
glab=(nv1-1)*(nv2-1)
glad=1
glE=(nv1*nv2+1)*(J-1)-(J-1)   
glT=(nv1*nv2+1)*J-1

QMB=SQB/glB
QMa=SQa/gla
QMb=SQb/glb
QMab=SQab/glab
QMad=SQad/glad
QME=SQE/glE
QMT=SQT/glT

FcB=QMB/QME
Fca=QMa/QME
Fcb=QMb/QME
Fcab=QMab/QME
Fcad=QMad/QME

pv.fs=c(1-pf(Fca,gla,glE), 1-pf(Fcb,glb,glE))

#Montando a tabela da ANAVA
anavaT<-data.frame("DF"=c(glB,gla, glb, glab, glad, glE, glT ),
"SS"=c(round(c(SQB,SQa,SQb,SQab,SQad,SQE,SQT),5)),
"MS"=c(round(c(QMB,QMa,QMb,QMab,QMad,QME,QMT),5)),
"Fc"=c(round(c(FcB,Fca,Fcb,Fcab,Fcad),4),'',''),
"Pr>Fc"=c(round(c(1-pf(FcB,glB,glE), pv.fs, 1-pf(Fcab,glab,glE), 1-pf(Fcad,glad,glE)),4),' ', ' '))
colnames(anavaT)[5]="Pr>Fc"
rownames(anavaT)=c("Block",fac.names[1],fac.names[2],paste(fac.names[1],'*',fac.names[2],sep=''),"Ad vs Factorial","Residuals","Total")
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(anavaT)
cat('------------------------------------------------------------------------\n\n')

#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------\nShapiro-Wilk normality test\n')
cat('p-value: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('WARNING: at 5% of significance, residuals can not be considered normal!
------------------------------------------------------------------------\n')}
else{cat('According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.
------------------------------------------------------------------------\n')}

#Contraste Ad vs Fatorial
cat('Contrast of the additional treatment with the factorial
------------------------------------------------------------------------\n')
x<-mean(respAd)
y<-mean(resp)

if(1-pf(Fcad,glad,glE)>sigF) { C1<-data.frame("Means"=c(x,y))
                                  rownames(C1)=c("Additional","Factorial")
                                  colnames(C1)<-c("Means")
                                  cat('According to the F test, the means of the two groups are statistical equal.\n')
                                  print(C1) }  else{
C2<-data.frame("Mean"=c(x,y),
" "=c(letters[1],letters[2]))
rownames(C2)=c("Additional","Factorial")
colnames(C2)<-c("Means"," ")
print(C2)
}
cat('------------------------------------------------------------------------\n')

#Para interacao nao significativa, fazer...
if(1-pf(Fcab,glab,glE)>sigF) {                            
cat('\nNot significant interaction: analyzing the simple effect
------------------------------------------------------------------------\n')
fatores<-data.frame('fator 1'=factor1,'fator 2' = factor2)

for(i in 1:2){

#Para os fatores QUALITATIVOS, teste de Tukey
if(quali[i]==TRUE && pv.fs[i]<=sigF) {
    cat(fac.names[i])
      if(mcomp=='tukey'){
    tukey(resp,fatores[,i],anavaT[6,1],anavaT[6,2],sigT)                                 
                    }
  if(mcomp=='duncan'){
    duncan(resp,fatores[,i],anavaT[6,1],anavaT[6,2],sigT)
                    }
  if(mcomp=='lsd'){
    lsd(resp,fatores[,i],anavaT[6,1],anavaT[6,2],sigT)
                    }
  if(mcomp=='lsdb'){
    lsdb(resp,fatores[,i],anavaT[6,1],anavaT[6,2],sigT)
                    }
  if(mcomp=='sk'){
    scottknott(resp,fatores[,i],anavaT[6,1],anavaT[6,2],sigT)
                    }
  if(mcomp=='snk'){
    snk(resp,fatores[,i],anavaT[6,1],anavaT[6,2],sigT)
                    }
  if(mcomp=='ccboot'){
    ccboot(resp,fatores[,i],anavaT[6,1],anavaT[6,2],sigT)
                    }
                   }
if(quali[i]==TRUE && pv.fs[i]>sigF) {
    cat(fac.names[i])
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && pv.fs[i]<=sigF){                              
    cat(fac.names[i])
    reg.poly(resp, fatores[,i], anavaT[6,1],anavaT[6,2], anavaT[i+1,1], anavaT[i+1,2])
}

if(quali[i]==FALSE && pv.fs[i]>sigF) {
    cat(fac.names[i])
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------')
                            }

cat('\n')
}

}

#Se a interacao for significativa, desdobrar a interacao
if(1-pf(Fcab,glab,glE)<=sigF){                                          
cat("\n\n\nSignificant interaction: analyzing the interaction
------------------------------------------------------------------------\n")

#Desdobramento de FATOR 1 dentro do niveis de FATOR 2
cat("\nAnalyzing ", fac.names[1], ' inside of each level of ', fac.names[2], '
------------------------------------------------------------------------\n')

des1<-aov(resp~ Block + Fator2/Fator1)

l1<-vector('list',nv2)
names(l1)<-names(summary(Fator2))
v<-numeric(0)
for(j in 1:nv2) {
        for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
        l1[[j]]<-v
        v<-numeric(0)
                }
des1.tab<-summary(des1,split=list('Fator2:Fator1'=l1))[[1]]

#Montando a tabela de ANAVA do des1
glf1=c(as.numeric(des1.tab[4:(nv2+3),1]))

SQf1=c(as.numeric(des1.tab[4:(nv2+3),2]))

QMf1=SQf1/glf1

Fcf1=QMf1/QME

rn<-numeric(0)
for(j in 1:nv2){ rn<-c(rn, paste(paste(fac.names[2],':',fac.names[1],sep=''),lf2[j]))}

anavad1<-data.frame("DF"=c(glB,glb, glf1, glad, glE, glT),
"SS"=c(round(c(SQB,SQb,SQf1,SQad,SQE,SQT),5)),
"MS"=c(round(c(QMB,QMb,QMf1,QMad,QME,QMT),5)),
"Fc"=c(round(c(FcB,Fcb,Fcf1,Fcad),4),'',''),
"Pr>Fc"=c(round(c(1-pf(FcB,glB,glE),1-pf(Fcb,glb,glE),1-pf(Fcf1,glf1,glE), 1-pf(Fcad,glad,glE)),4),' ', ' '))
colnames(anavad1)[5]="Pr>Fc"
rownames(anavad1)=c("Block",fac.names[2],rn,"Ad vs Factorial","Residuals","Total")
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(anavad1)
cat('------------------------------------------------------------------------\n\n')


ii<-0
for(i in 1:nv2) {
ii<-ii+1  
  if(1-pf(Fcf1,glf1,glE)[ii]<=sigF){                                       
    if(quali[1]==TRUE){
                      cat('\n\n',fac.names[1],' inside of the level ',lf2[i],' of ',fac.names[2],'
------------------------------------------------------------------------')
                        if(mcomp=='tukey'){
                          tukey(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],anavaT[6,1],anavaT[6,2],sigT)
                                          }                                  
                        if(mcomp=='duncan'){
                          duncan(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],anavaT[6,1],anavaT[6,2],sigT)            
                                           }                   
                        if(mcomp=='lsd'){
                          lsd(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],anavaT[6,1],anavaT[6,2],sigT)
                                        }
                        if(mcomp=='lsdb'){
                          lsdb(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],anavaT[6,1],anavaT[6,2],sigT)
                                         }
                        if(mcomp=='sk'){
                          scottknott(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],anavaT[6,1],anavaT[6,2],sigT)
                                       }
                        if(mcomp=='snk'){
                          snk(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],anavaT[6,1],anavaT[6,2],sigT)
                                        }
                        if(mcomp=='ccboot'){
                          ccboot(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],anavaT[6,1],anavaT[6,2],sigT)
                                        }
                      }
    else{  #regressao
    cat('\n\n',fac.names[1],' inside of the level ',lf2[i],' of ',fac.names[2],'
------------------------------------------------------------------------')
    reg.poly(resp[Fator2==lf2[i]], factor1[Fator2==lf2[i]], anavaT[6,1], anavaT[6,2], anavad1[i+2,1], anavad1[i+2,2])          
        }
                              }
    else{cat('\n\n',fac.names[1],' inside of the level ',lf2[i],' of ',fac.names[2],'\n')
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n')
        mean.table<-tapply.stat(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],mean)
        colnames(mean.table)<-c('  Levels','    Means')
        print(mean.table)
        cat('------------------------------------------------------------------------\n')
        }                          
                 }
cat('\n\n')

#Desdobramento de FATOR 2 dentro do niveis de FATOR 1
cat("\nAnalyzing ", fac.names[2], ' inside of each level of ', fac.names[1], '
------------------------------------------------------------------------\n')
    
des2<-aov(resp~ Block + Fator1/Fator2)

l2<-vector('list',nv1)
names(l2)<-names(summary(Fator1))
v<-numeric(0)
for(j in 1:nv1) {
        for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
        l2[[j]]<-v
        v<-numeric(0)
                }
des2.tab<-summary(des2,split=list('Fator1:Fator2'=l2))[[1]]
#Montando a tabela de ANAVA do des2
glf2=c(as.numeric(des2.tab[4:(nv1+3),1]))

SQf2=c(as.numeric(des2.tab[4:(nv1+3),2]))

QMf2=SQf2/glf2

Fcf2=QMf2/QME

rn<-numeric(0)
for(i in 1:nv1){ rn<-c(rn, paste(paste(fac.names[1],':',fac.names[2],sep=''),lf1[i]))}

anavad2<-data.frame("DF"=c(glB,gla, glf2, glad, glE, glT),
"SS"=c(round(c(SQB,SQa,SQf2,SQad,SQE,SQT),5)),
"MS"=c(round(c(QMB,QMa,QMf2,QMad,QME,QMT),5)),
"Fc"=c(round(c(FcB,Fca,Fcf2,Fcad),4),'',''),
"Pr>Fc"=c(round(c(1-pf(FcB,glB,glE),1-pf(Fca,gla,glE), 1-pf(Fcf2,glf2,glE), 1-pf(Fcad,glad,glE)),4),' ', ' '))
colnames(anavad2)[5]="Pr>Fc"
rownames(anavad2)=c("Block",fac.names[1],rn,"Ad vs Factorial","Residuals","Total")
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(anavad2)
cat('------------------------------------------------------------------------\n\n')


ii<-0
for(i in 1:nv1) {
  ii<-ii+1
  if(1-pf(Fcf2,glf2,glE)[ii]<=sigF){                               
    if(quali[2]==TRUE){
                      cat('\n\n',fac.names[2],' inside of the level ',lf1[i],' of ',fac.names[1],'
------------------------------------------------------------------------')
                      if(mcomp=='tukey'){
                          tukey(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],anavaT[6,1],anavaT[6,2],sigT)
                                          }
                        if(mcomp=='duncan'){
                          duncan(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],anavaT[6,1],anavaT[6,2],sigT)            
                                           }                   
                        if(mcomp=='lsd'){
                          lsd(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],anavaT[6,1],anavaT[6,2],sigT)
                                        }
                        if(mcomp=='lsdb'){
                          lsdb(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],anavaT[6,1],anavaT[6,2],sigT)
                                         }
                        if(mcomp=='sk'){
                          scottknott(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],anavaT[6,1],anavaT[6,2],sigT)
                                       }
                        if(mcomp=='snk'){
                          snk(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],anavaT[6,1],anavaT[6,2],sigT)
                                        }
                        if(mcomp=='ccboot'){
                          ccboot(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],anavaT[6,1],anavaT[6,2],sigT)
                                        }
                      }
    else{  #regressao
        cat('\n\n',fac.names[2],' inside of the level ',lf1[i],' of ',fac.names[1],'
------------------------------------------------------------------------')
        reg.poly(resp[Fator1==lf1[i]], factor2[Fator1==lf1[i]], anavaT[6,1],anavaT[6,2], anavad2[i+2,1], anavad2[i+2,2])
        }
                             }
    else{cat('\n\n',fac.names[2],' inside of the level ',lf1[i],' of ',fac.names[1],'\n')
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n')
        mean.table<-tapply.stat(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],mean)
        colnames(mean.table)<-c('  Levels','    Means')
        print(mean.table)
        cat('------------------------------------------------------------------------\n')
        }                          

                }
}




}
