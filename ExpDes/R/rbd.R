rbd <-
function(treat, block, resp, quali=TRUE, mcomp='tukey', sigT=0.05, sigF=0.05) {


Trat<-factor(treat)
Bloco<-factor(block)
anava<-aov(resp~Trat+Bloco)
tab<-summary(anava)

colnames(tab[[1]])<-c('DF','SS','MS','Fc','Pr>Fc')
tab[[1]]<-rbind(tab[[1]],c(apply(tab[[1]],2,sum)))
rownames(tab[[1]])<-c('Treatament','Block','Residuals','Total')
cv<-round(sqrt(tab[[1]][3,3])/mean(resp)*100, 2)
tab[[1]][4,3]=' '
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(tab[[1]])
cat('------------------------------------------------------------------------\nCV =',cv,'%\n')


#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------\nShapiro-Wilk normality test\n')
cat('p-value: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('WARNING: at 5% of significance, residuals can not be considered normal!
------------------------------------------------------------------------\n')}
else{cat('According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.
------------------------------------------------------------------------\n')}

if(tab[[1]][1,5]<sigF){ 

if(quali==TRUE) {
  
  if(mcomp=='tukey'){
    tukey(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
                    }
  if(mcomp=='duncan'){
    duncan(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)            
                    }                   
  if(mcomp=='lsd'){
    lsd(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
                    }
  if(mcomp=='lsdb'){
    lsdb(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
                    }
  if(mcomp=='sk'){
    scottknott(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
                    }
  if(mcomp=='snk'){
    snk(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
                    }
                }                   
else{   
    reg.poly(resp, treat, tab[[1]][3,1], tab[[1]][3,2], tab[[1]][1,1], tab[[1]][1,2])
}
                       }
else {
    cat('\nAccording to the F test, the means can not be considered distinct.\n')
mean.table<-tapply.stat(resp,treat,mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}

}
