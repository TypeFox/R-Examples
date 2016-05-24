tables <-
function(data){
d=data
        s=as.list(d)
        t1=lapply(s, table)
        t2=lapply(t1, prop.table)
        f1=function(i){round(t2[[i]]*100,2)}
        i=1:length(t2)
        t3=lapply(i, f1); names(t3)=names(t2)
        t4=lapply(t1, chisq.test)
        d2=d[,-1];ii=1:length(d2)
        f2=function(ii){table(d[,1],d2[,ii])}
        t5=lapply(ii,f2)
        n1=names(d); n1=n1[1]; n2=names(d[,-1])
        f3=function(ii){paste(n1,"vs",n2[ii])}
        n3=lapply(ii, f3); names(t5)=n3
        t6=lapply(t5, prop.table, 1)
        t2=t6
        t7=lapply(ii, f1); names(t7)=names(t6)
        t8=lapply(t5, chisq.test)
        t9=lapply(t5, fisher.test)
        
        l=list(t1,t3,t4,t5,t7,t8,t9)
        names(l)=c("Absolute frequency for all variables","Percentage frequency for all variables", "Chi-square test for all variables", "Contingency tables (absolute frequency)", "Contingency tables (percentage frequency)", "Chi-square test for contingency tables", "Fisher's exact test for contingency tables")
        return(l)
    }
