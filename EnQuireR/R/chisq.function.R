"chisq.function"=function(dataset,Y,X,method="proba",print=TRUE){
      res.chisq <- list()           # Creation of the R object (type list) aimed to store the results
      a=1
      b=0                          # a is an integer for iterations
    #varY=vector()
    #varX=vector()
    testfaits=vector()

      for (j in 1:length(Y)){
    #varY=c(varY,Y[j])
        for (i in 1:length(X)){
                if (Y[j]!=X[i]&!paste(Y[j],X[i],collapse="=")%in%testfaits){          # If the same variable is selected in both X and Y, the test crossing this variable with itself will not be performed
                    testfaits=c(testfaits,paste(Y[j],X[i],collapse="="),paste(X[i],Y[j],collapse="="))
                    #if(!X[i]%in%varX){
                    #varX=c(varX,X[i])
                    #}
          # Creation of the formula to be used by the xtabs function
                fmla=as.formula(paste("~",paste(c(Y[j],X[i]),collapse="+")))

                # Creation of the contingency table by the xtabs function
                tab=xtabs(formula=fmla,data=dataset)

                # Performs the Chi² test on the contingency table created above
                test=chisq.test(tab,correct=FALSE)

                # contrib is a matrix containing the contributions of each couple of levels to the Chi² statistic
                contrib=as.matrix(test$residuals^2)
            contrib1=as.data.frame(contrib[1:nrow(contrib),1:ncol(contrib)]) # Converting into a data frame
            contrib1=as.matrix(contrib1) #pour que coltable fonctionne dans R.8.1

            if (method=="mean"){
            # threshold is equal to the mean of contributions
                threshold=test$statistic/(nrow(contrib)*ncol(contrib))

          if (length(Y)*length(X)<11&&print==TRUE){
                # Using the function coltable developped for the SensoMineR package, we plot the contribution table
                # The contributions greater than 'threshold' are colored in blue. The ones that are smaller than 'threshold' are colored in pink
                coltable(contrib1, level.lower =threshold, level.upper = threshold, main.title = "")
                tpolice=par("cex") #Police titre
                title(paste(c(Y[j],X[i]),collapse=" X "),cex=tpolice)
                mtext(paste(c("Statistic  =",signif(test$statistic,4)),collapse=" "), side = 3, line=-0.1,cex=0.8,adj=0)
                mtext(paste(c("p.value  =",signif(test$p.value,digits=4)),collapse=" "), side = 3, line=-0.9,cex=0.8,adj=0)
          }
          else{ b=1
          }
          }
          #else{if (method=="proba"){
            nb.modalite <- length(levels(dataset[,Y[j]]))
            old.warn = options("warn")
                options(warn = -1)
                marge.li = xtabs(as.formula(paste(c("~",Y[j]),collapse="")),data=dataset)
                nom = tri = structure(vector(mode = "list", length = nb.modalite), names = levels(dataset[,Y[j]]))
                  Table <- xtabs(as.formula(paste(c("~",paste(c(Y[j],X[i]),collapse="+")),collapse="")),data=dataset)
                  marge.col = xtabs(as.formula(paste(c("~",X[i]),collapse="")),data=dataset)
                       for (l in 1:nlevels(dataset[,Y[j]])) {
                         for (k in 1:nlevels(dataset[,X[i]])) {
                        aux2 = Table[l,k]/marge.li[l]
                          aux3 = marge.col[k]/sum(marge.col)
                          if (aux2 > aux3) aux4 = phyper(Table[l,k]-1,marge.li[l],sum(marge.li)-marge.li[l],marge.col[k],lower.tail=FALSE)
                        else  aux4 = phyper(Table[l,k],marge.li[l],sum(marge.li)-marge.li[l],marge.col[k])
                                aux5 = (1-2*as.integer(aux2>aux3))*qnorm(aux4)
                                aux1 = Table[l,k]/marge.col[k]
                                tri[[l]] = rbind(tri[[l]],c(aux1,aux2,aux3,aux4,aux5))
                                nom[[l]] = rbind(nom[[l]],c(levels(dataset[,X[i]])[k],X[i]))

                      }
                }
           colormatr=matrix(data=NA,nrow=0,ncol=nlevels(dataset[,X[i]]))
           for (l in 1:nb.modalite){
           if (!is.null(tri[[l]])){
             oo = rev(order(tri[[l]][,5]))
             tri[[l]] = tri[[l]][oo,]
             nom[[l]] = nom[[l]][oo,]
             tri[[l]] = matrix(tri[[l]],ncol=5)
             rownames(tri[[l]]) = paste(nom[[l]][,1])
               colormatr=rbind(colormatr,t(tri[[l]][c(levels(dataset[,X[i]])),4]))
           rownames(colormatr)[l]=names(tri)[l]
           colnames(tri[[l]]) =  c("Cla/Mod","Mod/Cla","Global","p.value","V-test")
           }
           }
           colormatr=colormatr[c(levels(dataset[,Y[j]])),]
       if(method=="proba"){
           if (length(Y)*length(X)<11&&print==TRUE){
           coltable(contrib1,colormatr, level.lower =0.05, col.lower ="lightblue", col.upper = "mistyrose",level.upper = 0.05, main.title = "")
           tpolice=par("cex") #Police titre
           title(paste(c(Y[j],X[i]),collapse=" X "),cex=tpolice)
           mtext(paste(c("Statistic  =",signif(test$statistic,4)),collapse=" "), side = 3, line=-0.1,cex=0.8,adj=0)
           mtext(paste(c("p.value  =",signif(test$p.value,digits=4)),collapse=" "), side = 3, line=-0.9,cex=0.8,adj=0)
           }
           else{
            if (print==FALSE){
            }else{
            b=1
            }
           }
           }
        #}
          # Filling the res.chisq object with the results from this iteration
          res.chisq[[a]]=list(Test=test,Contrib=contrib,Proba.levels=colormatr)
            if (a==1) res.names=c(paste(c(Y[j],X[i]),collapse="."))
            else res.names=c(res.names,paste(c(Y[j],X[i]),collapse="."))
          a=a+1
          }
        }
    }
    names(res.chisq)=res.names
   if (b==1) print("Too many tables to be plotted, they will not be plotted but the results have been generated")
   #assign("varY",varY,envir=.GlobalEnv)
   #assign("varX",varX,envir=.GlobalEnv)
   return(res.chisq)
}
