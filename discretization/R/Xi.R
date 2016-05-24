Xi <-
function(data){
          data <- as.data.frame(data)
          p <- ncol(data)
          n <- nrow(data)
          Class <- dim(table(data[,p])) #class level
          data.ex <- data[,-p]#except class attribute
          if(p==2){data.ex=as.matrix(data.ex)} 
          mat <- unique(data.ex)   #matching instance
          if(dim(mat)[1]==dim(data)[1]) return(Xi=0)
          w <- which(duplicated(data.ex)==TRUE)#instance number of duplicated data 
                                          # - duplicated data is removed in mat by unique function.
          dup <- unique(data.ex[w,]) # duplicated instance
          if(is.numeric(dup)==TRUE){dup=as.matrix(dup)} 
          n.dup <- nrow(dup)
          num <- 0
             for(i in 1:n.dup){
                  same <- apply(as.numeric(dup[i,1:(p-1)])==t(data.ex),2,prod) #To find same instance for all attribute.
                                                                        #same==1 : It means same instances except class.
                  ww <- as.numeric(which(same==1)) # Same instance's instance number
                  ww.class <- as.numeric(data[ww,p])  # class of each same instance 

                  fac <- factor(ww.class) #class of each instance - for matching instances
                  same.cl <- table(fac) # number of class for same instance

                  firs <- ww[1] #ww[1] is representative matching instances
                  num <- as.numeric(ifelse(same==0,0,firs))+num  # dupliced data's instance number makes same.
                  fnum <- ifelse(num==0,rownames(data),num) #final numbering
                  fnum <- as.numeric(fnum)

                  freq <- table(fnum,data[,p]) # class frequency for matching instance
             }
            
             CED <- numeric()
             for(j in 1:dim(freq)[2]){
                  Ctable <- freq[freq[,j]!=0,]
                  if(class(Ctable)=="integer"){Ctable=as.matrix(Ctable,ncol=Class)}
                  Csum <- sum(Ctable)
                  rSum <- apply(Ctable,1,sum)
                  rMax <- apply(Ctable,1,max)
                  misClass <- sum(rSum-rMax)
                  CED[j] <- 1-(misClass/Csum)
             }
          m1 <- ifelse(all(0.5<=CED)==FALSE,0,1-min(CED[0.5<=CED]))
          m2 <- ifelse(all(CED<=0.5)==FALSE,0,max(CED[CED<=0.5]))
          Xi <- max(m1,m2)
          return(Xi)
          
}
