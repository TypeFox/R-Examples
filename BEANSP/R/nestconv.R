nestconv <-
function(type,nn,jj,data){
        y<-rep(0,nn)
        ul<-rep(0,nn)
        ur<-rep(0,nn)
        zl<-rep(0,nn)
        zr<-rep(0,nn)
        if(type==0){
            # interval format
            for(i in 1:nn){
                if(data[i,2]=="F"){y[i]=0}
                else {y[i]=1}
                ul[i]=data[i,3]
                ur[i]=data[i,4]
                zl[i]=data[i,5]
                zr[i]=data[i,6]
            }    
        }
        else{
            # date format
            for(i in 1:nn){
                if(data[i,5]=="F") {y[i]=0}
                else {y[i]=1}
                ### zl
                if(as.numeric(as.Date(data[i,2],"%m/%d/%Y"))==as.numeric(as.Date(data[i,4],"%m/%d/%Y")))
                {zl[i]=1}
                else
                {zl[i]=2+as.numeric(as.Date(data[i,3],"%m/%d/%Y"))-as.numeric(as.Date(data[i,2],"%m/%d/%Y"))}
                ### zr   
                zr[i]=1+as.numeric(as.Date(data[i,4],"%m/%d/%Y"))-as.numeric(as.Date(data[i,2],"%m/%d/%Y"))
                ### ul    
                if(y[i]==0) {ul[i]=1}
                else {ul[i]=jj-zr[i]+1}
                ### ur
                ur[i]=jj-zl[i]+1
            }
        }
        id<-seq(1,nn)
        for(i in 1:nn){
            ### update the bounds
            if(ur[i]>(jj-zl[i]+1)){
                ur[i]=(jj-zl[i]+1)
                if(ul[i]>(jj-zl[i]+1)){
                    ul[i]=jj-zl[i]+1
                }
            }
            if(y[i]==1){
                if(ul[i]<(jj-zr[i]+1)){
                    ul[i]=jj-zr[i]+1
                }
            }    
        }   
        temp<-cbind.data.frame(id,ul,ur,zl,zr,y)
        return(temp)
    }
