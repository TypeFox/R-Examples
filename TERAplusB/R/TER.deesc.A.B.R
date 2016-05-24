TER.deesc.A.B <-
function(Dose,Prob.Dose,A,B,C,D,E)
{
  A1.f<-function()
  {  n<-nrow(keep.all)
     ID<<-ID+1
     keep.temp<-data.frame(keep.all)
     parent.ID<-as.numeric(as.character(keep.temp$ID[n]))
     parent.time<-as.numeric(as.character(keep.temp$current.time[n]))
     parent.dose.id<-as.numeric(as.character(keep.temp$current.dose.id[n]))
     parent.dose<-as.character(keep.temp$current.dose[n])
     parent.state<-keep.temp$current.state[n]
     current.time<-parent.time+1
     current.dose.id<-parent.dose.id+1
     current.dose<-Dose[current.dose.id]
     current.state<-"A1"
     if(current.dose.id<=K)
     {  current.n<-A
        MTD<-NA
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
        A1.f() 
        parent.ID<-as.numeric(as.character(keep.all$ID[n]))
        parent.time<-as.numeric(as.character(keep.all$current.time[n]))
        parent.dose.id<-as.numeric(as.character(keep.all$current.dose.id[n]))
        parent.dose<-as.character(keep.all$current.dose[n])
        parent.state<-keep.all$current.state[n]
        current.time<-parent.time+1
        current.dose.id<-parent.dose.id+1
        current.dose<-Dose[current.dose.id]
        MTD<-NA
        current.state<-"A2"
        current.n<-A
        ID<<-ID+1
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
        if(current.dose.id<=K)
          A2.f() 
        current.state<-"A3"
        current.n<-A
        ID<<-ID+1
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
        if(current.dose.id<=K)
          A3.f() 
     } else
     {  current.n<-0
        current.time<-parent.time+1
        current.dose.id<-parent.dose.id
        current.dose<-Dose[current.dose.id]
        MTD<-Dose[current.dose.id]
        current.state<-"Final"
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
     } 
  }

  search.dose<-function(parent.ID,current.dose)
  {  flag<-TRUE
     temp.ID<-as.character(parent.ID)
     n<-0
     while(flag)
     {  sel.id<-which(keep.all$ID==temp.ID)
        if(keep.all$current.dose[sel.id]==current.dose)
           n<-n+as.numeric(as.character(keep.all$current.n[sel.id]))
        temp.ID<-as.character(keep.all$parent.ID[sel.id])
        flag<-!is.na(temp.ID)    
     }
     return(n)
  }

  A3.f<-function()
  {  n<-nrow(keep.all)
     parent.ID<-as.numeric(as.character(keep.all$ID[n]))
     parent.time<-as.numeric(as.character(keep.all$current.time[n]))
     parent.dose.id<-as.numeric(as.character(keep.all$current.dose.id[n]))
     parent.dose<-as.character(keep.all$current.dose[n])
     parent.state<-keep.all$current.state[n]
     current.time<-parent.time+1
     current.dose.id<-parent.dose.id-1
     current.dose<-ifelse(current.dose.id==0,"D0",Dose[current.dose.id])
     if(current.dose==keep.all$parent.dose[1])
     {  current.state<-"Final"
        current.n<-0
        ID<<-ID+1      
        MTD<-as.character(keep.all$parent.dose[1])
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
     } else
     {  nn<-search.dose(parent.ID,current.dose)
        if(nn<A+B)
        {  current.state<-"B1"
           current.n<-B
           ID<<-ID+1         
           MTD<-NA
           keep.temp<-as.matrix(keep.all)
           keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id, 
                            parent.dose,parent.state,current.time,current.dose.id, 
                            current.dose,current.state,current.n,MTD))
           keep.all<<-data.frame(keep.temp)
           if(current.dose.id<K)
              A1.f() 
           current.state<-"B2"
           current.n<-B
           ID<<-ID+1         
           MTD<-NA
           keep.temp<-as.matrix(keep.all)
           keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                            parent.dose,parent.state,current.time,current.dose.id,
                            current.dose,current.state,current.n,MTD))
           keep.all<<-data.frame(keep.temp)
           current.state<-"Final"
           current.n<-0
           ID<<-ID+1         
           MTD<-ifelse(current.dose.id==1,"D0",Dose[current.dose.id-1])
           keep.temp<-as.matrix(keep.all)
           keep.temp<-rbind(keep.temp,c(ID,ID-1,parent.time,parent.dose.id,
                            parent.dose,parent.state,current.time+1,current.dose.id,
                            current.dose,current.state,current.n,MTD))
           keep.all<<-data.frame(keep.temp)
        } else
        {  current.state<-"Final"
           current.n<-0
           ID<<-ID+1
           MTD<-Dose[current.dose.id]
           keep.temp<-as.matrix(keep.all)
           keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                            parent.dose,parent.state,current.time,current.dose.id,
                            current.dose,current.state,current.n,MTD))
           keep.all<<-data.frame(keep.temp)
        }
     }   
  }

  A2.f<-function()
  {  n<-nrow(keep.all)
     parent.ID<-as.numeric(as.character(keep.all$ID[n]))
     parent.time<-as.numeric(as.character(keep.all$current.time[n]))
     parent.dose.id<-as.numeric(as.character(keep.all$current.dose.id[n]))
     parent.dose<-as.character(keep.all$current.dose[n])
     parent.state<-keep.all$current.state[n]
     current.time<-parent.time+1
     current.dose.id<-parent.dose.id
     current.dose<-Dose[current.dose.id]
     nn<-search.dose(parent.ID,current.dose)
     if(nn<A+B)
     {  current.state<-"B1"
        current.n<-B
        MTD<-NA
        ID<<-ID+1
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
        if(current.dose.id<=K)
           A1.f() 
        current.state<-"B2"
        ID<<-ID+1
        current.n<-B
        MTD<-NA
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
        current.state<-"Final"
        ID<<-ID+1
        current.n<-0
        MTD<-ifelse(current.dose.id==1,"D0",Dose[current.dose.id-1])
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,ID-1,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time+1,current.dose.id,
                         current.dose,current.state,current.n,MTD))
        keep.all<<-data.frame(keep.temp)
     } else
     {  current.state<-"Final"
        current.n<-0
        ID<<-ID+1
        MTD<-Dose[current.dose.id]
        keep.temp<-as.matrix(keep.all)
        keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                         parent.dose,parent.state,current.time,current.dose.id,
                         current.dose,current.state,current.n,MTD))

        keep.all<<-data.frame(keep.temp)
     }
  }

  B1.f<-function()
  {  n<-nrow(keep.all)
     ID<<-ID+1
     parent.ID<-as.numeric(as.character(keep.all$ID[n]))
     parent.time<-as.numeric(as.character(keep.all$current.time[n]))
     parent.dose.id<-as.numeric(as.character(keep.all$current.dose.id[n]))
     parent.dose<-as.character(keep.all$current.dose[n])
     parent.state<-keep.all$current.state[n]
     current.time<-parent.time+1
     current.dose.id<-parent.dose.id+1
     current.dose<-Dose[current.dose.id]
     current.state<-"B1"
     current.n<-B
     MTD<-NA
     keep.temp<-as.matrix(keep.all)
     keep.temp<-rbind(keep.temp,c(ID,parent.ID,parent.time,parent.dose.id,
                      parent.dose,parent.state,current.time,current.dose.id,
                      current.dose,current.state,current.n,MTD))
     keep.all<<-data.frame(keep.temp)
     if(current.dose.id<=K)
       A1.f() 
  }

  P.ij<-function(n,pi,j)
  {  choose(n,j)*pi^j *(1-pi)^(n-j)
  }  

  calc.Prob<-function(x)
  {  Pi<-Prob.Dose[as.numeric(as.character(x[8]))] 
     n<-as.numeric(as.character(x[11]))
     state<-x[10]
     if(state=="A1")
     {  Prob<-sum(P.ij(n,Pi,0:(C-1)))
     } else if(state=="A2")
     {  Prob<-sum(P.ij(n,Pi,C:(D-1)))
     } else if(state=="A3")
     {  Prob<-sum(P.ij(n,Pi,D:n))
     } else if(state=="B1")
     {  Prob<-sum(P.ij(n,Pi,0:(E-D)))#Prob<-sum(P.ij(n,Pi,0:(C-1)))#
     } else if(state=="B2")
     {  Prob<-sum(P.ij(n,Pi,(E-D+1):n))#Prob<-sum(P.ij(n,Pi,C:n))#
     } else if(state=="Final")
     { Prob<-1
     }
     return(Prob)
  }

  search.id<-function(tot.result,t.ID)
  {  id.keep<-id.temp<-as.numeric(as.character(tot.result$parent.ID[t.ID]))
     flag<-TRUE
     while(flag)
     {  if(tot.result$parent.ID[id.temp]!=0& !is.na(tot.result$parent.ID[id.temp]))
        {  id.temp<-as.numeric(as.character(tot.result$parent.ID[id.temp]))
           id.keep<-c(id.keep,id.temp)
        } else
        {  flag<-FALSE
        }
     }
     id.keep<-c(t.ID,id.keep)
     return(id.keep)
  }
    
  search.id.new<-function(tot.result,t.ID)
  {  id.keep<-id.temp<-t.ID
     flag<-TRUE
     while(flag)
     {  if(tot.result$parent.ID[id.temp]!=0& !is.na(tot.result$parent.ID[id.temp]))
        {  id.temp<-as.numeric(as.character(tot.result$parent.ID[id.temp]))
           id.keep<-c(id.keep,id.temp)
        } else
        {  flag<-FALSE
        }
     }
     return(id.keep)
  }
        
  Prob.search<-function(tot.result,MTD)
  {  MTD.id<-which(tot.result$MTD==MTD)
     prob.save<-NULL
     for(i in MTD.id)
     {  temp.id<-search.id(tot.result,i)
        temp.prob<-prod(tot.result$Prob.tot[temp.id])
        prob.save<-c(prob.save,temp.prob)
     }
     return(sum(prob.save))
  }
  
  Prob.search1<-function(tot.result,MTD)
  {  MTD.id<-which(tot.result$MTD==MTD)
     prob.save<-NULL
     for(i in MTD.id)
     {  temp.id<-search.id(tot.result,i)
        temp.prob<-prod(tot.result$Prob.tot[temp.id])
        prob.save<-c(prob.save,temp.prob)
     }
     return(prob.save)
  }  
  
  n.search<-function(tot.result,MTD)
  {  MTD.id<-which(tot.result$MTD==MTD)
     n.save<-NULL
     for(i in MTD.id)
     {  temp.id<-search.id(tot.result,i)
        temp.n<-sum(as.numeric(as.character(tot.result$current.n[temp.id])))
        n.save<-c(n.save,temp.n)
     }
     return(mean(n.save))
  }  
  
  n.search1<-function(tot.result,MTD)
  {  MTD.id<-which(tot.result$MTD==MTD)
     n.save<-NULL
     for(i in MTD.id)
     {  temp.id<-search.id(tot.result,i)
        temp.n<-sum(as.numeric(as.character(tot.result$current.n[temp.id])))
        n.save<-c(n.save,temp.n)
     }
     return(n.save) 
  }
   
  time.search<-function(tot.result,MTD)
  {  MTD.id<-which(tot.result$MTD==MTD)
     time.save<-NULL
     for(i in MTD.id)
     {  temp.id<-search.id(tot.result,i)
        temp.time<-as.numeric(as.character(tot.result$current.time[temp.id[1]]))
        time.save<-c(time.save,temp.time)
     }
     return(mean(time.save))
  }  
  
  time.search1<-function(tot.result,MTD)
  {  MTD.id<-which(tot.result$MTD==MTD)
     time.save<-NULL
     for(i in MTD.id)
     {  temp.id<-search.id(tot.result,i)
        temp.time<-as.numeric(as.character(tot.result$current.time[temp.id[1]]))
        time.save<-c(time.save,temp.time)
     }
     return(time.save)
  }    
  
  keep.all<-NULL
  parent.time<-0
  parent.ID<-NA
  parent.dose.id<-0
  parent.dose<-"D0"
  parent.state<-NA
  current.time<-0
  current.dose.id<-0
  current.dose<-"D0"
  current.state<-NA
  current.n<-0
  MTD<-NA
  ID<<-0
  K<-length(Dose)
  keep.all<-rbind(keep.all,c(ID,parent.ID,parent.time,parent.dose.id,
                  parent.dose,parent.state,current.time,current.dose.id,
                  current.dose,current.state,current.n,MTD))
  colnames(keep.all)<-c("ID","parent.ID","parent.time","parent.dose.id",
                        "parent.dose","parent.state","current.time",
                        "current.dose.id","current.dose","current.state",
                        "current.n","MTD")
  keep.all<<-data.frame(keep.all)
  A1.f()
  temp<-as.matrix(keep.all)
  Prob.tot<- c(1,apply(temp[-1,],1,calc.Prob))
  temp.all<-data.frame(keep.all,Prob.tot=Prob.tot)
  keep.all<<-temp.all
  tot.result<-temp.all[-1,]
 
  A.prob<-Prob.search(tot.result,c("D0"))
  A.n<-n.search(tot.result,c("D0"))  
  A.time<-time.search(tot.result,c("D0"))-1  
  A.prob1<-Prob.search1(tot.result,c("D0"))
  A.n1<-n.search1(tot.result,c("D0"))  
  A.time1<-time.search1(tot.result,c("D0"))-1  
  
  tot.list<-list()
  tot.list[[1]]<-cbind(A.n1,A.prob1,A.time1) 
   E.n<-0
  E.time<-0 
  for(i in 1:K)
  { A.prob<-c(A.prob, Prob.search(tot.result,Dose[i]))
    A.n<-c(A.n, n.search(tot.result,Dose[i]))
    A.time<-c(A.time, time.search(tot.result,Dose[i])-1)
    A.prob1<-Prob.search1(tot.result,Dose[i])
    A.n1<-n.search1(tot.result,Dose[i])
    A.time1<-time.search1(tot.result,Dose[i])-1
    E.n<-E.n+sum(A.prob1*A.n1)
    E.time<-E.time+sum(A.prob1*A.time1)    
    tot.list[[i+1]]<-cbind(A.n1,A.prob1,A.time1)    
  }  
  Prob.result<-data.frame(Dose=c("D0",Dose),Prob.Dose=c("NA",Prob.Dose),Prob=A.prob,Time=A.time,N=A.n)
  E.toxrate<-sum(Prob.Dose*A.prob[-1]/(1-A.prob[1]))
  tot.list<-temp.all
  Tot.result<-list(tot.list=tot.list,Prob.result=Prob.result,E.toxrate=E.toxrate,E.n=E.n,E.time=E.time)
  return(Tot.result)
}

