
############################################################
# FUNCTIONS OF THIS FILE:
# discrepSA_LHS,discrepESE_LHS: low-discrepancy LHS
# maximinSA_LHS,maximinESE_LHS: Maximin LHS

# Optimization algorithms:
# SA  = Simulated Annealing 
# ESE = Enhanced Stochastic Evolutionary 

#####REFERENCES#####

#Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs:
#optimization of Latin Hypercube Samples and subprojection properties, Journal of Simulation, 7:276-289, 2013.

#D.Morris and J.Mitchell. Exploratory designs for computationnal experiments. Journal of 
#Statistical Planning and Inference, 43:381-402, 1995.

#R.JIN, W.CHEN, A.SUDJIANTO. An efficient algorithm for constructing optimal design
#of computer experiments. Journal of Statistical Planning and Inference, 134:268-287, 2005

###############################################################


#####discrepSA_LHS#####
#####L2 DISCREPANCY LHS VIA SIMULATED ANNEALING OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  m     : the design                                                 |
#        T0    : the initial temperature                                    |
#        c     : paramater regulating the temperature                       |
#        it    : the number of iterations                                   |
#    criterion : criterion to be optimized ("C2","W2" or "L2star")          |
#      profile : temperature down profile                                   |
#                "GEOM" or "GEOM_MORRIS" or "LINEAR". By default : "GEOM"   |
#       Imax   :parameter adjusting number of iterations without improvement|
#               (if you choose "GEOM_MORRIS" profile)                       |
#output        : a list containing all the input arguments plus:            |
#       low L2_discrepancy design                                           |
#       vector of criterion values along the iterations                     |
#       vector of temperature values along the iterations                   |
#       vector of acceptation probability values along the iterations       |
#depends :  discrepancyCriteria, discrepancyW2_EP, lhs_EP, discrepancyL2_EP |
#           discrepancyC2_EP, lhs_EP                                        |
#---------------------------------------------------------------------------|
   

discrepSA_LHS<-function(design,T0=10,c=0.95,it=2000,criterion="C2",profile="GEOM",Imax=100) 

{ 
  m <- design
  crit <- NULL ; temp <- NULL ; proba <- NULL
  
  if(criterion=="C2")
    
  {if(profile=="GEOM"){
    
    i<-0
    T<-T0
    v<-discrepancyCriteria(m,type='C2')[[1]]  
    crit <- v
    while(T>0 & i<it){
      
      G<-lhs_EP(m)
      g<-discrepancyC2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
      g<-sqrt(g)
      
      diff<-min(exp((v-g)/T),1)
      if (diff == 1){m<-G[[1]]
                     v<-g  
      }
      else         {
        Bernoulli<-rbinom(1,1,diff)
        if (Bernoulli==1){m=G[[1]]
                          v<-g    }}
      i<-i+1
      crit <- c(crit,v) ; temp <- c(temp,T) ; proba <- c(proba,diff)
      T<-(c^i)*(T0)
    }
    
    List.res <- list(design,T0,c,it,criterion,profile,Imax,m,crit,temp,proba)
    names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
  }
  
  if(profile=="LINEAR"){
  
  i<-0
  T<-T0
  v<-discrepancyCriteria(m,type='C2')[[1]]  
  crit <- v
  while(T>0 & i<it){
    
    G<-lhs_EP(m)
    g<-discrepancyC2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
    g<-sqrt(g)
    
    diff<-min(exp((v-g)/T),1)
    if (diff == 1){m<-G[[1]]
                   v<-g  
    }
    else         {
      Bernoulli<-rbinom(1,1,diff)
      if (Bernoulli==1){m=G[[1]]
                        v<-g    }}
    i<-i+1
    crit <- c(crit,v) ; temp <- c(temp,T) ; proba <- c(proba,diff)
    T<-T0*(1-i/it)
  }
  
  List.res <- list(design,T0,c,it,criterion,profile,Imax,m,crit,temp,proba)
  names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
}

if(profile=="MC"){
  
  i<-0
  T<-T0
  v<-discrepancyCriteria(m,type='C2')[[1]]  
  crit <- v
  
  n <- dim(design)[1]
  nbvar <- dim(design)[2]
  
  while(T>0 & i<it){
    
#    G<-lhs_EP(m)
    j <- sample(1:n,1)
    new <- runif(nbvar)
    G <- m
    G[j,] <- new
    
#    g<-discrepancyC2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
    g <- discrepancyCriteria(G,type='C2')[[1]] 
    
#    g<-sqrt(g)
    
    diff<-min(exp((v-g)/T),1)
    if (diff == 1){m<-G
                   v<-g  
    }
    else         {
      Bernoulli<-rbinom(1,1,diff)
      if (Bernoulli==1){m=G
                        v<-g    }}
    i<-i+1
    crit <- c(crit,v) ; temp <- c(temp,T) ; proba <- c(proba,diff)
    T<-(c^i)*(T0)
  }
  
  List.res <- list(design,T0,c,it,criterion,profile,Imax,m,crit,temp,proba)
  names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
}
 
    if(profile=="GEOM_MORRIS"){

       
       T<-T0
       Dbest<-m
       v<-discrepancyCriteria(m,type='C2')[[1]]
       ref<-v
       crit <- v
       for (i in 1:it)
         {
          flag<-0
          I<-1
          while(I<Imax){
                     
             G<-lhs_EP(m)
             g<-discrepancyC2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
             g=sqrt(g)
             diff<-min(exp((v-g)/T),1)
             if (diff == 1){m<-G[[1]]
                         v<-g   
                         flag=1
                       
                       if(v<ref){ref=v
                                  Dbest=m
                                  I=1}
                          else {I=I+1}
                         }
                        
             else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G[[1]]
                                    v=g 
                                    flag=1  
                         if (v<ref){ref=v
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
             }
             crit <- c(crit,ref) ; temp <- c(temp,T) ; proba <- c(proba,diff)
         } 
          
       if (flag==1) {T=T*c}
               else {break}
        } 

        List.res <- list(design,T0,c,it,criterion,profile,Imax,Dbest,crit,temp,proba)
        names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
    }
   } # end if(criterion=="C2")


   
   if(criterion=="L2star")
    {if(profile=="GEOM")
          {i<-0
           T<-T0
           v<-discrepancyCriteria(m,type='L2star')[[1]]	
           crit <- v
           while(T>0 & i<it){
           
                    G<-lhs_EP(m)
                    g<-discrepancyL2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
                    g<-sqrt(g)
                 
                    diff<-min(exp((v-g)/T),1)
                        if (diff == 1){m<-G[[1]]
                                  v<-g    
                                 }
                    else         {
                                 Bernoulli<-rbinom(1,1,diff)
                                 if (Bernoulli==1){m=G[[1]]
                                    v<-g    }}
                  i<-i+1
                  crit <- c(crit,v) ; temp <- c(temp,T) ; proba <- c(proba,diff)
                  T<-(c^i)*(T0)
           }
           List.res <- list(design,T0,c,it,criterion,profile,Imax,m,crit,temp,proba)
           names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
    }   
    if(profile=="LINEAR")
{i<-0
 T<-T0
 v<-discrepancyCriteria(m,type='L2star')[[1]]  
 crit <- v
 while(T>0 & i<it){
   
   G<-lhs_EP(m)
   g<-discrepancyL2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
   g<-sqrt(g)
   
   diff<-min(exp((v-g)/T),1)
   if (diff == 1){m<-G[[1]]
                  v<-g    
   }
   else         {
     Bernoulli<-rbinom(1,1,diff)
     if (Bernoulli==1){m=G[[1]]
                       v<-g    }}
   i<-i+1
   crit <- c(crit,v) ; temp <- c(temp,T) ; proba <- c(proba,diff)
   T<-T0*(1-i/it)
 }
 List.res <- list(design,T0,c,it,criterion,profile,Imax,m,crit,temp,proba)
 names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
}   
      if(profile=="GEOM_MORRIS")

       {
        T<-T0
        Dbest<-m
        v<-discrepancyCriteria(m,type='L2star')[[1]]  
        ref<-v
        crit <- ref
        for (i in 1:it)
         {
          flag<-0
          I<-1
          while(I<Imax){
                     
             G<-lhs_EP(m)
             g<-discrepancyL2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
             g=sqrt(g)
             diff<-min(exp((v-g)/T),1)
             if (diff == 1){m<-G[[1]]
                         v<-g   
                         flag=1
                       
                       if(v<ref){ref=v
                                  Dbest=m
                                  I=1}
                          else {I=I+1}
                         }
                        
             else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G[[1]]
                                    v=g 
                                    flag=1  
                         if (v<ref){ref=v
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
                 }
             crit <- c(crit,ref) ; temp <- c(temp,T) ; proba <- c(proba,diff)
         } 
          
       if (flag==1) {T=T*c}
               else {break}
        } 
        
        List.res <- list(design,T0,c,it,criterion,profile,Imax,Dbest,crit,temp,proba)
        names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
      }
   } # end if(criterion=="L2star")
  
   if(criterion=="W2")
    {if(profile=="GEOM")
          {i<-0
           T<-T0
           v<-discrepancyCriteria(m,type='W2')[[1]]	
           crit <- v
           while(T>0 & i<it){
           
                G<-lhs_EP(m)
                g<-discrepancyW2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
                g<-sqrt(g)
                 
                diff<-min(exp((v-g)/T),1)
                      if (diff == 1){m<-G[[1]]
                                  v<-g   
                                 }
                    else         {
                                 Bernoulli<-rbinom(1,1,diff)
                                 if (Bernoulli==1){m=G[[1]]
                                    v<-g    }}
                     i<-i+1
                crit <- c(crit,v) ; temp <- c(temp,T) ; proba <- c(proba,diff)
                T<-(c^i)*(T0)
                    }
           List.res <- list(design,T0,c,it,criterion,profile,Imax,m,crit,temp,proba)
           names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
           
           }
     if(profile=="LINEAR")
     {i<-0
      T<-T0
      v<-discrepancyCriteria(m,type='W2')[[1]]  
      crit <- v
      while(T>0 & i<it){
        
        G<-lhs_EP(m)
        g<-discrepancyW2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
        g<-sqrt(g)
        
        diff<-min(exp((v-g)/T),1)
        if (diff == 1){m<-G[[1]]
                       v<-g   
        }
        else         {
          Bernoulli<-rbinom(1,1,diff)
          if (Bernoulli==1){m=G[[1]]
                            v<-g    }}
        i<-i+1
        crit <- c(crit,v) ; temp <- c(temp,T) ; proba <- c(proba,diff)
        T<-T0*(1-i/it)
      }
      List.res <- list(design,T0,c,it,criterion,profile,Imax,m,crit,temp,proba)
      names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
     }
     
       if(profile=="GEOM_MORRIS")

       {
        T<-T0
        Dbest<-m
        v<-discrepancyCriteria(m,type='W2')[[1]]
        ref<-v
        crit <- ref
        for (i in 1:it)
         {
          flag<-0
          I<-1
          while(I<Imax){
                     
             G<-lhs_EP(m)
             g<-discrepancyW2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
             g=sqrt(g)
             diff<-min(exp((v-g)/T),1)
             if (diff == 1){m<-G[[1]]
                         v<-g   
                         flag=1
                       
                       if(v<ref){ref=v
                                  Dbest=m
                                  I=1}
                          else {I=I+1}
                        }
                        
             else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G[[1]]
                                    v=g 
                                    flag=1  
                         if (v<ref){ref=v
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
                 }
             crit <- c(crit,ref) ; temp <- c(temp,T) ; proba <- c(proba,diff)
         } 
          
       if (flag==1) {T=T*c}
               else {break}
        } 
        
        List.res <- list(design,T0,c,it,criterion,profile,Imax,Dbest,crit,temp,proba)
        names(List.res) <-c("InitialDesign","TO","c","it","criterion","profile","Imax","design","critValues","tempValues","probaValues")
       }
     } # end if(criterion=="W2")
  
     return(List.res)
} # end function

#####Function to extract elements from lists composed with lists#####
extract_list<-function(l){return(l[[1]])}



#####discrepESE_LHS#####
#####L2 DISCREPANCY LHS VIA ESE OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  design     : the design                                            |
#        T0    : the initial temperature                                    |
#        inner_it  : number of iterations for inner loop                    |
#        J     : number of new proposed LHS in inner loop                   |
#        it    : number of iterations for outer loop                        |
#        criterion: "C2", "W2" or "L2star"                               |
#output        : a list containing all the input arguments plus:            |
#       low L2_discrepancy design                                           |
#       vector of criterion values along the iterations                     |
#       vector of temperature values along the iterations                   |
#       vector of acceptation probability values along the iterations       |
#depends :  discrepancyL2_EP_ESE, discrepancyW2_EP_ESE, discrepancyL2_STAR  |
#           discrepancyC2_EP_ESE, discrepancyCriteria                       |
#---------------------------------------------------------------------------|


discrepESE_LHS<-function(design,T0=0.005*discrepancyCriteria(design,type='C2')[[1]],inner_it=100,J=50,it=2,criterion="C2")
{
  m<-design
  crit <- NULL ; temp <- NULL ; proba <- NULL
  
  if(criterion=="C2")
  {
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  dC2<-discrepancyCriteria(m,type='C2')[[1]]                   
  best<-dC2
  crit <- dC2

  for (q in 1:it)
   {
     
     BOLD<-Best
     bold<-best                                   # Best=new LHS built at every step        
                                                # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) #liste de liste
     
     g<-lapply(l,discrepancyC2_EP_ESE,k=modulo+1,p=discrepancyCriteria(m,type='C2')[[1]]^2)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]  
                                             
       Delta<-a-dC2
                                                
     if((Delta)<=(Temperature*runif(1)))                 # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                                # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           dC2<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
       crit <- c(crit,best)
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio
  
  temp <- c(temp,rep(Temperature,inner_it)) ; proba <- c(proba,rep(v1,inner_it))   

  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }

                                               # else, it is the exploratory step

    }
  List.res <- list(design,T0,inner_it,J,it,criterion,Best,crit,temp,proba)
  names(List.res) <-c("InitialDesign","TO","inner_it","J","it","criterion","design","critValues","tempValues","probaValues") 
  
  }


  if(criterion=="L2star")
  {
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  dL2<-discrepancyCriteria(m,type='L2star')[[1]]              
  best<-dL2              

  for (q in 1:it)
   {

     BOLD<-Best 
     bold<-best                                   # Best=new LHS built at every step        
                                                # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) #liste de liste
     
     g<-lapply(l,discrepancyL2_EP_ESE,k=modulo+1,p=dL2^2)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]
                          
                                               
       Delta<-a-dL2
                                                
     if((Delta)<=(Temperature*runif(1)))                 # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                                # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           dL2<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
       crit <- c(crit,best)
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio

  temp <- c(temp,rep(Temperature,inner_it)) ; proba <- c(proba,rep(v1,inner_it))

  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }

                                               # else, it is the exploratory step



    }
  List.res <- list(design,T0,inner_it,J,it,criterion,Best,crit,temp,proba)
  names(List.res) <-c("InitialDesign","TO","inner_it","J","it","criterion","design","critValues","tempValues","probaValues") 
  
  }

 if(criterion=="W2")
  {
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  dW2<-discrepancyCriteria(m,type='W2')[[1]]                  
  best<-dW2              

  for (q in 1:it)
   {

     BOLD<-Best
     bold<-best                                    # Best=new LHS built at every step        
                                                # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) #liste de liste
     
     g<-lapply(l,discrepancyW2_EP_ESE,k=modulo+1,p=discrepancyCriteria(m,type='W2')[[1]]^2)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]
                          
                                               
       Delta<-a-dW2
                                                
     if((Delta)<=(Temperature*runif(1)))                 # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                               # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           dW2<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
       crit <- c(crit,best)
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio

  temp <- c(temp,rep(Temperature,inner_it)) ; proba <- c(proba,rep(v1,inner_it))
     
  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }

                                               # else, it is the exploratory step

    }
  List.res <- list(design,T0,inner_it,J,it,criterion,Best,crit,temp,proba)
  names(List.res) <-c("InitialDesign","TO","inner_it","J","it","criterion","design","critValues","tempValues","probaValues") 
  
 }

return(List.res)

}


#####maximinSA_LHS#####
#####Maximin LHS VIA SIMULATED ANNEALING OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  m     : the design                                                 |
#        T0    : the initial temperature                                    |
#        c     : parameter regulating the temperature                       |
#        it    : the number of iterations                                   |
#        p     : power required in phiP criterion                           |
#      profile : temperature down profile                                   |
#                "GEOM" or "GEOM_MORRIS" or "LINEAR". By default : "GEOM"   |
#output        : a list containing all the input arguments plus:            |
#       a mindist optimized design                                          |
#       vector of criterion values along the iterations                     |
#       vector of temperature values along the iterations                   |
#       vector of acceptation probability values along the iterations       |
#depends :  phiP,lhs_EP                                                     |
#---------------------------------------------------------------------------|
   

maximinSA_LHS<-function(design,T0=10,c=0.95,it=2000,p=50,profile="GEOM",Imax=100) 
  { 
  crit <- NULL ; temp <- NULL ; proba <- NULL
  
  if(profile=="GEOM"){
  m<-design
  i<-0
  T<-T0
  fi_p<-phiP(m,p)
  crit <- fi_p
  
  while(T>0 & i<it)

         {           
               G<-lhs_EP(m)[[1]]
               fi_p_ep<-phiP(G,p)
               diff<-min(exp((fi_p-fi_p_ep)/T),1)
               if (diff == 1){m<-G
                       fi_p<-fi_p_ep
                      }
               else { Bernoulli<-rbinom(1,1,diff)
                      if (Bernoulli==1){m<-G
                                    fi_p<-fi_p_ep    }}

               i<-i+1
               crit <- c(crit,fi_p) ; temp <- c(temp,T) ; proba <- c(proba,diff)
               T<-(c^i)*(T0)
           }
  
  List.res <- list(design,T0,c,it,p,profile,Imax,m,crit,temp,proba)
  names(List.res) <-c("InitialDesign","TO","c","it","p","profile","Imax","design","critValues","tempValues","probaValues") 
  }
  
if(profile=="LINEAR"){
    m<-design
    i<-0
    T<-T0
    fi_p<-phiP(m,p)
    crit <- fi_p
    
    while(T>0 & i<it)
      
    {           
      G<-lhs_EP(m)[[1]]
      fi_p_ep<-phiP(G,p)
      diff<-min(exp((fi_p-fi_p_ep)/T),1)
      if (diff == 1){m<-G
                     fi_p<-fi_p_ep
      }
      else { Bernoulli<-rbinom(1,1,diff)
             if (Bernoulli==1){m<-G
                               fi_p<-fi_p_ep    }}
      
      i<-i+1
      crit <- c(crit,fi_p) ; temp <- c(temp,T) ; proba <- c(proba,diff)
      T<-T0*(1-i/it)
    }
    
    List.res <- list(design,T0,c,it,p,profile,Imax,m,crit,temp,proba)
    names(List.res) <-c("InitialDesign","TO","c","it","p","profile","Imax","design","critValues","tempValues","probaValues") 
  }
  
  if(profile=="GEOM_MORRIS") 
     
  {m<-design
   T<-T0
   Dbest<-m
   fi_p=phiP(m,p)
   ref<-fi_p
   crit <- ref
   for (i in 1:it)
   {
     flag=0
     I=1
     while(I<Imax){
                     
           G=lhs_EP(m)[[1]]
           fi_p_ep=phiP(G,p)
           diff=min(exp((fi_p-fi_p_ep)/T),1)
           if (diff == 1){m=G
                       fi_p=fi_p_ep    
                       flag=1
                       
                       if(fi_p<ref){ref=fi_p
                                  Dbest=m
                                  I=1}
                          else {I=I+1}
                       }
                        
           else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G
                                    fi_p=fi_p_ep 
                                    flag=1  
                       if (fi_p<ref){ref=fi_p
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
                 }
           crit <- c(crit,ref) ; temp <- c(temp,T) ; proba <- c(proba,diff)
         } 
          
       if (flag==1) {T=T*c}
               else {break}
   } 
   
   List.res <- list(design,T0,c,it,p,profile,Imax,Dbest,crit,temp,proba)
   names(List.res) <-c("InitialDesign","TO","c","it","p","profile","Imax","design","critValues","tempValues","probaValues") 
   
   }
  
  
  if(profile=="MC"){
    m<-design
    i<-0
    T<-T0
    fi_p<-phiP(m,p)
    crit <- fi_p
    n <- dim(design)[1]
    nbvar <- dim(design)[2]
    
    while(T>0 & i<it)
      
    {           
      
#      G<-lhs_EP(m)[[1]]
      j <- sample(1:n,1)
      new <- runif(nbvar)
      G <- m
      G[j,] <- new
      fi_p_ep<-phiP(G,p)
      diff<-min(exp((fi_p-fi_p_ep)/T),1)
      if (diff == 1){m<-G
                     fi_p<-fi_p_ep
      }
      else { Bernoulli<-rbinom(1,1,diff)
             if (Bernoulli==1){m<-G
                               fi_p<-fi_p_ep    }}
      
      i<-i+1
      crit <- c(crit,fi_p) ; temp <- c(temp,T) ; proba <- c(proba,diff)
      T<-(c^i)*(T0)
    }
    
    List.res <- list(design,T0,c,it,p,profile,Imax,m,crit,temp,proba)
    names(List.res) <-c("InitialDesign","TO","c","it","p","profile","Imax","design","critValues","tempValues","probaValues") 
  }
  
  
   return(List.res)
}

#####maximinESE_LHS#####
#####Maximin LHS VIA ESE OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  design: the design                                                 |
#        T0    : the initial temperature                                    |
#        inner_it  : number of iterations for inner loop                    |
#        J     : number of new proposed LHS in inner loop                   |
#        it    : number of iterations for outer loop                        |
#        p     : the power in phiP criterion                                |
#output        : a list containing all the input arguments plus:            |
#       a mindist optimized design                                          |
#       vector of criterion values along the iterations                     |
#       vector of temperature values along the iterations                   |
#       vector of acceptation probability values along the iterations       |
#depends :     phiP_EP_ESE, phiP                                            |
#---------------------------------------------------------------------------|

maximinESE_LHS<-function(design,T0=0.005*phiP(design,p=50),inner_it=100,J=50,it=1,p=50)
{
  m<-design
  crit <- NULL ; temp <- NULL ; proba <- NULL
  
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  fi_p<-phiP(m,p)                   
  best<-fi_p
  crit <- fi_p

  for (q in 1:it)
   {
     BOLD<-Best
     bold<-best                              # B=new LHS built at every step        
                                             # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) 
     
     g<-lapply(l,phiP_EP_ESE,k=modulo+1,p=p)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]
                                               
       Delta<-a-fi_p
                                                
     if((Delta)<=(Temperature*runif(1)))        # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                                # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           fi_p<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
       crit <- c(crit,best)
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio
     
  temp <- c(temp,rep(Temperature,inner_it)) ; proba <- c(proba,rep(v1,inner_it))

  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }
}
  
List.res <- list(design,T0,inner_it,J,it,p,Best,crit,temp,proba)
names(List.res) <-c("InitialDesign","TO","inner_it","J","it","p","design","critValues","tempValues","probaValues") 
     # else, it is the exploratory step

return(List.res)
}
