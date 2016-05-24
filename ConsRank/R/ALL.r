#          library('gregmisc')
#          library('MASS')
#

#FAST BRANCH AND BOUND ALGORITHM BY AMODIO, D'AMBROSIO AND SICILIANO (2015)
#Find multiple solutions for the consensus ranking problem.


#CR=FASTcons(X,Wk)
#X is a ranking data matrix (IMPORTANT: IT MUST BE A MATRIX N(JUDGES) X M(OBJECTS)
#Wk is a vector of weights
#CR$cons = Consensus Ranking
#CR$Tau = TauX rank correlation coefficient
#the code works under the dependencies of:
#
#          library('gregmisc')
#          library('MASS')
#          library('proxy')
#
#

#This code recall the following ocdes:
#BBconsensus.r      : principal code
#findconsensus.r    : first approximate solution
#kemenyd.r          : kemeny distance
#PenaltyBB2.r       : searching for penalties
#kemenydesign.r     : design matrix to compute kemeny distance
#combimpmatr.r      : combined input matrix
#scorematrix.r      : score matrix
#ReorderingBB.r     : particular reordering necessity




FASTcons = function(X, Wk=NULL, maxiter=50, FULL=FALSE, PS=FALSE)   {

  if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
  }

M = nrow(X)
N=ncol(X)

tic = proc.time()[3]
if (M==1) {
  CR = X
  Taux = 1
  } else {
  if (!is.null(Wk)) {
    
    if (is.numeric(Wk)) {
      Wk=matrix(Wk,ncol=1)
    }
    
    cij = combinpmatr(X,Wk)
    } else {
    cij = combinpmatr(X)
  }

    if ( sum(cij==0) == nrow(cij)^2 ) {
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 

  CR=matrix(0,maxiter,ncol(X))
  for (iter in 1:maxiter) {

    if (iter%%2==0) {
        R=matrix(sample(1:ncol(X),replace=T),1,ncol(X))
       } else {
       R=matrix(sample(1:ncol(X)),1,ncol(X))
       }

  consensus1 = BBconsensus(R,cij, FULL)
  cons=matrix(consensus1$cons,1,ncol(X))
  consensus = BBconsensus(cons,cij, FULL)
  #print(cons)
  #print(R)
  #flush.console()
  CR[iter,]=matrix(consensus$cons,1,ncol(X))
  if (PS==TRUE) {
    
    dsp1=paste("Iteration",iter,sep=" ")
    print(dsp1)
  }

  }
}
#d=kemenyd(X,consensus$cons)

Taux=matrix(0,nrow(CR),1)
for (k in 1:nrow(CR)) {
  Sij=scorematrix(matrix(CR[k,],1,ncol(X)))
  if (!is.null(Wk)){
      Taux[k,]=sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
     } else {
      Taux[k,]=sum(cij*Sij) / (  M*(N*(N-1)) )
      }
     }

CR=reordering(CR)
indice=which(Taux==max(Taux));
Taux=max(Taux)
CR=matrix(CR[indice,],ncol=N)
if (nrow(CR>1)){
  CR=unique(CR)
}
if (!is.null(dim(CR))) {
    Taux=matrix(rep(Taux,nrow(CR)))
    }

 colnames(CR)=colnames(X) 
 toc = proc.time()[3]
 eltime=toc-tic
 return(list(Consensus=CR, Tau=Taux, Eltime=eltime) )
}

#-------------------------------------------------------------------------------

#QUICK BRANCH AND BOUND ALGORITHM BY AMODIO, D'AMBROSIO AND SICILIANO
#Find a solution (or a solution really close to the final one). This
#program allows to find up to 4 multiple consensus rankings.

#CR=QuickCons(X,Wk)
#X is a ranking data matrix (IMPORTANT: IT MUST BE A MATRIX N(JUDGES) X M(OBJECTS)
#Wk is a vector of weights
#CR$cons = Consensus Ranking
#CR$Tau = TauX rank correlation coefficient


#This code recall the following ocdes:
#BBconsensus.r      : principal code
#findconsensus.r    : first approximate solution
#kemenyd.r          : kemeny distance
#PenaltyBB2.r       : searching for penalties
#kemenydesign.r     : design matrix to compute kemeny distance
#combimpmatr.r      : combined input matrix
#scorematrix.r      : score matrix
#ReorderingBB.r     : particular reordering necessity

#the code works under the dependencies of:
#
#          library('gregmisc')
#          library('MASS')
#          library('proxy')
#
#

QuickCons = function(X,Wk=NULL, FULL=FALSE,PS=FALSE)   {


if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
}

M = nrow(X)
N=ncol(X)
tic = proc.time()[3]

if (M==1) {
  consensus = X
  TauX = 1
  } else {
  if (!is.null(Wk)) {
    
    if (is.numeric(Wk)) {
      Wk=matrix(Wk,ncol=1)
    }
    
    cij = combinpmatr(X,Wk)
    } else {
    cij = combinpmatr(X)
    }
    
    if (sum(cij==0)==nrow(cij)^2){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 

  R=findconsensusBB(cij)
  R1=(N+1)-R
  consensusA = BBconsensus(R,cij,FULL,PS)$cons
  consensusB = BBconsensus(consensusA,cij,FULL,PS)$cons
  consensusC = BBconsensus(R1,cij,FULL,PS)$cons
  consensusD = BBconsensus(consensusC,cij,FULL,PS)$cons
  consensus = unique(reordering(rbind(consensusA,consensusB,consensusC,consensusD)))
  howcons = nrow(consensus)


  }
#d=kemenyd(X,consensus$cons)

Taux=matrix(0,nrow(consensus),1)
for (k in 1:nrow(consensus)) {

  #Sij=scorematrix(t(as.matrix(consensus[k,])))
  Sij=scorematrix(matrix(consensus[k,],1,N))

  if (!is.null(Wk)){
    Taux[k,1]=sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
    Taux[k,1]=sum(cij*Sij) / (  M*(N*(N-1)) )
    }

  }

if (howcons>1) {
  nco=which(Taux==max(Taux))
  if (length(nco)>1) {
    consensus=consensus[nco,]
    Taux=matrix(rep(max(Taux),nrow(consensus)),nrow(consensus),1)
    } else {
    Taux=max(Taux)
    #consensus = t(matrix(consensus[nco,]))
    consensus = matrix(consensus[nco,],1,N)
    }
   }
colnames(consensus)=colnames(X) 
 toc = proc.time()[3]
 eltime=toc-tic
 return(list(Consensus=reordering(consensus), Tau=Taux, Eltime=eltime) )
}

#-------------------------------------------------------------------------------


EMCons = function(X,Wk=NULL,PS=TRUE)  {
#Emond and Mason Branch and Bound algorithm to find median ranking
#X is a data matrix in which the rows are the judges and the columns indicates the objects
#Wk is the vector of weigths
if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
}



M = nrow(X)
N=ncol(X)
tic = proc.time()[3]
if (M==1) {
  consensus = X
  TauX = 1
  } else {
  if (!is.null(Wk)) {
    
    if (is.numeric(Wk)) {
      Wk=matrix(Wk,ncol=1)
    }
    
    cij = combinpmatr(X,Wk)
    } else {
    cij = combinpmatr(X)
    }
    
    if (sum(cij==0)==nrow(cij)^2){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
  R=findconsensusBB(cij)
  cons1=BBconsensus(R,cij,PS)
  consensus1=cons1$cons
  Po=cons1$pen
  consensus=BBconsensus2(consensus1,cij,Po,PS)
  }


if (nrow(consensus)==1) {

    Sij=scorematrix(consensus)

        if (!is.null(Wk)){
        TauX=sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
        } else {
        TauX=sum(cij*Sij) / (  M*(N*(N-1)) )
        }

} else {

    TauX=matrix(0,nrow(consensus),1)

    for (k in 1:nrow(consensus)) {

        Sij=scorematrix(t(matrix(consensus[k,])))

        if (!is.null(Wk)) {

            TauX[k,1] = sum(cij*Sij) / ( sum(Wk)*(N*(N-1)) )

        } else {

            TauX[k,1] = sum(cij*Sij) / (M*(N*(N-1)))

        }

    }

}
toc = proc.time()[3]
colnames(consensus)=colnames(X) 
#consensus=reordering(consensus)
eltime=toc-tic
return(list(Consensus=reordering(consensus), Tau=TauX, Eltime=eltime) )
}





#-------------------------------------------------------------------------------

### COMBINED INPUT MATRIX as defined by Emond and Mason


combinpmatr = function (X,Wk=NULL) {
if (is.null(Wk)) {
#X must be data matrix with n judges (on the rows) ranking m objects (on the columns)
  CI=matrix(0,ncol(X), ncol(X))
  for (i in 1:nrow(X)){
    sm=scorematrix(t(as.matrix(X[i,])))
    CI=CI+sm
  }
  } else {
    
    if (is.numeric(Wk)) {
      Wk=matrix(Wk,ncol=1)
    }
    
  CI=matrix(0,ncol(X), ncol(X))
  for (i in 1:nrow(X)){
    sm=scorematrix(t(as.matrix(X[i,])))*Wk[i]
    CI=CI+sm
  }
  }
CI
}

#-------------------------------------------------------------------------------

BBconsensus = function(RR,cij,FULL=FALSE,PS=FALSE) {

#Branch and Bound Algorithm to find the the consensus ranking PART I As modified by D'AMBROSIO (2008).
#Find the first approximation to the consensus ranking. Most of the time CR
#is a solution, maybe not unique
#Input:
#       RR -> First solution candidate to be the consensus ranking
#       cij -> Combined Input Matrix of the M individuals which judge n
#       objects
#Output:
#       CR -> Consensus Ranking
#
#
#References: Amodio et al.,2015; D'Ambrosio et al., 2015.


CR=RR;
sij=scorematrix(RR);
Po = sum(abs(cij))-sum(cij*sij)
a = t(matrix(sort(RR,decreasing = TRUE)))
ord = t(matrix(order(RR,decreasing = TRUE)))
R=RR;
addpenalty=matrix(0,length(a),1)

# exploration of the initial solution
for (k in 2:length(a)) {
  #print(k)
  b = 1:k
  R = ReorderingBB(R)
  KR=t(matrix(R[ord[b]]))
  KR=KR[-length(KR)]
  MO=max(KR)
  MI=min(KR)
  aa=1
  KO=1
  KR[length(KR)+1]=MO+1;
  R[ord[b]]=KR
  candidate=matrix(0,nrow(RR), ncol(RR))
  Pb = matrix(0, 1, 1)
  while (KO==1)  {
    #browser()
    candidate=rbind(candidate,R)
    #if (ncol(candidate>ncol(RR))) {

    #print(dim(candidate))
    #print(dim(R))
    #print(class(candidate))
    #print(class(R))

   # }
    if (aa==1){
      candidate=matrix(candidate[-1,],1,ncol(candidate))
      }

    Sij=scorematrix(matrix(candidate[aa,],1,ncol(R)))
  # print(Sij)
   #print(candidate)
   #flush.console()
    Pb=rbind(Pb,sum(abs(cij))-sum(cij*Sij))
    if (aa==1) {
      Pb=matrix(Pb[-1,],1,1)
    }
   # print(Pb)
    if (Pb[aa]==0) {

      CR = R
      Po = 0
      Pc = 0
      break
      }
      Pc=1
      if(FULL==TRUE){
        R[ord[b[length(b)]]] = R[ord[b[length(b)]]]-2 }else{
        R[ord[b[length(b)]]] = R[ord[b[length(b)]]]-1
        }
    if (MI-R[ord[b[length(b)]]] > 1) {
      KO = 0
      }
    aa=aa+1
    
  }
  
  if (PS==TRUE) {
    
    dsp2=paste("evaluated",nrow(candidate),"branches",sep=" ")
    print(dsp2)
    
  }

    if (Pc == 0) {
      break
      }

   minp=min(Pb)
   posp=which(Pb==min(Pb))

   if (minp<=Po) {
    Po=minp
    CR=t(matrix(candidate[posp[1],]))
    R=CR
    addpenalty[k,1]=PenaltyBB2(cij,R,ord[b])
    } else {
    R = CR
    addpenalty[k,1]=PenaltyBB2(cij,R,ord[b])
    }

    candidate = mat.or.vec(nrow(R), ncol(R))
    Pb = mat.or.vec(1, 1)

    }

if (Pc==0) {
  Po=0
  addpenalty = 0
  }
  else {
  Poo=sum(addpenalty)
  }

  SIJ = scorematrix(CR)
  Po=sum(addpenalty)

return(list(cons=CR,pen=Po))
}

#-------------------------------------------------------------------------------

## Place labels in a data matrix X of rankings (N judges by M objects)
#m is the number of objects
#label (optional) is the vector of the objects names
#labs = 1 or 2

labels = function(x, m, label = 1:m, labs ){
#source('reordering.r')
# if the class of the object is different from 'matrix' transform it in 'matrix'
  if(class(x) != 'matrix'){
  obs = length(x)
  XX = matrix(x, ncol = obs)
  } else {
  XX = x
  }

nj = nrow(XX)
nob = ncol(XX)

## if length of the object is higher than m, last number is the penalty
#if(nob > m){
## if the number of rows is 1 is a vector
#  if(nj == 1){
#    pens = x[m+1]
#    X = matrix(reordering(XX[1:m]), m, ncol = m)
#    } else {
#    pens = x[,m+1]
#    X = t(apply(x, 1, function(g) reordering(g, m)))
#    }
#} else {
X = XX
#}
if(labs ==1){
let = label
} else if(labs == 2){
let = LETTERS[label]
}

out = rep(0, nj)
for(i in 1:nj){

ord = rank(X[i,])
orders = tapply(let, ord, sort)

names1 = NULL
for(j in 1:length(orders)){
if(length(orders[[j]]) > 1){
nams = paste('(', paste(orders[[j]], sep = '', collapse = ''), ')', sep = '', collapse='')
} else {
nams = paste(orders[[j]], collapse = '')
}
names1 = c(names1, nams)
}
out[i] = paste(names1, collapse = ' ' )
}
out = matrix(out, nrow = nj)

#if(nob > m){
#dat = data.frame(data = out, pens = pens)
#} else {
dat = out
#}

return(dat)
}

#-------------------------------------------------------------------------------

### SCORE MATRIX OF RANK DATA ACCORDING EMOND AND MASON


scorematrix = function (X) {
  
if (is.numeric(X) & !is.matrix(X)){
  X=matrix(X,ncol=length(X))
}
  
#X must be a row vector containing a ranking of m objects
sm=matrix(0,ncol(X), ncol(X))
index=combinations(ncol(X),2)
for (i in 1:nrow(index)){
  if (!is.na(X[index[i,1]]) && !is.na(X[index[i,2]]) ) {
    if (X[index[i,1]]>X[index[i,2]]) {    #if object i is prerferred to object j
        sm[index[i,1],index[i,2]] = -1
        sm[index[i,2],index[i,1]] = 1
        } else if (X[index[i,1]]<X[index[i,2]]) {  #if object i is not prerferred to object j
        sm[index[i,1],index[i,2]] = 1
        sm[index[i,2],index[i,1]] = -1
        } else if (X[index[i,1]] == X[index[i,2]]) {  #if object i is in a tie with object j
        sm[index[i,1],index[i,2]] = 1
        sm[index[i,2],index[i,1]] = 1
     }
     }
     else {
     sm[index[i,1],index[i,2]] = 0
     sm[index[i,2],index[i,1]] = 0
     }
}
sm
}

#-------------------------------------------------------------------------------

tabulaterows = function(X,miss=FALSE) {
      
#given a sample of preference rankings, it counts the judges that have equal preferences
#and it tabulates the row of the data matrix
            
if (sum(is.na(X))>0) {
  miss=TRUE
  X[is.na(X)]=-10
}
            
coun = table(apply(X, 1, paste, collapse=","))
nam = names(coun)
spl = (strsplit(nam, ","))
kkn  = lapply(spl, as.numeric)
tab = t(as.data.frame(kkn))
cek =  cbind(tab,coun)
coun=as.matrix(coun)
rownames(coun)=NULL
rownames(tab)=NULL
if (miss==TRUE) {
  tab[tab==-10]=NA
}
      
          
return(list(X=tab, Wk=coun, tabfreq=cbind(tab,coun)))
}
          
          


#-------------------------------------------------------------------------------

branches = function(brR,cij,b,Po,ord,Pb,FULL=FALSE) {

candidate = findbranches(brR,ord,b,FULL)
Pb =matrix( rep(Pb,nrow(candidate)))

CR=mat.or.vec(nrow(candidate),ncol(candidate))
addpenalty=matrix(0,nrow(candidate),1)
QR=mat.or.vec(nrow(candidate),ncol(candidate))

for (gm in 1:nrow(candidate)) {

    CR[gm,]=candidate[gm,]
    addpenalty[gm,]=PenaltyBB2(cij,candidate[gm,],ord[b])

    if (Pb[gm]+addpenalty[gm] > Po) {

        CR[gm,]=-10.0e+15
        addpenalty[gm]=-10.0e+15

        }
    QR[gm,]=CR[gm,]
    }
Pbr=addpenalty+Pb
idp=Pbr<0

if (sum(idp)==0) {

    R=QR

    } else if (sum(idp==F)==nrow(QR)) {

    Pbr=NULL
    Pb=NULL
    R=NULL

    } else {
    Pbr=t(matrix(Pbr[idp==FALSE,],1))
    if (sum(idp==F)==1) {
        R=t(matrix(QR[idp==FALSE,]))
        } else {
        R=QR[idp==FALSE,]
        }
    }

return(list(cR=R,pcR=Pbr))
}

#-------------------------------------------------------------------------------

# given a ranking multiplied by 2, it returns a "normal" ranking
reordering = function (X) {

if (nrow(X)==1) {
    G=X
    OX = order(X)
    SX=matrix(sort(X),nrow=1)
    SX=SX-min(SX)+1
    DC=rbind(0,diff(t(SX)))


    for (i in 1:(ncol(X)-1)) {
        if (DC[i+1,1] >= 1) {
          SX[1,i+1]=SX[1,i]+1
          } else if (DC[i+1,1] == 0) {
          SX[1,i+1]=SX[1,i]
          }
      }

    G[1,OX]=SX

    } else {

    G=X
    for (j in 1:nrow(X)) {

        OX = order(X[j,])
        SX=matrix(sort(X[j,]),nrow=1)
        SX=SX-min(SX)+1
        DC=rbind(0,diff(t(SX)))

        for (i in 1:(ncol(X)-1)) {

            if (DC[i+1,1] >= 1) {
            SX[1,i+1]=SX[1,i]+1
            } else if (DC[i+1,1] == 0) {
            SX[1,i+1]=SX[1,i]
            }
        }
    G[j,OX]=SX
  }

 }

G
}

#-------------------------------------------------------------------------------

findbranches = function(R,ord,b,FULL=FALSE) {

  KR=t(matrix(R[ord[b]]))
  KR=KR[-length(KR)]
  MO=max(KR)
  MI=min(KR)
  aa=1
  KO=1
  KR[length(KR)+1]=MO+1;
  R[ord[b]]=KR
  candidate=mat.or.vec(nrow(R), ncol(R))

  while (KO==1)  {
      candidate=rbind(candidate,R)

      if (aa==1){
          candidate=matrix(candidate[-1,],1,ncol(candidate))
      }
      if (FULL==FALSE){
        R[ord[b[length(b)]]]=R[ord[b[length(b)]]]-1 }else{
          R[ord[b[length(b)]]]=R[ord[b[length(b)]]]-2
        }

      if (MI-R[ord[b[length(b)]]] > 1) {

          KO=0

      }

      aa=aa+1

  }

  Rt=candidate

}

#-------------------------------------------------------------------------------

## DETERMINATION OF PENALTIES FOR THE BRANCH AND BOUND ALGORITHM

PenaltyBB2 = function(cij,candidate,ord)   #indice must be order(CR)
{

Ds=t(mat.or.vec(1,(length(ord)-1)));
addpenalty=t(mat.or.vec(1,(length(ord)-1)));

for (k in 1:(length(ord)-1)) {

  Ds[k,1]=sign(candidate[ord[length(ord)]]-candidate[ord[k]]);

  if (Ds[k,1]==1) {


    if ( sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
       addpenalty[k,1]=cij[ord[length(ord)],ord[k]]-cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 & sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 & sign(cij[ord[k],ord[length(ord)]]) == 1) {
      addpenalty[k,1]=cij[ord[length(ord)],ord[k]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 & sign(cij[ord[k],ord[length(ord)]]) == 1) {
      addpenalty[k,1]=0
      }
    }
  else if (Ds[k,1]==-1) {
      if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
       addpenalty[k,1]=0
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
      addpenalty[k,1]=cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
      addpenalty[k,1]=cij[ord[k],ord[length(ord)]]-cij[ord[length(ord)],ord[k]]
      }
  }

 else if (Ds[k,1]==0) {
       if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
       addpenalty[k,1]=-cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
      addpenalty[k,1]=0
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
      addpenalty[k,1]=-cij[ord[length(ord)],ord[k]]
      }
 }
}

addpenalty=sum(addpenalty)

}

#-------------------------------------------------------------------------------

findconsensusBB = function(cij) {

X=mat.or.vec(1,ncol(cij)) + 1
N=ncol(X)
indici=combinations(N,2)
for (j in 1:nrow(indici)) {
  if ( sign(cij[indici[j,1],indici[j,2]]) == 1 & sign(cij[indici[j,2],indici[j,1]]) == -1 ) {
    X[indici[j,1]]=X[indici[j,1]]+1
    } else if ( sign(cij[indici[j,1],indici[j,2]]) == -1 & sign(cij[indici[j,2],indici[j,1]]) == 1 ) {
    X[indici[j,2]]=X[indici[j,2]] + 1
    } else if (sign(cij[indici[j,1],indici[j,2]]) == -1 & sign(cij[indici[j,2],indici[j,1]]) == -1 ) {
    X[indici[j,1]]= NA
    } else if (sign(cij[indici[j,1],indici[j,2]]) == 1 & sign(cij[indici[j,2],indici[j,1]]) == 1 ){
    X[indici[j,1]]=X[indici[j,1]]+1
    X[indici[j,2]]=X[indici[j,2]] + 1
    }
      
  }

X=(N+1)-X;
return(X)
}

#-------------------------------------------------------------------------------

#Compute the design matrix to compute Kemeny distance
kemenydesign = function(X) {
  
  if (is.numeric(X) & !is.matrix(X)) {
    X=matrix(X,ncol=length(X))
  }
  
  
  
M = ncol(X)
N = nrow(X)
indice=combinations(M,2)
KX=mat.or.vec(N,(M*(M-1)/2) )
for (j in 1:nrow(indice)) {
  KX[,j]=sign(X[,indice[j,1]] - X[,indice[j,2]])*-1
  }
 KX
}

#-------------------------------------------------------------------------------

##Kemeny Distance
kemenyd = function(X,Y=NULL) {
  
  if (is.numeric(X) & !is.matrix(X)) {
    X=matrix(X,ncol=length(X))
  }
  
if (is.null(Y)) {
  X = kemenydesign(X)
  d=dist(X,"manhattan")
  } else {
    
    if (is.numeric(Y) & !is.matrix(Y)) {
      Y=matrix(Y,ncol=length(Y))
    }

    
  X=kemenydesign(X)
  Y=kemenydesign(Y)
  d=dist(X,Y,"manhattan")
 }
d
}

#-------------------------------------------------------------------------------

ReorderingBB = function(RR) {

RR = RR+1
R = RR;
k = ncol(R)
neword = order(R)
indexing = mat.or.vec(1, ncol(R)-1)
for (j in (k-1):1) {
  indexing[j] = R[neword[j+1]]-R[neword[j]]
  }

if (sum(indexing==0)>0) {
  J = 1
  while (J<=ncol(indexing)) {
    if (indexing[J]==0) {
      R[neword[J+1]]=R[neword[J]]
      J=J+1
      }
      else if (indexing[J]>0) {
      R[neword[J+1]] = R[neword[J]]+2
      J=J+1
      }
    }
   }
  else  {
    J = 1
    while (J<= ncol(indexing)) {
      R[neword[J+1]] = R[neword[J]] + 2
      J=J+1
      }
    }
    R
  }
  
  #-------------------------------------------------------------------------------

  ## DETERMINATION OF PENALTIES FOR THE BRANCH AND BOUND ALGORITHM

Penalty = function(CR,cij,indice)   #indice must be order(CR)
{
  if (CR[indice[1,1]] < CR[indice[1,2]]) { #case 1, the first object is preferred
                                          #to the second object
     if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
     Po=0
     } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
     Po = cij[indice[1,2],indice[1,1]]-cij[indice[1,1],indice[1,2]]  #     cji-cij
     } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
     Po = cij[indice[1,2],indice[1,1]] #cji
     }
  } else if (CR[indice[1,1]] > CR[indice[1,2]]) { #case 2 the first object is not
                                                      #preferred to the second one
     if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
     Po= cij[indice[1,1],indice[1,2]]-cij[indice[1,2],indice[1,1]]  #cij-cji
     } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
     Po = 0
     } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
     Po = cij[indice[1,1],indice[1,2]] #cij
     }
  } else if (CR[indice[1,1]] == CR[indice[1,2]]) { #case 3 they are in a tie
     if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
     Po = -cij[indice[1,2],indice[1,1]] #-cj
     } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
     Po = -cij[indice[1,1],indice[1,2]] #-cij
     } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1 |
     sign(cij[indice[1,1],indice[1,2]]) == 0 & sign (cij[indice[1,2],indice[1,1]]) == 0 )  {
     Po = 0
     }

  }
  Po
 }
          
#---------------------------------------------------------------------------------
##Core code for the computation of the consensus ranking. Branch-and-bound 
##algorithm by Emond and Mason
          
BBconsensus2 = function(RR,cij,Po,PS=TRUE,FULL=FALSE) {
#Core function in computing consensus ranking
CR=RR;
a = t(matrix(sort(RR)))
ord = t(matrix(order(RR)))
r=ReorderingBB(RR)
BR.R=r #initialize ranking
BR.P=0 #initialize penalty
WCi=1
lambda=1
nobj=ncol(RR)
while (WCi == 1){
  if (PS==TRUE) {
      dsp1=paste("round",lambda,sep=" ")
      print(dsp1) 
    }
              
 for (k in 2:ncol(a)) { #primary loop: add the k^th better object ranked
                
    B = nrow(BR.R)    #branches size
                #print(B)
                #flush.console()
                
    b=1:k
          
                
    for (nb in 1:B) { #Secondary loop:  check the branches created by "nb"
                   
      BR.R[nb,] = ReorderingBB(t(matrix(BR.R[nb,])))
      rpbr=branches(matrix(BR.R[nb,],1,),cij,b,Po,ord,matrix(BR.P[nb]),FULL)
      R=rpbr$cR
      Pbr=rpbr$pcR        
              
                  
      if (is.null(R)) {
                  
                    #if (nb==1) {
                    #
                    #    JR=0
                    #    
                    #    next
                    #    
                    #} else {     
                    
        next
                    
                    #}
        } else {         #process good rankings
                    
        if (nb==1) {     #starting point
                      
                      #JR=nrow(R)
                      
            KR.R=R
            KR.P=Pbr
              
                      
                      
            } else {   #it is not the starting point
                      
            KR.R=rbind(KR.R,R)
            KR.P=rbind(KR.P,Pbr)
                  
                      
        }
                    
     }
                  
                  #JR = nrow(KR.R)  #update size of branches
                  
    } #end secondary loop
                
    if (is.null(R)) {
                  
      if (nb==B & nb!=1) { #If at the last iteration of the secondary loop all of thye rankings are not good 
                    
          rm(BR.R) #BR.R = NULL    #rm(BR.R)
          rm(BR.P) #BR.P = NULL    #rm(BR.P)
          BR.R = KR.R
          BR.P = KR.P
          KR.R = NULL
          KR.P = NULL
          #rm(JR) #JR = NULL
                    
          } else {
                    
          next
                    
      }
                  
    } else {
                  
    rm(BR.R) #BR.R = NULL  #rm(BR.R)
    rm(BR.P) #BR.P = NULL  #rm(BR.P)
    BR.R = KR.R
    BR.P = matrix(KR.P)
    #KR.R = NULL
    #KR.P = NULL
    #rm(JR)  #JR = NULL
                
    }
    
    if (PS==TRUE) {
      
      dsp2=paste("evaluating",B,"branches",sep=" ")
      print(dsp2)
      
    }
                
 } #end primary loop
              
              #AccPen = BR.P
              
              
 SSP=matrix(which(BR.P==min(BR.P)))
 MinP=min(BR.P)
 PenMin=Po-MinP
              
 if (PenMin==0) {
                
    #CR=t(matrix(BR.R[SSP,]))
   CR=matrix(BR.R[SSP,],length(SSP),nobj)
   WCi = 0
                
   } else {
                
   Po=MinP
   WCi=1
   lambda=lambda+1
   #nRR=t(matrix((BR.R[SSP[1],])))
   nRR=matrix((BR.R[SSP[1],]),1,nobj)
   rm(BR.R)
   rm(BR.P)
   BR.R=nRR
   BR.P=0
   a = t(matrix(sort(BR.R)))
   ord = t(matrix(order(BR.R)))
   rm(nRR)
                
  }
              
 }   #end while
            
            CR
}

#-------------------------------------------------------------
Tau_X = function(X,Y=NULL) {
  
  if(is.numeric(X) & !is.matrix(X)){
    X=matrix(X,ncol=length(X))
  }
  
  
  N=ncol(X)
  maxd=N*(N-1)
  if (is.null(Y)){
    d=kemenyd(X)
    tau = 1-(2*d/maxd)
  } else {
    if(is.numeric(Y) & !is.matrix(Y)){
      Y=matrix(Y,ncol=length(Y))
    }
    d=kemenyd(X,Y)
    tau = 1-(2*d/maxd)
  }
 tau 
}
  
#------------------------------------------------------------------------
polyplot = function(X=NULL,L=NULL,Wk=NULL,nobj=3){

  
  if (nobj==3){
    
    #rankings in the polytope
    ranks=rbind(
      c(1,2,3),
      c(1,2,2),
      c(1,3,2),
      c(1,2,1),
      c(2,3,1),
      c(2,2,1),
      c(3,2,1),
      c(2,1,1),
      c(3,1,2),
      c(2,1,2),
      c(2,1,3),
      c(1,1,2),
      c(1,1,1)
    )
    
    if (is.null(L)){
      rr=labels(ranks,3,labs=2)
    } else {
      rr=labels(ranks,3,L,labs=1)
    }
    
    #coordinates of polytope
    coord=rbind(
      c(  0.000000e+00,  4.082483e-01), #ABC
      c(  1.767767e-01,  3.061862e-01), #A(BC)
      c(  3.535534e-01,  2.041241e-01), #ACB
      c(  3.535534e-01,  1.817407e-16), #(AC)B
      c(  3.535534e-01, -2.041241e-01), #CAB
      c(  1.767767e-01, -3.061862e-01), #C(AB)
      c(  5.192593e-17, -4.082483e-01), #CBA
      c( -1.767767e-01, -3.061862e-01), #(BC)A
      c( -3.535534e-01, -2.041241e-01), #BCA
      c( -3.535534e-01, -1.038519e-16), #B(AC)
      c( -3.535534e-01,  2.041241e-01), #BAC
      c( -1.767767e-01,  3.061862e-01), #(AB)C
      c(  5.517130e-17,  4.543519e-17)  #(ABC)
    )
    
    
    plot(coord,ylim=c(-0.5,0.5),xlim=c(-0.5,0.5),axes=FALSE,ann=FALSE)
    lines(coord[1:11,])
    lines(c(coord[11,1],coord[1,1]),c(coord[11,2],coord[1,2]))
    lines(c(coord[2,1],coord[8,1]),c(coord[2,2],coord[8,2]),lty=2)
    lines(c(coord[4,1],coord[10,1]),c(coord[4,2],coord[10,2]),lty=2)
    lines(c(coord[6,1],coord[12,1]),c(coord[6,2],coord[12,2]),lty=2)
    
    t1=c(1,2,3,11,12,13)
    t2=c(4,10)
    t3=c(5,6,7,8,9)
    tcoord=coord #text coordinates
    tcoord[t1,2]=tcoord[t1,2]+0.1
    tcoord[t3,2]=tcoord[t3,2]-0.1
    tcoord[t2[1],1]=tcoord[t2[1],1]+0.1
    tcoord[t2[2],1]=tcoord[t2[2],1]-0.1
    #text(tcoord,rr)
    
    if (is.null(X)){
      indplot=matrix(1:13,ncol=1)
    } else {
      
      #o = outer(seq_len(nrow(X)), seq_len(nrow(ranks)), Vectorize(
      #  function(i, j) all(X2[i,]==ranks[j,])
      #))
      #ranksinplot=ranks[apply(o, 2, any),]
      
      o2 = outer(seq_len(nrow(ranks)), seq_len(nrow(X)), Vectorize(
        function(i, j) which(all(ranks[i,]==X[j,]))
      ))
      
      indexing=o2==1
      #print(indexing)
      indexing[indexing==TRUE]=1
      indexing[is.na(indexing)]=0
      indplot=which(rowSums(indexing)==1)
      #print(indplot)
    }
    
    # if (is.null(X)){
    #   X=ranks
    # }
    #X
    #put labels to rankings
    if (is.null(Wk)){
      
      points(coord[indplot,1],coord[indplot,2],pch=16,cex=0.8,col="blue")
      
    }else{
      
      if (is.numeric(Wk)) {
        
        Wk=matrix(Wk,ncol=1)
        
      }
      
      idwk=matrix(0,nrow(X),1)
      counter=0
      for (i in 1:13){
        for (j in 1:nrow(X)){
          check=sum(pos=X[j,]==ranks[i,])
          if (check==3){
            counter=counter+1
            idwk[counter]=j
            break}
        }
      }
      
      
      
      
      points(coord[indplot,1],coord[indplot,2],pch=16,cex=sqrt(100*((Wk[idwk]/sum(Wk))/pi)/2),col="blue")
      
    }

      
    text(tcoord[indplot,],rr[indplot,])
    
  }else{ ##4 objects
    
    
    #---------------------------------------------
    
    #-Exagon A first
    E1=rbind(
      c(0.5,  0.5,  1.4142135),  #   ...    %%'A B C D'     [1 2 3 4]   1
      c(1.0,  1.0,  0.70710677), #  ...    %%'A C B D'     [1 3 2 4]   2
      c(1.5,  0.5,  0.0), #                 %%'A C D B'     [1 4 2 3]   3
      c(1.5, -0.5,  0.0), #         ...    %%'A D C B'     [1 4 3 2]   4
      c(1.0, -1.0,  0.70710677), #  ...    %%'A D B C'     [1 3 4 2]   5
      c(0.5, -0.5,  1.4142135) #   ...    %%'A B D C'     [1 2 4 3]   6
    )
    
    
    MA=apply(E1,2,mean) #  center of exagon A   A(BCD)  [1 2 2 2]  7
    E1T=rbind(
      c(1.2500,    0.7500,    0.3536),#;...  %% 'A C {BD}'   [1 3 2 3]   8
      c(0.7500,   -0.7500,    1.0607),#;...  %% 'A {BD} C'   [1 2 3 2]   9
      c(0.5000,         0,    1.4142),#;...  %% 'A B {CD}'   [1 2 3 3]   10
      c(1.5000,         0,         0),#;...  %% 'A {CD} B'   [1 3 2 2]   11
      c(1.2500,   -0.7500,    0.3536),#;...  %% 'A D {BC}'   [1 3 3 2]   12
      c(0.7500,    0.7500,    1.0607)# ;...  %% 'A {BC} D'   [1 2 2 3]   13
    )
    
    
    
    ranksA=rbind(
      c(1,2,3,4),c(1,3,2,4),c(1,4,2,3),c(1,4,3,2),c(1,3,4,2),c(1,2,4,3),c(1,2,2,2),
      c(1,3,2,3),c(1,2,3,2),c(1,2,3,3),c(1,3,2,2),c(1,3,3,2),c(1,2,2,3)
    )
    
    #---------------------------------------------
    #-Exagon B first
    E2=rbind(
      c(-1.5, -0.5,  0.0),#        ...    %%'B D C A'     [4 1 3 2]   14
      c(-1.5,  0.5,  0.0),#        ...    %%'B C D A'     [4 1 2 3]   15 
      c(-1.0,  1.0,  0.70710677),# ...    %%'B C A D'     [3 1 2 4]   16 
      c(-0.5,  0.5,  1.4142135),#  ...    %%'B A C D'     [2 1 3 4]   17
      c(-0.5, -0.5,  1.4142135),#  ...    %%'B A D C'     [2 1 4 3]   18
      c(-1.0, -1.0,  0.70710677)# ...    %%'B D A C'     [3 1 4 2]   19
    )
    
    MB=apply(E2,2,mean) #center of exagon B  B(ACD)        [2 1 2 2]  20
    
    E2T=rbind(
      c(-0.5000,         0,    1.4142),#; ... %% 'B A {CD}'   [2 1 3 3]   21
      c(-1.5000,         0,         0),#;...  %% 'B {CD} A'   [3 1 2 2]   22
      c(-1.2500,    0.7500,    0.3536),#;...  %% 'B C {AD}'   [3 1 2 3]   23
      c(-0.7500,   -0.7500,    1.0607),#;...  %% 'B {AD} C'   [2 1 3 2]   24
      c(-1.2500,   -0.7500,    0.3536),#;...  %% 'B D {AC}'   [3 1 3 2]   25
      c(-0.7500,    0.7500,    1.0607)# ;...  %% 'B {AC} D'   [2 1 2 3]   26
    )
    
    
    ranksB=rbind(
      c(4,1,3,2),c(4,1,2,3),c(3,1,2,4),c(2,1,3,4),c(2,1,4,3),c(3,1,4,2),c(2,1,2,2),
      c(2,1,3,3),c(3,1,2,2),c(3,1,2,3),c(2,1,3,2),c(3,1,3,2),c(2,1,2,3)
    )
    
    
    
    #-------------------------------------------------
    #exagon C first
    
    E3=rbind(
      c(-1.0,  1.0, -0.70710677),# ...    %%'C B D A'     [4 2 1 3]   27
      c(-0.5,  0.5, -1.4142135),#  ...    %%'C D B A'     [4 3 1 2]   28
      c(0.5,  0.5, -1.4142135),#   ...    %%'C D A B'     [3 4 1 2]   29
      c(1.0,  1.0, -0.70710677),#  ...    %%'C A D B'     [2 4 1 3]   30
      c(0.5,  1.5,  0.0),#         ...    %%'C A B D'     [2 3 1 4]   31
      c(-0.5,  1.5,  0.0) #        ...    %%'C B A D'     [3 2 1 4]   32
    )
    
    MC=apply(E3,2,mean) #center of exagon C  C(ABD)        [2 2 1 2] 33
    
    E3T=rbind(
      c(0,    0.5000,   -1.4142),      #;...  %% 'C D {AB}'   [3 3 1 2]   34
      c(0,    1.5000,         0),      #;...  %% 'C {AB} D'   [2 2 1 3]   35
      c(0.7500,    1.2500,   -0.3536), #;...  %% 'C A {BD}'   [2 3 1 3]   36
      c(-0.7500,    0.7500,   -1.0607),#;...  %% 'C {BD} A'   [3 2 1 2]   37
      c(0.7500,    0.7500,   -1.0607), #;...  %% 'C {AD} B'   [2 3 1 2]   38
      c(-0.7500,    1.2500,   -0.3536) #;...  %% 'C B {AD}'   [3 2 1 3]   39
    )
    
    ranksC=rbind(c(4,2,1,3),c(4,3,1,2),c(3,4,1,2),c(2,4,1,3),c(2,3,1,4),c(3,2,1,4),c(2,2,1,2),
                 c(3,3,1,2),c(2,2,1,3),c(2,3,1,3),c(3,2,1,2),c(2,3,1,2),c(3,2,1,3)             
    )
    
    
    #--------------------------------------
    #exagon D first
    
    E4=rbind(
      c(-1.0, -1.0, -0.70710677),# ...    %%'D B C A'     [4 2 3 1]   40
      c(-0.5, -0.5, -1.4142135),#  ...    %%'D C B A'     [4 3 2 1]   41
      c(0.5, -0.5, -1.4142135),#   ...    %%'D C A B'     [3 4 2 1]   42
      c(1.0, -1.0, -0.70710677),#  ...    %%'D A C B'     [2 4 3 1]   43
      c(0.5, -1.5,  0.0),#         ...    %%'D A B C'     [2 3 4 1]   44
      c(-0.5, -1.5,  0.0) #        ...    %%'D B A C'     [3 2 4 1]   45
    )
    
    MD=apply(E4,2,mean) #center of exagon D  D(ABC)        [2 2 2 1]  46
    
    E4T=rbind(
      c(0,   -0.5000,   -1.4142),      #;...  %% 'D C {AB}'   [3 3 2 1]   47
      c(0,   -1.5000,         0),      #;...  %% 'D {AB} C'   [2 2 3 1]   48
      c(-0.7500,   -1.2500,   -0.3536),#;...  %% 'D B {AC}'   [3 2 3 1]   49
      c(0.7500,   -0.7500,   -1.0607), #;...  %% 'D {AC} B'   [2 3 2 1]   50
      c(-0.7500,   -0.7500,   -1.0607),#;...  %% 'D {BC} A'   [3 2 2 1]   51
      c(0.7500,   -1.2500,   -0.3536)  #;...  %% 'D A {BC}'   [2 3 3 1]   52
    )
    
    ranksD=rbind(c(4,2,3,1),c(4,3,2,1),c(3,4,2,1),c(2,4,3,1),c(2,3,4,1),c(3,2,4,1),c(2,2,2,1),
                 c(3,3,2,1),c(2,2,3,1),c(3,2,3,1),c(2,3,2,1),c(3,2,2,1),c(2,3,3,1)             
    )
    
    
    #squares----------------------------------------------------------
    ESQ=rbind(
      c(0,    0.5000,    1.4142), #      ; ... %% '{AB} C D'   [1 1 2 3]   53
      c(0,   -0.5000,    1.4142), #      ; ... %% '{AB} D C'   [1 1 3 2]   54
      c(-0.5000,         0,   -1.4142), # ;...  %% '{CD} B A'   [3 2 1 1]   55
      c(0.5000,         0,   -1.4142), #  ;...  %% '{CD} A B'   [2 3 1 1]   56
      c(1.2500,    0.7500,   -0.3536), #  ;...  %% '{AC} D B'   [1 3 1 2]   57
      c(0.7500,1.2500,0.3536),       #  ;...  %% '{AC} B D'   [1 2 1 3]   58
      c(-1.2500,   0.7500,   -0.3536), #  ;...  %% '{BC} D A'   [3 1 1 2]   59
      c(-0.7500,    1.2500,    0.3536), # ;...  %% '{BC} A D'   [2 1 1 3]   60
      c(1.2500,   -0.7500,   -0.3536), #  ;...  %% '{AD} C B'   [1 3 2 1]   61
      c(0.7500,   -1.2500,    0.3536), #  ;...  %% '{AD} B C'   [1 2 3 1]   62
      c(-1.2500,   -0.7500,   -0.3536), # ;...  %% '{BD} C A'   [3 1 2 1]   63
      c(-0.7500,   -1.2500,    0.3536) #  ;...  %% '{BD} A C'   [2 1 3 1]   64
    )
    
    ranksSQ=rbind(c(1,1,2,3),c(1,1,3,2),c(3,2,1,1),c(2,3,1,1),c(1,3,1,2),c(1,2,1,3),
                  c(3,1,1,2),c(2,1,1,3),c(1,3,2,1),c(1,2,3,1),c(3,1,2,1),c(2,1,3,1)
    )
    
    #------------------------------------------------------------------------
    
    M_AB_CD=apply(rbind(c(0,0.5,1.4142),c(-0.5,0,1.4142),c(0,-0.5,1.4142),
                        c(0.5,0,1.4142)),2,mean) # '{AB}{CD}'   [1 1 2 2]   65
    
    M_AC_BD=apply(rbind(c(1.25,0.75,0.3536),c(1.25,0.75,-0.3536),c(0.75,1.25,-0.3536),
                        c(0.75,1.25,0.3536)),2,mean) #'{AC}{BD}'   [1 2 1 2] 66
    
    M_BC_AD=apply(rbind(c(-0.75,1.25,-0.3536),c(-1.25,0.75,-0.3536),c(-1.25,0.75,0.3536),
                        c(-0.75,1.25,0.3536)),2,mean) # '{BC}{AD}'   [2 1 1 2] 67
    
    EMID=rbind(M_AB_CD, M_AC_BD, M_BC_AD, M_AB_CD*-1, M_AC_BD*-1, M_BC_AD*-1,MA*-1,MB*-1,MC*-1,MD*-1) 
    #M_AB_CD*-1 = '{CD}{AB}'   [2 2 1 1] 68
    #M_AC_BD*-1 = '{BD}{AC}'   [2 1 2 1] 69
    #M_BC_AD*-1;= '{AD}{BC}'   [1 2 2 1] 70
    #MA*-1;=      '{BCD}A'   [2 1 1 1] 71
    #MB*-1;=      '{ACD}B'   [1 2 1 1] 72
    #MC*-1;=      '{ABD}C'   [1 1 2 1] 73
    #MD*-1;=      '{ABC}D'   [1 1 1 2] 74
    #last =        {ABCD}    [1 1 1 1] 75 
    
    rankres=rbind(c(1,1,2,2),c(1,2,1,2),c(2,1,1,2),c(2,2,1,1),c(2,1,2,1),
                  c(1,2,2,1),c(2,1,1,1),c(1,2,1,1),c(1,1,2,1),c(1,1,1,2),c(1,1,1,1)
    )
    
    
    
    
    EE=rbind(E1,MA,E1T,E2,MB,E2T,E3,MC,E3T,E4,MD,E4T,ESQ,EMID)
    coord=rbind(EE,apply(EE,2,mean))
    ranks=rbind(ranksA,ranksB,ranksC,ranksD,ranksSQ,rankres)
    
    if (is.null(L)){
      rr=labels(ranks,4,labs=2)
    } else {
      rr=labels(ranks,4,L,labs=1)
    }
    
    
    
    plot3d(coord, type = 'p', xlab = '', ylab = '', zlab = '', add = T, 
           aspect = T, box = F, axes = F, col = 1)
    
    #exagon A first
    segments3d(E1[c(1,2),], lwd=1, col = 1)
    segments3d(E1[c(2,3),], lwd=1, col = 1)
    segments3d(E1[c(3,4),], lwd=1, col = 1)
    segments3d(E1[c(4,5),], lwd=1, col = 1)
    segments3d(E1[c(5,6),], lwd=1, col = 1)
    segments3d(E1[c(1,6),], lwd=1, col = 1)
    segments3d(E1T[c(1,2),], lwd=0.5, col = 'gray')
    segments3d(E1T[c(3,4),], lwd=0.5, col = 'gray')
    segments3d(E1T[c(5,6),], lwd=0.5, col = 'gray')
    #exagon B first
    segments3d(E2[c(1,2),], lwd=1, col = 1)
    segments3d(E2[c(2,3),], lwd=1, col = 1)
    segments3d(E2[c(3,4),], lwd=1, col = 1)
    segments3d(E2[c(4,5),], lwd=1, col = 1)
    segments3d(E2[c(5,6),], lwd=1, col = 1)
    segments3d(E2[c(1,6),], lwd=1, col = 1)
    segments3d(E2T[c(1,2),], lwd=0.5,col = 'gray')
    segments3d(E2T[c(3,4),], lwd=0.5,col = 'gray')
    segments3d(E2T[c(5,6),], lwd=0.5,col = 'gray')
    #exagon C first
    segments3d(E3[c(1,2),], lwd=1, col = 1)
    segments3d(E3[c(2,3),], lwd=1, col = 1)
    segments3d(E3[c(3,4),], lwd=1, col = 1)
    segments3d(E3[c(4,5),], lwd=1, col = 1)
    segments3d(E3[c(5,6),], lwd=1, col = 1)
    segments3d(E3[c(1,6),], lwd=1, col = 1)
    segments3d(E3T[c(1,2),], lwd=0.5,col = 'gray')
    segments3d(E3T[c(3,4),], lwd=0.5,col = 'gray')
    segments3d(E3T[c(5,6),], lwd=0.5,col = 'gray')
    #exagon D first
    segments3d(E4[c(1,2),], lwd=1, col = 1)
    segments3d(E4[c(2,3),], lwd=1, col = 1)
    segments3d(E4[c(3,4),], lwd=1, col = 1)
    segments3d(E4[c(4,5),], lwd=1, col = 1)
    segments3d(E4[c(5,6),], lwd=1, col = 1)
    segments3d(E4[c(1,6),], lwd=1, col = 1)
    segments3d(E4T[c(1,2),], lwd=0.5,col = 'gray')
    segments3d(E4T[c(3,4),], lwd=0.5,col = 'gray')
    segments3d(E4T[c(5,6),], lwd=0.5,col = 'gray')
    #squares
    segments3d(rbind(E1[1,],E2[4,]), lwd=1, col = 1)
    segments3d(rbind(E1[6,],E2[5,]), lwd=1, col = 1)
    segments3d(rbind(E3[2,],E4[2,]), lwd=1, col = 1)
    segments3d(rbind(E3[3,],E4[3,]), lwd=1, col = 1)
    segments3d(rbind(E1[5,],E4[5,]), lwd=1, col = 1)
    segments3d(rbind(E1[4,],E4[4,]), lwd=1, col = 1)
    segments3d(rbind(E2[1,],E4[1,]), lwd=1, col = 1)
    segments3d(rbind(E2[6,],E4[6,]), lwd=1, col = 1)
    segments3d(rbind(E1[2,],E3[5,]), lwd=1, col = 1)
    segments3d(rbind(E1[3,],E3[4,]), lwd=1, col = 1)
    segments3d(rbind(E2[2,],E3[1,]), lwd=1, col = 1)
    segments3d(rbind(E2[3,],E3[6,]), lwd=1, col = 1)
    #other exagons
    segments3d(rbind(ESQ[1,],ESQ[2,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[3,],ESQ[4,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[5,],ESQ[6,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[7,],ESQ[8,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[9,],ESQ[10,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[11,],ESQ[12,]), lwd=0.5, col = 'gray')
    segments3d(rbind(E1T[1,],E3T[3,]), lwd=0.5, col = 'gray')
    segments3d(rbind(E2T[3,],E3T[6,]), lwd=0.5, col = 'gray')
    segments3d(rbind(E1T[3,],E2T[1,]), lwd=0.5, col = 'gray')
    segments3d(rbind(E1T[5,],E4T[6,]), lwd=0.5, col = 'gray')
    segments3d(rbind(E2T[5,],E4T[3,]), lwd=0.5, col = 'gray')
    segments3d(rbind(E3T[1,],E4T[1,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[1,],E3T[2,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[6,],E2T[6,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[8,],E1T[6,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[2,],E4T[2,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[7,],E4T[5,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[5,],E4T[4,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[12,],E1T[2,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[10,],E2T[4,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[4,],E1T[4,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[9,],E3T[5,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[3,],E2T[2,]), lwd=0.5, col = 'gray')
    segments3d(rbind(ESQ[11,],E3T[4,]), lwd=0.5, col = 'gray')
    #-----------------------------
    
    if (is.null(X)){
      indplot=matrix(1:75,ncol=1)
    } else {
      

      #ranksinplot=ranks[apply(o, 2, any),]
      
      o2 = outer(seq_len(nrow(ranks)), seq_len(nrow(X)), Vectorize(
        function(i, j) which(all(ranks[i,]==X[j,]))
      ))
      
      indexing=o2==1
      #print(indexing)
      indexing[indexing==TRUE]=1
      indexing[is.na(indexing)]=0
      indplot=which(rowSums(indexing)==1)
      #print(indplot)
#       ranksinplot=ranks[indplot,]
#       
#       o = outer(seq_len(nrow(X)), seq_len(nrow(ranksinplot)), Vectorize(
#         function(i, j) all(X[i,]==ranksinplot[j,])
#       ))
#       
#       
#       
#       indexlabs=o==1
#       #print(indexing)
#       indexlabs[indexlabs==TRUE]=1
#       indexlabs[is.na(indexlabs)]=0
#       indlab=matrix(nrow=ncol(indexlabs),ncol=1)
#       for (j in 1:nrow(indlab)){
#         indlab[j,1]=which(indexlabs[j,]==1)
#       }

      
      
    }
    #points(coord[indplot,1],coord[indplot,2],pch=16)
    #points3d(coord[indplot,],col="blue",cex=sqrt(100*((Wk/sum(Wk))/pi)))

    if (is.null(Wk)){
      
      spheres3d(coord[indplot,],col="blue",radius=0.02)
      
    } else {
      idwk=matrix(0,nrow(X),1)
      counter=0
      for (i in 1:75){
        for (j in 1:nrow(X)){
          check=sum(pos=X[j,]==ranks[i,])
          if (check==4){
            counter=counter+1
            idwk[counter]=j
            break}
        }
      }
      
      spheres3d(coord[indplot,],col="blue",radius=sqrt(((Wk[idwk]/sum(Wk))/(25*pi))))
      
    }
    #text(tcoord[indplot,],rr[indplot,])
    text3d(coord[indplot,1]+0.1, coord[indplot,2]+0.1, 
           coord[indplot,3]+0.1,rr[indplot,],col=1,cex=0.7)
  }
}

#----------------------------------------------------------------------------------------

#------------------------------------

BBFULL = function(X,Wk=NULL,PS=TRUE)  {
  #Branch and Bound algorithm to find median ranking in the space of full rankings
  #X is a data matrix in which the rows are the judges and the columns indicates the objects
  #Wk is the vector of weigths
  if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
  }
  
  
  
  M = nrow(X)
  N=ncol(X)
  tic = proc.time()[3]
  if (M==1) {
    consensus = X
    TauX = 1
  } else {
    if (!is.null(Wk)) {
      
      if (is.numeric(Wk)) {
        Wk=matrix(Wk,ncol=1)
      }
      
      cij = combinpmatr(X,Wk)
    } else {
      cij = combinpmatr(X)
    }
    
    if (sum(cij==0)==nrow(cij)^2){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    R=findconsensusBB(cij)
    cons1=BBconsensus(R,cij,FULL=TRUE)
    consensus1=cons1$cons
    Po=cons1$pen
    consensus=BBconsensus2(consensus1,cij,Po,PS,FULL=TRUE)
  }
  
  
  if (nrow(consensus)==1) {
    
    Sij=scorematrix(consensus)
    
    if (!is.null(Wk)){
      TauX=sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      TauX=sum(cij*Sij) / (  M*(N*(N-1)) )
    }
    
  } else {
    
    TauX=matrix(0,nrow(consensus),1)
    
    for (k in 1:nrow(consensus)) {
      
      Sij=scorematrix(t(matrix(consensus[k,])))
      
      if (!is.null(Wk)) {
        
        TauX[k,1] = sum(cij*Sij) / ( sum(Wk)*(N*(N-1)) )
        
      } else {
        
        TauX[k,1] = sum(cij*Sij) / (M*(N*(N-1)))
        
      }
      
    }
    
  }
  toc = proc.time()[3]
  colnames(consensus)=colnames(X) 
  #consensus=reordering(consensus)
  eltime=toc-tic
  return(list(Consensus=reordering(consensus), Tau=TauX, Eltime=eltime) )
}
