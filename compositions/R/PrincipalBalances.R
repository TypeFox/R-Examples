gsi.PrinBal = function(x, method="PBhclust"){
  # This function is a wrapper around the existing methods of principal balances
  # All principal balance methods MUST start with the "PB" marker and 
  #   require an acomp data set; no check is done
  if(method=="PBhclust"){
    if(class(x)!="acomp") warning("gsi.PrinBal: data given in x is not acomp! Do not trust results")
    dd = as.dist(variation(x))
    hh = hclust(dd, method="ward")
    sgn = gsi.merge2signary(hh$merge)
  }else{
  print("gsi.PrinBal: starting computations of slow methods")
  }
  if(method=="PBangprox"){
    if(ncol(x)>10) warning("gsi.PrinBal: data given in x has many columns: computation is expected to be slow")
    # start with principal components and an empty vector
    PC = t(princomp(cpt(x))$loadings)[-ncol(x),]
    PC = unclass(PC)
    vc = rep(0,ncol(x)) 
    names(vc) = colnames(x)
    write.table(t(as.matrix(names(vc))),".APtable", 
                row.names=FALSE, col.names=FALSE, append=FALSE, quote=FALSE)
    # apply recursion
    lsn = gsi.APrecursiu(colnames(x), x, PC, vc)
    sgn = t(read.table(".APtable", header=TRUE))
    write.table("This file can be removed without harm",quote=FALSE,
                file = ".APtable", row.names=FALSE, col.names=FALSE, append=TRUE)
  }
  if(method=="PBmaxvar"){
    if(class(x)!="acomp") warning("gsi.PrinBal: data given in x is not acomp! Do not trust results")
    sgn = t(gsi.mvOPTIMAL(cpt(x)))
  }
  V = gsi.buildilrBase(sgn)
  return(V)
}


# function gsi.APdivideix: looks for combinations of parts and compute angles with PC's 
#   used in angular proximity principal balance
gsi.APdivideix = function(x, PC){
  # set of all directions
  vcs = rep(0, ncol(x))
  names(vcs) = colnames(x)
  matsub = sapply(1:(ncol(x)-1), function(m){
    num = combn(x=colnames(x),m)
    res = sapply(1:ncol(num), function(i){
      nn = length(num[,i]) 
      nd = ncol(x)-nn
      v = vcs-nn
      v[num[,i]] = nd
      return(v)
    })
    return(res)
  }
  )
  
  # destroy the list structure
  matsub = unlist(matsub)
  dim(matsub) = c(ncol(x), length(matsub)/ncol(x))
  matsub = t(matsub)
  matsub = unclass(matsub)
  # normalize and redimensionalize
  matsub = sapply(1:nrow(matsub), function(i){matsub[i,]/sqrt(sum(matsub[i,]^2))} )
  matsub = t(matsub)
  
  # angles
  angles = matsub %*% t(PC[,colnames(x)])
  tk = (angles == max(angles)) %*% rep(1,nrow(PC))
  
  # best balance, finding parts in the numerator and in the denominator
  balance = matsub[tk==1,]
  dp = colnames(x)[balance<0]
  np = colnames(x)[balance>0]
  
  return(list(np, dp))
}

# function gsi.APrecursiu: comutes angular proximity recursively
gsi.APrecursiu = function(nomvar,x,PC,vc){  
  dv = gsi.APdivideix(x[,nomvar],PC)
  vc2 = vc*0
  vc2[dv[[1]]] = +1 
  vc2[dv[[2]]] = -1
  dim(vc2) = c(1,length(vc2))
  #output signs as they are computed
  write.table(vc2,".APtable", row.names=FALSE, col.names=FALSE, append=TRUE)
  nm = lapply(dv, function(noms){
    if(length(noms)>1){
      gsi.APrecursiu(noms, x, PC, vc)
    }
  })
}






##### FUNCTION RSCODE ##################################
# computes the number of +,- in SBPcode[1:D-1,1:D]
# on return RSCODE[1:D-1,2] contains in column 1 num of +
#                                    in column 2 num of -
gsi.mvRSCODE = function(SBPcode) {
#.....  compute num +  -
   r = rep(0,nrow(SBPcode))
   s = r
for(i in 1:nrow(SBPcode)) {
   for(j in 1:ncol(SBPcode)) {
     if(SBPcode[i,j]==1) r[i] = r[i]+1
     if(SBPcode[i,j]== -1) s[i] = s[i]+1
    }
}
return(cbind(r,s))  } 
########## END RSCODE function ###########################

########## FUNCTION NEXTACTIVE ###########################
# Given a partial SBP, computes the next group of parts to
# be partitioned.
# jbase[1:npart-1,1:npart] contains the partial SBPcode up to
#           the iorder; higher orders are zeroes of other 
#           constants
# jstatus[1:npart-1,1:npart] is a working matrix containing
#           information about partial SBP up to order iorder-1
#           jstatus[i,j]=npart+1 means that j-part is not involved
#                in further partitions
#                For iorder=1 jstatus[1,] contain zeros
#           jstatus is modified on the output up to iorder+1
#                At output with iorder=npart-1 jstatus[npart,]
#                must be npart+1  
#  on return, NEXTACTIVE contains an updated version of jstatus
#           completed up to order iorder+1
#           Next gruop of parts to be partitioned is that in
#           jstatus[iorder+1,] has minimum and equal codes 
#           active[1:npart]:indicator vector pointing out
#           the parts to be partitioned in the next step of the
#           SBP can be extracted by
#             lev = min(jstatus[iorder+1,]
#             for(j in 1:npart){
#               if(jstatus[iorder+1,j] ==lev) active[j]=1
#               else  active[j]=0  }
gsi.mvNEXTACTIVE = function(npart,iorder,jstatus,jbase){

rs = gsi.mvRSCODE(jbase)
iorder1 = iorder +1

for(ipart in 1:npart){
   if(jbase[iorder,ipart] == 0 ){
      if(jstatus[iorder,ipart] != npart+1) jstatus[iorder1,ipart]=jstatus[iorder,ipart]+1
      else     jstatus[iorder1,ipart]=npart+1
   }
   if(jbase[iorder,ipart] == -1 & jstatus[iorder,ipart] != npart+1){
      if(rs[iorder,2] == 1)  jstatus[iorder1,ipart]=npart+1   
      else             jstatus[iorder1,ipart]=jstatus[iorder,ipart]+1
   }
   if(jbase[iorder,ipart] == 1  & jstatus[iorder,ipart] != npart+1){
      if(rs[iorder,1] == 1)   jstatus[iorder1,ipart]=npart+1   
      else            jstatus[iorder1,ipart]=jstatus[iorder,ipart]
   }
}  #ipart

return(jstatus)
}
######## END FUNCTION NEXTACTIVE #############


############# FUNCTION PROJclr ###############
gsi.mvPROJclr = function(npart,iorder,SBPcode,xclr) {

# projects the clr-data-set into the orthogonal
# subspace defined by the SBPcode up to the order iorder
# the projected clr-data-set is returned in xclrp
iorder1=iorder+1
npart1=npart-1
code = SBPcode

for(i in iorder1:npart1) code[i,]=rep(0,length(npart))
Vcontr = gsi.buildilrBase(t(code))
for(j in iorder1:npart1){
  for(i in 1:npart){
     Vcontr[i,j]=0
}}
# clr in subspace partial SBP
# parece que ya es clr (suma 0)
xclrsub =  xclr %*% Vcontr %*% t(Vcontr)
xclrp=xclr-xclrsub 

return(xclrp)  }
############# END FUNCTION PROJclr ###############



############# FUNCTION SVDSUGG ###############
gsi.mvSVDSUGG = function(npart,xclr,active){
# carries out a SVD of xclr[,npart]
# xclr is a sample clr (assumed centred)
# using the first column of loadings in V a partition
#  of parts indicated by active[1:npart]
# the suggested partition is stored in partsugg[1:npart]

xsvd = svd(xclr)
V =xsvd$v
partsugg = V[,1]
for(i in 1:npart){
  if(active[i]==0)  partsugg[i]=0
  else {
     if(partsugg[i]>=0)  partsugg[i]=1
     else  partsugg[i]=-1    }
}

return(partsugg) }
############# END FUNCTION SVDSUGG ###############

############# FUNCTION REFINESBP #################
gsi.mvREFINESBP = function(npart,iorder,SBPcode,xclr) {
# tries to change SBPcode[iorder,] to optimize
# explained variance of the balance.
# It changes one part in the balance SBPcode[iorder,]
#   from +1 to -1 and viceversa. 
#   it stops when no change is carried out or when
#   the maximum number of iterations (50) has been attained
# The new (changed or unchanged) SBPcode[1:npart1,1:npart]=code
#   is returned. Previous orders of SBP are not changed.
print("order partition")
print(iorder)
rs = gsi.mvRSCODE(SBPcode) 
r=rs[iorder,1]
s=rs[iorder,2]
prevr = r
prevs = s
jcode = rep(0,npart)
jcode = SBPcode[iorder,1:npart]
inivar = gsi.mvVARBAL(xclr,jcode)
prevvar = inivar
newvar = rep(0,npart)
prevcode = jcode

nchange=1
iter=0

while(nchange == 1 & iter<50){
   nchange=0
   for(ich in 1:npart){
     activechange=0
     newcode=prevcode
     if(prevcode[ich] == 1 & prevr>1 ) {
       newcode[ich]=-1
       newr=prevr-1
       news=prevs+1
       activechange=1
     }
     if(prevcode[ich] == -1 & prevs>1 ) {
       newcode[ich]=1
       newr=prevr+1
       news=prevs-1 
       activechange=1
     }
     if(activechange == 1){
        newvar[ich] = gsi.mvVARBAL(xclr,newcode)
     }
   }  # for ich
   print("variance comparison")
   print(prevvar)
   print(newvar)
   newmaxvar=max(newvar)
   newich=which.max(newvar)
   if(newmaxvar > prevvar){
     newcode=prevcode
     if(prevcode[newich] == 1) {
       newcode[newich]=-1
       newr=prevr-1
       news=prevs+1
     }
     if(prevcode[newich] == -1 ) {
       newcode[newich]=1
       newr=prevr+1
       news=prevs-1 
     }
     prevcode=newcode
     prevr=newr
     prevs=news
     prevvar=newmaxvar
     newvar=rep(0,npart)
     nchange=1
     iter=iter+1
   }   # if newmaxvar

}  # while nchange

code=SBPcode
code[iorder,1:npart]=prevcode

return(code)  }

############# END FUNCTION REFINESBP  #################

############ FUNCTION VARBAL #################
### funtion computing variance of one balance
### using the logs of sample
### balance coded in jcode as +1,-1, 0 's
### as a vector
gsi.mvVARBAL = function(xclr, jcode){
  ## attention: this function requires R-ization
rs = gsi.mvRSCODE(t(as.matrix(jcode)))
r  =  rs[1,1]
s  =  rs[1,2]

#..... compute balances
bal = rep(0,nrow(xclr))
for(i in 1:nrow(xclr)) {
   for(j in 1:ncol(xclr)) {
      if(jcode[j]== 1) bal[i] = bal[i] + sqrt(s/(r*(r+s)))*xclr[i,j]
      if(jcode[j]== -1) bal[i] = bal[i] - sqrt(r/(s*(r+s)))*xclr[i,j]
}}
return(var(bal))
}
############ END FUNCTION VARBAL #################



############## FUNCTION gsi.mvOPTIMAL #############

gsi.mvOPTIMAL = function(xclr){

#
# carries out the loop in order partition searching
# for the optimal SBP corresponding to Mvar (maximum variance)
# principal balances.
# Input: xclr(1:nrow, 1:npart)
#      should be the CENTRED clr of the compositional data set 
# Output: on return, code of of the optimal SBP found 
#         SBPcode(1:npart , 1:npart)
########################################################

npart = ncol(xclr)
npart1 = npart -1

#  define SBPcode and xclrorth (working clr)
SBPcode = matrix(0,nrow=npart1,ncol=npart)
xclrorth=xclr
#  define active(1:npart) working indicator
active=rep(1,length=npart)
#  define jstatus(1:npart , 1:npart) working SBP-like matrix
jstatus= matrix(0,nrow=npart,ncol=npart)

# Main loop

for(iorder in 1:npart1) {

#    SVD suggesting a partition
  partsugg = gsi.mvSVDSUGG(npart,xclrorth,active)
#    suggestion taken as optimum
  SBPcode[iorder,1:npart]=partsugg
#    refine suggested partition
  partrefin =gsi.mvREFINESBP(npart,iorder,SBPcode,xclr)
  SBPcode = partrefin
#    orthogonal projection of clr-data
  if(iorder < npart1){
   xclrorth = gsi.mvPROJclr(npart,iorder,SBPcode,xclrorth)
#    select next step for partition
   if(iorder < npart1){
    jstatus = gsi.mvNEXTACTIVE(npart,iorder,jstatus,SBPcode)
       lev = min(jstatus[iorder+1,])
       for(j in 1:npart){
             if(jstatus[iorder+1,j] ==lev) active[j]=1
              else  active[j]=0  }
       }  # j
    } # if(iorder 
}  # iorder

colnames(SBPcode)=colnames(xclr)
return(SBPcode)  }
##############END FUNCTION gsi.mvOPTIMAL ####################





