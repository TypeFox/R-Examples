
# very simple function, used for demonstration purposes...
mpause=function(text1="",text2="")
{
 if(text1!="") print(text1)
 if(text2=="") text2="-- Pause: press enter to continue --"
 print(text2)
 readLines(n=1)
 return(TRUE)
}

#----------------------------------------------------------------------------------------------------
#-- Auxiliary private functions, in principle you should not need to use these:
#----------------------------------------------------------------------------------------------------
# 1) binary variables should not be scaled.
# 2) in some "rare" cases, by sampling, some variables may be constant. In such cases, 
#    scale leads to NA values. 
# This is why these 2 functions were made:
#-- scale only the inputs:
scaleinputs<-function(data,outindex)
{ 
 #print(" > SCALE INPUTS............")
 #print("  >... scale...")
 #print(data[1,])
 #print(paste("first scale out index:",outindex))
 #print(names(data))

 L<-NCOL(data)
 sx<-vector(length=L)
 cx<-vector(length=L)

#print(summary(data))
 for(i in 1:L)
   {
#cat("i:",i,"\n")
     if(is.factor(data[,i][1]) || i==outindex ) # do not scale a factor or output
     { sx[i]<- 0; cx[i]<- 0} # do not scale
     else 
     { M<-mean(data[,i])
       SD<-sd(data[,i])
       if(SD==0 || (SD==1 && M==0)) # do not scale if data constant or already scaled
              { sx[i]<-0; cx[i]<-0 }
       else { # finally: do scale!!!
            data[,i]<-(data[,i]-M)/SD 
            sx[i]<-SD
            cx[i]<-M
          }
     }
   }
 return(list(data=data,cx=cx,sx=sx))
}

# transform scaled inputs into their normal values
scaleinputs2<-function(data,cx,sx)
{ 
 #print(paste("2nd scale L Data:",ncol(data)," L sx:",length(sx)))
 L<-length(sx)
 for(i in 1:L) if(sx[i]!=0)  data[,i]<-(data[,i]-cx[i])/sx[i]
 return(data)
}


# xtransform and not "transform" from base
xtransform=function(x,transform,A=0,B=0,attributes=NULL)
{
 if(transform=="none") return (TRUE)
 if(transform=="log" || transform=="logpositive") x<-(log(x+1))  
 else if(transform=="scale" && B!=0) # use the scaling function
      x<-(x-A)/B
 # add other transformations here ...
 else if(transform=="bool")
 { # A = outindex
  if(A==0) A=NCOL(x)
  NOUT=names(x)[A]
  ATT=setdiff((1:NCOL(x)),A)
  I=as.data.frame(x[,ATT]>0)
  for(i in ATT) I[,i]=as.numeric(I[,i])
  x=cbind(I,x[,A])
  names(x)[A]=NOUT
  #print(" >> X:")
  #print(x)
 }
 else if(transform=="TF")
 {
# current problem: x is smaller than data, addapt???
   # A= outindex, B= vector computed previously 
   if(A==0) A=NCOL(x)
   if(is.null(attributes)) ATT=setdiff((1:NCOL(x)),A)
   else { 
          ATT=setdiff(attributes,A)
        }
   if(length(B)==1)  # B==0, no values computed!!!
   {
    AUX=vector(length=length(ATT))
    N=nrow(x)
    for(i in ATT) 
    {  AUX[i]=log10(N/sum(x[,i]>0))
       if(AUX[i]==Inf) AUX[i]=0 # to avoid Inf/NA errors...
       x[,i]=log10(x[,i]+1) # TF

       x[,i]=x[,i]*AUX[i] # TF-IDF

       temp=sum(x[,i]^2) # IDF
       if(temp==0) temp=0.000001 else temp=sqrt(temp)
       x[,i]=x[,i]/temp
       AUX[i]=AUX[i]/temp
    }
    #print(" >> X:")
    #print(x)
    x=list(x=x,B=AUX)
   }
   else # B is a vector, already computed the values...
   {
    if(is.data.frame(B)) B=as.numeric(B)
    for(i in ATT) x[,i]=log10(x[,i]+1)*B[i] # TF-IDF
   }
 }
    #print("----\n>> X:")
    #print(x)
 return (x)
}

invtransform=function(x,transform,A=0,B=0)
{
 if(transform=="log") return (exp(x)-1) 
 else if(transform=="scale" && B!=0) # 
      return (A+B*x)
 else if(transform=="positive") # 
     {
      x[x<0]<-0 # transform all negative values into zero
      return (x)
     }
 else if(transform=="logpositive")
     {
      x<-(exp(x)-1) 
      x[x<0]<-0 # transform all negative values into zero
      return (x)
     }
 # add other transformations here ... 
 else return (x)
}

# ====== MISSING DATA ======================================================
#---------------------------------------------------------
# substitutes missing values (NA) for a given attribute/feature with Value
# Value can be one number or a vector. if vector, length should be equal to the number of missing values
# if Attribute is NULL performs imputation over all data.frame!!!
#---------------------------------------------------------
imputation=function(imethod="value",D,Attribute=NULL,Missing=NA,Value=1)
{
 return(switch(imethod,value=impvalue(D,Attribute,Missing=Missing,Value=Value),hotdeck=hotdeck(D,Attribute,Missing,K=Value),D))
}

impvalue=function(D,Attribute=NULL,Missing=NA,Value)
{
 LV=length(Value)
 if(is.null(Attribute)) N=NCOL(D)
 else 
 { N=1
   if(class(Attribute)=="character") 
    { ND<-names(D); A<-which(ND==Attribute)}
   else A<-Attribute
 }
 
 for(j in 1:N)
 {
  if(is.null(Attribute)) A=j
  AUX<-D[,A]
  if(is.na(Missing)) I<-which(is.na(AUX))
  else I<-which(AUX==Missing)
  if(length(Value)==1) 
 	for(i in 1:length(I)) 
        {
		if(LV==1) AUX[I[i]]<-Value
                else AUX[I[i]]=Value[i]
        }
  else for(i in 1:length(I)) AUX[I[i]]<-Value[i]
  D[,A]<-AUX
 }
 return(D)
}

# auxiliary function, please do not use this!!!
missingatts=function(D,Missing=NA)
{
 C=NCOL(D); RES=NULL
 for(i in 1:C)
  {
    if(is.na(Missing)) S=sum(is.na(D[,i]))
    else S=sum(D[,i]==Missing) 
    if(!is.na(S) && S>0) RES=c(RES,i)
  }
 return(RES)
}

# uses a nearest neighbour (default 1-NN) to replace the missing value with similar cases
# if Attribute is NULL performs hotdeck over all data.frame
# D is a data.frame (data)
hotdeck<-function(D,Attribute=NULL,Missing=NA,K=1)
{
  #D=d; Attribute=3; Missing=NA; K=1
  #D=d; Attribute=NULL;Missing=NA;K=1
  if(is.null(Attribute)) N=NCOL(D)
  else 
    { N=1
      if(class(Attribute)=="character") 
        { ND<-names(D); A<-which(ND==Attribute)
        }
      else A<-Attribute
    }
  MA=missingatts(D,Missing) 
  ALL=setdiff(1:NCOL(D),MA)
  for(i in 1:N)
  { 
   if(is.null(Attribute)) A=i
   if(is.na(Missing)) I<-which(is.na(D[,A])) else I<-which(D[,A]==Missing) 
   if(length(I)>0)
   { 
    ALL1=c(ALL,A)
    DTR<-D[-I,ALL1]; DTS<-D[I,ALL1];YDTS=NCOL(DTS)
    if(is.factor(D[,A][1])){L=levels(D[,A][1]);DTS[,YDTS]=factor(rep(L[1],length(I)),levels=levels(D[,A])); task="class";}
    else { DTS[,YDTS]=0; task="reg";}
    if(!is.null(names(D))) x<-as.formula(paste(names(D)[A]," ~ .",sep="")) else x<-as.formula("y ~ .")
    M<-fit(x,data=DTR,model="kknn",search=K,task=task)
    P<-predict(M,DTS)
    D=impvalue(D,A,Missing=Missing,Value=P)
   }
  }
  return(D)
}
#===========================================================================


#---------------------------------------------------------
# reduces the number of labels for a given factor
# all labels not included in new reduced levels are transformed into "_OTHER"
#---------------------------------------------------------
delevels<-function(x,levels,label=NULL)
{
 L<-levels(x); LL<-length(levels);
 for(k in 1:LL)
  { I<-which(L==levels[k])
    if(is.null(label)) L[I]<-"_OTHER"
    else L[I]<-label
  }
 levels(x)<-L
 return(x)
}

# -------- addition of 2 factors with the same levels into one factor
# improve this later?
addfactor<-function(f1,f2)
{
 if(!is.null(f1)) {L<-levels(f1[1]); L1<-length(f1);}
 else {L<-levels(f2[1]); L1<-0;}
 L2<-length(f2);
 R<-factor(rep(L[1],length=(length(f1)+length(f2))),levels=L) 
 if(L1>0) R[1:L1]<-f1[1:L1]
 R[(L1+1):length(R)]<-f2[1:L2]
 return(R)
}

#---------------------------------------------------------
# transforms a factor into a numeric variable
#---------------------------------------------------------
factor2numeric<-function(x,levels,numbers)
{
 L<-levels(x); LL<-length(levels);
 res<-vector(length=length(x))
 for(k in 1:LL)
  { I<-which(x==levels[k])
    res[I]<-numbers[k]
  }
 return(res)
}

# ----- get the mode (most common class index) of a factor -------
mostcommon<-function(x){return(which.max(table(x)[]))}

# ----- get the average/median class of an ordered factor -------
middleclass<-function(x,method="mean")
{return(round(switch(method,mean=mean(as.numeric(x)),median(as.numeric(x)))))}

# value: single value or vector of values
filter_equal<-function(D,Attribute,Value,reverse=FALSE)
{
  if(class(Attribute)=="character") 
  { ND<-names(D); A<-which(ND==Attribute)}
  else A<-Attribute
  
  res=NULL; ires=NULL;
  for(i in 1:length(Value))
  {
    if(is.na(Value))
    {
       if(reverse) I=which(!is.na(D[,A]))
       else I=which(is.na(D[,A]))
    }
    else
    {
       if(reverse) I=which(D[,A]!=Value[i])
       else I=which(D[,A]==Value[i])
    }
    if(length(I)>0) ires=c(ires,I) 
  }
  if(length(ires)>0) return (D[ires,])
  else return (NULL)
}

# f is a factor or matrix with factor variables
one_of_c=function(f)
{ 
  
  LF=length(f)
  L=levels(f)
  LN=length(L)
  m=matrix(0,ncol=LN,nrow=LF)
  for(i in 1:LN) m[(which(f==L[i])),i]=1
  return (m)
}


# similar to median except that when an impar number of elements returns the first of the middle and not the average:
medianfirst=function(x)
{
 X=sort.int(x,index.return=TRUE)
 LX=length(X$ix);if(LX%%2==1) mid=(LX-1)/2+1 else mid=LX/2
 return (list(mid=X$ix[mid],val=X$x[mid]))
}


#--- END OF AUXILIAR FUNCTIONS --------------------------------------------

# transform a TS intro a data.frame matrix of #W inputs and y output variables
# W - vector with the sliding time windows
# start - default is one
# end - default is the length of the series
CasesSeries=function(t,W,start=1,end=length(t))
{
 LW=length(W)
 LL=W[LW]
 JL=(end-start+1)-LL
 I=matrix(ncol=(LW+1),nrow=JL)
 S=start-1
 for(j in 1:JL)
 {
  for(i in 1:LW)
  	I[j,i]=t[(S+LL-W[LW-i+1]+j)]
  I[j,(LW+1)]=t[(S+LL+j)]
 }
 D=data.frame(I)
 N=names(D)
 LN=length(N)
 for(i in 1:(LN-1)) N[LN-i]<-paste("lag",W[i],sep="")
 N[LN]="y"
 names(D)<-N
 return (D)
}

lforecast=function(M,data,start,horizon)
{
 #cat("start:",start,"horizon:",horizon,"\n")
 #print(data)

 Y=NCOL(data)
 NW=names(data)
 LW=length(NW)-1
 W=vector(length=LW)
 for(i in 1:LW)
 {
   W[i]=strsplit(NW[LW-i+1],"lag")[[1]][2] 
 }
 W=as.numeric(W)
 
 ML=W[length(W)]
 s1=start
 N=NROW(data)
 if(horizon==-1) horizon=(N-s1)

 F=vector(length=horizon+ML) 
 for(i in 1:ML) F[i]=data[(s1-(ML-i+1)),Y]

 for(i in 1:horizon)
 {
   x=data[1,] 
   for(j in 1:LW)
    {
      x[1,j]=F[(ML+i)-W[(LW-j+1)]]
    }
   F[ML+i]=predict(M,x)
   #cat("i:",i,"d:",as.numeric(x[1,1:LW]),"f:",F[ML+i],"\n")
 }
 return (F[(ML+1):(ML+horizon)])
}

# given d data, returns the number of levels for each attribute
datalevels=function(d,L=7,Lfactor=FALSE)
{
 NC=ncol(d)
 res=rep(L,NC)
 for(i in 1:NC)
              {
                NL=length(levels(d[1,i]))
                if(NL>0){if(Lfactor) res[i]=NL else res[i]=min(NL,L)}
              }
 return(res) 
}

factorize<-function(x,limits,labels)
{ warning("Deprecated function, please use instead: cut(x,limits,labels)")
 N<-length(labels)
 res<-vector(length=length(x))
 for(i in 1:N)
   { 
     Ind<-which(x>=limits[i] & x<limits[i+1])
     res[Ind]<-labels[i]
   }
 return(factor(res))
}

# ----------NOT USED CURRENTLY (stored only for backup purposes): --------
# time series: get a reasonable value for first sazonality period or DEFAULT
# load a .ts file into a vector time series
if(FALSE){

#---------------------------------------------------------
# equal to cut !!!
# returns the variable x (vector made of a continuous/numeric values) into
# a vector of discrete values (factor)
# 
#
# Parameters:
# x - vector of numeric values
# levels - a vector with the limits for each vector. 
# labels - a vector of length(levels)-1 with the labels for each 2 limits
#     
# Each label is assigned with level[i] <= label < level[i+1]
#
# example: 
# N<-1:10
# F<-factorize(N,c(0,5,11),c("less-than-five","higher-than-five"))
#---------------------------------------------------------

readts=function(filename,header=FALSE)
{
 TSNAME<-paste(name,".ts",sep="")
 t<-read.table(TSNAME,sep=";",header=header)
 return (t[,1])
}

guess_k=function(x,Max=300) # get a reasonable value for first sazonality period or DEFAULT
{
 DEFAULT=10
 L=length(x)
 ULIM=1.96*1/sqrt(L)
 Max=min(L,Max)
 A=acf(x,lag.max=Max)$acf
 
 I=which(A<ULIM)
 if(length(I)>1)
 {  Min=I[1]
    A=A[Min:length(A)]
    I=which.max(A)
  #cat(" I",I," A",A[i],"Def:",DEFAULT,"\n")
    if(A[I]>ULIM) return (Min+I-2)
    else return (DEFAULT)
 }
 else return (DEFAULT) 
}


# nominal : currently this function is not needed anymore, as nnet, ksvm and lm use
#           already an internal similar process. Yet I will keep this function here
#           in case I need it in the future!
#---------------------------------------------------------
# Performs a 1-of-C coding on the nominal (factor with 3 or more discrete values)
# DMData - matrix or dataframe with the DM dataset
# thresh - maximum number of labels (C) for the 1-of-C encoding?
# positive - label for the positive value (e.g. 1)
# negative - label for the negative value (e.g. -1 or 0)
# warning - only works with dataframes with 2 or more attributes... need to correct this... 
#---------------------------------------------------------
nominal <-function(DMData, thresh = 20,positive=1,negative=-1,exclude=-1) 
{
   AttrNo <- NCOL(DMData)
   result <- DMData
   k<-1
   for (i in 1:AttrNo) 
   {
       if (i!=exclude && is.factor(DMData[ , i]))
       {
        L<-length(levels(DMData[,i]))
        #print(paste("i: ",i," k:", k, " L:", L))
        if (L>2) # is nominal!!!
        {
          result <- substituteNV(result, k, thresh,positive,negative)
          if(L<thresh) k<-k+L
          else k<-k+1
        }
        else { 
               BIN<-as.numeric(DMData[,i])
               BIN[BIN==1]<-negative  
               BIN[BIN==2]<-positive
               result[k]<-BIN
               k<-k+1
             }
       } 
       else k<-k+1
   } 
   return (result)       
}
# ------------------------------------------------------------------
# internal R function used by nominal: do not use this
# by Milan Legat and Martin Gruber 2005 (very fast):
# by Paulo Cortez 2006, some minor corrections and adaptations
# ------------------------------------------------------------------
"substituteNV" <-
function(DMData, colNo, thresh = 20, positive=1,negative=-1) {
         
             NAMES<-names(DMData)
             
             result <- DMData[ , ]
             rowCount <- nrow(DMData)
             colF <- as.factor(DMData[ , colNo])
             
             levelsNo <- nlevels(colF)
             origName <- names(DMData)[colNo]

             if (levelsNo < thresh) 
              {
                subs <- data.frame(matrix(negative, nrow = rowCount, ncol = levelsNo))
                for (i in 1:levelsNo) names(subs)[i] <- paste(origName, '_', levels(colF)[i], sep = '')
                for (i in 1:levelsNo) 
                   {
                    replaceWhat <- which(DMData[ , colNo] == levels(colF)[i])
                    if(length(replaceWhat>0)) subs[replaceWhat, i] <- positive 
                   } 
              } # if (levelsNo < thresh) {
             else 
              {
                subs <- data.frame(matrix(0, nrow = rowCount, ncol = 1))
                names(subs) <- paste(origName, '_', sep = '')
                for (i in 1:levelsNo) 
                  {
                    replaceWhat <- which(DMData[ , colNo] == levels(colF)[i])
                    if(length(replaceWhat>0))subs[replaceWhat, ] <- i
                  } 
              } 
             
             Xcol<-NCOL(result)
             if(colNo==1) { 
                            result <- cbind(subs, result[,(2:Xcol)])                          
                            NAMES  <- c(names(subs),NAMES[(2:Xcol)])
                          }
             else if(colNo==Xcol) { result <- cbind(result[,(1:(colNo-1))], subs)                          
                                    NAMES  <- c(NAMES[(1:(colNo-1))],names(subs))
                                  }
             else 
               { 
                result <- cbind(result[,(1:(colNo-1))], subs, result[,((colNo+1):Xcol)])
                NAMES  <- c(NAMES[(1:(colNo-1))],names(subs),NAMES[((colNo+1):Xcol)])
               }
             names(result)<-NAMES
             return (result)
             
} # substituteNV <- function(column) {

}
