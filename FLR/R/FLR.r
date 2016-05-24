library(combinat)


##################################
#                                #
#        Data Preprocessing      #
#                                #
##################################
prepData<-function(data){
dataT<-data
if (class(data)=='list'){

data<-dataT[[1]]
}

if(names(data)[length(data)]!="class"){
stop("Last attribute MUST be the class attribute and be labeled as 'class'!")
}#if

data$class<-as.factor(data$class)
ch<-c()

# Removes instances with missing class.
for(i in 1:length(data[[1]])){
	if (data[i,length(data)]!="?"){ 
	
	ch<-c(ch,i)}#if
}#for i


data<-data[ch,]


# Converts nominal attributes into numeric.
for(i in 1:(length(data)-1)){
data[,i]<-as.numeric(data[,i])
}#for i

# data1 is data.frame without missing class instances
data1<-data.frame(data[[1]])

for (i in 2:length(data)){
	data1<-data.frame(data1,data[[i]])
}#for i

names(data1)<-names(data)

data1<-data.frame(data1,categ=0)

#Adds categ attribute.
for(i in 1:length(data1[[1]])){
	for(j in 1:length(levels(data1$class))){
		if (data1$class[i]==levels(data1$class)[j]){
				data1$categ[i]<-j

		}#if
	}#for j
}#for i

data1$class<-as.factor(as.matrix(data1$class)) #This removes '?' from levels(data1$class)
if (levels(data1[, length(data1)-1])[1]=="?"){
data1[,length(data1)]<-data1[, length(data1)]-1
}
if (class(dataT)=='list'){
dataT[[1]]<-data1}
else{
	dataT<-data1}
return(dataT)
}


#Linear arrangement of spatial data.

spatdt<-function(data,idx,mat,pre_order=0,snd=0){
if (length(pre_order)==length(mat)){
ord<-c()
for(i in 1:length(pre_order)){
ord[i]<-levels(data[,idx])[pre_order[i]]
}

winrt<-as.factor(ord)
wincost<-get.cost2(winrt,mat)
}#if pre_order==length(mat)
else if (pre_order==0){
z<-levels(data[,idx])

if (length(z)!=length(mat)){
stop("Mat length does not match the levels of the spatial attribute you have chosen (idx).")
}#if

zzz<-c()
zz<-permn(z)
k<-1

#When snd!=0 
if(snd!=0){

 for(i in 1:length(zz)){
	if(zz[[i]][1]==z[snd]){
		zzz[k]<-zz[i]
		k<-k+1}#if zz[[i]][1]==z[snd]
		}#for
	}#if snd!=0

if(snd==0){
 zzz<-zz
 }
 

zzz<-as.data.frame(zzz)
names(zzz)<-1:length(zzz)
cost<-get.cost(zzz,mat)
winrt<-zzz[,winner.route(cost)]
wincost<-cost[winner.route(cost)]
}#else if pre_order==0
else{
 stop("pre_order must be either 0 or a vector of dimensions equal to the length of mat")
}
dat<-as.matrix(data[,idx])

for(i in 1:length(dat)){
	for(j in 1:length(winrt)){
		if(dat[i]==winrt[j]){
			dat[i]<-as.numeric(j)
		}#if
	}#for j
}#for i

data[,idx]<-dat
data<-list("data"=data,"win_route"=winrt,"cost"=wincost)
return(data)
}

# Normalizes Data in the range of [0,1]
normData<-function(data1){
dataT<-data1
if (class(data1)=='list'){

data1<-dataT[[1]]
}

bounds<-set_bounds(data1)
for(i in 1:length(data1[[1]])){
for(j in 1:(length(data1)-2)){
data1[i,j]<-(data1[i,j]-bounds[j,1])/(bounds[j,2]-bounds[j,1])
	}
}
if (class(dataT)=='list'){
dataT[[1]]<-data1}
else{ dataT<-data1 }
return(dataT)
}

# Dernormalize Data
denormData<-function(data1,bounds){
for(i in 1:length(data1[[1]])){
for(j in 1:(length(data1)-2)){
data1[i,j]<-(data1[i,j]*(bounds[j,2]-bounds[j,1]))+bounds[j,1]
	}
}

return(data1)
}

# Denormalize Fuzzy lattices

denormDatal<-function(fuzlat,bounds){
for(i in 1:length(fuzlat[[1]][[1]])){
for(j in 1:((length(fuzlat[[1]])/2)-1)){
fuzlat[[1]][i,2*j-1]<-(fuzlat[[1]][i,2*j-1]*(bounds[j,2]-bounds[j,1]))+bounds[j,1]
fuzlat[[1]][i,2*j]<-(fuzlat[[1]][i,2*j]*(bounds[j,2]-bounds[j,1]))+bounds[j,1]
}
}

return(fuzlat)
}


# Creates bounds file.
set_bounds<-function(data1){

res<-data.frame(min=rep(0,times=length(data1)-2),max=rep(0,times=length(data1)-2))
for (i in 1:(length(data1)-2)){
res[[1]][i]<-min(data1[[i]])
res[[2]][i]<-max(data1[[i]])
}

return(res)
}

#Flags instances to be used as training(0) or testing(1) data with the ratio depending on variable gg.
sepFlag<-function(gg,data1){
testSample<-sort(sample(1:length(data1[[1]]),as.integer(length(data1[[1]])*gg),replace=F))
data2<-data.frame(data1,flag=0)
data2[testSample,length(data2)]<-1
return(data2)
}

#Returns testing data.
testD<-function(data2){
data2<-data2[data2$flag==1,1:length(data2)-1]
return(data2)
}

#Returns training data.
trainD<-function(data2){
data2<-data2[data2$flag==0,1:length(data2)-1]
return(data2)
}


# Constructs a Fuzzy Lattice from an instance.
fuzzyLatticec<-function(dF,dR,bounds){
k<-1

for (i in 1:(length(dR)-2)){
#length(dR) is dR.classIndex()
	if(i!=length(dR)){
		if(dR[i]!="?"){
			dF[2*k-1]<-max(dR[i],bounds$min[k])
			dF[2*k]<-min(dR[i],bounds$max[k])
			k<-k+1
			}#if
	else {
		dF[2*k-1]<-bounds$min[k]
		dF[2*k]<-bounds$max[k]
		k<-k+1
		}#else
	}#if
}#for
	dF$categ<-dR$categ
	dF$class<-dR$class
	dF<-as.data.frame(dF)
	return(dF)
	}	


# Calculates the valuation function of the Fuzzylattice.V2
valuation<-function(fuzlat,x0,l,param){
k<-0
for (i in 1:((length(fuzlat)-2)/2)){
k<-k+ufun(theta(as.numeric(fuzlat[2*i-1]),x0,param),x0,l,param)+ufun(as.numeric(fuzlat[2*i]),x0,l,param)
}
return(k)
}	


# Implements theta function.
theta<-function(x,x0,param){
res<-0
if (param=='linear'){
res<-(1-x)
}#if
if (param=='sigmoid'){
res<-2*x0-x
}#if
return(res)
}

	
# Implements u function.
ufun<-function(x,x0,l,param){
res<-0
if(param=='linear'){
res<-x
}#if
if(param=='sigmoid'){
res<-1/(1+exp(-l*(x-x0)))
}#if
return(res)
}



# Implements the Join Function.
join<-function(inpBuf,num){

jn<-inpBuf

for(i in 1:((length(inpBuf)-2)/2)){
jn[2*i-1]<-min(inpBuf[2*i-1],num[2*i-1])
jn[2*i]<-max(inpBuf[2*i],num[2*i])
}#for

jn$categ<-inpBuf$categ
jn$class<-inpBuf$class
return(jn)
}

# Creates namesList parameter

createNlist<-function(trainData){
ramen<-names(trainData)
ramenMax<-0
ramenMin<-0
for(i in 1:length(ramen)-2){
ramenMax[i]<-paste(ramen[i],"Max",sep="")
ramenMin[i]<-paste(ramen[i],"Min",sep="")
}#for

ramenList<-c(ramenMin[1],ramenMax[1])
for(i in 2:(length(ramen)-2)){
ramenList[2*i-1]<-ramenMin[i]
ramenList[2*i]<-ramenMax[i]
}#for

ramenList[length(ramenList)+1]<-"class"
ramenList[length(ramenList)+1]<-"categ"


nameslist<-list(0)
length(nameslist)<-length(ramenList)
names(nameslist)<-ramenList
return(nameslist)
}

# Create namesframe parameter
createNframe<-function(trainData){
ramen<-names(trainData)
ramenMax<-0
ramenMin<-0
for(i in 1:length(ramen)-2){
ramenMax[i]<-paste(ramen[i],"Max",sep="")
ramenMin[i]<-paste(ramen[i],"Min",sep="")
}#for

ramenList<-c(ramenMin[1],ramenMax[1])
for(i in 2:(length(ramen)-2)){
ramenList[2*i-1]<-ramenMin[i]
ramenList[2*i]<-ramenMax[i]
}#for

ramenList[length(ramenList)+1]<-"class"
ramenList[length(ramenList)+1]<-"categ"


namesframe<-data.frame(0)
length(namesframe)<-length(ramenList)
names(namesframe)<-ramenList
return(namesframe)
}

##################################
#                                #
#         Training phase         #
#                                #
##################################

trainNow<-function(trainData,param,rhoa=0.5,l=6,x0=0.5,EPSILON=10^(-6)){
boundss<-set_bounds(trainData)
parameters=data.frame("param"=param,"rhoa"=rhoa,"x0"=x0,"l"=l,"EPSILON"=EPSILON)
# Set the first instance to be the first Rule in the model.
learnedCode<-createNframe(trainData)
lat<-createNlist(trainData)
lrnC<-fuzzyLatticec(lat,trainData[1,],boundss)
learnedCode<-as.data.frame(lrnC)


searching<-1

# Training iteration.
for (i in 2:length(trainData[[1]])){ 	# for all instances

	inst<-trainData[i,]
	flag<-0
	for(w in 1:(length(trainData)-2)){
		if(w!=length(trainData)&inst[w]=="?"){
			flag<-flag+1
		}#if
	}#for w



	if(flag!=(length(trainData)-1)){
		inputBuffer<-fuzzyLatticec(lat,inst,boundss)
		sigma<-rep(0,times=length(learnedCode[[1]]))
		names(inputBuffer)<-names(learnedCode)
	
	
		for(j in 1:length(learnedCode[[1]])){
		
			num<-learnedCode[j,]
			den<-join(inputBuffer,num)
			numden<-valuation(num,x0,l,param)/valuation(den,x0,l,param)
			sigma[j]<-numden	
		}#for j
	searching<-1
		while (searching==1){
		winner<-1
		winnerf<-sigma[1]
		if (length(learnedCode[[1]])>1){
		for( j in 2:length(learnedCode[[1]])){
			if( winnerf<sigma[j] ){
			winner<-j
			winnerf<-sigma[j]
				}#if winnerf<sigma[j]
			}#for j
		}#if length(learnedCode[[1]])>1
		num<-inputBuffer;
		winnerBox<-learnedCode[winner,]
		den<-join(winnerBox,num)
		numden<-valuation(num,x0,l,param)/valuation(den,x0,l,param);
		
		
		if ((inputBuffer$categ==winnerBox$categ)&(rhoa<numden)){
		
		learnedCode[winner,]<-join(inputBuffer,winnerBox)
            searching <-0
		}#if (inputBuffer$categ==winnerBox$categ)&(rhoa<numden)
		else{
		sigma[winner]<-0
		rhoa<-rhoa+EPSILON
		searching<-0
		for(j in 1:length(learnedCode[[1]])){
			if (sigma[j]!=0)
			{searching<-1}#if sigma[j]!=0
			}#for j
		
		
		if (searching==0){
			
              learnedCode[(length(learnedCode[[1]])+1),]<-inputBuffer;
              
		}#if searching==0
		
		}#else
		
		}#while
	}#if flag!=(length(trainData)-1)
}#for i

learnedCode<-list("rules"=learnedCode,"parameters"=parameters)
return(learnedCode)
}


##################################
#                                #
#         Testing phase          #
#                                #
##################################

testNow<-function(testData,learnedCode){
testData1<-testData
parameters<-learnedCode[[2]]
param<-parameters$param
x0<-parameters$x0
l<-parameters$l
learnedCode<-learnedCode[[1]]


boundss<-set_bounds(testData)
testData$categ<-0
nameslist<-createNlist(testData)

for(i in 1:length(testData[[1]])){
inst<-testData[i,]
inputBuffer<-fuzzyLatticec(nameslist,inst,boundss)
sigma<-rep(0,times=length(learnedCode[[1]]));

for(j in 1:length(learnedCode[[1]])){
num<-learnedCode[j,];
den<-join(inputBuffer,num);
sigma[j]<-valuation(num,x0,l,param)/valuation(den,x0,l,param);
}#for j

winner<-1;
winnerf<-sigma[1];

for(j in 2:length(learnedCode[[1]])){
	if(winnerf<sigma[j]){
	winner<-j;
	winnerf<-sigma[j];
	}#if winnerf<sigma[j]
}#for j

currentBox<-learnedCode[winner,];
testData$categ[i]<-currentBox$categ
testData$class[i]<-currentBox$class
}#for i
acc<-accIs(testData1,testData)
confusion_mat<-table(testData1$class,testData$class)
testData<-list("Original test data"=testData1,"Predicted test data"=testData,"accuracy"=acc,"confusion matrix"=confusion_mat )
return(testData)
}



#Returns a vector that contains the number of rules created for each class.
indexCalc<-function(learnedCode){

indexx<-rep(0,times=max(learnedCode[,length(learnedCode)]))
for(i in 1:length(learnedCode[[1]])){
indexx[learnedCode[i,]$categ]<-indexx[learnedCode[i,]$categ]+1
}
return(indexx)
}


#Calculates the accuracy of the predictions made.
accIs<-function(testData,testDataB){
minato<-0;
for(i in 1:length(testData[[1]])){

if (testData$categ[i]==testDataB$categ[i]){
	minato<-minato+1;
	}#if
}#for i
acc<-(minato/length(testData[[1]]))
return(acc)
}

#Used in spatdt function.
get.pos<-function(instz){
p<-levels(instz)
ord<-c(1:length(instz))
for(i in 1:length(instz)){
	for(j in 1:length(p)){
		if( instz[i]==p[j]){
			ord[i]<-j

		}	#if
	}#for j
}#for i
return(ord)
}


#Used in spatdt function.
get.cost<-function(zzz,mat){
cost<-rep(0,times=length(zzz))
for(i in 1:length(zzz)){
		instz<-zzz[,i]
		ord<-get.pos(instz)
	for(j in 1:(length(ord)-1)){
		cost[i]<-cost[i]+mat[ord[j],ord[j+1]]
	}#for j	
}#for i
return(cost)
}


#Used in spatdt function(when order is definer by user).
get.cost2<-function(pre_order,mat){

cost<-0
	

		ord<-get.pos(pre_order)
	for(j in 1:(length(ord)-1)){
		cost<-cost+mat[ord[j],ord[j+1]]
	}#for j	

return(cost)
}


#Used in spatdt function.
winner.route<-function(cost){
for(i in 1:length(cost)){

	if(cost[i]==min(cost)){
		winner1<-i
	}#if
}#for i
return(winner1)
}