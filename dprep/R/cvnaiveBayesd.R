cvnaiveBayesd <-
function(data, repet,method=c("ew","ef","1R","chiMerge"))
{
#This function estimates the misclassification error rate
#for the naive Bayes classifier using  10-fold vrossvalidation
# The dataset is discretized using ether the equal width
# or the ChiMerge method 
# It requieres the naiveBayes function from the e1071 library 
#inputs:
# data: the name of the dataset
# repet: the number of rep
 #------------------------------------------
#require(e1071)
 n=dim(data)[1]
p= dim(data)[2]
if(method=="ew")
{data=disc.ew(data,1:(p-1),out="num")}
if(method=="ef")
{data=disc.ef(data,1:(p-1),10, out="num")}
if(method=="1R")
{data=disc.1r(data,1:(p-1),out="num")}
 if(method=="chiMerge")
 {data=chiMerge(data,1:(p-1),out="num")}
for(i in 1:p)
 {data[,i]=as.factor(data[,i])}
 nombres<-colnames(data)
f1=as.formula(paste(nombres[p],".",sep="~"))
#print(f1)
ECV <- rep(0, repet)
for(i in 1:repet) {
salida <- matrix(0, 1, 10)
azar <- data[rank(runif(n)),  ]
azar[, p] <- as.factor(azar[, p])
parti <- floor(n/10)
for(j in 1:10) {
cc= ((j - 1) * parti + 1):(j * parti
 )
if(j == 10) {
cc <- ((j - 1) * parti + 1):
n
 }
 datap <- azar[cc,  ]
 datat <- azar[ - cc,  ]
tempo =e1071::naiveBayes(f1,data=datat)
 tempo1=predict(tempo,datap[,-p])
 salida[j] <- sum(tempo1 != as.numeric(datap[, p]))
}
ECV[i] <- sum(salida)/n
 }
cat("The error estimations in each repetition are:\n"
 )
 print(ECV)
 ECV1 <- mean(ECV)
 cat("The mean error estimation by cross-validation using all the repetitions is: \n"
 )
 ECV1
 }
