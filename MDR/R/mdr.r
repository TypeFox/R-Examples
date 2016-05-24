############ Code for Rmdr function; modify this to incorporate 3-way split

require(lattice)

compare<-function(mat,vec,k) {  #this function is for internal use only; it matches an individuals data vector with a target
b<-1
match<- 1:dim(mat)[1]
while (b <= (k+1)) {
match<- match[mat[match,b]==as.numeric(vec[b])]
b<-b+1
} #close while loop
return(length(match))
}


###### function to use MDR, calculate balanced accuracy, save top x models

mdr<-function(split,comb,x,ratio,equal="HR",genotype=c(0,1,2)) { #split is the split of the data, comb is the matrix of combinations to consider (row=comb, col=loci)
          # x is how many models to save, ratio is threshold to compare each cell, equal is how to define cells where cases=ctrls ("HR or LR")
          # includes balanced accuracy and missing data (where missing is encoded with NA)
if (is.vector(comb)) comb<-as.matrix(comb)  #need to make sure its read as a matrix
    x.new<-ifelse(x>dim(comb)[1], dim(comb)[1], x) #if x is bigger than the number of combinations, reduce to the number of combinations
    N<-dim(split)[1]
    n<-dim(split)[2]
    n.comb<-dim(comb)[1]
    k<-dim(comb)[2]  

    results <- matrix(0, n.comb, k+1)  #matrix to store results (classification)

    results[, 1:k] <- comb  #first nbr columns (model size) replaced with all possible combinations
    g<-length(genotype)

    geno <- list(genotype) # count three genotypes
    case <- cbind(rep(1, g^k), expand.grid(rep(geno, k))) #creates design matrix with all combinations (contingency for cases)
    ctrl <- cbind(rep(0, g^k), expand.grid(rep(geno, k))) #this acts as contingency table for controls
    counts <- matrix(0, dim(case)[1], 3)  #store number of cases, controls, ratio for all combinations

       snps <- split[, 2:n]                               #snps for training 
       hr.lr<-matrix(0,n.comb,dim(case)[1])   # matrix to store hr/lr combination (1=high, 0=low)

        for (i in 1:n.comb) {  #loop over all genotype combinations

            if (is.vector(comb)) {model<-comb[i]} else {model<-comb[i,]}

part<-split[, c(1, model + 1)]
miss<-is.na(part) #which entries are missing?
miss2<-apply(miss,1,sum) #which rows have missings? (sums > 0?)
if (sum(miss2)>0  ) {part<-part[-which(miss2>0),]} #eliminate the rows which have missings (only if there is missing data)
  
            counts[,1]<-apply(case,1,compare,mat=part,k=k)   #which rows match (cases with given combination)?
	                               # gives number of cases for each genotype combination 
  
            counts[,2]<-apply(ctrl,1,compare,mat=part,k=k)  #which rows match (controls with given combination)?
	                                # gives number of controls for each genotype combination

            counts[, 3] <- counts[, 1]/counts[, 2]   #observed ratio
	if (equal=="HR") {hr.lr[i,]<-(counts[, 3] >= ratio)} else { if (equal=="LR") { 
			hr.lr[i,]<-(counts[, 3] > ratio)} else {print("Improper value for argument 'equal' ")}}  #define which combinations are hr and lr 
		tp<- sum(counts[which(hr.lr[i,]==1), 1]) #true positives
		fp<- sum(counts[which(hr.lr[i,]==1), 2]) #false positives
		tn<- sum(counts[which(hr.lr[i,]==0), 2]) #true negatives
		fn<- sum(counts[which(hr.lr[i,]==0), 1]) #false negatives

            results[i, k+1] <- 100*( tp/(tp+fn) + tn/(tn+fp) ) / 2
                   										 # maximize balanced accuracy           
        }
    sorted<-order(results[,k+1], decreasing=TRUE) #want maximum
    best.models<-results[sorted[1:x.new],]
    best.hrlr<-hr.lr[sorted[1:x.new],]
    if (x.new>1) {models<-best.models[,1:k]} else {models<-best.models[1:k]}
    if (x.new>1) {acc<-best.models[,k+1]} else {acc<-best.models[k+1]}

    return(list("models"=models, "balanced accuracy"=acc, "high risk/low risk"=best.hrlr)  )

}   #end mdr function

mdr.hr<-function(split,model,hr,genotype=c(0,1,2)) { #split is the split of the data, model combination of loci to consider (row=comb, col=loci), 
          # hr is an indicator for which cmb are high risk
    N<-dim(split)[1]
    n<-dim(split)[2]
    k<-length(model)
    g<-length(genotype)

    geno <- list(genotype) # count three genotypes
    case <- cbind(rep(1, g^k), expand.grid(rep(geno, k))) #creates design matrix with all combinations (contingency for cases)
    ctrl <- cbind(rep(0, g^k), expand.grid(rep(geno, k))) #this acts as contingency table for controls
    counts <- matrix(0, g^k, 2)  #store number of cases, controls

       snps <- split[, 2:n]                               #snps for training 

part<-split[, c(1, model + 1)]
miss<-is.na(part) #which entries are missing?
miss2<-apply(miss,1,sum) #which rows have missings? (sums > 0?)
if (sum(miss2)>0  ) {part<-part[-which(miss2>0),]} #eliminate the rows which have missings (only if there is missing data)
  
            counts[,1]<-apply(case,1,compare,mat=part,k=k)   #which rows match (cases with given combination)?
	                               # gives number of cases for each genotype combination 
  
            counts[,2]<-apply(ctrl,1,compare,mat=part,k=k)  #which rows match (controls with given combination)?
	                                # gives number of controls for each genotype combination
		
		tp<- sum(counts[which(hr==1), 1]) #true positives
		fp<- sum(counts[which(hr==1), 2]) #false positives
		tn<- sum(counts[which(hr==0), 2]) #true negatives
		fn<- sum(counts[which(hr==0), 1]) #false negatives

            results<- 100*( tp/(tp+fn) + tn/(tn+fp) ) / 2      								      
    return(list("balanced accuracy"=results)  )
}   #end mdr.hr function


mdr.cv<-function (data, K, cv, ratio = NULL, equal="HR",genotype=c(0,1,2)) #K is highest level of interaction to consider, cv is number of cv intervals
{
    if (is.null(ratio)) 
	 ratio <- length(which(data[,1] == 1))/length(which(data[,1] == 0))   #set threshold to be ratio of cases to controls in whole data      
    N <- dim(data)[1]   #number of individuals
    p <- dim(data)[2]   # number of snps + 1
    n <- N %/% cv  # number of individuals per split (minus remainder)
    r<- N-n*cv   #remainder
 
        data <- data[sample(1:N), ]  #randomly order the data for the cv split

###### store results
models<-vector("list",K)
acc<-matrix(0,K,3) #save ca, pa, cvc

######## Perform MDR for first split for 1-nbr levels of interaction (model building)
for (i in 1:K) {

best.mod<-matrix(0,cv,i)
ca<-rep(0,cv)
pa<-rep(0,cv)

######### start cross validation  
for (k in 1:cv) {



#need to split the data
if (k==cv) { build<- data[-((1+(k-1)*n):(k*n+r)),] } else { build<- data[-((1+(k-1)*n):(k*n)),]
                                                            } #make sure we put the remainder individuals in a prediction interval
	geno <- combn(p - 1, i)  #list of all possible combinations
        res<-mdr(build,t(geno),1,ratio,equal,genotype) #pick only the best model
        best.mod[k,]<-res$models
        ca[k]<-res$bal    #put best models and accuracies into the list to store (by cv interval, then model size)
if (k==cv) {predict<- data[(1+(k-1)*n):(k*n+r),]} else { predict<- data[(1+(k-1)*n):(k*n),] }

        res<-mdr.hr(predict,t(best.mod[k,]),res$high,genotype)  #get prediction accuracies only for this list
        pa[k]<-res$bal  #put best models into the list to store 

} #close cv loop

#get avg ca, pa, and cvc for this level and pick best model
unique.results<-as.matrix(unique(best.mod[,1:i]))
cvc<-rep(0,nrow(unique.results))

for (j in 1:length(cvc)) {   #count up cross-val consistency
for (k in 1:cv) {
if (sum(best.mod[k,]==unique.results[j,])==i) {
    cvc[j]<-cvc[j]+1 }
} }

if (i>1) { best<-matrix(unique.results[which(cvc==max(cvc)),], ncol=i) } else { best<-as.matrix(unique.results[which(cvc==max(cvc)),]) } #make sure we have correct dimensions
  #what about ties in cvc?

if (nrow(best)==1) {  sum.ca <- 0
  for(k in 1:cv){  #loop across cv intervals
    
     if (sum(best.mod[k,]==best)==i) {    
      sum.ca <- sum.ca + ca[k]  #recursion to add up classification accuracy
    } }  #close if and for loop
  
avg.ca <- sum.ca/max(cvc) } else { b<-dim(best)[1]

sum.ca <- rep(0,b)
  
for(k in 1:cv){  #loop across cv intervals
for(l in 1:b) { #loop over high cvc models    
     if (sum(best.mod[k,]==best[l,])==i) {    
      sum.ca[l] <- sum.ca[l] + ca[k]  #recursion to add up classification accuracy
    } }  #close if and for loop
}
avg.ca <- sum.ca/max(cvc) 
best<-best[which(avg.ca==max(avg.ca)),]
avg.ca<-max(avg.ca)
}

models[[i]]<-best  #chose model with max cvc and max classification accuracy

  ### best.model.PA   #same idea, find average prediction accuracy of best model
    sum.pa <- 0
  for(k in 1:cv){  #loop across cv intervals
    if (sum(best.mod[k,]==best)==i) {    
      sum.pa <- sum.pa + pa[k]  #recursion to add up classification error
    }
  }

  avg.pa <- sum.pa/max(cvc)  #average prediction error

acc[i,1]<-avg.ca
acc[i,2]<-avg.pa
acc[i,3]<-max(cvc)
colnames(acc)<-list("classification accuracy", "prediction accuracy", "cross-validation consistency")

} #close K loop


final.acc<-which(acc[,2]==max(acc[,2]))	#final model maximizes p accuracy across levels
final.cvc<-which(acc[,3]==max(acc[,3]))	#final model maximizes cvc across levels

if (length(final.cvc)>1) { if (length(which(final.cvc==final.acc))>0) {
final.cvc<-final.cvc[which(final.cvc==final.acc)] } } #in case some models have the same CVC, use the one with maximum acc

final<-min(c(final.acc, final.cvc))  #pick index which is smallest of the two, if not equal (parsimony)

final.model<-models[[final]]
hrlr<-mdr(data,final.model,1,ratio) #get hr/lr of final model

obj<-list("final model"=final.model, "final model accuracy"= acc[final,2], "top models"=models, "top model accuracies"=acc, "high-risk/low-risk"=hrlr$high, "genotypes"=genotype, "validation method"="CV")
class(obj)<-"mdr" #define mdr class
return(obj)

}  #end cv function


mdr.3WS<-function (data, K, x=NULL, proportion=NULL, ratio = NULL, equal="HR",genotype=c(0,1,2)) #x is threshold to send to testing set, K is highest level of interaction to consider, split is the ratio training:testing:validation (default is 2:2:1)
{
    if (is.null(ratio)) 
        ratio <- length(which(data[,1] == 1))/length(which(data[,1] == 0))   #set threshold to be ratio of cases to controls in whole data
    N <- dim(data)[1]   #number of individuals
    p <- dim(data)[2]   # number of snps + 1
    if (is.null(proportion)) 
        proportion <- c(2,2,1)   #set default proportions to 2:2:1
    if (is.null(x))     
        x<-p-1 #default is the number of snps
        data <- data[sample(1:N), ] #randomly order the data for the three-way split
    prop<- sum(proportion) #number of proportions for the split
    n <- N%/%prop  # number of individuals per proportion (minus remainder)
#need to split the data
   test<-data[1:(proportion[2]*n),] #second split
   val<-data[(proportion[2]*n+1):((proportion[2]+proportion[3])*n),] #third split
   train<-data[((proportion[2]+proportion[3])*n+1):N,]  #first split; want training set to be largest so add remainder

################## Store results
models<-vector("list",K) #store final model for each level of K
accuracy<-matrix(0,K,3) #save training, testing, and prediction accuracy for each level of K

######## Perform MDR for first split for 1-nbr levels of interaction (training)

results1<-vector("list",K)  #create matrices to store results for each level

for (i in 1:K) {
	geno <- combn(p - 1, i)  #list of all possible combinations
        x.new<-ifelse(x>dim(geno)[2], dim(geno)[2], x) #if x is bigger than the number of combinations, reduce to the number of combinations
        res<-mdr(train,t(geno),x.new,ratio,equal,genotype)
        results1[[i]]<-res }  #put best models into the list to store

######## Perform MDR for second split for 1-nbr levels of interaction (testing)

results2<-vector("list",K)  #create matrices to store results for each level

for (i in 1:K) {
	comb<-results1[[i]]$models  #extract models of interest
        res<-mdr(test,comb,1,ratio,equal,genotype) #only save the best model for each level
        results2[[i]]<-res }  #put best models into the list to store

####### Perform MDR for third split and return model + prediction error
results3<-vector("list",K)  #create matrices to store results for each level
acc<-NULL
for (i in 1:K) {
	comb<-t(as.matrix(results2[[i]]$models))  #extract models of interest; make sure it reads it as matrix and not vector
        res<-mdr(val,comb,1,ratio,equal,genotype) #only save the best model for each level
        acc<-c(acc,res$bal)
        results3[[i]]<-res }

final<-which(acc==max(acc))	#final model maximizes accuracy across levels (minimizes prediction error)
final<-min(final) #choose model of smallest size (parsimony) if more than one model has the minimum prediction error

final.model<-results3[[final]]$models
BA<-results3[[final]]$bal
hrlr<-mdr(data,t(final.model),1,ratio) #get hr/lr of final model

###get final model for each level K and accuracy
for (i in 1:K) {
models[[i]]<-results3[[i]]$models
index<-j<-0
while (index<1) {
j<-j+1
if (sum(as.matrix(results1[[i]]$models)[j,]==results3[[i]]$models)==i) index<-j
}

accuracy[i,1]<-results1[[i]]$'balanced accuracy'[index] #training accuracy
accuracy[i,2]<-results2[[i]]$'balanced accuracy' #testing accuracy
accuracy[i,3]<-results3[[i]]$'balanced accuracy' #prediction accuracy
}

colnames(accuracy)<-list("training accuracy", "testing accuracy", "prediction accuracy")
obj<-list("final model"=final.model, "final model accuracy"=BA, "top models"=models, "top model accuracies"=accuracy, "high-risk/low-risk"=hrlr$high, "genotypes"=genotype, "validation method"="3WS")
class(obj)<-"mdr" #define mdr class
return(obj)

}  #end 3way function

###############################
# post-hoc functions
###############################

mdr.ca.adj<-function(data,model,hr,prev,genotype=c(0,1,2)) { #model combination of loci to consider (row=comb, col=loci) read as vector, 
          # hr is an indicator for which cmb are high risk, prev is prevalence estimate
    N<-dim(data)[1]
    n<-dim(data)[2]
    k<-length(model)
    g<-length(genotype)

p.case <- length(which(data[,1] == 1))/ N   #estimated probability of being a case in the sample

    geno <- list(genotype) # count three genotypes
    case <- cbind(rep(1, g^k), expand.grid(rep(geno, k))) #creates design matrix with all combinations (contingency for cases)
    ctrl <- cbind(rep(0, g^k), expand.grid(rep(geno, k))) #this acts as contingency table for controls
    counts <- matrix(0, g^k, 2)  #store number of cases, controls

       snps <- data[, 2:n]                               #snps for training 

part<-data[, c(1, model + 1)]
miss<-is.na(part) #which entries are missing?
miss2<-apply(miss,1,sum) #which rows have missings? (sums > 0?)
if (sum(miss2)>0  ) {part<-part[-which(miss2>0),]} #eliminate the rows which have missings (only if there is missing data)
  
            counts[,1]<-apply(case,1,compare,mat=part,k=k)   #which rows match (cases with given combination)?
	                               # gives number of cases for each genotype combination 
  
            counts[,2]<-apply(ctrl,1,compare,mat=part,k=k)  #which rows match (controls with given combination)?
	                                # gives number of controls for each genotype combination
		
		tp<- sum(counts[which(hr==1), 1]) #true positives
		fp<- sum(counts[which(hr==1), 2]) #false positives
		tn<- sum(counts[which(hr==0), 2]) #true negatives
		fn<- sum(counts[which(hr==0), 1]) #false negatives

            results<- 100-(100/N)*( ((1-prev)/(1-p.case))*fp + (prev/p.case)*fn )      								      
    return(list("adjusted classification accuracy"=results, "adjusted classification error"=100-results)  )
}   #end mdr.ca.adj function; this does not implement CV

boot.error<-function(data, prev, model, hr, b, genotype=c(0,1,2)) {# inputs are case/control dataset, prevalence estimate, and the final mdr model loci (as a vector), and high-risk/low-risk, b is number of bootstrap samples

data<-data[, c(1, model + 1)] #reduce data to only the model of interest
miss<-is.na(data) #which entries are missing?
miss2<-apply(miss,1,sum) #which rows have missings? (sums > 0?)
if (sum(miss2)>0  ) {data<-data[-which(miss2>0),]} #eliminate the rows which have missings (only if there is missing data)

n<-dim(data)[1]

cases<-data[which(data[,1]==1),]
ctrls<-data[which(data[,1]==0),]
n.case<- as.integer(round(prev*n,0))
n.ctrl<- n-n.case
boot.error<-rep(0,b) #vector to save error estimates for each bootstrap sample

for (i in 1:b) { #loop across bootstrap samples
index.cases<-sample(seq(1:(dim(cases)[1])),n.case,replace=TRUE)
index.ctrls<-sample(seq(1:(dim(ctrls)[1])),n.ctrl,replace=TRUE)
boot.cases<-cases[index.cases,]  #this is bootstrap resample
boot.ctrls<-ctrls[index.ctrls,]

k<-length(model)
g<-length(genotype)

geno <- list(genotype) # define possible genotypes
x.case <- cbind(rep(1, g^k), expand.grid(rep(geno, k))) #creates design matrix with all combinations (contingency for cases)
x.ctrl <- cbind(rep(0, g^k), expand.grid(rep(geno, k))) #this acts as contingency table for controls
n.comb <- dim(x.case)[1]  # number of genotype combinations

counts<-matrix(0,n.comb,3) #store number of cases, ctrls, and sample prev for each combination
snps.case<-boot.cases[,-1]
snps.ctrl<-boot.ctrls[,-1]

counts[, 1] <- apply(x.case, 1, compare,mat=boot.cases,k=k)  # gives number of cases for each genotype combination 
counts[, 2] <- apply(x.ctrl, 1, compare,mat=boot.ctrls,k=k)  # gives number of controls for each genotype combination

		tp<- sum(counts[which(hr==1), 1]) #true positives
		fp<- sum(counts[which(hr==1), 2]) #false positives
		tn<- sum(counts[which(hr==0), 2]) #true negatives
		fn<- sum(counts[which(hr==0), 1]) #false negatives

error <- 100 * (fp + fn)/n
boot.error[i]<-error  } #end loop across bootstrap samples

est<-mean(boot.error) #estimate is the mean across b samples
return(list("classification error estimate"=est, "classification accuracy estimate"=100-est))
} #end function

permute.mdr<-function(accuracy, loci, N.permute, method = c("CV", "3WS","none"), 
    data, cv, K, x = NULL, proportion = NULL, ratio = NULL, equal = "HR", 
    genotype = c(0, 1, 2), LRT = FALSE) 
{
    if (is.null(ratio)) 
        ratio <- length(which(data[, 1] == 1))/length(which(data[, 
            1] == 0))
    mdr.lrt <- function(data, loci, x, ratio, equal, genotype) {
        X <- matrix(0, dim(data), length(loci))
        for (i in 1:length(loci)) {
            hr <- mdr(data, loci[i], x, ratio, equal, genotype)
            hr <- hr$"high risk/low risk"
            f <- as.factor(data[, (loci[i] + 1)])
            s <- nlevels(f)
            d <- diag(s)[f, ]
            X[, i] <- d %*% hr
        }
        form1 <- paste(paste("X[,", 1:length(loci), "]", sep = ""), 
            collapse = "*")
        com <- combn(1:length(loci), length(loci) - 1)
        term <- NULL
        for (j in 1:dim(com)[2]) {
            term <- c(term, paste(paste("X[,", com[, j], "]", 
                sep = ""), collapse = "*"))
        }
        form2 <- paste(term, collapse = "+")
        lrfit.full <- glm(as.formula(paste("data[,1]~", form1, 
            sep = "")), family = binomial(link = logit))
        lrfit.reduced <- glm(as.formula(paste("data[,1]~", form2, 
            sep = "")), family = binomial(link = logit))
        x2 <- lrfit.reduced$deviance - lrfit.full$deviance
        return(x2)
    }
    if (LRT == TRUE) {
        x2 <- mdr.lrt(data = data, loci = loci, x = 1, ratio = ratio, 
            equal = equal, genotype = genotype)
        chi <- rep(0, N.permute)
    }
    acc <- rep(0, N.permute)
    if (method == "CV") {
        for (i in 1:N.permute) {
            permute <- sample(data[, 1], length(data[, 1]), replace = FALSE)
            new <- cbind(permute, data[, -1])
            fit <- mdr.cv(data = new, K, cv, ratio, equal, genotype)
            acc[i] <- fit$"final model accuracy"
            if (LRT == TRUE) {
                chi[i] <- mdr.lrt(data = new, loci = fit$"final model", 
                  x = 1, ratio = ratio, equal = equal, genotype = genotype)
            }
        }
    }
    if (method == "3WS") {
        for (i in 1:N.permute) {
            permute <- sample(data[, 1], length(data[, 1]), replace = FALSE)
            new <- cbind(permute, data[, -1])
            fit <- mdr.3WS(data = new, K, x, proportion, ratio, 
                equal, genotype)
            acc[i] <- fit$"final model accuracy"
            if (LRT == TRUE) {
                chi[i] <- mdr.lrt(data = new, loci = fit$"final model", 
                  x = 1, ratio = ratio, equal = equal, genotype = genotype)
            }
        }
    }
    if (method == "none") {
        for (i in 1:N.permute) {
            permute <- sample(data[, 1], length(data[, 1]), replace = FALSE)
            new <- cbind(permute, data[, -1])
            fit <- mdr(split = new, comb=loci, x=1, ratio=ratio,equal=equal,genotype=genotype)
            acc[i] <- fit$"balanced accuracy"
            if (LRT == TRUE) {
                chi[i] <- mdr.lrt(data = new, loci = fit$"models", 
                  x = 1, ratio = ratio, equal = equal, genotype = genotype)
            }
        }
    }	
    acc <- sort(acc)
    pval <- sum(acc > accuracy)/N.permute
    result <- list(`Permutation P-value` = pval, `Permutation Distribution` = acc)
    if (LRT == TRUE) {
        chi <- sort(chi)
        lrt.p <- sum(chi > x2)/N.permute
        result <- list(`Permutation P-value` = pval, `Permutation Distribution` = acc, 
            `LRT P-value` = lrt.p, `LRT Distribution` = chi)
    }
    return(result)
}

################### define methods for the 'mdr' class

summary.mdr<-function(object, ...) { #previously fit mdr object

k<-length(object$'top models')  #define sizes of interaction considered
level<-1:k  #define the sequences of levels
best.models<-matrix(NA,k,k) #matrix to store the best models for each level

for (i in 1:k) {
mod<-object $'top models'[[i]]
if (length(mod)==k) {
best.models[i,]<-mod} else {
add<-k-length(mod)
best.models[i,]<-c(mod, rep(NA,add))  #store the best model for level i
}
}

model<-object$'final model'
best<-NULL
for (i in 1:k) {
if (i==length(model)) {best<-c(best,"*") } else {best<-c(best,"") } 
}

if (object$'validation method'=="3WS") {

	train<-round(object $'top model accuracies'[,1], digits=2 )#save training, testing, and validation accuracies
	test<-round(object $'top model accuracies'[,2], digits=2 )
	val<-round(object $'top model accuracies'[,3], digits=2 )

	tab<-cbind(level,best.models,train,test,val)
	rownames(tab)<-best
	colnames(tab)<-c("Level", "   Best Models", paste(rep("", k-1), sep=""), "   Training Accuracy", "   Testing Accuracy", "   Validation Accuracy")
	print(tab, na.print="")
	cat("", "\n")
	cat("'*' indicates overall best model")
	}    else {  if   (object$'validation method'=="CV") { 

		ca<-round(object $'top model accuracies'[,1], digits=2 )#save classification and prediction accuracy and CVC
		pa<-round(object $'top model accuracies'[,2], digits=2 )
		cvc<-round(object $'top model accuracies'[,3], digits=2 )

		tab<-cbind(level,best.models,ca,pa,cvc)
		rownames(tab)<-best
		colnames(tab)<-c("Level", "   Best Models", paste(rep("", k-1), sep=""), "   Classification Accuracy", "   Prediction Accuracy", "   Cross-Validation Consistency")
		print(tab, na.print="")
		cat("", "\n")
		cat("'*' indicates overall best model")
		}  else {  print("invalid mdr object")
		} }#close nested ifs
 } #close function 

predict.mdr<-function(object,new.data,...){ #previously fit mdr object, and a new dataset with the same original variables in the same order (minus the response)
model<-object$'final model'
hr<-object$high
k<-length(model)
geno<-list(object$genotypes) #need to define the genotype combinations for the model
comb<-as.matrix(expand.grid(rep(geno, k)))
new.data<-as.matrix(new.data[, model]) #reduce data to only the model of interest
match<-matrix(0,dim(new.data)[1],dim(comb)[1])#define matrix to store genotype combination
for (i in 1:dim(new.data)[1]) {
for (j in 1:dim(comb)[1]) {
match[i,j]<- sum(new.data[i,]==comb[j,])==k  #create indicator for genotype combination j for individual i
}#close loop over combinations
}#close for loop over individuals
predict<-match%*%hr #indicator for genotype combination of individual i is high risk
return(predict)
}#close function

plot.mdr<-function(x,data,main="",xlab="",ylab="Count", table=FALSE,...){ #previously fit MDR object and associated dataset; table arg gives a summary table of the counts

mod<-x$'final model'
geno<-list(as.factor(x$genotypes)) #read as categorical for appropriate labeling
g<-length(x$genotypes)
k<-length(mod)
m<-k%/%2
hr<-x$high

r<- g^m #number of rows for plotting window
c<- g^(k-m) #number of columns for plotting window

case<-data[which(data[,1]==1),mod+1] #split data into cases and controls
ctrl<-data[which(data[,1]==0),mod+1]

tab1<-table(case)  #create contingency table of multifactorial genotypes
tab0<-table(ctrl)

comb<-expand.grid(rep(geno,k))
dat<-c(as.vector(tab1),as.vector(tab0))
dat<-data.frame(dat)
dat<-cbind(dat,c(rep("Case",g^k),rep("Control",g^k)))
comb2<-rbind(comb,comb)
hr2<-rbind(hr,hr)
dat<-cbind(dat,comb2,hr)

snps<-NULL
for (i in 1:k) {
snps<-c(snps,paste("SNP",as.character(mod[i]),sep="")) }

colnames(dat)<-c("count","status",snps,"high-risk")
if (table) {print(dat)}

genotypes<-snps[1]
if (k>1) {for (i in 2:k) {
genotypes<-paste(genotypes,snps[i],sep="*") } }

form<-formula(paste("count~status |",genotypes))

panel.color<-function(x,y,...) { 
color<-c("black","white")
if(y[1]>y[2])panel.fill(col="gray")
panel.barchart(x,y,col=color,horizontal=FALSE)
} 

barchart(form, data=dat,layout=c(c,r),
        ylab=ylab, xlab=xlab,
        strip=strip.custom(strip.names=c(TRUE,TRUE),strip.levels=c(TRUE,TRUE)),
        panel=panel.color, key=list(text=list(c("High-Risk","Low-Risk"),columns=2,title="Legend"),rectangles=list(col=c("gray","white"),columns=2)), main=main)

} #close function