
semiMarkov<-function(data,cov=NULL,states,mtrans,cov_tra=NULL,cens=NULL,
		dist_init=NULL,proba_init=NULL,coef_init=NULL,alpha_ci=0.95,Wald_par=NULL,
		eqfun=NULL,ineqLB=NULL,ineqUB=NULL,control = list()){

#length of cov_tra
if(missing(cov_tra)){length_cov<-0}
else{length_cov<-length(cov_tra)}
#what are the states in data set
states_data<-sort(unique(c(unique(data[,2]),unique(data[,3]))))
#checking conditions
if (missing(data))
	stop("Argument 'data' is missing with no default")
if (missing(states))
	stop("Argument 'states' is missing with no default")
if (missing(mtrans))
	stop("Argument 'mtrans' is missing with no default")
if (missing(cov) && missing(cov_tra)==FALSE)
	stop("To indicate the transitions for covariates you need to define the covariates matrix")
if(!missing(cov)){
if(!is.data.frame(cov))
	stop("Argument 'cov' must be a data.frame")
if(dim(data)[1]!=dim(cov)[1])
	stop("Argument 'data' and 'cov' need to have the same number of rows")
}
if(!is.data.frame(data))
	stop("Argument 'data' must be a data.frame")
 if (nrow(mtrans) != ncol(mtrans)) 
        stop("Argument 'mtrans' must be a quadratic  matrix.")
if (!is.list(control)){warning("Argument 'control' must be a list. The default parameters for optimization are considered.")
						control <- list()}



##transitions on diagonal
k<-0
for(i in 1:dim(mtrans)[1]){
if(mtrans[i,i]!=FALSE)k<-k+1}
if (k> 0) 
       stop("Transitions into the same state are not allowed")
if (nrow(mtrans) != length(states)) 
        stop("The row number of 'mtrans' must be equal to the number of states.")
 if (length(states) != length(unique(states))) 
        stop("The state names must be unique.")

if(!all(noquote(states)%in%states_data))stop("States in vector 'states' should be the same as states in data")
#we change the names of states for numbers
for(i in 1:dim(data)[1]){
for(j in 1:length(states)){
if(data[i,2]==states[j])data[i,2]<-j
if(data[i,3]==states[j])data[i,3]<-j
}
}
for(i in 1:length(states)){states[i]<-i}

#transitions specified in mtrans must be the same as in data
s_0<-0
for(i in 1:dim(mtrans)[1]){
if( all(mtrans[i,]==F))s_0<-s_0+1
for(j in 1:dim(mtrans)[2]){
if(i!=j && mtrans[i,j]==FALSE && paste(i,j)%in%unique(paste(data[,2],data[,3])))stop("All observed transitions in data must be specified in matrix 'mtrans'")
}}

#Defining number of states

s<-length(states)

if(s_0>0){ss<-s-s_0}
else{ss<-s}

### censorship and table.state
if(!missing(cens) && !cens%in%data[,3])stop("Wrong format of the argument 'cens'")

if(!missing(cens)){
for(i in 1:dim(data)[1]){
if(data[i,3]==cens){
if(i<dim(data)[1] && data[i+1,1]==data[i,1]){

stop("The censure is not the last state of individuals")}

}
}
}
#trans.h - all states left by the process in data
#trans.j - all states entered by the process in data
trans.h<-data[data[,2]!=data[,3],2]
trans.j<-data[data[,2]!=data[,3],3]

#if censorship is not defined by user
if(missing(cens)){
trans.cens<-data[data[,2]==data[,3],3]}
else{
# else we take censorship as indicated by user
trans.cens<-data[data[,3]==cens,3]}

#defining table of states and counting the number of each transition
table.state<-matrix(ncol=s, nrow=s) 
for(i in 1:s){
for (j in 1:s){
if(i!=j){
table.state[i,j]<-1*(sum(1*(trans.h==states[i] & trans.j==states[j])))}
else{table.state[i,j]<-1*(sum(1*(trans.cens==states[i])))}
}
}

#if there is a censorship defined by user we change it into transition h->h
if(!missing(cens)){
for(i in 1:dim(data)[1]){
if(data[i,3]==cens){

data[i,3]<-data[i,2]

}
}}

#number of censured patients
Ncens<-length(which(data[,2]==data[,3]))

#if the particular transitions for covariates are indicated we need to check if they are defined for all the covariates in data.frame 'cov'
t<-0
if(length_cov>=1){
if(dim(cov)[2]!=length_cov)
	stop("The number of transitions indicated must be the same as the number of covariates of the model")
for(i in 1:length_cov){
for(j in 1:length(cov_tra[[i]])){
#we define transitions for every covariate and we check if they are possible in 'mtrans'
m<-as.numeric(substring(cov_tra[[i]][j],first=1,last=1))
k<-as.numeric(substring(cov_tra[[i]][j],first=2,last=2))
if(mtrans[m,k]==FALSE){
t<-t+1}}}

if(t>0)stop("The indicated transitions in 'cov_tra' must specified as possible in the matrix 'mtrans'")
}

#max time
last<-max(data[,4])

time.ind<-NA

#number of the covariates r
if(missing(cov)){r<-0

cov=as.matrix(rep(0,dim(data)[1]))}
else{
if(length_cov>=1){r<-c()
for(i in 1:length_cov){
r<-c(r,length(cov_tra[[i]]))
cov_tra[[i]]<-sort(cov_tra[[i]])
}
cov<-as.matrix(cov)

}
else{cov<-as.matrix(cov)
r<-dim(cov)[2]

}
#we check if covariates are dependent from time
time.ind<-rep(TRUE,dim(cov)[2])
for(j in 1:dim(cov)[2]){
l<-0
for(i in 2:dim(cov)[1]){
if(data[i,1]==data[i-1,1]){
if(cov[i,j]!=cov[i-1,j]){
time.ind[j]<-FALSE
i<-dim(cov)[1]+1
}}}}}



#Construction of the variable that identifies transitions to a different state
#h!=j -> 1
#h=j -> 0
data[,5]<-1*(data[,2]!=data[,3])

#Construction of the variable that specifies the transitions 
data[,6]<-as.character(paste(data[,2], data[,3], sep=""))

#All observed transitions h!=j
trans.hj<-c()
for(i in 1:s){
for(j in 1:s){
if(mtrans[i,j]!=FALSE)trans.hj[length(trans.hj)+1]<-as.character(paste(states[i],states[j],sep=""))
}}

#All observed transitions h=j
trans.hh<-c()
trans.hh<-sort(unique(data[data[,5]==0,6]))
trans.hh2<-trans.hh
if(length(trans.hh)>0){
for(i in 1:length(trans.hh)){
a<-0
for(j in 1:length(trans.hh)){
 
if(mtrans[as.numeric(substring(trans.hh[i],first=1,last=1)),j]==FALSE){
a<-a+1
}

if(a==length(trans.hh)){
trans.hh2<-trans.hh[-which(as.numeric(substring(trans.hh,first=1,last=1))==i)]}}}}

trans.hh<-trans.hh2

#Number of transitions to estimate
p<-length(which(as.vector(mtrans)!=FALSE))


###Initialisation of the parameters

#Initial Transition Matrix P
#t vector of number of essential transitions in a row
#npos vector for all possible transitions
t<-rep(0,s)
pos<-list()
npos<-c()
matrix.P<-matrix(ncol=s, nrow=s) 

#defining transtions' distributions
dist<-c()
for(i in 1:s){
k<-0
for (j in 1:s){
if(mtrans[i,j]!=FALSE){
dist<-c(dist,mtrans[i,j])
if(sum(1*(trans.h==states[i]))>0){
matrix.P[i, j]<-1*(i!=j)*(sum(1*(trans.h==states[i] & trans.j==states[j]))/sum(1*(trans.h==states[i])))}
else{matrix.P[i,j]<-0}

k<-k+1
pos[[length(pos)+1]]<-c(states[i],states[j])
npos[length(npos)+1]<-as.character(paste(states[i],states[j],sep=""))
}
else{matrix.P[i,j]<-0}

}
t[i]<-k
}
t<-t-1


## codes for distributions

for(i in 1:length(dist)){
if(!dist[i]%in%c("E","Exp","Exponential","W","Weibull","EW","EWeibull","Exponentiated Weibull","TRUE"))stop("Wrong format of matrix 'mtrans'")
if(dist[i]%in%c("E","Exp","Exponential"))dist[i]<-1
if(dist[i]%in%c("W","Weibull","TRUE"))dist[i]<-2
if(dist[i]%in%c("EW","EWeibull","Exponentiated Weibull"))dist[i]<-3}
dist<-as.numeric(dist)

pos.temp<-list()

#deleting the last transitions in the row from the calculations of P
if(s>2){
for(i in 1:(length(pos)-1)){

if(pos[[i]][[1]]==pos[[i+1]][[1]])pos.temp[length(pos.temp)+1]<-pos[i]
if(i>1){
if(pos[[i]][[1]]==pos[[i+1]][[1]])pos.temp[length(pos.temp)+1]<-pos[i]# && pos[[i]][[1]]!=pos[[i-1]][[1]]
}}
pos<-unique(pos.temp)

npos<-c(length(pos))
for(i in 1:length(pos)){

npos[i]<-as.character(paste(pos[[i]][[1]],pos[[i]][[2]],sep=""))
}
}

##auxilary matrice of logical values
mtrans.log<-matrix(FALSE,ncol=s,nrow=s)
for(i in 1:dim(mtrans)[1]){
for(j in 1:dim(mtrans)[2]){
if(mtrans[i,j]!=FALSE)mtrans.log[i,j]<-TRUE
}}

#number of probabilites in the parameters matrix
proba_new<-c()
proba<-proba_init

if(length(states)>2){
	nprob<-length(npos)

	for(i in 1:s){
		if(1*sum(mtrans.log[i,])>0){
			if(!missing(proba_init)){
				p_in_row<-1*sum(mtrans.log[i,])

				proba_new<-c(proba_new,proba[1:p_in_row])
				proba_new<-proba_new[1:(length(proba_new)-1)]

				proba<-proba[(p_in_row+1):length(proba)]}

			}
		}
	}
else{nprob<-0
npos<-c()}



if(!missing(proba_init)){
proba_init<-proba_new}

#Matrix of the Parametres
# e parametres sigma
# trans.e for which transitions
# w parametres sigma
# trans.w for which transitions
# ew parameters theta
# trans.ew for which transitions
###########
e<-0
trans.e<-trans.hj
w<-0
trans.w<-c()
ew<-0
trans.ew<-c()
for(i in 1:length(dist)){
if(dist[i]==1)e<-e+1
if(dist[i]==2){e<-e+1
w<-w+1
trans.w<-c(trans.w,trans.hj[i]) }
if(dist[i]==3){e<-e+1
w<-w+1
trans.w<-c(trans.w,trans.hj[i])
ew<-ew+1
trans.ew<-c(trans.ew,trans.hj[i])}}
ndist<-e+w+ew

# sum(t)  initial probabilities of Markov chain
# r coeff of regression Beta
nbeta<-0
if(length(r)>1 || r>0){
if(length_cov>=1){
label_beta<-c()
nbeta<-0
trans_beta<-c()
val_beta<-c()
for(i in 1:length(r)){
label_beta<-c(label_beta,sort(rep(paste("Beta", i, sep=""), r[i])))
nbeta<-nbeta+r[i]
trans_beta<-c(trans_beta,cov_tra[[i]])
val_beta<-rep(0,length(trans_beta))
}

}
else{label_beta<-sort(rep(paste("Beta", 1:r, sep=""), p))
nbeta<-p*r
trans_beta<-rep(trans.hj,r)
val_beta<-rep(0,p*r)
}}
else{
label_beta<-c()
trans_beta<-c()
val_beta<-c()}

#The initial parameters
#checking if the vectors are OK


if(!missing(dist_init)){
	if(length(dist_init)!=ndist || dist_init[1:(ndist)]<=0 || !is.numeric(dist_init))	
		stop("Wrong format of the vector 'dist_init'")
}
if(!missing(proba_init)){
	if(length(proba_init)!=nprob || proba_init<=0 || proba_init>=1 || !is.numeric(proba_init))
		stop("Wrong format of the vector 'proba_init'")
}
if(!missing(coef_init)){
	if(length(coef_init)!=nbeta || !is.numeric(coef_init))
		stop("Wrong format of the vector 'coef_init'")
}

#Parameters for the Wald test
if(!missing(Wald_par)){
	if(length(Wald_par)!=ndist+nbeta || !is.numeric(Wald_par))	
		stop("Wrong format of the vector 'Wald_par'")
}

#conditions for the constraints for solnp
if (!missing(eqfun) && !is.list(eqfun)) {warning("Argument 'eqfun' must be a list of vectors of length 4.")
											eqfun <- NULL}
	else if(!missing(eqfun)){
	for(i in 1:length(eqfun)){
	if(eqfun[[i]][1]=="coef" && nbeta==0){warning("There are no covariates in the model. The default constraints for optimization are considered.")
											eqfun <- NULL
											break}	
											
	if(length(states)>2){		

	if(!is.vector(eqfun[i]) || length(eqfun[[i]])!=4 || as.numeric(eqfun[[i]][2])>max(ndist,nbeta,nprob) || as.numeric(eqfun[[i]][3])>max(ndist,nbeta,nprob) || !(eqfun[[i]][1]%in%c("coef","proba","dist"))){
	warning("1Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break}
											else{
	
	if(eqfun[[i]][1]=="dist" && (!(as.numeric(eqfun[[i]][2])%in%c(1:ndist)) | !(as.numeric(eqfun[[i]][3])%in%(1:ndist))) ){warning("2Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break}
	if(eqfun[[i]][1]=="proba" && (!(as.numeric(eqfun[[i]][2])%in%c(1:nprob)) | !(as.numeric(eqfun[[i]][3])%in%(1:nprob))) ){warning("4Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break}										
											
	if(eqfun[[i]][1]=="coef" && (!(as.numeric(eqfun[[i]][2])%in%c(1:nbeta)) | !(as.numeric(eqfun[[i]][3])%in%(1:nbeta)))){warning("3Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break}		
	
	}}else{
	if(eqfun[[i]][1]=="proba"){warning("5Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break	
	}else{
	if(!is.vector(eqfun[i]) || length(eqfun[[i]])!=4 || as.numeric(eqfun[[i]][2])>max(ndist,nbeta) || as.numeric(eqfun[[i]][3])>max(ndist,nbeta) || !(eqfun[[i]][1]%in%c("coef","dist"))){
	warning("6Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break}
											else{
	
	if(eqfun[[i]][1]=="dist" && (!(as.numeric(eqfun[[i]][2])%in%c(1:ndist)) | !(as.numeric(eqfun[[i]][3])%in%(1:ndist))) ){warning("7Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break}
		if(eqfun[[i]][1]=="coef" && (!(as.numeric(eqfun[[i]][2])%in%c(1:nbeta)) | !(as.numeric(eqfun[[i]][3])%in%(1:nbeta)))){warning("8Wrong format of 'eqfun'.The default constraints for optimization are considered.")
											eqfun <- NULL
											break}		
	
	}}
	}
		
	}
	}


if (!missing(ineqLB) && !is.list(ineqLB)) {warning("Argument 'ineqLB' must be a list of vectors of length 2. The default constraints for optimization are considered.")
											ineqLB <- NULL}
	else if(!missing(ineqLB)){
	for(i in 1:length(ineqLB)){
	if(length(states)>2){
	if(!is.vector(ineqLB[i]) || length(ineqLB[[i]])!=3 || as.numeric(ineqLB[[i]][2])>max(ndist,nbeta,nprob) || !(ineqLB[[i]][1]%in%c("coef","dist","proba" ))){
	warning("Wrong format of 'ineqLB'. The default constraints for optimization are considered.")
											i <- length(ineqLB)
											ineqLB <- NULL}
											}
	else{
		if(ineqLB[[i]][1]=="proba" ){
	warning("Wrong format of 'ineqLB'. The default constraints for optimization are considered.")
											i <- length(ineqLB)
											ineqLB <- NULL}
		if(!is.vector(ineqLB[i]) || length(ineqLB[[i]])!=3 || as.numeric(ineqLB[[i]][2])>max(ndist,nbeta) || !(ineqLB[[i]][1]%in%c("coef","dist"))){
	warning("Wrong format of 'ineqLB'. The default constraints for optimization are considered.")
											i <- length(ineqLB)
											ineqLB <- NULL}
	}
	}
	}

	if (!missing(ineqUB) && !is.list(ineqUB)) {warning("Argument 'ineqUB' must be a list of vectors of length 3. The default constraints for optimization are considered.")
											ineqUB <- NULL}
	else if(!missing(ineqUB)){
	for(i in 1:length(ineqUB)){
	if(length(states)>2){
	if(!is.vector(ineqUB[i]) || length(ineqUB[[i]])!=3 || as.numeric(ineqUB[[i]][2])>max(ndist,nbeta,nprob) || !(ineqUB[[i]][1]%in%c("coef","dist","proba" ))){
	warning("Wrong format of 'ineqLB'. The default constraints for optimization are considered.")
											i <- length(ineqUB)
											ineqUB <- NULL}
											}
	else{
		if(ineqUB[[i]][1]=="proba"){
	warning("Wrong format of 'ineqLB'. The default constraints for optimization are considered.")
											i <- length(ineqUB)
											ineqUB <- NULL}
		if(!is.vector(ineqUB[i]) || length(ineqUB[[i]])!=3 || as.numeric(ineqUB[[i]][2])>max(ndist,nbeta) || !(ineqUB[[i]][1]%in%c("coef","dist"))){
	warning("Wrong format of 'ineqLB'. The default constraints for optimization are considered.")
											i <- length(ineqUB)
											ineqUB <- NULL}
	}
	}
	}
	
	
temp2<-data.frame( transit=data[,6],covariable=cov)

#data frame of parameters with columns Label, Transition and Value
if(length(r)==1){
if(r==0){l<-0
} else if(length(r)>1 || r>0){l<-p}
parameters<-data.frame(label=c(rep("sigma",e),rep("nu", w), rep("theta",ew),
 rep("P",nprob),label_beta ),
transit=c(trans.e,trans.w,trans.ew,npos,trans_beta),
 value=c(rep(1, e),rep(1,w),rep(1,ew),rep(0,nprob),val_beta))
for(i in 1:s){
for(j in 1:s){
parameters$value[parameters$label=="P" & parameters$transit==as.character(paste(states[i], states[j], sep=""))]<-matrix.P[i,j]
}}
parameters<-parameters[parameters$value<=1,]

if(!missing(dist_init))parameters[1:ndist,3]<-dist_init
if(!missing(proba_init))parameters[(ndist+1):(ndist+nprob),3]<-proba_init
if(!missing(coef_init))parameters[(ndist+nprob+1):(ndist+nprob+nbeta),3]<-coef_init

}else{
l<-p
parameters<-data.frame(label=c(rep("sigma",e),rep("nu", w), rep("theta",ew),
 rep("P",nprob),label_beta ),
transit=c(trans.e,trans.w,trans.ew,npos,trans_beta),
 value=c(rep(1, e),rep(1,w),rep(1,ew),rep(0,nprob),val_beta))
for(i in 1:s){
for(j in 1:s){
parameters$value[parameters$label=="P" & parameters$transit==as.character(paste(states[i], states[j], sep=""))]<-matrix.P[i,j]
}}
parameters<-parameters[parameters$value<=1,]

if(!missing(dist_init))parameters[1:(ndist),3]<-dist_init
if(!missing(coef_init))parameters[(ndist+nprob+1):(ndist+nprob+nbeta),3]<-coef_init
}
data1<-data[data[,2]!=data[,3],]
cov1<-as.matrix(cov[data[,2]!=data[,3],])

data2<-data[data[,2]==data[,3],]  
cov2<-as.matrix(cov[data[,2]==data[,3],])

indic1<-(rep(data1[,6], length(trans.hj))==sort(rep(trans.hj, length(data1[,6]))))
indic2<-(rep(data2[,6], length(trans.hh))==sort(rep(trans.hh, length(data2[,6]))))

cov_pos<-c()
if(missing(cov_tra)){}else{
position<-c()
for(i in 1:length_cov){
for(j in 1:length(cov_tra[[i]])){
position_temp<-which(parameters[,2]==cov_tra[[i]][j] )

position<-c(position,position_temp[position_temp<p+1])


}
cov_pos<-c(cov_pos,position)}}
if(length(cov_pos)>0){
cov_pos<-sort(position)}


################
#likelihood
################

if(missing(cov_tra)){
V<-function(x){

logdensity<-log(.dens(nprob,mtrans,dist,data1[,4],r,cov1, x,parameters)[indic1])
logmarginal<-log(.marg(nprob,mtrans,dist,data2[,4],r,cov2,x,parameters)[indic2])


return(-1*(sum(logdensity) + sum(logmarginal)))
}}
else{
V<-function(x){

  logdensity<-log(.dens_trans(nprob,mtrans,dist,data1[,4],dim(cov)[2],cov1, x,cov_pos,parameters)[indic1])
logmarginal<-log(.marg_trans(nprob,mtrans,dist,data2[,4],dim(cov)[2],cov2,x,cov_pos,parameters)[indic2])

  return(-1*(sum(logdensity) + sum(logmarginal)))
}
}

#number of states with probabilities to estimate
###########


if(length(r)>1 || r>0)ncov<-length(parameters[(ndist+nprob+1):dim(parameters)[1],1])
else ncov<-0
##########
#optimisation
###########

if(length(states)>2){
#constraints for optimisation
	hin<-function(par){
	
	vec_contr<-c()
  #temporary parameter for the lower bound
	a<-par[ndist+1]
 
  #temporary paramter for the upper bound
		for(i in (ndist+1):(ndist+nprob)){
    #we add the proba as long as it belongs to the same row of transitions
		if(i<ndist+nprob && substring(parameters[i,2],first=1,last=1)==substring(parameters[i+1,2],first=1,last=1))
			{
      a<-par[i+1]+a

			}else{
        #if the next parameter belong to other row of transitions we finish summing and we put constraints 
        #on temporary parameter a and b
        
			#a_last<-1-a
      #a<-a+a_last
					vec_contr<-c(vec_contr,a)
			a<-par[i+1]
							}

		}
	h<-c(vec_contr)

	
	
		
	return(h)
	}

	
	ineq_low <- c(rep(0.001,ndist),rep(0.001,nprob),rep(-1000,ncov))
	if(!is.null(ineqLB)){
	for(i in 1:length(ineqLB)){
	if(ineqLB[[i]][1]=="dist" )ineq_low[as.numeric(ineqLB[[i]][2])] <- as.numeric(ineqLB[[i]][3])
	if(ineqLB[[i]][1]=="proba" )ineq_low[as.numeric(ineqLB[[i]][2])+ndist] <- as.numeric(ineqLB[[i]][3])
	if(ineqLB[[i]][1]=="coef" )ineq_low[as.numeric(ineqLB[[i]][2])+ndist+nprob] <- as.numeric(ineqLB[[i]][3])
	if(ineqLB[[i]][1]=="dist" && ineq_low[as.numeric(ineqLB[[i]][2])]<0.001)ineq_low[as.numeric(ineqLB[[i]][2])]<-0.001
	if(ineqLB[[i]][1]=="proba" && ineq_low[as.numeric(ineqLB[[i]][2])+ndist]<0.001)ineq_low[as.numeric(ineqLB[[i]][2])+ndist]<-0.001
	if( ineqLB[[i]][1]=="coef" && ineq_low[as.numeric(ineqLB[[i]][2])+ndist+nprob]<(-1000))ineq_low[as.numeric(ineqLB[[i]][2])+ndist+nprob]<--1000
	}
	}
	
	#UB
	ineq_up <- c(rep(1000,ndist),rep(0.999,nprob),rep(1000,ncov))
	
	if(!is.null(ineqUB)){
	for(i in 1:length(ineqUB)){
	if(ineqUB[[i]][1]=="dist" )ineq_up[as.numeric(ineqUB[[i]][2])] <- as.numeric(ineqUB[[i]][3])
	if(ineqUB[[i]][1]=="proba" )ineq_up[as.numeric(ineqUB[[i]][2])+ndist] <- as.numeric(ineqUB[[i]][3])
	if(ineqUB[[i]][1]=="coef" )ineq_up[as.numeric(ineqUB[[i]][2])+ndist+nprob] <- as.numeric(ineqUB[[i]][3])
	if(ineqUB[[i]][1]=="dist" && ineq_up[as.numeric(ineqUB[[i]][2])]>1000)ineq_up[as.numeric(ineqUB[[i]][2])]<-1000
	if( ineqUB[[i]][1]=="proba" && ineq_up[as.numeric(ineqUB[[i]][2])+ndist]>0.999)ineq_up[as.numeric(ineqUB[[i]][2])+ndist]<-0.999
	if( ineqUB[[i]][1]=="coef" && ineq_up[as.numeric(ineqUB[[i]][2])+ndist+nprob]>1000)ineq_up[as.numeric(ineqUB[[i]][2])+ndist+nprob]<-1000
	}
}
	
	if(!is.null(eqfun)){
	eqfun1 <- function(par){
	
	h <- c()
	for(i in 1:length(eqfun)){
	if(eqfun[[i]][1]=="dist"){
	h <- c(h,par[as.numeric(eqfun[[i]][2])]-par[as.numeric(eqfun[[i]][3])]*as.numeric(eqfun[[i]][4]))}
	else if(eqfun[[i]][1]=="proba"){
	h <- c(h,par[as.numeric(eqfun[[i]][2])+ndist]-par[as.numeric(eqfun[[i]][3])+ndist]*as.numeric(eqfun[[i]][4]))
	}else
	{
	h <- c(h,par[as.numeric(eqfun[[i]][2])+ndist+nprob]-par[as.numeric(eqfun[[i]][3])+ndist+nprob]*as.numeric(eqfun[[i]][4]))
	}
		}
	return(h)
	}
	
	}else{eqfun1<-NULL}
	
	

if(length(states)>3){

res <- solnp( pars = parameters[,3], fun=V,eqfun=eqfun1,eqB=rep(0,length(eqfun)),ineqfun=hin,ineqUB=c(rep(1,length(hin(parameters[,3])))),ineqLB=rep(0,length(hin(parameters[,3]))),
                     LB=ineq_low, UB=ineq_up,control = control)

}
else{


 res <- solnp( pars = parameters[,3], fun=V,eqfun=eqfun1,eqB=rep(0,length(eqfun)),
			LB=ineq_low, UB=ineq_up,control = control)

  }

  solution<-res$par
objective<-res$value[-1]
iterations<-res$outer.iterations

}else{
# 2 states 
ineq_low <- c(rep(0.001,ndist),rep(-1000,ncov))
	if(!is.null(ineqLB)){
	for(i in 1:length(ineqLB)){
	if(ineqLB[[i]][1]=="dist")ineq_low[as.numeric(ineqLB[[i]][2])] <- as.numeric(ineqLB[[i]][3])
	if( ineqLB[[i]][1]=="coef")ineq_low[as.numeric(ineqLB[[i]][2])+ndist] <- as.numeric(ineqLB[[i]][3])
	if(ineqLB[[i]][1]=="dist" && ineq_low[as.numeric(ineqLB[[i]][2])]<0.001)ineq_low[as.numeric(ineqLB[[i]][2])]<-0.001
	if( ineqLB[[i]][1]=="coef" && ineq_low[as.numeric(ineqLB[[i]][2])+ndist]<(-1000))ineq_low[as.numeric(ineqLB[[i]][2])+ndist]<--1000
	}
	}
	
	#UB
	ineq_up <- c(rep(1000,ndist),rep(1000,ncov))
	
	if(!is.null(ineqUB)){
	for(i in 1:length(ineqUB)){
	ineq_up[as.numeric(ineqUB[[i]][2])] <- as.numeric(ineqUB[[i]][3])
	if(ineqUB[[i]][1]=="dist" && ineq_up[as.numeric(ineqUB[[i]][2])]>1000)ineq_up[as.numeric(ineqUB[[i]][2])]<-1000
	if(ineqUB[[i]][1]=="coef" && ineq_up[as.numeric(ineqUB[[i]][2])+ndist]>1000)ineq_up[as.numeric(ineqUB[[i]][2])+ndist]<-1000
	}
}
if(!is.null(eqfun)){
	eqfun1 <- function(par){
	
	h <- c()
	for(i in 1:length(eqfun)){

	if(eqfun[[i]][1]=="dist"){
	h <- c(h,par[as.numeric(eqfun[[i]][2])]-par[as.numeric(eqfun[[i]][3])]*as.numeric(eqfun[[i]][4]))}
	else{
	h <- c(h,par[as.numeric(eqfun[[i]][2])+ndist+nprob]-par[as.numeric(eqfun[[i]][3])+ndist+nprob]*as.numeric(eqfun[[i]][4]))
	}
		}
	return(h)
	}
	
	}else{eqfun1<-NULL}


 res <- solnp( pars = parameters[,3], fun=V,eqfun=eqfun1,eqB=rep(0,length(eqfun)),
			LB=ineq_low, UB=ineq_up,control = control)
solution<-res$par
objective<-res$value[-1]
iterations<-res$outer.iterations
}
if(res$convergence>0)stop("Problem in optimization. Verify the initial parameters or the constraints")
############
#definition of the class "semiMarkov"
############

hessian<-round(diag(ginv(hessian(V,solution))),5)
#Transitions
transitions<-paste(substring(parameters[,2],first=1,last=1),"->",substring(parameters[,2],first=2,last=2))

#Initial vectors for the statistics
sd<-c(1:dim(parameters)[1])
lower<-c(1:dim(parameters)[1])
upper<-c(1:dim(parameters)[1])
StatTest<-c(1:dim(parameters)[1])
pvalue<-c(1:dim(parameters)[1])
Wald_H0<-c(1:dim(parameters)[1])

#Standard Deviation for parameters of waiting time distribution

for(i in 1:dim(parameters)[1]){
if(hessian[i]>=0){
sd[i]<-formatC(round(sqrt(hessian[i]),digits=2),2,format="f")

lower[i]<-formatC(round(solution[i]-qnorm(.975)*sqrt(hessian[i]),digits=2),2,format="f")
upper[i]<-formatC(round(solution[i]+qnorm(.975)*sqrt(hessian[i]),digits=2),2,format="f")
#StatTest[i]<-round(solution[i]^2 / sqrt(hessian[i])^2 ,digits=2)
if(i<=e+w+ew){
if(missing(Wald_par)){
StatTest[i]<-formatC(round((solution[i]-1)^2 / hessian[i] ,digits=2),2,format="f")
Wald_H0[i]<-formatC(round(1,2),2,format="f")}
else {
StatTest[i]<-formatC(round((solution[i]-Wald_par[i])^2 / hessian[i] ,digits=2),2,format="f")
Wald_H0[i]<-formatC(round(Wald_par[i],2),2,format="f")}
}

else if(i>e+w+ew+nprob){
if(missing(Wald_par)){
StatTest[i]<-formatC(round((solution[i]-0)^2 / hessian[i],digits=2),2,format="f")
Wald_H0[i]<-formatC(round(0,2),2,format="f")}
else{StatTest[i]<-formatC(round((solution[i]-Wald_par[i-nprob])^2 / hessian[i],digits=2),2,format="f")
Wald_H0[i]<-formatC(round(Wald_par[i-nprob],2),2,format="f")
} 
}else{
StatTest[i]<-formatC(round((solution[i]-0)^2 / hessian[i],digits=2),2,format="f")}

if(StatTest[i]!='-'){

pvalue[i]<-formatC(round(1-pchisq(as.numeric(StatTest[i]),1),digits=4),4,format="f")

if(pvalue[i]=="0.0000")pvalue[i]<-"<0.0001"}
}
else if(parameters$label[i]=="P" && hessian[i]>=0){
sd[i]<-'-'
lower[i]<-'-'
upper[i]<-'-'
StatTest[i]<-'-'
pvalue[i]<-'-'
Wald_H0[i]<-'-'

}else{sd[i]<-NaN
lower[i]<-NaN
upper[i]<-NaN
StatTest[i]<-NaN
pvalue[i]<-NaN
Wald_H0[i]<-NaN}}




if(w==0 && ew==0){
table.dist1<-data.frame(Type="dist",Index=c(1:e),Transition=as.character(transitions[1:e]),Estimation=as.character(round(solution[1:e],digits=3)),SD=as.character(sd[1:e]),Lower_CI=as.character(lower[1:e]),Upper_CI=as.character(upper[1:e]),Wald_H0 =as.character(Wald_H0[1:e]),Wald_test=as.character(StatTest[1:e]), p_value=as.character(pvalue[1:e]))
table.dist<-list(Sigma=table.dist1)

}
if(w>0 && ew==0){
table.dist1<-data.frame(Type="dist",Index=c(1:e),Transition=as.character(transitions[1:e]),Sigma=as.character(round(solution[1:e],digits=3)),SD=as.character(sd[1:e]),Lower_CI=as.character(lower[1:e]),Upper_CI=as.character(upper[1:e]),Wald_H0 =as.character(Wald_H0[1:e]),Wald_test=as.character(StatTest[1:e]),p_value=as.character(pvalue[1:e]))
table.dist2<-data.frame(Type="dist",Index=c((e+1):(e+w)),Transition=as.character(transitions[(e+1):(e+w)]),Nu=as.character(round(solution[(e+1):(e+w)],digits=3)),SD=as.character(sd[(e+1):(e+w)]),Lower_CI=as.character(lower[(e+1):(e+w)]),Upper_CI=as.character(upper[(e+1):(e+w)]),Wald_H0 =as.character(Wald_H0[(e+1):(e+w)]),Wald_test=as.character(StatTest[(e+1):(e+w)]),p_value=as.character(pvalue[(e+1):(e+w)]))

table.dist<-list(Sigma=table.dist1,Nu=table.dist2)
}
if(w>0 && ew>0){
table.dist1<-data.frame(Type="dist",Index=c(1:e),Transition=as.character(transitions[1:e]),Sigma=as.character(round(solution[1:e],digits=3)),SD=as.character(sd[1:e]),Lower_CI=as.character(lower[1:e]),Upper_CI=as.character(upper[1:e]),Wald_H0 =as.character(Wald_H0[1:e]),Wald_test=as.character(StatTest[1:e]),p_value=as.character(pvalue[1:e]))
table.dist2<-data.frame(Type="dist",Index=c((e+1):(e+w)),Transition=as.character(transitions[(e+1):(e+w)]),Nu=as.character(round(solution[(e+1):(e+w)],digits=3)),SD=as.character(sd[(e+1):(e+w)]),Lower_CI=as.character(lower[(e+1):(e+w)]),Upper_CI=as.character(upper[(e+1):(e+w)]),Wald_H0 =as.character(Wald_H0[(e+1):(e+w)]),Wald_test=as.character(StatTest[(e+1):(e+w)]),p_value=as.character(pvalue[(e+1):(e+w)]))

table.dist3<-data.frame(Type="dist",Index=c((e+w+1):(e+w+ew)),Transition=as.character(transitions[(e+w+1):(e+w+ew)]),Theta=as.character(round(solution[(e+w+1):(e+w+ew)],digits=3)),SD=as.character(sd[(e+w+1):(e+w+ew)]),Lower_CI=as.character(lower[(e+w+1):(e+w+ew)]),Upper_CI=as.character(upper[(e+w+1):(e+w+ew)]),Wald_H0 =as.character(Wald_H0[(e+w+1):(e+w+ew)]),Wald_test=as.character(StatTest[(e+w+1):(e+w+ew)]),p_value=as.character(pvalue[(e+w+1):(e+w+ew)]))
table.dist<-list(Sigma=table.dist1,Nu=table.dist2,Theta=table.dist3)
}
trans<-matrix(ncol=length(states),nrow=length(states))
colnames(trans)<-states
rownames(trans)<-states
for(i in 1:length(states)){
for(j in 1:length(states)){
if(mtrans[i,j]==FALSE)trans[i,j]<-"-"
if(mtrans[i,j]=="E" || mtrans[i,j]=="Exp" || mtrans[i,j]=="Exponential")trans[i,j]<-"Exponential"
if(mtrans[i,j]=="W" || mtrans[i,j]==TRUE || mtrans[i,j]=="Weibull")trans[i,j]<-"Weibull"
if(mtrans[i,j]=="EW" || mtrans[i,j]=="EWeibull" || mtrans[i,j]=="Exponentiated Weibull")trans[i,j]<-"Exponentiated Weibull"
}}

semiMarkovobject<-list(minus2loglik=objective[-1],
	solution=data.frame(label=parameters[,1],transition=parameters[,2],solution=solution),
	opt.iter=iterations,opt.message=res$message,call=match.call(),nstates=length(states),Transition_matrix=trans,
	table.dist=table.dist,param.init=parameters,time.ind=time.ind)


##################
#we look for the probabilities of transitions that were not optimised as the may calculated 1 - sum of the proba in a row
#################

if(nprob>0){
transitionsP<-paste(substring(parameters[(e+w+ew+1):(e+w+ew++nprob),2],first=1,last=1),"->",substring(parameters[(e+w+ew++1):(e+w+ew++nprob),2],first=2,last=2))
trP<-paste(substring(parameters[(e+w+ew++1):(e+w+ew++nprob),2],first=1,last=1),substring(parameters[(e+w+ew++1):(e+w+ew++nprob),2],first=2,last=2),sep="")

which<-c()
for(i in 1:length(trP)){
a<-which(trans.hj==trP[i])
which<-c(which,a)}


trP<-trans.hj[-which]
transitionsP<-c(transitionsP,paste(substring(trP[1:length(trP)],first=1,last=1),"->",substring(trP[1:length(trP)],first=2,last=2)))


################
# Statistics for the proba 1-pr
################

## variance
SD_last<-c()
pozycje<-which(substring(parameters[,1],first=1,last=1)=="P")
macierz<-round(ginv(hessian(V,solution))[pozycje,pozycje],4)


#how many states not absorbent
s_trans<-0
which_rows<-c()
for(i in 1:s){
if(sum(1*(trans.h==states[i]))>0){
s_trans<-s_trans+1
which_rows<-c(which_rows,i)}
}

for(i in 1:s_trans){
licz<-length(which(substring(parameters[(e+w+ew+1):(e+w+ew+nprob),2],first=1,last=1)==states[which_rows[i]]))

if(is.matrix(macierz) && licz>0){
if(all(diag(macierz[(1:licz),(1:licz)])>=0)){
SD_last<-c(SD_last,as.character(formatC(round(sqrt(sum(as.matrix(macierz[(1:licz),(1:licz)]))),2),2,format="f")))
    }
else{SD_last<-c(SD_last,NaN)}
macierz<-macierz[-(1:licz),-(1:licz)]
}else if(licz>0){SD_last<-c(SD_last,macierz)}
else(SD_last <- c(SD_last, '-'))


}

#SD_last <- SD_last[1:length(unique(substring(parameters$transit[(e+w+ew+1):(e+w+ew+nprob)],first=1,last=1)))]
#SD_last<-round(sqrt(SD_last),2)
#if(length(trP)>length(SD_last)){
#  SD_addit <- length(trP)#-length(SD_last)
#}else{
#  SD_addit <- 0
#}
############
# table proba
#############

table.probabilities<-data.frame(Index=c(1:nprob,rep('-',ss)),Transition=transitionsP,Probability=c(solution[(e+w+ew+1):(e+w+ew+nprob)],rep(0,length(trP))),SD=c(sd[(e+w+ew+1):(e+w+ew+nprob)],SD_last),
				Lower_CI=rep('-',nprob+length(trP)),Upper_CI=rep('-',nprob+length(trP)),Wald_H0=rep('-',nprob+length(trP)),Wald_test=rep('-',nprob+length(trP)),
				p_value=rep('-',nprob+length(trP)))


j<-1
#pr<-table.probabilities[1,3]

if(nprob>1){
  table.probabilities<-table.probabilities[ order(table.probabilities$Transition), ]
  j <- 0
  pr <- 0
  k <- 1
for(i in 1:dim(table.probabilities)[1]){
if(i<dim(table.probabilities)[1]){
  if(substring(table.probabilities[i,2],first=1,last=1)==substring(table.probabilities[i+1,2],first=1,last=1)){
pr<-pr+table.probabilities[i,3]
j <- j+1
}else{
table.probabilities[i,3]<-1-pr
table.probabilities$SD[i]<-noquote(SD_last[k])
table.probabilities$Lower_CI[i]<-'-'
table.probabilities$Upper_CI[i]<-'-'
table.probabilities$Wald_test[i]<-'-' 
table.probabilities$p_value[i]<-'-'
pr<-0#as.numeric(table.probabilities[i+1,3])
k <- k+1
j<-0}
}else if(i==dim(table.probabilities)[1]){

table.probabilities[i,3]<-1-pr
table.probabilities$SD[i]<-noquote(SD_last[k])
table.probabilities$Index[i]<-'-'
table.probabilities$Lower_CI[i]<-'-'
table.probabilities$Upper_CI[i]<-'-'
table.probabilities$Wald_H0[i]<-'-'
table.probabilities$p_value[i]<-'-'}}

table.probabilities2 <- table.probabilities
table.probabilities<-data.frame(Type=rep("proba",dim(table.probabilities2)[1]),Index=table.probabilities2$Index,Transition=table.probabilities2[,2],Probability=table.probabilities2[,3],SD=table.probabilities2[,4],
				Lower_CI=table.probabilities2[,5],Upper_CI=table.probabilities2[,6],Wald_H0=table.probabilities2[,7],Wald_test=table.probabilities2[,8],p_value=table.probabilities2[,9])

}
else{
table.probabilities<-table.probabilities[ order(table.probabilities$Transition), ]
table.probabilities2<-table.probabilities[ order(table.probabilities$Transition), ]


for(i in 1:(nrow(table.probabilities2)-1)){

if(table.probabilities$Probability[i]!=0){
table.probabilities2$Probability[i+1]<-1-table.probabilities2$Probability[i]}
}

for(i in 1:nrow(table.probabilities2)){
if(table.probabilities2$Probability[i]==0){
table.probabilities2$Probability[i]<-1}}
table.probabilities<-data.frame(Type=rep("proba",dim(table.probabilities2)[1]),Index=table.probabilities2$Index, Transition=table.probabilities2[,2],Probability=table.probabilities2[,3],SD=table.probabilities2[,4],
				Lower_CI=table.probabilities2[,5],Upper_CI=table.probabilities2[,6],Wald_H0=table.probabilities2[,7],Wald_test=table.probabilities2[,8],p_value=table.probabilities2[,9])

}


semiMarkovobject$table.proba<-table.probabilities


SD<-c(sd,table.probabilities$SD[(nprob+1):dim(table.probabilities)[1]])}
else{
table.probabilities<-data.frame(Type=c(),Index=c(),Transition=c(),Probability=c(),SD=c(),
				Lower_CI=c(),Upper_CI=c(),Wald_H0=c(),Wald_test=c(),p_value=c())
semiMarkovobject$table.probabilities<-table.probabilities
SD<-sd
}

if(length(r)>1 || r>0){


if(nprob>0){
	semiMarkovobject$table.coef<-data.frame(Type=rep("coef",nbeta),Index=1:nbeta,Transition=paste(substring(parameters[(e+w+ew+nprob+1):nrow(parameters),2],first=1,last=1),"->",substring(parameters[(e+w+ew+nprob+1):nrow(parameters),2],first=2,last=2)),
		Covariates=parameters[(e+w+ew++nprob+1):(dim(parameters)[1]),1],
		Estimation=solution[(e+w+ew+nprob+1):dim(parameters)[1]],SD=sd[(e+w+ew+nprob+1):(dim(parameters)[1])],
		Lower_CI=lower[(e+w+ew+nprob+1):(dim(parameters)[1])],Upper_CI=upper[(e+w+ew+nprob+1):(dim(parameters)[1])],
		Wald_H0 = Wald_H0[(e+w+ew+nprob+1):dim(parameters)[1]],
		Wald_test=StatTest[(e+w+ew+nprob+1):dim(parameters)[1]],
		p_value=pvalue[(e+w+ew+nprob+1):(dim(parameters)[1])])
}else{
semiMarkovobject$table.coef<-data.frame(Type=rep("coef",nbeta),Index=1:nbeta,Transition=paste(substring(parameters[(e+w+ew+nprob+1):nrow(parameters),2],first=1,last=1),"->",substring(parameters[(e+w+ew+nprob+1):nrow(parameters),2],first=2,last=2)),
		Covariates=parameters[(e+w+ew+nprob+1):(dim(parameters)[1]),1],
		Estimation=solution[(e+w+ew+nprob+1):dim(parameters)[1]],SD=sd[(e+w+ew+nprob+1):(dim(parameters)[1])],
		Lower_CI=lower[(e+w+ew+nprob+1):(dim(parameters)[1])],Upper_CI=upper[(e+w+ew+nprob+1):(dim(parameters)[1])],
		Wald_H0 = Wald_H0[(e+w+ew+nprob+1):dim(parameters)[1]],
		Wald_test=StatTest[(e+w+ew+nprob+1):dim(parameters)[1]],
		p_value=pvalue[(e+w+ew+nprob+1):(dim(parameters)[1])])

}
}else{
semiMarkovobject$table.coef<-c()}


# Possible transitions matrix
#############
colnames(mtrans)<-1:dim(mtrans)[2]
rownames(mtrans)<-1:dim(mtrans)[1]
semiMarkovobject$transition.matrix<-mtrans
semiMarkovobject$nproba<-nprob
#semiMarkovobject$pos.trans<-data.frame(from=paste(substring(trans.hj,first=1,last=1)),
#	to=paste(substring(trans.hj,first=2,last=2)))
semiMarkovobject$trans.hj<-trans.hj
semiMarkovobject$last<-last
###############
#table.parameters

if(length(r)>1 || r>0){
if(nprob>0){
 
  
table.param<-data.frame(Type=c(rep("dist",ndist),rep("proba",p),rep("coef",nbeta)),Index=c(1:ndist,as.character(table.probabilities$Index),1:nbeta),Label=c(as.character(parameters[1:(e+w+ew),1]),rep("P",p),as.character(parameters[(e+w+ew+nprob+1):dim(parameters)[1],1])),
	Transition=c(transitions[1:(e+w+ew)],as.character(table.probabilities$Transition),transitions[(e+w+ew+nprob+1):length(transitions)]),
	Estimation=round(c(solution[1:(e+w+ew)],table.probabilities$Probability,solution[(e+w+ew+nprob+1):length(solution)]),digits=3),
	SD=c(as.character(SD[1:(e+w+ew)]),as.character(table.probabilities$SD),as.character(SD[(e+w+ew+nprob+1):length(solution)])),
	Lower_CI=c(as.character(lower[1:(e+w+ew)]),as.character(table.probabilities$Lower_CI),as.character(lower[(e+w+ew+nprob+1):length(solution)])),
	Upper_CI=c(as.character(upper[1:(e+w+ew)]),as.character(table.probabilities$Upper_CI),as.character(upper[(e+w+ew+nprob+1):length(solution)])),
	Wald_H0=c(as.character(Wald_H0[1:(e+w+ew)]),as.character(table.probabilities$Wald_H0),as.character(Wald_H0[(e+w+ew+nprob+1):length(solution)])),	
	Wald_statistic=c(as.character(StatTest[1:(e+w+ew)]),as.character(table.probabilities$Wald_test),as.character(StatTest[(e+w+ew+nprob+1):length(solution)])),
	Wald_p.value=c(as.character(pvalue[1:(e+w+ew)]),as.character(table.probabilities$p_value),as.character(pvalue[(e+w+ew+nprob+1):length(solution)])))
}else{

table.param<-data.frame(Type=c(rep("dist",ndist),rep("coef",nbeta)),Index=c(1:ndist,1:nbeta),Label=c(as.character(parameters[1:(e+w+ew),1]),as.character(parameters[(e+w+ew+1):dim(parameters)[1],1])),
	Transition=transitions,
	Estimation=round(c(solution[1:length(solution)]),digits=3),
	SD=c(as.character(SD[1:length(solution)])),
	Lower_CI=c(as.character(lower[1:(e+w+ew)]),as.character(lower[(e+w+ew+1):length(solution)])),
	Upper_CI=c(as.character(upper[1:(e+w+ew)]),as.character(upper[(e+w+ew+1):length(solution)])),
	Wald_H0 = c(as.character(Wald_H0[1:(e+w+ew)]),as.character(Wald_H0[(e+w+ew+nprob+1):length(solution)])),
	Wald_statistic=c(as.character(StatTest[1:(e+w+ew)]),as.character(StatTest[(e+w+ew+nprob+1):length(solution)])),
	Wald_p.value=c(as.character(pvalue[1:(e+w+ew)]),as.character(pvalue[(e+w+ew+nprob+1):length(solution)])))

}}
else{
if(nprob>0){

table.param<-data.frame(Type=c(rep("dist",ndist),rep("proba",p)),Index=c(1:ndist,as.character(table.probabilities$Index)),Label=c(as.character(parameters[1:(e+w+ew),1]),rep("P",p)),
	Transition=c(transitions[1:(e+w+ew)],as.character(table.probabilities$Transition)),
	Estimation=round(c(solution[1:(e+w+ew)],table.probabilities$Probability),digits=3),
	SD=c(as.character(SD[1:(e+w+ew)]),as.character(table.probabilities$SD)),
	Lower_CI=c(as.character(lower[1:(e+w+ew)]),as.character(table.probabilities$Lower_CI)),
	Upper_CI=c(as.character(upper[1:(e+w+ew)]),as.character(table.probabilities$Upper_CI)),
	Wald_H0 = c(as.character(Wald_H0[1:(e+w+ew)]),as.character(table.probabilities$Wald_H0)),
	Wald_statistic=c(as.character(StatTest[1:(e+w+ew)]),as.character(table.probabilities$Wald_test)),
	Wald_p.value=c(as.character(pvalue[1:(e+w+ew)]),as.character(table.probabilities$p_value)))
	
}
else{
table.param<-data.frame(Type=c(rep("dist",ndist)),Index=c(1:ndist),Label=as.character(parameters[1:(e+w+ew),1]),
	Transition=transitions[1:(e+w+ew)],
	Estimation=round(c(solution[1:length(solution)]),digits=3),
	SD=c(as.character(SD[1:length(solution)])),
		Lower_CI=c(as.character(lower[1:(e+w+ew)]),as.character(table.probabilities$Lower_CI)),
	Upper_CI=c(as.character(upper[1:(e+w+ew)]),as.character(table.probabilities$Upper_CI)),
	Wald_H0 = c(as.character(Wald_H0[1:(e+w+ew)]),as.character(table.probabilities$Wald_H0)),
	Wald_statistic=c(as.character(StatTest[1:(e+w+ew)]),as.character(table.probabilities$Wald_test)),
	Wald_p.value=c(as.character(pvalue[1:(e+w+ew)]),as.character(table.probabilities$p_value)))
}

}
if(length(which(table.param$SD=="NaN"))>0)warning("NaN for standard deviation estimation caused by negative values in the hessian matrix. Consider modifying the optional constraints on parameters.")
semiMarkovobject$states<-states
dimnames(table.state)<-list(states,states)
semiMarkovobject$table.state<-table.state
semiMarkovobject$Ncens<-Ncens
semiMarkovobject$table.param<-table.param
class(semiMarkovobject)<-"semiMarkov"
semiMarkovobject



}

#############################################
print.semiMarkov<-function(x,CI=TRUE,Wald.test=TRUE,...){
if (!inherits(x, "semiMarkov")) 
        stop("'x' must be of class 'semiMarkov'")
if(!is.logical(CI))
	stop("'CI' must be 'logical'")
if(!is.logical(Wald.test))
	stop("'Wald.test' must be 'logical'")

cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
cat("Multi-State Semi-Markov Model for ", x$nstates," states.\n\n",sep="")
cat("Transition matrix\n")
print(x$Transition_matrix)
cat("\n-2 * log-likelihood: ", x$minus2loglik, "\n\n")
cat("Estimated parameters:\n")
if(CI==TRUE && Wald.test==TRUE){
print(x$table.param)}
else if(CI==TRUE && Wald.test==FALSE){
print(x$table.param[,1:6])
}else if(CI==FALSE && Wald.test==TRUE){
print(x$table.param[,c(1,2,3,4,7,8)])
}else{
print(x$table.param[,1:4])}

error<-which(x$table.param[,4]=="NaN")

if(length(error)>=1)print(x$opt.message)
}

#######################################################
summary.semiMarkov<-function(object,all=TRUE,transitions=NULL,...){
 if (!inherits(object, "semiMarkov")) 
        stop("'object' must be of class 'semiMarkov'")
res<-list()
res$param.init<-object$param.init
dimnames(object$table.state)<-list(object$states,object$states)
res$table.state<-object$table.state
res$Ncens<-object$Ncens
if(all==FALSE && is.vector(transitions)){
if(!all(transitions %in% object$trans.hj))
	stop("'transitions' must be a vector of possible transitions")
tables<-list()
for(i in 1:length(transitions)){
name<-paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2),"\n",object$Transition_matrix[substring(transitions[i],first=1,last=1),substring(transitions[i],first=2,last=2)],"distribution\n")
distribution<-paste(object$Transition_matrix[substring(transitions[i],first=1,last=1),substring(transitions[i],first=2,last=2)],"distribution")
table<-data.frame(Parameter=object$table.param[which(object$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2))),1],
	Estimation=object$table.param[which(object$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2))),3],
	SD=object$table.param[which(object$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2))),4],
	Lower_CI=object$table.param[which(object$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2))),5],
	Upper_CI=object$table.param[which(object$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2))),6],
	Wald_H0 = object$table.param[which$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2)),7],
	Wald_statistic=object$table.param[which$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2)),6],
	Wald_p.value=object$table.param[which(object$table.param[,2]==paste(substring(transitions[i],first=1,last=1),"->",substring(transitions[i],first=2,last=2))),9])

tables$name<-list(tables,table)
}
}else {
if(all==FALSE && missing(transitions))
	warning("The transitions need to be specified. All the possible transitions are displayed.")
if(all==TRUE && length(transitions)>0)
	warning("Option for the transitions was chosen.")
if(!all(all %in% c(TRUE,FALSE)))
	stop("'all' must be of logical value")
tables<-list()
for(i in 1:nrow(object$table.dist[[1]])){
name<-paste(as.character(object$table.dist[[1]][i,1]),"; ",object$Transition_matrix[as.numeric(substring(object$table.dist[[1]][i,1],first=1,last=1)),as.numeric(substring(object$table.dist[[1]][i,1],first=6,last=6))],"distribution")
table<-data.frame(Parameter=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),1],
	Estimation=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),3],
	SD=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),4],
	Lower_CI=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),5],
	Upper_CI=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),6],
	Wald_H0=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),7],
	Wald_statistic=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),8],
	Wald_p.value=object$table.param[which(object$table.param[,2]==object$table.dist[[1]][i,1]),9]
)

tables[[name]]<-table
}
}
res$table.param<-tables

class(res)<-"summary.semiMarkov"
res
}


########################################################################################################################

################################
##distributions
#################################



#############
#Weibull
#####################


.densityW<-function(t,cova,cova.mat,sigma,nu,beta=0)
{
if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(c(beta),ncol=1))
}else{z<-1}

result<-nu*((1/sigma)^nu)*(t^(nu-1))*z*exp(-z*(t/sigma)^nu)
change<-which(result<1e-320)
result[change]<-1e-320
changeInf<-which(is.nan(result))
result[changeInf]<-1e+307
result
}

.densityW2<-function(t,cova,cova.mat,sigma,nu,beta=0)
{
if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(c(beta),ncol=1))
}else{z<-1}
small<-0

nu*((1/sigma)^nu)*(t^(nu-1))*z*exp(-z*(t/sigma)^nu)
}
.survivalW<-function(t,cova,cova.mat,sigma,nu,beta=0)
{if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(c(beta),ncol=1))
}else{z<-1}
 
 result<-exp((-(t/sigma)^nu)*z)
change<-which(result<1e-320)
result[change]<-1e-320
changeInf<-which(is.nan(result))
result[changeInf]<-1e+307
result

}


.hazardW<-function(t,cova, cova.mat,sigma,nu,beta)
{
if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(c(beta),ncol=1))
}else{z<-1}
 result<-nu*((1/sigma)^nu)*(t^(nu-1))*z
change<-which(result<1e-320)
result[change]<-1e-320
changeInf<-which(is.nan(result))
result[changeInf]<-1e+307
result
}


#############
# Exponentiated Weibull 
#################################

.densityEW<-function(t,cova,cova.mat,sigma,nu,theta,beta=0)
{
if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(c(beta),ncol=1))
}else{z<-1}

a<-exp(-(t/sigma)^nu)
result<-((nu*theta*((1-a)^(theta-1))*a*((t/sigma)^(nu-1))*z)/sigma)*(1-(1-a)^theta)^(z-1)
change<-which(result<1e-3200)
result[change]<-1e-320
changeInf<-which(is.nan(result))
result[changeInf]<-1e+307
result

}

.survivalEW<-function(t,cova,cova.mat,sigma,nu,theta,beta=0)
{if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(beta,ncol=1))
}else{z<-1}

a<-exp(-(t/sigma)^nu)
result<-((1-((1-a)^theta))^z)
change<-which(result<1e-320)
result[change]<-1e-320

changeInf<-which(is.nan(result))
result[changeInf]<-1e+307
result

}

.hazardEW<-function(t,cova,cova.mat,sigma,nu,theta,beta=0)
{
if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(beta,ncol=1))
}else{z<-1}
a<-exp(-(t/sigma)^nu)

for(i in 1:length(a)){
if(a[i]<1e-320)a[i]<-1e-320
}

result<-((nu*theta*((1-a)^(theta-1))*a*((t/sigma)^(nu-1)))/(sigma*(1-(1-a)^theta)))*z
change<-which(result<1e-320)
result[change]<-1e-320
changeInf<-which(is.nan(result))
result[changeInf]<-1e+307

result

}

#############
#  Exponential
##########################
.densityE<-function(t,cova,cova.mat,sigma,beta=0)
{
if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(beta,ncol=1))
}else{z<-1}
change<-which(z==Inf)
z[change]<-1e+307
result<-(1/sigma)*z*(exp(-t/sigma)^z)
change<-which(result<1e-320)
result[change]<-1e-320
changeInf<-which((1/sigma)*z==Inf)
result[changeInf]<-1e+307
result

}
.densityE2<-function(t,cova,cova.mat,sigma,beta=0)
{
if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(beta,ncol=1))
}else{z<-1}
change<-which(z==Inf)
z[change]<-1e+307
result<-(1/sigma)*z*(exp(-t/sigma)^z)

change<-which(result<1e-320)
result[change]<-1e-320
changeInf<-which((1/sigma)*z==Inf)
result[changeInf]<-1e+307
result

}
.survivalE<-function(t,cova,cova.mat,sigma,beta=0)
{if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(beta,ncol=1))
}else{z<-1}
change<-which(z==Inf)
z[change]<-1e+307
result<-exp(-t/sigma)^z

change<-which(result<1e-320)
result[change]<-1e-320
changeInf<-which((1/sigma)*z==Inf)
result[changeInf]<-1e+307
result

}
.survivalE2<-function(t,cova,cova.mat,sigma,beta=0)
{if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(beta,ncol=1))
}else{z<-1}
result<-exp(-t/sigma)^z

change<-which(result==0)
result[change]<-1e-320

result

}
.hazardE<-function(t,cova,cova.mat,sigma,beta=0)
{if(cova!=FALSE){z<-exp(as.matrix(cova.mat)%*%matrix(beta,ncol=1))
}else{z<-1}
result<-(1/sigma)^z
change<-which(result==0)
result[change]<-1e-320

result

}


#############################
# likelihood
############################

##########################
#Censored part of the Likelihood
##########################


.marg<-function(s,trans,d,t,covariates, cova,x,parameters)
{

##number of transitions
p<-length(which(trans!=FALSE))
## ndist
ndist<-length(which(parameters[,1]%in%c("sigma","nu","theta")))
# number of parameters e, w ,ew
ew<-length(which(parameters[,1]=="theta"))
w<-length(which(parameters[,1]=="nu"))
e<-length(which(parameters[,1]=="sigma"))

#marginal contribution to likelihood
marg_row<-c()
#auxillary indices
k2<-0
k3<-0
k<-0
m<-0
#loop for the rows
for(i in 1:(dim(trans)[1])){

#number of transitions in a row i
tr_in_row<-length(which(trans[i,]!=FALSE))

#marginal contribution for a row
marg_row_el<-rep(0,length(t))

#cumulated probability in a row
pr<-0


#Check if the transition occurs as a parameter for the probability

if(tr_in_row>1 && ndist+i+m<=dim(parameters)[1] && substring(parameters[i+k,2],first=1,last=1)==substring(parameters[ndist+i+m,2],first=1,last=1)&& parameters[ndist+i+m,1]=="P")
{

for(j in 0:(tr_in_row-1)){
 
  #vector of covariates for a transition
beta_init<-c()
for(l in 0:(dim(cova)[2]-1)){beta_init[length(beta_init+1)+1]<-c(x[ndist+s+k+i+j+p*l])
}


#marginal contribution for a row if the transition is not the last in the row
if(j!=(tr_in_row-1)){
if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalE(t,covariates,cova, sigma=x[i+j+k],beta=beta_init)))
}
if(d[i+j+k]==2){
#Weibull

k2<-k2+1
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))

}
if(d[i+j+k]==3){
#Exponentiated Weibull
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalEW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))
}

pr<-pr+x[ndist+i+j+m]
if(pr>1){
pr<-0.999
}

}else{
  
if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+((1-pr)*.survivalE(t,covariates,cova,  sigma=x[i+j+k],beta=beta_init)))

pr<-0

}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
marg_row_el<-c(marg_row_el+((1-pr)*.survivalW(t,covariates,cova,  sigma=x[i+j+k],nu=x[e+i+j+k],beta=beta_init)))

pr<-0
}
if(d[i+j+k]==3){
##Exponentiated Weibull
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el+((1-pr)*.survivalEW(t,covariates,cova,  sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))
pr<-0
}


f<-j
g<-j-1
}}
k<-k+f

m<-m+g


marg_row<-c(marg_row,marg_row_el)}
else if(tr_in_row>0){

##If the probability is not a parameter, so it always equals to 1

for(j in 0:(tr_in_row-1)){
beta_init<-c()
for(l in 0:(dim(cova)[2]-1)){beta_init[length(beta_init+1)+1]<-c(x[ndist+s+k+i+j+p*l])
}

if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+.survivalE(t,covariates,cova,  sigma=x[i+j+k],beta=beta_init))}
if(d[i+j+k]==2){
##Weibull
marg_row_el<-c(marg_row_el+.survivalW(t,covariates,cova,  sigma=x[i+j+k],nu=x[e+k2],beta=beta_init))
}
if(d[i+j+k]==3){
##Exponentiated Weibull

marg_row_el<-c(marg_row_el+.survivalEW(t,covariates,cova,  sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init))

}

f<-j
g<-j-1
}

k<-k+f
m<-m+g

marg_row<-c(marg_row,marg_row_el)}else{m<-m-1
k <- k-1}

}

matrix(marg_row,ncol=1)

}
#############################
#Density part of the Likelihood
#############################



.dens<-function(s,trans,d,t,covariates,cova, x,parameters)
{
##number of transitions
p<-length(which(trans!=FALSE))

## ndist
ndist<-length(which(parameters[,1]%in%c("sigma","nu","theta")))

# number of parameters e, w ,ew
ew<-length(which(parameters[,1]=="theta"))
w<-length(which(parameters[,1]=="nu"))
e<-length(which(parameters[,1]=="sigma"))

dens_row<-c()
k2<-0
k3<-0
k<-0
m<-0
#i for the rows
for(i in 1:(dim(trans)[1])){
#how many transitions in a row
tr_in_row<-length(which(trans[i,]!=FALSE))

pr<-0


#Check if the transition occurs as a parameter for the probability
if(tr_in_row>1 && ndist+i+m<=dim(parameters)[1] && substring(parameters[i+k,2],first=1,last=1)==substring(parameters[ndist+i+m,2],first=1,last=1)&& parameters[ndist+i+m,1]=="P")
{
#j for the columns
for(j in 0:(tr_in_row-1)){
beta_init<-c()
for(l in 0:(dim(cova)[2]-1)){beta_init[length(beta_init+1)+1]<-c(x[ndist+s+k+i+j+p*l])
}

if(j!=(tr_in_row-1)){
if(d[i+j+k]==1){
##Exponential
dens_row_el<-c((x[ndist+i+j+m]*.densityE(t,covariates,cova, sigma=x[i+j+k],beta=beta_init)))

}
if(d[i+j+k]==2){
#Weibull
k2<-k2+1
dens_row_el<-c((x[ndist+i+j+m]*.densityW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))
}
if(d[i+j+k]==3){
#Exponentiated Weibull
k2<-k2+1
k3<-k3+1
dens_row_el<-c((x[ndist+i+j+m]*.densityEW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))


}

pr<-pr+x[ndist+i+j+m]
if(pr>1){
pr<-0.999}
dens_row<-c(dens_row,dens_row_el)

}else{

if(d[i+j+k]==1){
##Exponential
dens_row_el<-c(((1-pr)*.densityE(t,covariates,cova, sigma=x[i+j+k],beta=beta_init)))
pr<-0
}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
dens_row_el<-c(((1-pr)*.densityW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))
pr<-0

}
if(d[i+j+k]==3){
##Exponentiated Weibull
k3<-k3+1
k2<-k2+1
dens_row_el<-c(((1-pr)*.densityEW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))
pr<-0
}


dens_row<-c(dens_row,dens_row_el)

f<-j
g<-j-1
}}
k<-k+f
m<-m+g}else if(tr_in_row>0){
##If the probability is not a parameter, so it always equals to 1
if(tr_in_row>0){
for(j in 0:(tr_in_row-1)){
beta_init<-c()
for(l in 0:(dim(cova)[2]-1)){beta_init[length(beta_init+1)+1]<-c(x[ndist+s+k+i+j+p*l])
}

if(d[i+j+k]==1){
##Exponential
dens_row_el<-c(.densityE(t,covariates,cova, sigma=x[i+j+k],beta=beta_init))
}
if(d[i+j+k]==2){
##Weibull
  k2<-k2+1
dens_row_el<-c(.densityW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],beta=beta_init))


}
if(d[i+j+k]==3){
##Exponentiated Weibull
  k2<-k2+1
  k3<-k3+1
dens_row_el<-c(.densityEW(t,covariates,cova, sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init))
}


dens_row<-c(dens_row,dens_row_el)

f<-j
g<-j-1
}}
k<-k+f
m<-m+g

}else{
  m<-m-1
k <- k-1}


}

matrix(dens_row,ncol=1)

}






###########################################################################################################################################

###################################
#likelihood.transitions
###################################

##########################
#Censored part of the Likelihood
##########################



.marg_trans<-function(s,trans,d,t,covariates,cova, x,cov_pos,parameters)
{

##number of transitions
p<-length(which(trans!=FALSE))

#Total number of covariates
dim<-dim(cova)[1]

## ndist
ndist<-length(which(parameters[,1]%in%c("sigma","nu","theta")))

# number of parameters e, w ,ew
ew<-length(which(parameters[,1]=="theta"))
w<-length(which(parameters[,1]=="nu"))
e<-length(which(parameters[,1]=="sigma"))

#marginal contribution to likelihood
marg_row<-c()

#auxillary indices
k2<-0
k3<-0
k<-0
n<-0
m<-0
#loop for the rows
for(i in 1:(dim(trans)[1])){

#number of transitions in a row i
tr_in_row<-length(which(trans[i,]!=FALSE))

#marginal contribution for a row
marg_row_el<-rep(0,length(t))

#cumulated probability in a row
pr<-0


#Check if the transition occurs as a parameter for the probability

if(tr_in_row>1 && ndist+i+m<=dim(parameters)[1] && substring(parameters[i+k,2],first=1,last=1)==substring(parameters[ndist+i+m,2],first=1,last=1)&& parameters[ndist+i+m,1]=="P")
{

for(j in 0:(tr_in_row-1)){
#vector of covariates for a transition
#Check if we consider this transition for the covariates
if(all((i+j+k) %in% cov_pos)){
transition<-parameters[i+j+k,2]

which_beta<-intersect(which(parameters[,2]==transition),which(substring(parameters[,1],first=1,last=1)=="B"))
beta_init<-x[which_beta]
which_cova<-as.numeric(substring(parameters[which_beta,1],first=5,last=5))

cova_temp<-cova[,which_cova]
covariates_temp<-length(beta_init)
}else{
#If we do not consider this covariates we do not have the coefficient Beta
beta_init<-NA
cova_temp<-rep(0,dim)
covariates_temp<-0}

#marginal contribution for a row if the transition is not the last in the row
if(j!=(tr_in_row-1)){

if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalE(t,covariates_temp,cova_temp,  sigma=x[i+j+k],beta=beta_init)))
}
if(d[i+j+k]==2){
#Weibull
k2<-k2+1
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))

}
if(d[i+j+k]==3){
#Exponentiated Weibull
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalEW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))

}

pr<-pr+x[ndist+i+j+m]
if(pr>1){
pr<-0.999}

}else{

if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+((1-pr)*.survivalE(t,covariates_temp,cova_temp,   sigma=x[i+j+k],beta=beta_init)))}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
marg_row_el<-c(marg_row_el+((1-pr)*.survivalW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))
}
if(d[i+j+k]==3){
##Exponentiated Weibull
k3<-k3+1
k2<-k2+1
marg_row_el<-c(marg_row_el+((1-pr)*.survivalEW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))

}


f<-j
g<-j-1
}}
k<-k+f
m<-m+g
marg_row<-c(marg_row,marg_row_el)}
else if(tr_in_row>0){

##If the probability is not a parameter, so it always equals to 1

for(j in 0:(tr_in_row-1)){
#vector of covariates for a transition 
beta_init<-c()

#Check if we consider this transition for the covariates

if(all((i+j+k) %in% cov_pos)){
transition<-parameters[i+j+k,2]
which_beta<-intersect(which(parameters[,2]==transition),which(substring(parameters[,1],first=1,last=1)=="B"))
beta_init<-x[which_beta]
which_cova<-as.numeric(substring(parameters[which_beta,1],first=5,last=5))

cova_temp<-cova[,which_cova]
covariates_temp<-length(beta_init)}else{
#If we do not consider this covariates we do not have the coefficient Beta
beta_init<-NA
cova_temp<-rep(0,dim)
covariates_temp<-0}
if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+.survivalE(t,covariates_temp,cova_temp,sigma=x[i+j+k],beta=beta_init))}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
marg_row_el<-c(marg_row_el+.survivalW(t,covariates_temp,cova_temp,sigma=x[i+j+k],nu=x[e+k2],beta=beta_init))
}
if(d[i+j+k]==3){
##Exponentiated Weibull
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el+.survivalEW(t,covariates_temp,cova_temp,sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init))

}

f<-j
g<-j-1
}

k<-k+f
m<-m+g

marg_row<-c(marg_row,marg_row_el)}else{k <- k-1
m <- m-1}
}
matrix(marg_row,ncol=1)
}



.marg2_trans<-function(s,trans,d,t,covariates,cova, x,cov_pos,parameters)
{

##number of transitions
p<-length(which(trans!=FALSE))

#Total number of covariates
dim<-dim(cova)[1]

## ndist
ndist<-length(which(parameters[,1]%in%c("sigma","nu","theta")))

# number of parameters e, w ,ew
ew<-length(which(parameters[,1]=="theta"))
w<-length(which(parameters[,1]=="nu"))
e<-length(which(parameters[,1]=="sigma"))

#marginal contribution to likelihood
marg_row<-c()

#auxillary indices
k2<-0
k3<-0
k<-0
n<-0
m<-0
#loop for the rows
for(i in 1:(dim(trans)[1])){

#number of transitions in a row i
tr_in_row<-length(which(trans[i,]!=FALSE))

#marginal contribution for a row
marg_row_el<-rep(0,length(t))

#cumulated probability in a row
pr<-0


#Check if the transition occurs as a parameter for the probability

if(tr_in_row>1 && ndist+i+m<=dim(parameters)[1] && substring(parameters[i+k,2],first=1,last=1)==substring(parameters[ndist+i+m,2],first=1,last=1)&& parameters[ndist+i+m,1]=="P")
{

for(j in 0:(tr_in_row-1)){
#vector of covariates for a transition
#Check if we consider this transition for the covariates
if(all((i+j+k) %in% cov_pos)){
transition<-parameters[i+j+k,2]

which_beta<-intersect(which(parameters[,2]==transition),which(substring(parameters[,1],first=1,last=1)=="B"))
beta_init<-x[which_beta]
which_cova<-as.numeric(substring(parameters[which_beta,1],first=5,last=5))

cova_temp<-cova[,which_cova]
covariates_temp<-length(beta_init)
}else{
#If we do not consider this covariates we do not have the coefficient Beta
beta_init<-NA
cova_temp<-rep(0,dim)
covariates_temp<-0}

#marginal contribution for a row if the transition is not the last in the row
if(j!=(tr_in_row-1)){

if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalE(t,covariates_temp,cova_temp,  sigma=x[i+j+k],beta=beta_init)))
}
if(d[i+j+k]==2){
#Weibull
k2<-k2+1
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))
}
if(d[i+j+k]==3){
#Exponentiated Weibull
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el+(x[ndist+i+j+m]*.survivalEW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))

}

pr<-pr+x[ndist+i+j+m]
if(pr>1){
pr<-0.999}

}else{

if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+((1-pr)*.survivalE(t,covariates_temp,cova_temp,   sigma=x[i+j+k],beta=beta_init)))

}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
marg_row_el<-c(marg_row_el+((1-pr)*.survivalW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))

}
if(d[i+j+k]==3){
##Exponentiated Weibull
k3<-k3+1
k2<-k2+1
marg_row_el<-c(marg_row_el+((1-pr)*.survivalEW(t,covariates_temp,cova_temp,  sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))

}


f<-j
g<-j-1
}}
k<-k+f
m<-m+g
marg_row<-c(marg_row,marg_row_el)}
else if(tr_in_row>0){

##If the probability is not a parameter, so it always equals to 1

for(j in 0:(tr_in_row-1)){
#vector of covariates for a transition 
beta_init<-c()

#Check if we consider this transition for the covariates

if(all((i+j+k) %in% cov_pos)){
transition<-parameters[i+j+k,2]
which_beta<-intersect(which(parameters[,2]==transition),which(substring(parameters[,1],first=1,last=1)=="B"))
beta_init<-x[which_beta]
which_cova<-as.numeric(substring(parameters[which_beta,1],first=5,last=5))

cova_temp<-cova[,which_cova]
covariates_temp<-length(beta_init)}else{
#If we do not consider this covariates we do not have the coefficient Beta
beta_init<-NA
cova_temp<-rep(0,dim)
covariates_temp<-0}
if(d[i+j+k]==1){
##Exponential
marg_row_el<-c(marg_row_el+.survivalE(t,covariates_temp,cova_temp,sigma=x[i+j+k],beta=beta_init))
}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
marg_row_el<-c(marg_row_el+.survivalW(t,covariates_temp,cova_temp,sigma=x[i+j+k],nu=x[e+k2],beta=beta_init))

}
if(d[i+j+k]==3){
##Exponentiated Weibull
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el+.survivalEW(t,covariates_temp,cova_temp,sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init))

}

f<-j
g<-j-1
}

k<-k+f
m<-m+g

marg_row<-c(marg_row,marg_row_el)}else{k <- k-1
m <- m-1}
}

matrix(marg_row,ncol=1)

}


#############################
#Density part of the Likelihood
#############################



.dens_trans<-function(s,trans,d,t,covariates,cova, x,cov_pos,parameters)
{

##number of transitions
p<-length(which(trans!=FALSE))

#number of coavariates
dim<-dim(cova)[1]

## ndist
ndist<-length(which(parameters[,1]%in%c("sigma","nu","theta")))

# number of parameters e, w ,ew
ew<-length(which(parameters[,1]=="theta"))
w<-length(which(parameters[,1]=="nu"))
e<-length(which(parameters[,1]=="sigma"))

dens_row<-c()
k2<-0
k3<-0
n<-0
k<-0
m<-0
for(i in 1:(dim(trans)[1])){

tr_in_row<-length(which(trans[i,]!=FALSE))

pr<-0


#Check if the transition occurs as a parameter for the probability

if(tr_in_row>1 && ndist+i+m<=dim(parameters)[1] && substring(parameters[i+k,2],first=1,last=1)==substring(parameters[ndist+i+m,2],first=1,last=1)&& parameters[ndist+i+m,1]=="P")
{

for(j in 0:(tr_in_row-1)){
#vector of covariates for a transition 

beta_init<-c()

#Check if we consider this transition for the covariates

if(all((i+j+k) %in% cov_pos)){
transition<-parameters[i+j+k,2]
which_beta<-intersect(which(parameters[,2]==transition),which(substring(parameters[,1],first=1,last=1)=="B"))
beta_init<-x[which_beta]
which_cova<-as.numeric(substring(parameters[which_beta,1],first=5,last=5))

cova_temp<-cova[,which_cova]
covariates_temp<-length(beta_init)


}else{

#If we do not consider this covariates we do not have the coefficient Beta

beta_init<-NA
cova_temp<-rep(0,dim)
covariates_temp<-0}

if(j!=(tr_in_row-1)){

if(d[i+j+k]==1){
##Exponential
dens_row_el<-c((x[ndist+i+j+m]*.densityE(t,covariates_temp,cova_temp,  sigma=x[i+j+k],beta=beta_init)))
}
if(d[i+j+k]==2){
#Weibull
k2<-k2+1
dens_row_el<-c((x[ndist+i+j+m]*.densityW(t,covariates_temp,cova_temp, sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))

}
if(d[i+j+k]==3){
#Exponentiated Weibull
k2<-k2+1
k3<-k3+1
dens_row_el<-c((x[ndist+i+j+m]*.densityEW(t,covariates_temp,cova_temp, sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))

}

pr<-pr+x[ndist+i+j+m]
if(pr>1){
pr<-0.999}
dens_row<-c(dens_row,dens_row_el)

}else{

if(d[i+j+k]==1){
##Exponential
dens_row_el<-c(((1-pr)*.densityE(t,covariates_temp,cova_temp, sigma=x[i+j+k],beta=beta_init)))

}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
dens_row_el<-c(((1-pr)*.densityW(t,covariates_temp,cova_temp, sigma=x[i+j+k],nu=x[e+k2],beta=beta_init)))
}
if(d[i+j+k]==3){
##Exponentiated Weibull
k2<-k2+1
k3<-k3+1
dens_row_el<-c(((1-pr)*.densityEW(t,covariates_temp,cova_temp, sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init)))

}

dens_row<-c(dens_row,dens_row_el)

f<-j
g<-j-1
}}
k<-k+f
m<-m+g}else if(tr_in_row>0){

##If the probability is not a parameter, so it always equals to 1
if(tr_in_row>0){
for(j in 0:(tr_in_row-1)){
beta_init<-c()

#Check if we consider this transition for the covariates

if(all((i+j+k) %in% cov_pos)){
transition<-parameters[i+j+k,2]
which_beta<-intersect(which(parameters[,2]==transition),which(substring(parameters[,1],first=1,last=1)=="B"))
beta_init<-x[which_beta]
which_cova<-as.numeric(substring(parameters[which_beta,1],first=5,last=5))
cova_temp<-cova[,which_cova]
covariates_temp<-length(beta_init)


}else{

#If we do not consider this covariates we do not have the coefficient Beta

beta_init<-NA
cova_temp<-rep(0,dim)
covariates_temp<-0}

if(d[i+j+k]==1){
##Exponential
dens_row_el<-c(.densityE(t,covariates_temp,cova_temp, sigma=x[i+j+k],beta=beta_init))

}
if(d[i+j+k]==2){
##Weibull
k2<-k2+1
dens_row_el<-c(.densityW(t,covariates_temp,cova_temp, sigma=x[i+j+k],nu=x[e+k2],beta=beta_init))
}
if(d[i+j+k]==3){
##Exponentiated Weibull
k2<-k2+1
k3<-k3+1
dens_row_el<-c(.densityEW(t,covariates_temp,cova_temp, sigma=x[i+j+k],nu=x[e+k2],theta=x[e+w+k3],beta=beta_init))
}


dens_row<-c(dens_row,dens_row_el)

f<-j
g<-j-1
}}
k<-k+f
m<-m+g

}else{k <- k-1
m <- m-1}}


matrix(dens_row,ncol=1)

}






###########################################################################################################################

###########################
#param.init
###########################

param.init<-function(data=NULL,cov=NULL,states,mtrans, cov_tra = NULL, cens = NULL,
			dist_init=NULL, proba_init=NULL, coef_init=NULL){

#length of cov_tra
if(missing(cov_tra)){length_cov<-0}
else{length_cov<-length(cov_tra)}

#checking conditions
if (missing(states))
	stop("Argument 'states' is missing with no default")
if (missing(mtrans))
	stop("Argument 'mtrans' is missing with no default")
 if (nrow(mtrans) != ncol(mtrans)) 
        stop("Argument 'mtrans' must be a quadratic  matrix.")
## is mtrans a quadratic matrix
if(!missing(mtrans) && nrow(mtrans)!=ncol(mtrans))stop("'mtrans' must be a quadratic matrix")

if(!missing(proba_init) && !missing(mtrans) && length(proba_init)!=length(which(mtrans!="FALSE")))
	stop("The length of 'proba_init' must be the same as the number of possible transitions in the matrix 'mtrans'")

##transitions on diagonal
k<-0
for(i in 1:dim(mtrans)[1]){
if(mtrans[i,i]!=FALSE)k<-k+1}
if (!missing(mtrans) && k> 0) 
       stop("Transitions into the same state are not allowed")
if (nrow(mtrans) != length(states)) 
        stop("The row number of 'mtrans' must be equal to the number of states.")
if (length(states) != length(unique(states))) 
        stop("The state names must be unique.")

#do we operate on real data or not
if (missing(data))
###################################
# no real data
##################################
{
if(missing(proba_init))stop("Argument 'proba_init' is missing with no default")

#number of the covariates r
if(length_cov>=1){r<-c()
for(i in 1:length_cov){
r<-c(r,length(cov_tra[[i]]))
cov_tra[[i]]<-sort(cov_tra[[i]])
}
cov<-as.matrix(cov)
}else{
if(!missing(coef_init)){
r<-c()
transitions<-length(which(mtrans!=FALSE))
if(length(coef_init)%/%transitions-floor(length(coef_init)%/%transitions)!=0){stop("Wrong format of the argument 'coef_init'")}
else{
r<-length(coef_init)%/%transitions

}}else{r<-0}}
s<-length(states)
#All observed transitions i!=j
trans.hj<-c()
for(i in 1:s){
for(j in 1:s){
if(mtrans[i,j]!=FALSE)trans.hj[length(trans.hj)+1]<-as.character(paste(states[i],states[j],sep=""))
}}
#All observed transitions i=j
trans.hh<-c()
for(i in 1:s){
if(missing(cens)){
trans.hh<-c(trans.hh,as.character(paste(states[i],states[i],sep="")))
}else{
trans.hh<-c(trans.hh,as.character(paste(states[i],cens,sep="")))
}}
a<-0
for(i in 1:length(trans.hh)){

for(j in 1:length(trans.hh)){
if(length(trans.hh)>0 && mtrans[as.numeric(substring(trans.hh[i],first=1,last=1)),j]==FALSE){
a<-a+1
}
if(a==length(trans.hh)){
trans.hh2<-trans.hh[-which(as.numeric(substring(trans.hh,first=1,last=1))==i)]}}}
trans.hh<-trans.hh2
##auxilary matrice of logical values
mtrans.log<-matrix(FALSE,ncol=s,nrow=s)
for(i in 1:dim(mtrans)[1]){
for(j in 1:dim(mtrans)[2]){
if(mtrans[i,j]!=FALSE)mtrans.log[i,j]<-TRUE
}}
#Number of transitions to estimate
p<-sum(1*mtrans.log)
#Initialisation of the parameters
#Initial Transition Matrix P
#t vector of number of essential transitions in a row
#npos vector fo all possible transitions
t<-rep(0,s)
pos<-list()
npos<-c()
dist<-c()
matrix.P<-matrix(0,ncol=s, nrow=s)
table.state<-matrix(ncol=s, nrow=s) 
m<-1
npos2<-matrix(0,nrow=s,ncol=s)
for(i in 1:s){
for (j in 1:s){
if(mtrans[i,j]!=FALSE){
dist<-c(dist,mtrans[i,j])
matrix.P[i, j]<-proba_init[m]
npos2[i,j]<-as.character(paste(i,j,sep=""))
m<-m+1}
else{matrix.P[i,j]<-0}
}}
##auxilary matrice of logical values
mtrans.log<-matrix(FALSE,ncol=s,nrow=s)
for(i in 1:dim(mtrans)[1]){
for(j in 1:dim(mtrans)[2]){
if(mtrans[i,j]!=FALSE)mtrans.log[i,j]<-TRUE
}}
#number of probabilites in the parameters matrix
nprob<-0
npos<-c()
proba_new<-c()
proba<-proba_init
if(length(states)>2){
for(i in 1:s){
if(1*sum(mtrans.log[i,])>0){
p_in_row<-1*sum(mtrans.log[i,])-1
proba_new<-c(proba_new,proba[1:p_in_row])

proba_new<-proba_new[1:(length(proba_new))]
proba<-proba[(p_in_row+1):length(proba)]
nprob<-nprob+1*sum(mtrans.log[i,])-1
npos<-c(npos,npos2[i,which(npos2[i,]>0)])
npos<-npos[1:(length(npos)-1)]
}
}
}
proba_init2<-proba_init
proba_init<-proba_new
## codes for distributions
for(i in 1:length(dist)){
if(dist[i]%in%c("E","Exp","Exponential"))dist[i]<-1
if(dist[i]%in%c("W","Weibull"))dist[i]<-2
if(dist[i]%in%c("EW","EWeibull","Exponentiated Weibull"))dist[i]<-3}
dist<-as.numeric(dist)
#Matrix of the Parametres
# e parametres sigma
# trans.e for which transitions
# w parametres sigma
# trans.w for which transitions
# ew parameters theta
# trans.ew for which transitions
e<-0
trans.e<-trans.hj
w<-0
trans.w<-c()
ew<-0
trans.ew<-c()
for(i in 1:length(dist)){
if(dist[i]==1)e<-e+1
if(dist[i]==2){e<-e+1
w<-w+1
trans.w<-c(trans.w,trans.hj[i]) }
if(dist[i]==3){e<-e+1
w<-w+1
trans.w<-c(trans.w,trans.hj[i])
ew<-ew+1
trans.ew<-c(trans.ew,trans.hj[i])}}
ndist<-e+w+ew

#Matrix of the Parametres
# p parametres nu
# p parametres sigma
# p parametres teta
# sum(t)  initial probabilities of Markov chain (tyle ile trzeba)
# r coeff of regression Beta
nbeta<-0
if(length(r)>1 || r>0){
if(length_cov>=1){
label_beta<-c()
nbeta<-0
trans_beta<-c()
val_beta<-c()
for(i in 1:length(r)){
label_beta<-c(label_beta,sort(rep(paste("Beta", i, sep=""), r[i])))
nbeta<-nbeta+r[i]
trans_beta<-c(trans_beta,cov_tra[[i]])
val_beta<-rep(0,length(trans_beta))
}
}
else{label_beta<-sort(rep(paste("Beta", 1:r, sep=""), p))
nbeta<-p*r
trans_beta<-rep(trans.hj,r)
val_beta<-rep(0,p*r)
}}
else{
label_beta<-c()
trans_beta<-c()
val_beta<-c()}
#The initial parameters
#checking if the vectors are OK

if(!missing(dist_init)){
	if(length(dist_init)!=ndist || dist_init[1:ndist]<=0 || !is.numeric(dist_init))	
		stop("Wrong format of the vector 'dist_init'")
}
if(length(proba_init)!=nprob || proba_init<0 || proba_init>=1 || !is.numeric(proba_init))
		stop("Wrong format of the vector 'proba_init'")
if(!missing(coef_init)){
	if(length(coef_init)!=nbeta || !is.numeric(coef_init))stop("Wrong format of the vector 'coef_init'")}
if(r==0){l<-0
} else if(length(r)>1 || r>0){l<-p}
parameters<-data.frame(Label=c(rep("sigma",e),rep("nu", w), rep("theta",ew),
 rep("P",nprob),label_beta ),
Transition=c(trans.e,trans.w,trans.ew,npos,trans_beta),
 Value=c(rep(1, e),rep(1,w),rep(1,ew),rep(0,nprob),val_beta))
for(i in 1:s){
for(j in 1:s){
parameters$Value[parameters$Label=="P" & parameters$Transition==as.character(paste(states[i],states[j], sep=""))]<-matrix.P[i,j]
}}
if(!missing(dist_init))parameters[1:ndist,3]<-dist_init
parameters[(ndist+1):(ndist+nprob),3]<-proba_init
if(!missing(coef_init))parameters[(ndist+nprob+1):(ndist+nprob+nbeta),3]<-coef_init
dimnames(matrix.P)<-list(states,states)
message<-paste("Initial values for the Multi-State Semi-Markov Model for ", s," states")
colnames(mtrans)<-states
rownames(mtrans)<-states
for(i in 1:length(states)){
for(j in 1:length(states)){
if(mtrans[i,j]==FALSE)mtrans[i,j]<-"-"
if(mtrans[i,j]=="E")mtrans[i,j]<-"Exponential"
if(mtrans[i,j]=="W")mtrans[i,j]<-"Weibull"
if(mtrans[i,j]=="EW")mtrans[i,j]<-"Exponentiated Weibull"
}}
dimnames(table.state)<-list(states,states)

if(missing(proba_init)){
proba.init<-parameters[1:e,]
proba_init2<-parameters[(ndist+1):(ndist+nprob),3]

licznik<-1
last_proba<-0
for(i in 1:(e-1)){
if(substring(parameters$Transition[i],first=1,last=1)==substring(parameters$Transition[i+1],first=1,last=1)){
proba.init[i,3]<-proba_init2[licznik]
last_proba<-last_proba+proba_init2[licznik]
licznik<-licznik+1}
else if(i<e-1){
proba.init[i,3]<-1-last_proba
last_proba<-0}
if(i==e-1){

proba.init[i+1,3]<-1-last_proba
}
}
proba.init$Label<-"P"
if(nbeta>0){
coef.init <- parameters[(ndist+nprob+1):(ndist+nprob+nbeta),]
rownames(coef.init) <- c(1:dim(coef.init)[1])
param.init.object<-list(nstates=length(states), matrix.P=matrix.P, last=0, Transition_matrix=mtrans,dist.init=parameters[1:(ndist),],proba.init=proba.init,coef.init=coef.init)
}else{
param.init.object<-list(nstates=length(states), matrix.P=matrix.P, last=0, Transition_matrix=mtrans,dist.init=parameters[1:(ndist),],proba.init=proba.init)

}}
else{
proba.init<-parameters[1:e,]
proba.init$Label<-"P"
for(i in 1:e){
proba.init$Value[i]<-proba_init2[i]
}
if(nbeta>0){
coef.init <- parameters[(ndist+nprob+1):(ndist+nprob+nbeta),]
rownames(coef.init) <- c(1:dim(coef.init)[1])
param.init.object<-list(nstates=length(states), matrix.P=matrix.P, last=0, Transition_matrix=mtrans,dist.init=parameters[1:(ndist),],proba.init=proba.init,coef.init=coef.init)
}else{
param.init.object<-list(nstates=length(states), matrix.P=matrix.P, last=0, Transition_matrix=mtrans,dist.init=parameters[1:(ndist),],proba.init=proba.init)

}}
}
################################
# with data
################
else{
if (missing(cov) && missing(cov_tra)==FALSE && missing(proba_init))
	stop("To indicate the transitions for covariates you need to define the covariates matrix")
if(!missing(cov)){
if(dim(data)[1]!=dim(cov)[1])
	stop("Argument 'data' and 'cov' need to have the same number of rows")
if(!is.data.frame(cov))
	stop("Argument 'cov' must be a data.frame")}
if(!is.data.frame(data))
	stop("Argument 'data' must be a data.frame")
if(!missing(cens) && !cens%in%data[,3])stop("Wrong format of the argument 'cens'")

if(!missing(proba_init) && !missing(mtrans) && length(proba_init)!=length(which(mtrans!=FALSE)))
	stop("The length of 'proba_init' must be the same as the number of possible transitions in the matrix 'mtrans'")
#number of states
s<-length(states)
if(!missing(cens)){
for(i in 1:dim(data)[1]){
if(data[i,3]==cens){
if(i<dim(data)[1] && data[i+1,1]==data[i,1]){
stop("The censure is not the last state of individuals")}
}
}
}
trans.h<-data[data[,2]!=data[,3],2]
trans.j<-data[data[,2]!=data[,3],3]
if(missing(cens)){
trans.cens<-data[data[,2]==data[,3],3]}
else{
trans.cens<-data[data[,3]==cens,3]}

table.state<-matrix(ncol=s, nrow=s) 
for(i in 1:s){
for (j in 1:s){
if(i!=j){
table.state[i,j]<-1*(sum(1*(trans.h==states[i] & trans.j==states[j])))}
else{table.state[i,j]<-1*(sum(1*(trans.cens==states[i])))}
}
}
if(!missing(cens)){
for(i in 1:dim(data)[1]){
if(data[i,3]==cens){
data[i,3]<-data[i,2]
}
}}
#number of censured patients
Ncens<-length(which(data[,2]==data[,3]))
time.ind<-NA
#number of the covariates r
if(missing(cov)){r<-0
cov=as.matrix(rep(0,dim(data)[1]))}
else{
if(length_cov>=1){r<-c()
for(i in 1:length_cov){
r<-c(r,length(cov_tra[[i]]))
cov_tra[[i]]<-sort(cov_tra[[i]])
}
cov<-as.matrix(cov)
}
else{cov<-as.matrix(cov)
r<-dim(cov)[2]
}

#we check if covariates are dependent from time
time.ind<-rep(TRUE,dim(cov)[2])
for(j in 1:dim(cov)[2]){
l<-0
for(i in 2:dim(cov)[1]){
if(data[i,1]==data[i-1,1]){
if(cov[i,j]!=cov[i-1,j]){
time.ind[j]<-FALSE
i<-dim(cov)[1]+1
}}}}}

#Construction of the variable that identyfies transitions to a different state
#i!=j -> 1
#i=j -> 0
data[,5]<-1*(data[,2]!=data[,3])
#Construction of the variable that specifies the transitions 
data[,6]<-as.character(paste(data[,2], data[,3], sep=""))
#All observed transitions i!=j
trans.hj<-c()
trans.hh<-c()
for(i in 1:s){
for(j in 1:s){
if(mtrans[i,j]!=FALSE)trans.hj[length(trans.hj)+1]<-as.character(paste(states[i],states[j],sep=""))
}}
#All observed transitions i=j
for(i in 1:s){
trans.hh<-sort(unique(data[data[,5]==0,6]))
}
trans.hh2<-trans.hh
for(i in 1:length(trans.hh)){
a<-0
for(j in 1:length(trans.hh)){
 
if(length(trans.hh)>0 && mtrans[as.numeric(substring(trans.hh[i],first=1,last=1)),j]==FALSE){
a<-a+1
}
if(a==length(trans.hh)){
trans.hh2<-trans.hh[-which(as.numeric(substring(trans.hh,first=1,last=1))==i)]}}}
trans.hh<-trans.hh2
##auxilary matrice of logical values
mtrans.log<-matrix(FALSE,ncol=s,nrow=s)
for(i in 1:dim(mtrans)[1]){
for(j in 1:dim(mtrans)[2]){
if(mtrans[i,j]!=FALSE)mtrans.log[i,j]<-TRUE
}}
#Number of transitions to estimate
p<-sum(1*mtrans.log)
#Initialisation of the parameters
#Initial Transition Matrix P
#t vector of number of essential transitions in a row
#npos vector fo all possible transitions
t<-rep(0,s)
pos<-list()
npos<-c()
dist<-c()
matrix.P<-matrix(ncol=s, nrow=s)
#table.state<-matrix(ncol=s, nrow=s) 
for(i in 1:s){
k<-0
for (j in 1:s){
if(mtrans[i,j]!=FALSE){
dist<-c(dist,mtrans[i,j])
if(sum(1*(trans.h==states[i]))>0){
matrix.P[i, j]<-1*(i!=j)*(sum(1*(trans.h==states[i] & trans.j==states[j]))/sum(1*(trans.h==states[i])))}
else{matrix.P[i,j]<-0}
k<-k+1
pos[[length(pos)+1]]<-c(states[i],states[j])
npos[length(npos)+1]<-as.character(paste(states[i],states[j],sep=""))
}
else{matrix.P[i,j]<-0}
#if(i!=j){
#table.state[i,j]<-1*(sum(1*(trans.h==states[i] & trans.j==states[j])))}
#else{table.state[i,j]<-1*(sum(1*(trans.cens==states[i])))}
}
t[i]<-k
}
t<-t-1
pos.temp<-list()
## codes for distributions
for(i in 1:length(dist)){
if(dist[i]%in%c("E","Exp","Exponential"))dist[i]<-1
if(dist[i]%in%c("W","Weibull"))dist[i]<-2
if(dist[i]%in%c("EW","EWeibull","Exponentiated Weibull"))dist[i]<-3}
dist<-as.numeric(dist)
#deleting the last transitions in the row from the calculations of P

if(s>2){
for(i in 1:(length(pos)-1)){
if(pos[[i]][[1]]==pos[[i+1]][[1]])pos.temp[length(pos.temp)+1]<-pos[i]
if(i>1){
if(pos[[i]][[1]]!=pos[[i+1]][[1]] && pos[[i]][[1]]!=pos[[i-1]][[1]])pos.temp[length(pos.temp)+1]<-pos[i]
}}
pos<-pos.temp

npos<-c(length(pos))
for(i in 1:length(pos)){

npos[i]<-as.character(paste(pos[[i]][[1]],pos[[i]][[2]],sep=""))
}
#number of probabilites in the parameters matrix
proba_new<-c()
proba<-proba_init

if(length(states)>2){
nprob<-length(npos)
for(i in 1:s){
if(1*sum(mtrans.log[i,])>0){
if(!missing(proba_init)){
p_in_row<-1*sum(mtrans.log[i,])
proba_new<-c(proba_new,proba[1:p_in_row])
proba_new<-proba_new[1:(length(proba_new)-1)]

proba<-proba[(p_in_row+1):length(proba)]}
}
}
}}
else{

  nprob<-0
npos<-c()}
if(!missing(proba_init)){
proba_init2<-proba_init
proba_init<-proba_new}
#Matrix of the Parametres
# e parametres sigma
# trans.e for which transitions
# w parametres sigma
# trans.w for which transitions
# ew parameters theta
# trans.ew for which transitions
e<-0
trans.e<-trans.hj
w<-0
trans.w<-c()
ew<-0
trans.ew<-c()
for(i in 1:length(dist)){
if(dist[i]==1)e<-e+1
if(dist[i]==2){e<-e+1
w<-w+1
trans.w<-c(trans.w,trans.hj[i]) }
if(dist[i]==3){e<-e+1
w<-w+1
trans.w<-c(trans.w,trans.hj[i])
ew<-ew+1
trans.ew<-c(trans.ew,trans.hj[i])}}
ndist<-e+w+ew


#Matrix of the Parametres
# p parametres nu
# p parametres sigma
# p parametres teta
# sum(t)  initial probabilities of Markov chain (tyle ile trzeba)
# r coeff of regression Beta
nbeta<-0
if(length(r)>1 || r>0){
if(length_cov>=1){
label_beta<-c()
nbeta<-0
trans_beta<-c()
val_beta<-c()
for(i in 1:length(r)){
label_beta<-c(label_beta,sort(rep(paste("Beta", i, sep=""), r[i])))
nbeta<-nbeta+r[i]
trans_beta<-c(trans_beta,cov_tra[[i]])
val_beta<-rep(0,length(trans_beta))
}

}
else{label_beta<-sort(rep(paste("Beta", 1:r, sep=""), p))
nbeta<-p*r
trans_beta<-rep(trans.hj,r)
val_beta<-rep(0,p*r)
}}
else{
label_beta<-c()
trans_beta<-c()
val_beta<-c()}

temp2<-data.frame( Transition=data[,6],covariates=cov)
#The initial parameters
#checking if the vectors are OK

if(!missing(dist_init)){
	if(length(dist_init)!=ndist || dist_init[1:ndist]<=0 || !is.numeric(dist_init))	
		stop("Wrong format of the vector 'dist_init'")
}

if(!missing(proba_init)){
	if(length(proba_init)!=nprob || proba_init<0 || proba_init>=1 || !is.numeric(proba_init))
		stop("Wrong format of the vector 'proba_init'")
}

if(!missing(coef_init)){
	if(length(coef_init)!=nbeta || !is.numeric(coef_init))stop("Wrong format of the vector 'coef_init'")}

if(r==0){l<-0
} else if(length(r)>1 || r>0){l<-p}
parameters<-data.frame(Label=c(rep("sigma",e),rep("nu", w), rep("theta",ew),
 rep("P",nprob),label_beta ),
Transition=c(trans.e,trans.w,trans.ew,npos,trans_beta),
 Value=c(rep(1, e),rep(1,w),rep(1,ew),rep(0,nprob),val_beta))
for(i in 1:s){
for(j in 1:s){
parameters$Value[parameters$Label=="P" & parameters$Transition==as.character(paste(states[i],states[j], sep=""))]<-matrix.P[i,j]
}}
if(!missing(dist_init))parameters[1:ndist,3]<-dist_init
if(!missing(proba_init))parameters[(ndist+1):(ndist+nprob),3]<-proba_init
if(!missing(coef_init))parameters[(ndist+nprob+1):(ndist+nprob+nbeta),3]<-coef_init

dimnames(matrix.P)<-list(states,states)
colnames(mtrans)<-states
rownames(mtrans)<-states
for(i in 1:length(states)){
for(j in 1:length(states)){
if(mtrans[i,j]==FALSE)mtrans[i,j]<-"-"
if(mtrans[i,j]=="E")mtrans[i,j]<-"Exponential"
if(mtrans[i,j]=="W")mtrans[i,j]<-"Weibull"
if(mtrans[i,j]=="EW")mtrans[i,j]<-"Exponentiated Weibull"
}}

message<-paste("Initial values for the Multi-State Semi-Markov Model for ", s," states")
dimnames(table.state)<-list(states,states)
if(missing(proba_init)){
proba.init<-parameters[1:e,]
proba_init2<-parameters[(ndist+1):(ndist+nprob),3]
licznik<-1
last_proba<-0
for(i in 1:(e-1)){
if(substring(parameters$Transition[i],first=1,last=1)==substring(parameters$Transition[i+1],first=1,last=1)){
proba.init[i,3]<-proba_init2[licznik]
last_proba<-last_proba+proba_init2[licznik]
licznik<-licznik+1}
else if(i<e-1){
proba.init[i,3]<-1-last_proba
last_proba<-0}
if(i==e-1){
proba.init[i+1,3]<-1-last_proba
}
}
proba.init$Label<-"P"
if(nbeta>0){
coef.init = parameters[(ndist+nprob+1):(ndist+nprob+nbeta),]
rownames(coef.init) <- c(1:dim(coef.init)[1])
param.init.object<-list(nstates=length(states), table.state=table.state, Ncens=Ncens, matrix.P=matrix.P,last=max(data[,4]),Transition_matrix=mtrans,dist.init=parameters[1:ndist,],proba.init=proba.init,coef.init=coef.init)
}else{
param.init.object<-list(nstates=length(states), table.state=table.state,  Ncens=Ncens, matrix.P=matrix.P,last=max(data[,4]),Transition_matrix=mtrans,dist.init=parameters[1:ndist,],proba.init=proba.init)
}}
else{
proba.init<-parameters[1:e,]
proba.init$Label<-"P"
for(i in 1:e){
proba.init$Value[i]<-proba_init2[i]
}
if(nbeta>0){
coef.init = parameters[(ndist+nprob+1):(ndist+nprob+nbeta),]
rownames(coef.init) <- c(1:dim(coef.init)[1])
param.init.object<-list(nstates=length(states), table.state=table.state, Ncens=Ncens, matrix.P=matrix.P,last=max(data[,4]),Transition_matrix=mtrans,dist.init=parameters[1:ndist,],proba.init=proba.init,coef.init=coef.init)
}else{

param.init.object<-list(nstates=length(states), table.state=table.state,  Ncens=Ncens, matrix.P=matrix.P,last=max(data[,4]),Transition_matrix=mtrans,dist.init=parameters[1:ndist,],proba.init=proba.init)
}}

}
class(param.init.object)<-"param.init"
param.init.object
}

############################################################################################################################

#############################
# hazard
#############################

hazard<-function(object,
	type="alpha",time=NULL,cov=NULL,s=0,t="last",Length=1000){

if(missing(object))stop("Argument 'object' is missing with no default")

 if (!inherits(object, "semiMarkov") && !inherits(object,"param.init")) 
        stop("Argument 'object' must be 'semiMarkov' or 'param.init' object")
if(type!="alpha" && type!="lambda")
	stop("Argument 'type' either 'alpha' or 'lambda'")
if(!missing(time) && !is.vector(time))stop("Argument 'time' must be a vector")
if(t=="last" ){
if(object$last==0)stop("The argument 't' must be given")
t<-object$last}
else if(!is.numeric(t)){
stop("'t' must be a number")
}
Nstates<-object$nstates
#by<-Length/(t-s)
k<-1
if(inherits(object,"semiMarkov")){
object1<-object[[2]]}
else{object1<-as.data.frame(rbind(object$dist.init,object$proba.init))}

v1<-c()
v<-list()
names<-c("v1")

#the names of arguments
call<-match.call()
if(as.character(call)[1]=="alpha"){
Names<-noquote(as.character(call))[-1]
}else{Names<-noquote(as.character(call))[2:(k+1)]}

#how many transitions
tr<-0
#possible transitions
vec<-c()
#possible states
states<-c()

for(i in 1:dim(object1)[1]){
if(object1[i,1]=="sigma"){tr<-tr+1
vec<-c(vec,as.character(object1[i,2]))
states<-c(states,substring(object1[i,2],first=1,last=1))}}
vec<-unique(vec)

#how many covariates
beta<-c()
transition_beta<-c()
#which distributions
dist<-c()
for(i in 1:dim(object$Transition_matrix)[1]){
for(j in 1:dim(object$Transition_matrix)[2]){
if(object$Transition_matrix[i,j]!="-")dist[length(dist)+1]<-object$Transition_matrix[i,j]
}}
for(i in 1:dim(object1)[1]){
if(substring(object1[i,1],first=1,last=1)=="B"){beta<-c(beta,object1[i,1])
transition_beta<-c(transition_beta,object1[i,2])}}
beta_unique<-unique(beta)
beta<-length(beta_unique)

########
#vector Time
########
if(!is.vector(time)){
#time<-c(s:(s+Length))/by
time<-seq(0, t, length = Length+1)
time<-time[-1]
}
else{s<-time[1]
t<-time[-1]
Length<-length(time)}

#############
#no covariates
#############
if(beta==0){
################
##time sojourn hazard rate
###############
if(type=="alpha"){
## ndist
ndist<-length(which(object$param.init[,1]%in%c("sigma","nu","theta")))
# number of parameters e, w ,ew
ew<-length(which(object$param.init[,1]=="theta"))
w<-length(which(object$param.init[,1]=="nu"))
e<-length(which(object$param.init[,1]=="sigma"))
trans_unique<-unique(object1[1:(ndist),2])
distribution<-c()

v_temp<-matrix(rep(0,tr*Length),nrow=Length)
colnames(v_temp)<-vec
k2<-0
k3<-0
for(i in 1:tr){
#Exponential distribution
if(dist[i]=="Exponential"){
distribution[i]<-"Exponential"
x<- .hazardE(time,FALSE,0, object1[i,3])
x<-as.vector(x)}
#Weibull distribution
else if(dist[i]=="Weibull"){
distribution[i]<-"Weibull"
k2<-k2+1
x<- .hazardW(time,FALSE,0, object1[i,3], object1[e+k2,3])
x<-as.vector(x)
#Exponentiated Weibull distribution
}else{
distribution[i]<-"Exponentiated Weibull"
k2<-k2+1
k3<-k3+1
x<- .hazardEW(time,FALSE,0, object1[i,3], object1[k2+e,3],object1[k3+e+w,3])
x<-as.vector(x)
}

v_temp[,i]<-x
}
table<-as.matrix(summary(v_temp[,1]))


if(length(vec)>1){
for(i in 2:length(vec)){

table<-cbind(table,summary(v_temp[,i]))
}}
colnames(table)<-vec
v<-list(vector=v_temp,Time=time,Covariates=NA,Summary=table,Transition_matrix=object$Transition_matrix,call=Names,Type=type,cova=cov)
###############################
##SM hazard rate
################################
}else{
r<-1
## ndist
ndist<-length(which(object$param.init[,1]%in%c("sigma","nu","theta")))
# number of parameters e, w ,ew
ew<-length(which(object$param.init[,1]=="theta"))
w<-length(which(object$param.init[,1]=="nu"))
e<-length(which(object$param.init[,1]=="sigma"))

v_temp<-matrix(rep(0,tr*Length),nrow=Length)
colnames(v_temp)<-vec	
s<-length(unique(states))
y0<-0
k2<-0
k3<-0
k<-0
m<-0
for(i in 1:s){
tr_in_row<-0
for(z in 1:tr){
if(substring(object1[z,2],first=1,last=1)==as.character(i)){
tr_in_row<-tr_in_row+1}}
marg_row_el<-rep(0,length(time))
pr<-0
for(j in 0:(tr_in_row-1)){
if(j!=(tr_in_row-1)){
if(dist[i+j+k]=="Exponential"){
#########################################
##	Exponential distribution
########################################
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalE(time,FALSE,0, sigma=object1[i+j+k,3]))
pr<-pr+object1[ndist+i+j+m,3]}
else if(dist[i+j+k]=="Weibull"){
#########################################
##	Weibull distribution
########################################
k2<-k2+1
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalW(time,FALSE,0, sigma=object1[i+j+k,3],nu=object1[e+k2,3]))
pr<-pr+object1[ndist+i+j+m,3]
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
#########################################
##	Exponentiated Weibull distribution
########################################
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalEW(time,FALSE,0, sigma=object1[i+j+k,3],nu=object1[e+k2,3],theta=object1[e+w+k3,3]))
pr<-pr+object1[ndist+i+j+m,3]
}
}else{
if(dist[i+j+k]=="Exponential"){
#########################################
##	Exponential distribution
########################################
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalE(time,FALSE,0,  sigma=object1[i+j+k,3]))}
else if(dist[i+j+k]=="Weibull"){
#########################################
##	Weibull distribution
########################################
k2<-k2+1
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalW(time,FALSE,0,  sigma=object1[i+j+k,3],nu=object1[e+k2,3]))
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
#########################################
##	Exponentiated Weibull distribution
########################################
k2<-k2+1
k3<-k3+1
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalEW(time,FALSE,0,  sigma=object1[i+j+k,3],nu=object1[e+k2,3],theta=object1[e+w+k3,3]))
}
f<-j
g<-j-1
}
}
k<-k+f
m<-m+g
y0<-y0+marg_row_el
}
k<-0
m<-0
k2<-0
k3<-0
for(i in 1:s){
tr_in_row<-0
for(z in 1:tr){
if(substring(object1[z,2],first=1,last=1)==as.character(i)){
tr_in_row<-tr_in_row+1}}
pr0<-0
for(j in 0:(tr_in_row-1)){
if(j!=(tr_in_row-1)){
if(dist[i+j+k]=="Exponential"){
#########################################
##	Exponential distribution
########################################
x0<- .densityE(time, FALSE,0,  object1[i+j+k,3])*object1[ndist+i+j+m,3]
pr0<-pr0+object1[ndist+i+j+m,3]}
else if(dist[i+j+k]=="Weibull"){
#########################################
##	Weibull distribution
#########################################
k2<-k2+1
x0<- .densityW(time, FALSE,0,  object1[i+j+k,3], object1[e+k2,3])*object1[ndist+i+j+m,3]
pr0<-pr0+object1[ndist+i+j+m,3]
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
#########################################
##	Exponentiated Weibull distribution
########################################
k2<-k2+1
k3<-k3+1
x0<- .densityEW(time, FALSE,0,  object1[i+j+k,3], object1[e+k2,3],object1[e+w+k3,3])*object1[ndist+i+j+m,3]
pr0<-pr0+object1[ndist+i+j+m,3]
}
z0<-x0/y0
z0<-as.vector(z0)
v_temp[,r]<-z0
r<-r+1
}else{
if(dist[i+j+k]=="Exponential"){
#########################################
##	Exponential distribution
########################################
x0<- .densityE(time, FALSE,0,   object1[i+j+k,3])*(1-pr0)}
else if(dist[i+j+k]=="Weibull"){
#########################################
##	Weibull distribution
########################################
k2<-k2+1
x0<- .densityW(time, FALSE,0,   object1[i+j+k,3], object1[e+k2,3])*(1-pr0)
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
#########################################
##	Exponentiated Weibull distribution
########################################
k2<-k2+1
k3<-k3+1
x0<- .densityEW(time, FALSE,0,   object1[i+j+k,3], object1[e+k2,3],object1[e+w+k3,3])*(1-pr0)
}
				
z0<-x0/y0
z0<-as.vector(z0)
f<-j
g<-j-1
v_temp[,r]<-z0
r<-r+1
}
}
k<-k+f
m<-m+g
}
table<-as.matrix(summary(v_temp[,1]))
if(length(vec)>1){
for(i in 2:length(vec)){
table<-cbind(table,summary(v_temp[,i]))
}}
colnames(table)<-vec
v<-list(vector=v_temp,Time=time,Covariates=NA,Summary=table,Transition_matrix=object$Transition_matrix,call=Names,Type=type,cova=cov)
}
##########################
## with covariates
###########################
}else{
##time sojourn hazard rate
if(type=="alpha"){
#the value for the vector of the covariates
Cov<-cov
if(!missing(cov)){
if(length(cov)==beta && is.numeric(cov)){
cova<-rep(cov[1],Length)
COVA<-as.data.frame(cova)
if(beta>1){
for(i in 2:(beta)){
COVA<-cbind(COVA,rep(cov[i],Length))
}}}
else{
#we check the length of vectors covariates an Time
if(beta>1){
if(!is.matrix(cov) || !is.vector(cov[[1]]))stop("Wrong format of 'cov'")
l<-dim(cov)[1]
COVA<-as.data.frame(cov)
}else{
if(!is.vector(cov) || length(cov)!=length(time))stop("Wrong format of argument 'cov'")
COVA<-as.data.frame(cov)
}}
}else{
cova<-rep(0,Length)
COVA<-as.data.frame(cova)
if(beta>1){
for(i in 2:(beta)){
COVA<-cbind(COVA,rep(0,Length))}
}}
## ndist
ndist<-length(which(object$param.init[,1]%in%c("sigma","nu","theta")))
# number of parameters e, w ,ew
ew<-length(which(object$param.init[,1]=="theta"))
w<-length(which(object$param.init[,1]=="nu"))
e<-length(which(object$param.init[,1]=="sigma"))
#how many parameters P
p<-0
for(i in 1:dim(object1)[1]){
if(object1[i,1]=="P"){p<-p+1}}
 #how many transitions for beta
transitions<-object1[1:(ndist),2]
trans_unique<-unique(transitions)
v_temp<-matrix(rep(0,length(trans_unique)*Length),nrow=Length)
colnames(v_temp)<-trans_unique
cov<-1
k2<-0
k3<-0
for(t in 1:length(trans_unique)){
#which covariates
tr.h<- as.numeric(substring(trans_unique[t],first=1,last=1))	
tr.j<- as.numeric(substring(trans_unique[t],first=2,last=2))
#which and how many covariates for this transition
which2<-which(object1[t,2]==object1[(ndist+p+1):dim(object1)[1],2])+ndist+p
which1<-c()
for(i in 1:length(which2)){
which1<-c(which1,as.numeric(substring(object1[which2[i],1],first=5,last=nchar(as.character(object1[which2[i],1])))))}
how_many<-length(which(object1[(ndist+p+1):dim(object1)[1],2]==object1[t,2]))
#which sigma
sigma<-object1[t,3]
#which P
P<-0
for(i in (ndist+1):(ndist+p)){
if(as.numeric(substring(object1[i,2],first=1,last=1))==tr.h && as.numeric(substring(object1[i,2],first=2,last=2))==tr.j){
P<-object1[i,3]}}
if(P==0){
P1<-0
for(i in (ndist+1):(ndist+p)){
if(as.numeric(substring(object1[i,2],first=1,last=1))==tr.h){P<-P+object1[i,3]}
}
P<-1-P1
}
if(dist[t]=="Exponential"){
#########################################
##	Exponential distribution
########################################
#covariate
if(!is.na(which1[1])){
x<- as.vector(.hazardE(time,how_many, COVA[,which1], sigma, object1[which2,3]))
cov<-cov+1}
else{
x<- as.vector(.hazardE(time,FALSE,0, sigma))
}
}
else if(dist[t]=="Weibull"){
#########################################
##	Weibull distribution
########################################
k2<-k2+1
nu<-object1[k2+e,3]
#covariate
if(!is.na(which1[1])){
x<- as.vector(.hazardW(time,how_many, COVA[,which1], sigma, nu,object1[which2,3]))
cov<-cov+1}
else{
x<- as.vector(.hazardW(time,FALSE, 0, sigma, nu))
}
}
else if(dist[t]=="Exponentiated Weibull"){
#########################################
##	Exponentiated Weibull distribution
########################################
k2<-k2+1
nu<-object1[k2+e,3]
#which theta
k3<-k3+1
theta<-object1[e+w+k3,3]
#covariate
if(!is.na(which1[1])){
x<- as.vector(.hazardEW(time,how_many, COVA[,which1], sigma, nu,theta,object1[which2,3]))
cov<-cov+1
}else{
x<- as.vector(.hazardEW(time,FALSE, 0, sigma, nu,theta))
}
}
v_temp[,t]<-x}
table<-as.matrix(summary(v_temp[,1]))
if(length(trans_unique)>1){
for(i in 2:length(trans_unique)){
table<-cbind(table,summary(v_temp[,i]))
}}
colnames(table)<-trans_unique
v<-list(vector=v_temp,Time=time,Covariates=COVA,Summary=table,Transition_matrix=object$Transition_matrix,call=Names,Type=type,cova=Cov)
######################
##SM hazard rate
######################
}else{
Cov<-cov
#the value for the vector of the covariates
if(!missing(cov)){
if(length(cov)==beta && is.numeric(cov)){
cova<-rep(cov[1],Length)
COVA<-as.data.frame(cova)
if(beta>1){
for(i in 2:(beta)){
COVA<-cbind(COVA,rep(cov[i],Length))
}}}
else{
if(beta>1){
if(!is.matrix(cov) || !is.vector(cov[[1]]))stop("Wrong format of 'cov'")
l<-dim(cov)[1]
if(length(time)!=l)stop("Formats of 'time' and 'cov' are not equal")
COVA<-as.data.frame(cov)
}else{
if(!is.matrix(cov) || length(cov)!=length(time))stop("Wrong format of argument 'cov'")
COVA<-as.data.frame(cov)
}}
}else{
cova<-rep(0,Length)
COVA<-as.data.frame(cova)
if(beta>1){
for(i in 2:(beta)){
COVA<-cbind(COVA,rep(0,Length))}
}}
#how many transitions
tr<-0
#possible transitions
vec<-c()
#possible states
states<-c()
for(i in 1:dim(object1)[1]){
if(object1[i,1]=="sigma"){tr<-tr+1
vec<-c(vec,as.character(object1[i,2]))
states<-c(states,substring(object1[i,2],first=1,last=1))}}
## ndist
ndist<-length(which(object$param.init[,1]%in%c("sigma","nu","theta")))
# number of parameters e, w ,ew
ew<-length(which(object$param.init[,1]=="theta"))
w<-length(which(object$param.init[,1]=="nu"))
e<-length(which(object$param.init[,1]=="sigma"))
#how many parameters P
p<-0
for(i in 1:dim(object1)[1]){
if(object1[i,1]=="P"){p<-p+1}}
#preparation of final table
trans_unique<-unique(object1[1:(ndist),2])
v_temp<-matrix(rep(0,length(trans_unique)*Length),nrow=Length)
colnames(v_temp)<-trans_unique
cov<-1
s<-length(unique(states))
y0<-0
 k2<-0
k3<-0
k<-0
m<-0
r<-1
for(i in 1:s){
tr_in_row<-0
for(z in 1:tr){
if(substring(object1[z,2],first=1,last=1)==as.character(i)){
tr_in_row<-tr_in_row+1}}
pr<-0
marg_row_el<-rep(0,length(time))
for(j in 0:(tr_in_row-1)){
#which and how many covariates for this transition
which2<-which(object1[i+j+k,2]==object1[(ndist+p+1):dim(object1)[1],2])+ndist+p
which1<-c()
for(t in 1:length(which2)){
which1<-c(which1,as.numeric(substring(object1[which2[t],1],first=5,last=nchar(as.character(object1[which2[t],1])))))}
how_many<-length(which(object1[(ndist+p+1):dim(object1)[1],2]==object1[i+j+k,2]))
#which transition
tr.h<- as.numeric(substring(trans_unique[i+j+k],first=1,last=1))	
tr.j<- as.numeric(substring(trans_unique[i+j+k],first=2,last=2))
if(j!=(tr_in_row-1)){
if(substring(object1[i+j+k,2],first=1,last=1)==tr.h && substring(object1[i+j+k,2],first=2,last=2)==tr.j)
{
if(dist[i+j+k]=="Exponential"){

##	Exponential distribution

#covariate
if(!is.na(which1[1])){
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalE(time,how_many,COVA[,which1], sigma=object1[i+j+k,3],object1[which2,3]))
}else{
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalE(time,FALSE,0, sigma=object1[i+j+k,3]))
}
pr<-pr+object1[ndist+i+j+m,3]
}
else if(dist[i+j+k]=="Weibull"){
k2<-k2+1

##	Weibull distribution

if(!is.na(which1[1])){
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalW(time,how_many,COVA[,which1], sigma=object1[i+j+k,3],nu=object1[e+k2,3],object1[which2,3]))
}else{
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalW(time,FALSE,0, sigma=object1[i+j+k,3],nu=object1[e+k2,3]))
}
pr<-pr+object1[ndist+i+j+m,3]
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
k2<-k2+1
k3<-k3+1

##	Exponentiated Weibull distribution

if(!is.na(which1[1])){
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalEW(time,how_many,COVA[,which1], sigma=object1[i+j+k,3],nu=object1[e+k2,3],theta=object1[e+w+k3,3],object1[which2,3]))
}else{
marg_row_el<-c(marg_row_el)+c(object1[ndist+i+j+m,3]*.survivalEW(time,FALSE,0, sigma=object1[i+j+k,3],nu=object1[e+k2,3],theta=object1[e+w+k3,3]))
}
pr<-pr+object1[ndist+i+j+m,3]
}		
}
}else{
if(dist[i+j+k]=="Exponential"){

##	Exponential distribution

#covariate?
if(!is.na(which1[1])){
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalE(time,how_many,COVA[,which1], sigma=object1[i+j+k,3],object1[which2,3]))
}else{
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalE(time,FALSE,0, sigma=object1[i+j+k,3]))
}
}
else if(dist[i+j+k]=="Weibull"){
k2<-k2+1

##	Weibull distribution

if(!is.na(which1[1])){
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalW(time,how_many,COVA[,which1], sigma=object1[i+j+k,3],nu=object1[e+k2,3],object1[which2,3]))
}else{
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalW(time,FALSE,0, sigma=object1[i+j+k,3],nu=object1[e+k2,3]))
}
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
k2<-k2+1
k3<-k3+1

##	Exponentiated Weibull distribution

if(!is.na(which1[1])){
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalEW(time,how_many,COVA[,which1], sigma=object1[i+j+k,3],nu=object1[e+k2,3],theta=object1[e+w+k3,3],object1[which2,3]))
}else{
marg_row_el<-c(marg_row_el)+c((1-pr)*.survivalEW(time,FALSE,0, sigma=object1[i+j+k,3],nu=object1[e+k2,3],theta=object1[e+w+k3,3]))
}
}	
}
}
k<-k+j
m<-m+j-1
y0<-y0+marg_row_el
}
k<-0
m<-0
k2<-0
k3<-0
for(i in 1:s){
tr_in_row<-0
for(z in 1:tr){
if(substring(object1[z,2],first=1,last=1)==as.character(i)){
tr_in_row<-tr_in_row+1}}
pr<-0
for(j in 0:(tr_in_row-1)){

#which and how many covariates for this transition
which2<-which(object1[i+j+k,2]==object1[(ndist+p+1):dim(object1)[1],2])+ndist+p
which1<-c()
for(t in 1:length(which2)){
which1<-c(which1,as.numeric(substring(object1[which2[t],1],first=5,last=nchar(as.character(object1[which2[t],1])))))}
how_many<-length(which(object1[(ndist+p+1):dim(object1)[1],2]==object1[i+j+k,2]))
#which transition
tr.h<- as.numeric(substring(trans_unique[i+j+k],first=1,last=1))	
tr.j<- as.numeric(substring(trans_unique[i+j+k],first=2,last=2))
if(j!=(tr_in_row-1)){
if(substring(object1[i+j+k,2],first=1,last=1)==tr.h && substring(object1[i+j+k,2],first=2,last=2)==tr.j)
				{
if(dist[i+j+k]=="Exponential"){

##	Exponential distribution

#covariate?
if(!is.na(which1[1])){
x0<- .densityE(time,how_many,COVA[,which1],  object1[i+j+k,3],object1[which2,3])*object1[ndist+i+j+m,3]
}else{
x0<- .densityE(time,FALSE,0,  object1[i+j+k,3])*object1[ndist+i+j+m,3]
}
pr<-pr+object1[ndist+i+j+m,3]
}
else if(dist[i+j+k]=="Weibull"){
k2<-k2+1

##	Weibull distribution

if(!is.na(which1[1])){
x0<- .densityW(time,how_many,COVA[,which1],  object1[i+j+k,3], object1[e+k2,3],object1[which2,3])*object1[ndist+i+j+m,3]
}else{
x0<- .densityW(time,FALSE,0,  object1[i+j+k,3], object1[e+k2,3])*object1[ndist+i+j+m,3]
}
pr<-pr+object1[ndist+i+j+m,3]
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
k2<-k2+1
k3<-k3+1

##	Exponentiated Weibull distribution

if(!is.na(which1[1])){
x0<- .densityEW(time,how_many,COVA[,which1],  object1[i+j+k,3], object1[e+k2,3],object1[e+w+k3,3],object1[which2,3])*object1[ndist+i+j+m,3]
}else{
x0<- .densityEW(time,FALSE,0,  object1[i+j+k,3], object1[e+k2,3],object1[e+w+k3,3])*object1[ndist+i+j+m,3]
}
pr<-pr+object1[ndist+i+j+m,3]
}		
}
}else{
if(dist[i+j+k]=="Exponential"){

##	Exponential distribution

#covariate?
if(!is.na(which1[1])){
x0<- .densityE(time,how_many,COVA[,which1],  object1[i+j+k,3],object1[which2,3])*(1-pr)
}else{
x0<- .densityE(time,FALSE,0,  object1[i+j+k,3])*(1-pr)
}
}
else if(dist[i+j+k]=="Weibull"){
k2<-k2+1

##	Weibull distribution

if(!is.na(which1[1])){
x0<- .densityW(time,how_many,COVA[,which1],  object1[i+j+k,3], object1[e+k2,3],object1[which2,3])*(1-pr)
}else{
x0<- .densityW(time,FALSE,0,  object1[i+j+k,3], object1[e+k2,3])*(1-pr)
}
}
else if(dist[i+j+k]=="Exponentiated Weibull"){
k2<-k2+1
k3<-k3+1

##	Exponentiated Weibull distribution
########################################
if(!is.na(which1[1])){
x0<- .densityEW(time,how_many,COVA[,which1],  object1[i+j+k,3], object1[e+k2,3],object1[e+w+i+j+k,3],object1[which2,3])*(1-pr)
}else{
x0<- .densityEW(time,FALSE,0,  object1[i+j+k,3], object1[e+k2,3],object1[e+w+k3,3])*(1-pr)
}
}	
}
z0<-x0/y0
z0<-as.vector(z0)

v_temp[,r]<-z0
r<-r+1
}
k<-k+j
m<-m+j-1
}
table<-as.matrix(summary(v_temp[,1]))
if(length(trans_unique)>1){
for(i in 2:length(trans_unique)){
table<-cbind(table,summary(v_temp[,i]))}}
}
colnames(table)<-trans_unique
v<-list(vector=v_temp,Time=time,Covariates=COVA,Summary=table,Transition_matrix=object$Transition_matrix,call=Names,Type=type,cova=Cov)
}
class(v)<-"hazard"
v
}



summary.hazard<-function(object,...){
if(object$Type=="alpha"){

cat(object$call," : Hazard rates of waiting times\n")
}else{
cat(object$call," : Hazard rates of the semi-Markov process\n")
}
cat("\nTransition_matrix\n")
print(object$Transition_matrix)
cat("\nSummary statistics\n")
print(object$Summary)
cat("\nVector time\n")
cat("s = ",min(object$Time))
cat("\nt = ",max(object$Time))
cat("\nlength = ",length(object$Time))
cat("\n\nCovariates : ")
if(is.data.frame(object$Covariates) || !is.na(object$Covariates)){cat("Yes\n")
}
else{cat("No\n")}


}


##############################
#plot.hazard
##############################

plot.hazard<-function(x,x2=NULL,x3=NULL,x4=NULL,x5=NULL,x6=NULL,x7=NULL,
x8=NULL,x9=NULL,x10=NULL,transitions=NULL,names=NULL,legend=TRUE,legend.pos=NULL,cex=NULL,colors=NULL,xlab="Time",ylab="Hazard function",
lwd=3,type="p",...){

if(missing(x))stop("Argument 'x' is missing with no default")

 if (!inherits(x, "hazard")) 
        stop("expected object to be a result of the 'hazard' function")
if(legend==FALSE && length(legend.pos>0))
	warning("No legend displayed")

k<-length(x$call)
x1<-x$vector
if(k>1){
for(i in 2:k){
x1<-cbind(x1,x$vector)}}
#boundaries for Time
#they need to be equal for all the arguments
t<-max(x$Time)
s<-min(x$Time)
length<-length(x$Time)
if(missing(names)){
N<-0
Names<-x$call
cov<-x$cova
}
else{
N<-1
Names<-names}
#type of hazard (alpha/lambda)
Type<-x$Type
if(length(x2[[1]])>1){
 if (!inherits(x2, "hazard")) 
        stop("expected object [[2]] to be a result of the 'hazard' function")
t1<-max(x2$Time)
s1<-min(x2$Time)
length1<-length(x2$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x2$call)
cov<-c(cov,x2$cova)}
Type1<-x2$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x2$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x2$call)
k<-k+k_temp
x1<-cbind(x1,x2$vector)
if(k_temp>1){
for(i in 2:k_temp){x1<-cbind(x1,x2$vector)}}
}
if(length(x3)>1){
if (!inherits(x3, "hazard")) 
        stop("expected object [[3]]  to be a result of the 'hazard' function")
t1<-max(x3$Time)
s1<-min(x3$Time)
length1<-length(x3$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x3$call)
cov<-c(cov,x3$cova)}
Type1<-x3$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x3$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x3$call)
k<-k+k_temp
x1<-cbind(x1,x3$vector)
if(k_temp>1){ 
for(i in 2:k_temp){x1<-cbind(x1,x3$vector)}}
}
if(length(x4)>1){
 if (!inherits(x4, "hazard")) 
        stop("expected object [[4]] to be a result of the 'hazard' function")
t1<-max(x4$Time)
s1<-min(x4$Time)
length1<-length(x4$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x4$call)
cov<-c(cov,x4$cova)}
Type1<-x4$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x4$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x4$call)
k<-k+k_temp
x1<-cbind(x1,x4$vector)
if(k_temp>1){
for(i in 2:k_temp){x1<-cbind(x1,x4$vector)}}
}
if(length(x5)>1){
 if (!inherits(x5, "hazard")) 
        stop("expected object [[5]] to be a result of the 'hazard' function")
t1<-max(x5$Time)
s1<-min(x5$Time)
length1<-length(x5$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x5$call)
cov<-c(cov,x5$cova)}
Type1<-x5$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x5$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x5$call)
k<-k+k_temp
x1<-cbind(x1,x5$vector)
if(k_temp>1){
for(i in 2:k_temp){x1<-cbind(x1,x5$vector)}}
}
if(length(x6)>1){
if (!inherits(x6, "hazard")) 
        stop("expected object [[6]] to be a result of the 'hazard' function")
t1<-max(x6$Time)
s1<-min(x6$Time)
length1<-length(x6$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x6$call)
cov<-c(cov,x6$cova)}
Type1<-x6$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x6$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x6$call)
k<-k+k_temp
x1<-cbind(x1,x6$vector)
if(k_temp>1){
 for(i in 2:k_temp){x1<-cbind(x1,x6$vector)}}
}
if(length(x7)>1){
 if (!inherits(x7, "hazard")) 
        stop("expected object [[7]] to be a result of the 'hazard' function")
t1<-max(x7$Time)
s1<-min(x7$Time)
length1<-length(x7$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x7$call)
cov<-c(cov,x7$cova)}
Type1<-x7$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x7$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x7$call)
k<-k+k_temp
x1<-cbind(x1,x7$vector)
if(k_temp>1){
for(i in 2:k_temp){x1<-cbind(x1,x7$vector)}}
}
if(length(x8)>1){
 if (!inherits(x8, "hazard")) 
        stop("expected object [[8]] to be a result of the 'hazard' function")
t1<-max(x8$Time)
s1<-min(x8$Time)
length1<-length(x8$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x8$call)
cov<-c(cov,x8$cova)}
Type1<-x8$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x8$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x8$call)
k<-k+k_temp
x1<-cbind(x1,x8$vector)
if(k_temp>1){
for(i in 2:k_temp){x1<-cbind(x1,x8$vector)}}
}
if(length(x9)>1){
 if (!inherits(x9, "hazard")) 
        stop("expected object [[9]] to be a result of the 'hazard' function")
t1<-max(x9$Time)
s1<-min(x9$Time)
length1<-length(x9$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x9$call)
cov<-c(cov,x9$cova)}
Type1<-x9$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x9$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x9$call)
k<-k+k_temp
x1<-cbind(x1,x9$vector)
if(k_temp>1){
for(i in 2:k_temp){x1<-cbind(x1,x9$vector)}}
}
if(length(x10)>1){
 if (!inherits(x10, "hazard")) 
        stop("expected object [[10]] to be a result of the 'hazard' function")
t1<-max(x10$Time)
s1<-min(x10$Time)
length1<-length(x10$Time)
if(t!=t1 || s!=s1 || length!=length1)
	stop("The compatibility of Time is required")
if(missing(names)){
Names<-c(Names,x10$call)
cov<-c(cov,x10$cova)}
Type1<-x10$Type
if(Type!=Type1)
	stop("All arguments need to represent the same type of risk")
if(dim(x10$vector)[2]!=dim(x$vector)[2])
	stop("All arguments need to have the same number of transitions considered")
k_temp<-length(x10$call)
k<-k+k_temp
x1<-cbind(x1,x10$vector)
if(k_temp>1){
for(i in 2:k_temp){x1<-cbind(x1,x10$vector)}}
}
if(length(names)>0 && length(names)!=k)stop("The length of vector 'names' does not correspond to the number of arguments")
if(length(colors)>0 && length(colors)!=k)stop("The number of colors must be the same as number of arguments")
name<-colnames(x1)[1]
len<-length(which(colnames(x1)==name))
if(missing(transitions)){
tr<-dim(x$vector)[2]
names<-colnames(x$vector)
}else{
if(!all(transitions%in%colnames(x$vector)))
	stop("Wrong format of argument 'transitions'")
tr<-length(transitions)
which<-which(colnames(x$vector)%in%transitions)
temp<-which

if(len>1){
for(i in 1:(len-1)){
which<-c(which,temp+i*dim(x$vector)[2])}
}
x1<-data.frame(x1[,unique(which)])
names<-transitions
}
if(tr>1 && ceiling(tr/2)<5){
par(mfrow=c(ceiling(tr/2),2), oma = c(0, 0, 4, 0))}
else{
#if(tr==1){par( oma = c(0, 0, 4, 0))}
if(ceiling(tr/2)>=5){
par(mfrow=c(1,1), oma = c(0, 0, 4, 0),ask=TRUE)
}}
col<-c(1,2)
l<-1
for(i in 1:tr){
col<-c(1)
plot_colors<-c(col)
Time<-as.data.frame(x$Time)
y<-x1[,i]
if(k>1){
for(w in 1:(k-1)){
Time[,dim(Time)[2]+1]<-x$Time
y<-cbind(y,x1[,i+w*tr])
plot_colors<-c(plot_colors,w+1)
}}
max<-max(y)
if(missing(colors)){
matplot(Time, y,xlab=xlab,ylab=ylab,lwd=lwd,xlim=c(s,t),type=type,pch=".",col=plot_colors)}
else{
matplot(Time, y,xlab=xlab,ylab=ylab,lwd=lwd,xlim=c(s,t),type=type,pch=".",col=colors)}
if(tr>1){	title(paste("Transition",names[i],collapse=" "))}
else{
if(Type=="alpha"){
title(paste("Sojourn time hazard rate \n Transition",names[i],collapse=" "),font = 2, line = 1, cex = 1.2,  outer = F) 
}else{
title(paste("Semi-Markov process hazard rate \n Transition",names[i],collapse=" "),font = 2, line = 1, cex = 1.2,  outer = F) 
 }
}
if(legend==TRUE){
#if(N==0){
Legend<-c()
for(i in 1:length(Names)){
if(is.vector(cov)){
Legend[i]<-paste(Names[i],", cov=",cov[i])}
else{
Legend[i]<-paste(Names[i])
}
}
if(missing(cex))cex=0.5
#}
	if(missing(legend.pos)){
		if(missing(colors)){
		legend("topleft", "(x,y)",legend = Legend, 
   		col=plot_colors, lwd=2, cex=.8, horiz = FALSE)}
		else{legend("topleft", "(x,y)",legend = Legend, 
   		col=colors, lwd=2, cex=cex, horiz = FALSE)}
	}else if(is.vector(legend.pos) && length(legend.pos)==2*tr){
	pos<-legend.pos[l]
	max<-legend.pos[l+1]
		if(missing(colors)){
		legend(pos,max,legend = Legend, 
   		col=plot_colors, lwd=2, cex=cex, horiz = FALSE)
		}else{
		legend(pos,max,legend = Legend, 
   		col=colors, lwd=2, cex=cex, horiz = FALSE)}
	}else{
		stop("Wrong format of the argument legend.pos")}
		
}
l<-l+2
if(tr>1){
if(Type=="alpha"){
mtext("Sojourn time hazard rate",font = 2, line = 1, cex = 1.2,  outer = TRUE) 
}else{
mtext("Semi-Markov process hazard rate",font = 2, line = 1, cex = 1.2,  outer = TRUE) }}
}
}


#################################
# print.hazard
#################################
print.hazard<-function(x,whole = FALSE, ...){
if(x$Type=="alpha"){
cat(x$call," : Hazard rates of waiting times\n")
}else{
cat(x$call," : Hazard rates of the semi-Markov process\n")
}
if(whole==TRUE){
cat("\nTransition_matrix\n")
print(x$Transition_matrix)
cat("\nHazard rate values \n")
print(as.data.frame(x$vector))
cat("\n")
Time<-as.data.frame(x$Time)
colnames(Time)<-"Time"
print(Time)
cat("\n")
if(!is.na(x$Covariates) || length(x$Covariates)>1){
print(as.data.frame(x$Covariates))}
}
else if(whole == FALSE){
cat("\nTransition_matrix\n")
print(x$Transition_matrix)
cat("\nHazard rates values \n")
print(head(as.data.frame(x$vector)))
cat("\n")
Time<-head(as.data.frame(x$Time))
colnames(Time)<-"Time"
print(Time)
cat("\n")

if(!is.na(x$Covariates) || length(x$Covariates)>1){
cova<-head(as.data.frame(x$Covariates))
if(dim(cova)[2]>1){
colnames<-c()
for(i in 1:dim(cova)[2]){
colnames[i]<-paste("cov",i)
}
colnames(cova)<-colnames}
print(cova)}
cat("\nSummary statistics\n")
print(x$Summary)

}
else{stop("Argument 'whole' must be logical")}

}

#############################################

##############
# table.state
##############

table.state<-function(data,states=NULL,mtrans=NULL,cens=NULL){

#checking conditions
if (missing(data))
	stop("Argument 'data' is missing with no default")
if(!missing(states) && (!all(states %in% unique(c(data[,2],data[,3]))) || !all(unique(c(data[,2],data[,3])) %in% states)))
	warning("The states defined differ from those in 'data'")
if(!is.data.frame(data))
	stop("Argument 'data' must be a data.frame")
if (!missing(mtrans) && nrow(mtrans) != ncol(mtrans)) 
      stop("Argument 'mtrans' must be a quadratic  matrix.")
## is mtrans a quadratic matrix
if(!missing(mtrans) && nrow(mtrans)!=ncol(mtrans))stop("'mtrans' must be a quadratic matrix")
##transitions on diagonal
if(!missing(mtrans)){
k<-0
for(i in 1:dim(mtrans)[1]){
if(mtrans[i,i]!=FALSE)k<-k+1}
}
if (!missing(mtrans) && k> 0) 
       stop("Transitions into the same state are not allowed")
if (!missing(mtrans) && !missing(states) && nrow(mtrans) != length(states)) 
       stop("The row number of 'mtrans' must be equal to the number of states.")
 if (length(states) != length(unique(states))) 
        stop("The state names must be unique.")

if(!missing(cens)){
for(i in 1:dim(data)[1]){
if(data[i,3]==cens){
if(i<dim(data)[1] && data[i+1,1]==data[i,1]){

stop("The censure is not the last state of individuals")}
}
}
}
if (missing(states)){states1<-sort(c(unique(data[,2]),data[,3]))
states2<-sort(unique(data[,3]))
states<-sort(unique(c(states1,states2)))
k<-0
if(!missing(cens)){
k<-which(states==cens)
if(k==0)warning("The censorship is not observed in data")
states<-states[-which(states==cens)]
}
}

s<-length(states)
trans.h<-data[data[,2]!=data[,3],2]
trans.j<-data[data[,2]!=data[,3],3]
if(missing(cens)){
trans.cens<-data[data[,2]==data[,3],3]}
else{
trans.cens<-data[data[,3]==cens,3]}
table.state<-matrix(ncol=s, nrow=s) 
for(i in 1:s){
for (j in 1:s){
if(i!=j){
table.state[i,j]<-1*(sum(1*(trans.h==states[i] & trans.j==states[j])))}
else{table.state[i,j]<-1*(sum(1*(trans.cens==states[i])))}
}
}
if(!missing(cens)){
for(i in 1:dim(data)[1]){
if(data[i,3]==cens){

data[i,3]<-data[i,2]
}
}}
#number of censured patients
Ncens<-length(which(data[,2]==data[,3]))
dimnames(table.state)<-list(states,states)
table.state<-list(table.state=table.state,Ncens=Ncens)
table.state
}
