################################################################################################
#  R file which peforms 10K bootstrap analysis !!!
#  We ask R to create a new seedand script file for use with JAGS in each bootstrap iteration
#  order based searches 320 x 32 cpu --> with for each, I can do it time 12!!!
################################################################################################
setwd("~/compute/PigsData_Analysis/Jags")
load("~/compute/PigsData_Analysis/RData/vienna_pigs.RData") 

library(abn);
library(coda);
orig.data<-pigs.vienna[,-11]; ## get the Pigs data - drop batch variable
max.par<-3;#parent limit for original data
start<-seq(1,10240,by=320);
stop<-seq(320,10240,by=320);

## now have the boot.data in identical format to original to now repeat exact search.
ban<-matrix(rep(0,dim(orig.data)[2]^2),ncol=dim(orig.data)[2]);
## the ban matrix must have names set
colnames(ban)<-rownames(ban)<-names(orig.data);

retain<-matrix(rep(0,dim(orig.data)[2]^2),ncol=dim(orig.data)[2]);# retain nothing
## again must set names
colnames(retain)<-rownames(retain)<-names(orig.data);

## setup distribution list for each node
mydists<-list(
  PC="binomial",
  PT="binomial",
  MS="binomial",
  HS="binomial",
  TAIL="binomial",
  Abscess="binomial",
  Pyaemia="binomial",
  EPcat="binomial",
  PDcat="binomial",
  plbinary="binomial"
);


dags<-list();
#########################################
for(i in start[index]:stop[index]){ #MASTER LOOP - each interation creates a bootstrap sample and finds mostprobable model
#########################################
   #create bootstrap data
   #1. create parameter file with unique random seed 
   init.file<-paste("init_",i,sep="");#tempfile(pattern=paste("_",index,"_",sep=""),tmpdir=getwd());#file to hold jags seed
   cat(paste("\".RNG.name\" <-\"base::Mersenne-Twister\"","\n",sep=""),file=init.file,append=FALSE);
   cat(paste("\".RNG.seed\" <- ",i,"\n",sep=""),file=init.file,append=TRUE);#note i is unique
   #2. create script file with unique output file name
   run.file<-paste("script_",i,sep="");
   #tempfile(pattern=paste("_",index,"_",sep=""),tmpdir=getwd());
   #file to hold jags seed

#this is needed verbatim     
cat("model in pigs_model.bug
data in  pigs_post_params.R
compile, nchains(1)
",file=run.file);
cat(paste("parameters in ",init.file,"\n",sep=""),file=run.file,append=TRUE);
cat("initialize
update 100000
monitor PC, thin(10)
monitor PT, thin(10)
monitor MS, thin(10)
monitor HS, thin(10)
monitor TAIL, thin(10)
monitor Abscess, thin(10)
monitor Pyaemia, thin(10)
monitor EPcat, thin(10)
monitor PDcat, thin(10)
monitor plbinary, thin(10)
update 250000, by(1000)
",file=run.file,append=TRUE);
   out.file<-paste("out_",i,sep="");
   cat(paste("coda *, stem(\"",out.file,"\")\n",sep=""),file=run.file,append=TRUE);
 
   #3. run the MCMC sampler
   system(paste("jags ",run.file,sep="")); # usr/sepp/bin/jags

   #4. read in mcmc data and convert to format suitable for mostprobable
   boot.data<-read.coda(paste(out.file,"chain1.txt",sep=""),
                        paste(out.file,"index.txt",sep=""));
   boot.data<-as.data.frame(boot.data);
   for(j in 1:dim(orig.data)[2]){if(is.factor(orig.data[,j])){
       boot.data[,j]<-as.factor(boot.data[,j]);
       levels(boot.data[,j])<-levels(orig.data[,j]);}}

   #5. run the MostProb search on the bootstrap data
   boot1.cache<-buildscorecache(data.df=boot.data,data.dists=mydists, 
                max.parents=max.par, dag.banned=ban,dag.retained=retain);
   dags[[i]]<-mostprobable(score.cache=boot1.cache);
   unlink(c(init.file,run.file,out.file,paste(out.file,"chain1.txt",sep=""),
                            paste(out.file,"index.txt",sep="")));#tidy up

datadir <- '~/compute/PigsData_Analysis/RData/Pigs_BootstrappingResults/' # Pigs_BootstrappingResults folder created later and so edit later, in order to be a bit more tidy and clear.

save(dags,file=paste(datadir,"mp10Kboot",index,".RData",sep=""));

}



