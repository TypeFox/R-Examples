# NOTE: This example needs demo04cpu.r to be run before, this will create demo04cpu.RData.
path <- paste(find.package("TDMR"), "demo01cpu/",sep="/");
tdm <- list(  filenameEnvT="demo04cpu.RData" );   # file with environment envT 
           
load(paste(path,tdm$filenameEnvT,sep="/"));
newdata=read.csv2(file=paste(path,"data/cpu.csv", sep=""), dec=".")[1:15,];     # take only the first 15 records
z=predict(envT,newdata);
print(z);

