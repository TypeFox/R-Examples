# Test Step 1 of ISOpure using the saved data

# load library 
library(ISOpureR);

# load the data from that path
data.path <-  paste0(file.path(system.file(package = "ISOpureR"), 'extdata', 'Beer'));  
load(file.path(data.path , 'beer.tumordata.1000.transcripts.30.patients.RData'));
load(file.path(data.path , 'beer.normaldata.1000.transcripts.RData'));

# the normaldata and tumourdata should be matrices
beer.tumordata <- as.matrix(beer.tumordata);
beer.normaldata <- as.matrix(beer.normaldata);

# run step 1 of ISOpureR
set.seed(123);
ISOpureS1model <- ISOpure.step1.CPE(beer.tumordata, beer.normaldata);