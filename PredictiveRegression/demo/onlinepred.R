# This program calls the IID predictor, the Gauss predictor, or the MVA predictor
# in the on-line mode.
# The estimated waiting time for a mid-range laptop in 2011:
# IID: about 20 sec; Gauss: below 10 sec; MVA: about 20 sec.

answer1 <- readline("Choose 'IID' (default), 'Gauss', or 'MVA' (the first letter is sufficient) ");
ans1 <- substr(answer1,1,1);
if ((ans1 == "G") || (ans1 == "g")) {
  subtitle = "for the Gauss predictor"
} else if ((ans1 == "M") || (ans1 == "m")) {
  subtitle = "for the MVA predictor"
} else {
  ans1 <- "I";
  subtitle = "for the IID predictor"}

answer2 <- readline("Choose the 'median-accuracy' (default) or 'error' plot (the first letter) ");
ans2 <- substr(answer2,1,1);
if ((ans2 == "e") || (ans2 == "E")) {
  title <- c("Error plot at 5%, 1%, and 0.5%\n",subtitle);
} else {
  ans2 <- "m";
  title <- c("Median-accuracy plot at 5%, 1%, and 0.5%\n",subtitle);}

start <- Sys.time()   # timing: start

# generate the artificial data set: start
N <- 600;  # number of observations
K <- 100;  # number of explanatory variables
XX <- array(rnorm(N*K),dim=c(N,K));  # the explanatory variables
beta <- array(0,dim=c(K,1));         # initializing the coefficients
beta[1:K,1] <- c(1,-1);
beta[1:10,1] <- c(10,-10);
yy <- 100 + XX %*% beta + array(rnorm(N),dim=c(N,1));	# the response variables
dataset <- cbind(XX,yy);
dim(dataset) <- c(N,K+1);
# generate the artificial data set: end

N <- dim(dataset)[1];           # number of observations
K <- dim(dataset)[2]-1;         # number of explanatory variables
K0 <- 10;                       # the first K0 attributes are regarded as more important
epsilons <- c(0.05,0.01,0.005); # the significance levels to be considered

# initializing err, up, and low with the values suitable for n=1 in the main loop
up <- array(Inf,c(N,length(epsilons)));
low <- array(-Inf,c(N,length(epsilons)));
err <- array(0,c(N,length(epsilons)));

for(n in 2:N) {     # main loop: start

  Kdagger <- ifelse(n<K+3,K0,K);      # taking into account
                                      # only the important explanatory variables
  if ((ans1 == "G") || (ans1 == "g")) # except for the Gauss predictor
    Kdagger <- K;                     # (see the paper, Remark 2)

  train <- array(0,dim=c(n-1,Kdagger+1));           # preparing the training set: start
  train[1:(n-1),1:Kdagger] <- dataset[1:(n-1),1:Kdagger];
  train[1:(n-1),Kdagger+1] <- dataset[1:(n-1),K+1]; # preparing the training set: end
  test <- array(0,dim=c(1,Kdagger));                # preparing the test set: start
  test <- dataset[n,1:Kdagger,drop=FALSE];
  test_response <- dataset[n,K+1];                  # preparing the test set: end

  if (ans1 == "I") {
    output <- iidpred(train,test,epsilons,0.01); } # using the IID predictor
  if ((ans1 == "G") || (ans1 == "g")) {
   output <- gausspred(train,test,epsilons); }     # using the Gauss predictor
  if ((ans1 == "M") || (ans1 == "m")) {
    output <- mvapred(train,test,epsilons,0.01); } # using the MVA predictor

  low[n,] <- output[[1]];
  up[n,] <- output[[2]];
  for (eps_index in 1:length(epsilons))
    err[n,eps_index] <- ifelse((low[n,eps_index]<=test_response)&&(test_response<=up[n,eps_index]),0,1);
}                   # main loop: end

width <- up-low;  # dimension: c(N,length(epsilons))
Err <- Med <- array(0,c(N,length(epsilons)));
for(n in 1:N) {
for (eps_index in 1:length(epsilons)) {
  Med[n,eps_index] <- median(width[1:n,eps_index]);
  Err[n,eps_index] <- sum(err[1:n,eps_index]);
}}

if (ans2 == "m") {
  Med <- pmin(Med,160);
  plot(1:N,Med[1:N,1],'l',ylim=c(0,150),main=title,
    xlab="observations",ylab="median accuracy");
  for (eps_index in 2:length(epsilons)) lines(1:N,Med[1:N,eps_index]);}

if ((ans2 == "e") || (ans2 == "E")) {
  plot(1:N,Err[1:N,1],'l',ylim=c(0,40),main=title,
    xlab="observations",ylab="cumulative errors");
  for (eps_index in 2:length(epsilons)) lines(1:N,Err[1:N,eps_index]);}

end <- Sys.time()   # timing: end
cat("Time elapsed:", difftime(end,start,units="secs"), "seconds\n")
