# Description: 	Modeling code for the example of ancient Chinese wars.  See page 125-127.
# Claudio Cioffi-Revilla and David Lai, 2001,
# "Chinese Warfare and Politics in the Ancient East Asian International System",
# hdl:1902.1/02016 http://id.thedata.org/hdl%3A1902.1%2F02016 Social Science Program
# Henry A. Murray Research Archive
# Center for International Relations, Department of Political Science, University of Colorado, Boulder, USA
# Download from http://dvn.iq.harvard.edu/dvn/dv/mra/faces/study/StudyPage.jsp?studyId=110&studyListingIndex=0_9aa809933d3c24eb76216fa477f8

library(nnet); library(mice)


#wars <- load("DA_cwp.rda")
wars <- read.table("data/DA_cwp.tab",header=TRUE)
wars <- cbind(as.matrix(wars[,4:15]))
imp.wars <- mice(wars,m=5)
complete.wars <- array(c(as.matrix(complete(imp.wars,1)), as.matrix(complete(imp.wars,2)), as.matrix(complete(imp.wars,2)),
                       as.matrix(complete(imp.wars,4)), as.matrix(complete(imp.wars,5))),dim=c(104,12,5))

wars <- complete.wars[,,1]   #  REPEAT FOR ARRAY DIMENSIONS 2-5
table(wars[,8:9])
wars <- cbind(wars,"SCOPE"=log( 10*wars[,8] + wars[,9]))
wars <- cbind(wars, "DURATION" = 1 + wars[,2]-wars[,1])
dimnames(wars) <- list(NULL,c("ONSET","TERM","EXTENT","ETHNIC","DIVERSE","ALLIANCE","DYADS","POL.LEV","COMPLEX","POLAR","BALANCE","TEMPOR","SCOPE","DURATION"))
wars <- data.frame(wars)

wars.lm <- lm(SCOPE ~ -1 + EXTENT + DIVERSE + ALLIANCE + DYADS + TEMPOR + DURATION, data=wars)
summary(wars.lm)

data(wars)
X <- cbind(as.numeric(wars$.EXTENT, wars$.DIVERSE,wars$.ALLIANCE,wars$.DYADS,wars$.TEMPOR,wars$.DURATION))
y <- as.numeric(wars$.SCOPE) 
n <- nrow(X); 
k <- ncol(X)
nu <- 5
num.sims <- 100000

war.samples <- matrix(NA,nrow=num.sims,(ncol=k+n+1))
beta <- rep(1,ncol(X)); sigma.sq <- 3; Omega <- 3*diag(n)

#for (i in 1:num.sims)  {
#    b <- solve(t(X) %*% X) %*% t(X) %*% y
#    b.star <- solve(t(X) %*% Omega %*% X) %*% t(X) %*% Omega %*% y
#    s.sq.star <- t(y-X%*%b) %*% solve(Omega) %*% (y-X%*%b)
#    u <- y - X %*% beta

#    beta <- as.vector( rmultinorm(1, b.star, sigma.sq *solve(t(X) %*% solve(Omega) %*% X) ) )
#    sigma.sq <- 1/rgamma(1, shape=(n-1)/2, rate=s.sq.star/2 )
#    for (j in 1:n)  Omega[j,j] <- 1/rgamma(1, shape=(nu+1)/2, rate=((sigma.sq^(-1))*u2 + nu)/2 )

#    war.samples[i,] <- c(beta,sigma.sq,diag(Omega))
#    if(i %% 100 == 0) print(i)
#}

#start <- 5001; stop <- num.sims
#round( cbind( apply(war.samples[start:num.sims,1:(k+1)],2,mean),
#       apply(war.samples[start:num.sims,1:(k+1)],2,sd),
#       apply(war.samples[start:num.sims,1:(k+1)],2,mean) - 1.96*apply(war.samples[start:num.sims,1:(k+1)],2,sd),
#       apply(war.samples[start:num.sims,1:(k+1)],2,mean) + 1.96*apply(war.samples[start:num.sims,1:(k+1)],2,sd) ) ,4)

