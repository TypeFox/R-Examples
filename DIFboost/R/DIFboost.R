DIFboost <-
function(Y, X, mstop = 400, trace = TRUE, cutoff = 0.9, B = 500, mc.cores = 1, q=0.6*I){

  vardiffs <- abs(apply(X,2,var)-1)
  
  if(sum(vardiffs)>1e-6)
    stop("X has to be standardized")
  
  if(sum(!apply(Y,2,unique) %in% c(0,1))>0)
    stop("Y may only contain 0 or 1")
  
  if(!is.data.frame(X))
    stop("X has to be a data.frame")
  
  if(!is.data.frame(Y))
    stop("Y has to be a data.frame")
  
# number of objects/persons
P <- nrow(Y)

# number of items 
I <- ncol(Y)

# number of covariates
l <- ncol(X)

# Convert data.frames into matrices
names.y <- names(Y)
names.x <- names(X)

Y <- matrix(as.double(as.matrix(Y)),ncol=I,byrow=FALSE)
X <- matrix(as.double(as.matrix(X)),ncol=l,byrow=FALSE)

# which objects to exclude with 0 or I right answers
exclu1 <- rowSums(Y)!=I & rowSums(Y)!=0
exclu <- rep(exclu1, each=I)

## create columns of design matrix with covariates
xp <- matrix(0,ncol=l*I,nrow=P*I)

suppressWarnings(
  for (qq in 1:P){
    xp[((qq-1)*I+1):(qq*I),] <- matrix(c(X[qq,],rep(0,l*I)),byrow=T,ncol=l*I,nrow=I)
  })

# create reponse vector
XP <- c()
for(t in 1:P){XP <- c(XP,Y[t,])}

# new number of object after excluding 
ex <- sum(!exclu1)
P <- P - ex

## create design matrix
help1 <- c(1,rep(0,P-2))

help2 <- c(rep(help1,I),0)

suppressWarnings(
  help3 <- matrix(help2, ncol = P-1,nrow=P*I, byrow=TRUE))

# columns for person parameters
help3[((P-1)*I+1):(P*I),] <- 0

help4 <- rep(diag(I),P)

# columns for item parameters
help5 <- matrix(help4,ncol=I,nrow=P*I,byrow = TRUE)

# columns for itemspecific parameters
xp <- xp[exclu,]

# design matrix
design.matrix <- cbind(help3,-help5,-xp)

# factorize response vector
XP <- as.factor(XP[exclu])

# create data set for boosting
boost.dat <- as.data.frame(cbind(design.matrix,XP))

boost.dat$XP <- as.factor(boost.dat$XP)

names(boost.dat) <- c(paste("V",1:ncol(design.matrix),sep=""),"XP")

# calculate ML estimates for offset in boosting
m.glm <- glm.fit(y=XP,x=design.matrix[,1:(P+I-1)],family=binomial(),intercept=FALSE)
offset <- m.glm$linear.predictors/2

### create formula for boosting
act <- P + I

# formula for grouped itemspecific parameters (gamma)
form2 <-paste("bols(",paste("V",act:(act+l-1),sep="",collapse=","),",intercept=FALSE,df=1)",sep="")
act <- act + l
for(u in 2:I){
form.help <- paste("V",act:(act+l-1),sep="",collapse=",")
act <- act + l
form.help2 <- paste("bols(",form.help,",intercept=FALSE,df=1)",sep="")
form2 <- paste(form2,form.help2,sep="+")
}

# boosting formula
form3 <- as.formula(paste("XP~",form2,sep=""))

# Boosting 
suppressWarnings(
 m1 <- gamboost(form3,data=boost.dat,family=Binomial(),offset=offset,control=boost_control(mstop=mstop,center=FALSE,trace=trace))
)

# create folds for stability selection
folds <- matrix(0,ncol=B,nrow=P*I)
index <- rep(1:P,each=I)
for(o in 1:B){
smp<-  sample(1:P,floor(P/2))
 folds[index %in% smp,o] <- 1
}

# perform stability selection
m2 <- stabsel(m1,q=q ,cutoff=cutoff, folds = folds, mc.cores=mc.cores)

# extract selected base learners/items
selected <- as.numeric(which(apply(m2$phat,1,max)>cutoff))

# extract selected columns from design matrix
sel.vars <- sort(matrix(seq(1,l*I),ncol=l,nrow=I,byrow=TRUE)[selected,])

####################
# which item will be the reference/anchor item
ref.item <- max(which(!(1:I) %in% selected))

help3 <- cbind(help3,c(rep(0,((P-1)*(I))),rep(1,I)))
help5 <- help5[,-ref.item]

# create new design matrix for refit
design.matrix <- cbind(help3,-help5,-xp[,sel.vars])

XP <- as.numeric(XP)-1

# create data set for refit
refit.data <- as.data.frame(cbind(design.matrix,XP))
names(refit.data) <- c(paste("V",1:ncol(design.matrix),sep=""),"XP")

# create formula for refit
form1 <- as.formula("~ 0")

form2 <- as.formula(paste("~ ",paste("V",1:ncol(design.matrix),sep="",collapse="+"),sep=""))

# refit with DIF items
m3 <- penalized(response=XP,unpenalized=form1,penalized=form2,lambda1=0,lambda2=0.0001,data=refit.data,model="logistic")

# extract relevant estimates
coefs <- coef(m3)

theta.hat <- head(coefs,P)
beta.hat <- coefs[(P+1):(P + I - 1)]
gamma.hat <- tail(coefs, length(sel.vars))

# create matrix with all gamma estimates, zeros for non DIF items
dif.mat <- matrix(0,nrow=l,ncol=I)
dif.mat[,selected] <- gamma.hat 
colnames(dif.mat) <- names.y
rownames(dif.mat) <- names.x

lin.pred <- design.matrix%*%coefs
 
returns <- list(model = m1, dif.mat = dif.mat, coefficients = coefs, theta = theta.hat, 
                beta = beta.hat, gamma = gamma.hat, P = P, I = I, names.x = names.x, names.y = names.y,
                design.matrix = design.matrix, PFER = m2$PFER, lin.pred = lin.pred, 
                DIF.Items = selected, ref.item = ref.item, phat = m2$phat, cutoff = cutoff)
  
class(returns) <- "DIFboost"

return(returns)
}
