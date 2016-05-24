pedroni99 <-
function(Y, X, kk=0, type.stat=1, ka=2){
# calculates Pedroni(1999) panel cointegration statistics, bivariate case
# the function strictly follows the text in :
#  Pedroni, Peter, 1999. "Critical Values for Cointegration Tests in Heterogeneous Panels with Multiple Regressors," Oxford Bulletin of Economics and Statistics, Department of Economics, University of Oxford, vol. 61(0), pages 653-70, Special I.
 
# Y, X - variables, in matrix T*N form, no NA values (balanced panels required) 
# (time in columns, each column is individual series, y and x must be equal size)
# kk - parameter (number of lags) for the Newey-West (1994) procedure, scalar (same value for all series) or vector 
# If different for each individual series, must be a vector; by default is 1/3 of the series lenght.
# type.stat - type of the main model - 1 - no intercept and trend, 2 - intercept, 3 - intercept and time trend
# ka - parameter for lag selection in parametric statistics, scalar (same value for all series) or vector

# extract residuals functions
ff<- function(Y1,X1){ NN=ncol(X1); sapply(1:NN, function (l) {lm(Y1[,l]~X1[,l]-1)$residuals})}
ff1<-function(Y1,X1){ NN=ncol(X1); sapply(1:NN, function (l) {lm(Y1[,l]~X1[,l])$residuals})}
ff2<-function(Y1,X1){ NN=ncol(X1); trend=1:nrow(X1); sapply(1:NN, function (l) {lm(Y1[,l]~X1[,l]+trend)$residuals})}

# Newey-West longterm variance 
nw<-function(xx,ki){
tt=length(xx); (1/tt)*sum(sapply(1:ki, function(s) {(1-s/(ki+1))*sum(xx[(s+1):tt]*xx[1:(tt-s)])} )) }

# ADF type function, needed for the parametric statistics
adfl<-function (ee, lags) {
nn<-length(ee)
z<-ee[(lags+1):nn]
zl<-ee[lags:(nn-1)]

zd<-matrix(cbind(rep(z,lags)),ncol=lags)
ii<-embed(1:nn,lags)
ii<-ii[-(nrow(ii)),]
zd<-zd-ee[ii]
zd<-zd[,-1]
z<-ee[(lags+1):nn]
zl<-ee[lags:(nn-1)]
return(lm(z ~ zl + zd -1)$residuals)
}

# data check
Y<-as.matrix(Y)
X<-as.matrix(X)
if (any((dim(Y) != dim(X)))) {stop("Y and X are not compatible.") };
na.fail(Y)
na.fail(X)

TD = nrow(X)
N = ncol(X)

# check parameter for Newey-West (1994) procedure
if (is.vector(kk) && length(kk)==N)  {
 k=kk } else if (kk>0) {
 k = rep(round(kk),N)
 }else {i=round(4*(TD/100)^(2/9)); k = rep(i,N) }

# check lags for parametric statistics
if (ka < 2 ) {ka = 2; warning("Parameter 'ka' was changed to 2.")}
ka<-as.vector(ka)
if (length(ka) != N)  {ka<-rep(ka[1],N) }

stats<-matrix(nrow=7,ncol=2)
rownames(stats)<-c("nipanel", "rhopanel", "tpanelnonpar", "tpanelpar", "rhogroup", "tgroupnonpar", "tgrouppar")
colnames(stats)<-c("empirical","standardized")

# set tabular values for the statistics - see Pedroni (1999), p. 666, Table 2
statsm<-cbind(c(6.982, -6.388, -1.662, -1.662, -9.889, -1.992, -1.992),
c(11.754, -9.495, -2.177, -2.177, -12.938, -2.453,  -2.453),
c(21.162, -14.011, -2.648, -2.648, -17.359, -2.872, -2.872))
rownames(statsm)<-c("nipanel", "rhopanel", "tpanel", "tpanelp", "rhogroup", "tgroup", "tgroupp")
colnames(statsm)<-c("none","intercept","trend")
statsv<-cbind(c(81.145, 64.288, 1.559, 1.559, 41.943, 0.649, 0.649),
c(104.546, 57.61, 0.964, 0.964, 51.49, 0.618, 0.618),
c(160.249, 64.219, 0.69, 0.69, 66.387, 0.555, 0.555))
rownames(statsv)<-c("nipanel", "rhopanel", "tpanel", "tpanelp", "rhogroup", "tgroup", "tgroupp")
colnames(statsv)<-c("none","intercept","trend")

# obtain the residuals from the main regression
e<-matrix(ncol=N,nrow=TD)
if (type.stat==2) {
e<-ff1(Y,X)
} else if (type.stat==3) {
e<-ff2(Y,X); 
} else {e<-ff(Y,X); type.stat=1}

De<-diff(e)

estar<-e #### 
Destar<-diff(estar)

# obtain the residuals from the regression on differences 
DX<-diff(X)
DY<-diff(Y)
eta<-matrix(ncol=ncol(DX),nrow=nrow(DX))
eta<-ff(DY,DX)

L11hat2<-sapply(1:N,function (i) {
(1/nrow(eta))*sum(eta[,i]^2) + 2*nw(eta[,i],k[i])})


mu<-matrix(ncol=ncol(DX),nrow=nrow(DX))
mu<-ff(e[2:TD,],e[1:(TD-1),])

lambdahat<-sapply(1:N,function (i) {nw(mu[,i],k[i])})

mustar<-matrix(ncol=ncol(DX),nrow=nrow(DX))
mustar<-sapply(1:N, function(i){adfl(e[,i],ka[i])})

shatstar2<-sapply(1:N,function (i) {(1/nrow(mustar))*sum(mustar[,i]^2)})
stildestar2<-(1/N)*sum(shatstar2)



shat2<-sapply(1:N,function (i) {(1/nrow(mu))*sum(mu[,i]^2)})

sigmahat2<-shat2+2*lambdahat
sigmatilde2<-(1/N)*sum(L11hat2^(-2)*sigmahat2)

nipa<-sum(sapply(1:N, function(i) {sum((L11hat2[i]^(-2))*(e[1:(TD-1),i]^2))}))

lel<-sum(sapply(1:N, function(i) {
(L11hat2[i]^(-2)) * sum(sapply(1:(nrow(De)), function(ttt) {
(e[(ttt),i]*De[ttt,i] - lambdahat[i]) } ) )  } ))
# the index (ttt) is the same for both e/De variables, because of the diff() reindexing

nipanel<-(TD^2)*(N^(3/2))*nipa^(-1)

stats[1,1]<-nipanel

rhopanel<-TD*(N^(1/2)) * (nipa^(-1)) * lel
stats[2,1]<-rhopanel

tpanelnonpar<- ((sigmatilde2 * nipa)^(-1/2)) * lel
stats[3,1]<-tpanelnonpar

tpanelpar<- ((stildestar2*sum(sapply(1:N, function(i) {
sum((L11hat2[i]^(-2)) * estar[1:(nrow(estar)-1),i]^2)})))^(-1/2)) * sum(sapply(1:N, function(i) {
 sum(sapply(1:(nrow(Destar)), function(ttt) {
(L11hat2[i]^(-2)) * (estar[ttt,i]*Destar[ttt,i]) } ) ) })) 
stats[4,1]<-tpanelpar

rhogroup<- TD*(N^(-1/2)) * sum(sapply(1:N, function(i) {
 ((sum(e[1:(nrow(e)-1),i]^2))^(-1)) * sum(sapply(1:(nrow(De)), function(ttt) {
 (e[ttt,i]*De[ttt,i] - lambdahat[i]) } ) )  }))  
stats[5,1]<-rhogroup

tgroupnonpar<- (N^(-1/2)) * sum(sapply(1:N, function(i) { ((sigmahat2[i] * sum(e[1:(nrow(e)-1),i]^2))^(-1/2)) * sum(sapply(1:(nrow(De)), function(ttt) { (e[(ttt),i]*De[ttt,i] - lambdahat[i]) } ) ) })) 
#tgroupnonpar<-1
stats[6,1]<-tgroupnonpar

tgrouppar<- (N^(-1/2)) * sum(sapply(1:N, function(i) {
 (sum(shat2[i]*estar[1:(nrow(estar)-1),i]^2))^(-1/2) *
 sum(estar[1:(nrow(estar)-1),i]*Destar[1:(nrow(estar)-1),i]) })) 
stats[7,1]<-tgrouppar

stats[,2]<-sapply(1:7, function (i) {(stats[i,1]-statsm[i,type.stat]*sqrt(N))/sqrt(statsv[i,type.stat])})

list(CALL=match.call(), METHOD="Pedroni(1999) panel tests for cointegration",  STATISTIC=stats)
}
