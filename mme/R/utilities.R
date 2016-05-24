
add <- function(x) Reduce("+", x)


############################################################################
#' Function to generate matrices and the initial values
#'
#' Based on the input data, this function generates some matrices that are
#' required in subsequent calculations and the initial values obtained from the function \code{\link[mme]{initial.values}}.
#'
#' @param fi input data set with (d x t) rows and 4+k+sum(pp) columns. The first four columns of the data set are, in this order: the area indicator (integer), the time indicator (integer), the sample size of each domain and the population size of each domain. The following k columns are the categories of the response variable. Then, the auxiliary variables: first, the auxiliary variables of the first category, second, the auxiliary variables of the second category, and so on. Examples of input data sets are in \code{\link[mme]{simdata}}, \code{\link[mme]{simdata2}} and \code{\link[mme]{simdata3}}. 
#' @param k number of categories of the response variable.
#' @param pp vector with the number of auxiliary variables per category.
#' @param mod a number specifying the type of models: 1=multinomial mixed model with one independent 
#' random effect in each category of the response variable
#' (Model 1), 2=multinomial mixed model with two independent random effects in each category of the 
#' response variable: one domain random effect and another independent time and domain random effect (Model 2) and
#'  3= multinomial model with two independent random effects in each category of the response variable: 
#' one domain random effect and another correlated time and domain random effect (Model 3).
#' @return A list containing the following components.
#' \item{n}{vector with the area sample sizes.}
#' \item{N}{vector with the area population sizes.}
#' \item{Z}{design matrix  of random effects.}
#' \item{Xk}{list of matrices with the auxiliary variables per category. The dimension of the list is the number of domains}
#' \item{X}{list of matrices with the auxiliary variables. The dimension of the list is the number of categories of the response variable minus one.}
#' \item{y}{matrix with the response variable. The rows are the domains and the columns are the categories of the response variable.}
#' \item{initial}{a list with the initial values for the fixed and the random effects obtained from \code{\link[mme]{initial.values}} .}
#' @seealso \code{\link[mme]{initial.values}}, \code{\link[mme]{wmatrix}},
#' \code{\link[mme]{phi.mult}}, \code{\link[mme]{prmu}},
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @keywords models
#' @export
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata2) #Data
#' mod=2
#'
#' ##Needed matrix and initial values
#' datar=data.mme(simdata2,k,pp,mod)

data.mme<-function(fi,k,pp,mod){
if (class(mmedata(fi,k,pp))!="mmedata") stop("Error: the input dataset is not mmedata class")	
d=length(levels(factor(fi[,1])))
t=length(levels(factor(fi[,2])))	
p=sum(pp)
MM=fi[,4]
MM=data.matrix(MM)
fi=data.matrix(fi)
D=nrow(fi)

M=fi[,3]
M=data.matrix(M)
y=fi[,5:(k+3)]
y=data.matrix(y)

#Auxiliary matrix per category

aux=list()
inn=k+5
f=k+4+pp[1]
for (i in 1:(k-1)){
	aux[[i]]=as.matrix(fi[,(inn:f)])
	colnames(aux[[i]])=colnames(fi)[(inn:f)]
	inn=f+1
	f=inn+pp[i+1]-1}


Xk=list()
X=list()
Z=list()
############# X
for(i in 1:D){
	X[[i]]=c(1,aux[[1]][i,],rep(0,sum(pp[2:(k-1)]+1)))
	for (j in 2:(k-1)){
		n0=sum(pp+1)-sum(pp[1:(j-1)]+1)-1-pp[j]
		if (n0<0) {n0=0}
		bb=c(rep(0,sum(pp[1:(j-1)]+1)),1,aux[[j]][i,],rep(0, n0))
		X[[i]]=rbind(X[[i]],bb)
		}}
############## Xk
for (i in 1:(k-1)){
	Xk[[i]]=cbind(1,aux[[i]])}


#################Z
for(i in 1:D){
	n0=(k-1)*(D-i)
	if (n0<0) {n0=0}
	Z[[i]]=cbind(matrix(0,nrow=k-1,ncol=(k-1)*(i-1)),diag(k-1),matrix(0,nrow=k-1,ncol=n0))}
#resp=cbind(y,(M-rowSums(y)))
resp=fi[,5:(k+4)]
datar=list(n=M,N=MM,Z=Z,Xk=Xk,X=X,y=resp)
initial=initial.values(d,pp,datar,mod)
result=list(n=M,N=MM,Z=Z,Xk=Xk,X=X,y=resp,d=d,t=t,initial=initial)
return(result)
}



############################################################################
#' Model variance-covariance matrix of the multinomial mixed models
#'
#' This function calculates the variance-covariance matrix of the multinomial mixed models.
#' Three types of multinomial mixed model are considered. The first model (Model 1), with one
#' random effect in each category of the response variable; Model 2, introducing
#' independent time effect; Model 3, introducing correlated time effect.
#'
#' @param M vector with area sample sizes.
#' @param pr matrix with the estimated probabilities for the categories of the
#' response variable obtained from \code{\link[mme]{prmu}} or \code{\link[mme]{prmu.time}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{phi.mult}}, \code{\link[mme]{prmu}}, \code{\link[mme]{prmu.time}}
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}
#' @return W a list with the model variance-covariance matrices for each domain.
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling,13,153-178.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=2 #type of model
#' data(simdata2)
#' datar=data.mme(simdata2,k,pp,mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#
#' ##The model variance-covariance matrix
#' varcov=wmatrix(datar$n,mean$estimated.probabilities)

wmatrix<-function(M,pr){
		sigmap=list()
		D=nrow(pr)
		k=ncol(pr)
			for (i in 1:D){
			sigmap[[i]]=M[i]*(diag(pr[i,(1:(k-1))])-(pr[i,(1:(k-1))])%*%t(pr[i,(1:(k-1))]))}
		return(W=sigmap)}


####################################################################
#' Standard deviation and p-values of the estimated model parameters
#'
#' This function calculates the standard deviations and the p-values
#' of the estimated model parameters. The standard deviations are obtained from the asymptotic Fisher information matrix in the fitting
#' algorithms \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}}, \code{\link[mme]{modelfit3}},
#' depending of the current multinomial mixed model.
#'
#' @param a vector with the estimated parameters obtained from \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}}.
#' @param F inverse of the Fisher Information Matrix obtained from \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}}.
#' @keywords models
#' @export
#' @seealso
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}},
#'  \code{\link[mme]{modelfit3}}.
#' @return A list containing the following components.
#' \item{Std.dev}{vector with the standard deviations of the parameters. The parameters are sorted per category. }
#' \item{p.value}{vector with the p-values of the parameters for testing H0:a=0.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata) #data
#' mod=1  #Type of model
#' datar=data.mme(simdata,k,pp,mod)
#'
#' #Model fit
#' result=modelfit1(pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N)
#' beta=result[[8]][,1] #fixed effects
#' Fisher=result[[3]] #Fisher information matrix
#'
#' ##Standard deviation and p-values
#' res=ci(beta,Fisher)

ci<-function(a,F){
qqq=sqrt(diag(F))
zz=as.vector(a/qqq)
ppp=as.vector(2*pnorm(abs(zz),lower.tail=FALSE))
output=list(Std.dev=qqq,p.value=ppp)}


####################################################################
#' Add items from a list
#'
#' This function adds items from a list of dimension d*t, where d is the number of areas
#' and t is the number of times periods.
#'
#' @param B_d a list in each area.
#' @param t number of time periods.
#' @param d number of areas.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{Fbetaf.ct}},
#' \code{\link[mme]{modelfit2}}, \code{\link[mme]{modelfit3}}
#' @return B_d a list of dimension d.
#' @examples
#' k=3  #number of categories for the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata2) # data
#' mod=2
#' datar=data.mme(simdata2,k,pp,mod)
#'
#' ##Add the time periods
#' l=addtolist(datar$X,datar$t,datar$d)


addtolist<-function(B_d,t,d){
n=nrow(B_d[[1]])
nn=ncol(B_d[[1]])
md=rep(t,d)
u=cumsum(md)
B2_d=list()
B2_d[[1]]=add(B_d[1:u[1]])
for (i in 2:d){
	B2_d[[i]]=add(B_d[(u[i-1]+1):u[i]])}
return(B2_d)
rm(B_d)
}

####################################################################
#' Add rows from a matrix
#'
#' This function adds rows from a matrix of dimension d*t*(k-1) times d*(k-1). 
#'
#' @param C2 a matrix of dimension d*t*(k-1) times d*(k-1).
#' @param d number of areas.
#' @param t number of time periods.
#' @param k number of categories of the response variable.
#'
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{Fbetaf.ct}},
#' \code{\link[mme]{modelfit2}},\code{\link[mme]{modelfit3}}
#' @return C22 a matrix of dimension d*(k-1) times d*(k-1).
#' @examples
#' k=3 #number of categories of the response variable
#' d=15 # number of areas
#' t=2  # number of time periods
#' mat=matrix(1,d*t*(k-1),d*(k-1)) # a matrix
#'
#' ##Add items in the matrix
#' mat2=addtomatrix(mat,d,t,k)


addtomatrix<-function(C2,d,t,k){
j=1
s=seq(1,((k-1)*d),by=(k-1))
C22=matrix(0,ncol(C2),ncol(C2))
for (i in s){
	for (l in 1:t){
		C22[i:(i+k-2),]=C22[i:(i+k-2),]+C2[j:(j+k-2),]
		j=j+k-1}}

return(C22)
rm(C2)
}
####################################################################
#' Bias and MSE using parametric bootstrap
#'
#' This function calculates the bias and the mse for the multinomial mixed effects models
#' using parametric bootstrap. Three types of multinomial mixed models are considered, with one independent domain random effect in each category of the response variable  (Model 1),
#' with two random effects: the first, with a domain random effect and with independent time and domain random effect (Model 2) and the second, with a domain random effect and with correlated time and domain random effect (Model 3).
#' See details of the parametric bootstrap procedure in Gonzalez-Manteiga et al. (2008) and in Lopez-Vizcaino et al. (2013)
#' for the adaptation to these three models. This function uses the output of \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}},
#' depending of the current multinomial mixed model.
#'
#' @param pp vector with the number of the auxiliary variables per category.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix of random effects obtained from \code{\link[mme]{data.mme}}.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @param resul output of the function \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}}.
#' @param B number of bootstrap replications.
#' @param mod a number specifying the type of models: 1=multinomial mixed model with one independent random effect in each category of the response variable
#' (Model 1), 2=multinomial mixed model with two independent random effects in each category of the response variable: one domain random effect and another independent time and domain random effect (Model 2) and
#'  3= multinomial model with two independent random effects in each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3).
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult}}, \code{\link[mme]{phi.mult.it}},
#' \code{\link[mme]{phi.mult.ct}}, \code{\link[mme]{prmu}}, \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{phi.direct.it}}, \code{\link[mme]{phi.direct.ct}}, \code{\link[mme]{sPhikf}}, \code{\link[mme]{sPhikf.it}},
#' \code{\link[mme]{sPhikf.ct}}, \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}},
#' \code{\link[mme]{modelfit3}}, \code{\link[mme]{omega}},
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{Fbetaf.ct}},
#' \code{\link[mme]{ci}}.
#' @return a list containing the following components.
#' \item{bias.pboot}{BIAS of the parametric bootstrap estimator of the mean of the response variable}
#' \item{mse.pboot}{MSE of the parametric bootstrap estimator of the mean of the response variable}
#' \item{rmse.pboot}{RMSE of the parametric bootstrap estimator of the mean of the response variable}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13 ,153-178.
#' @references
#' Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @references
#' Gonzalez-Manteiga, W, Lombardia, MJ, Molina, I, Morales, D, Santamaria, L (2008). Estimation of the mean
#' squared error of predictors of small area linear parameters under a logistic mixed model, Computational Statistics
#' and Data Analysis, 51, 2720-2733.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)
#' mod=1  # Type of model
#' datar=data.mme(simdata,k,pp,mod)
#' ##Model fit
#' result=modelfit1(pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],datar$n,datar$N)
#'
#' B=1    #Bootstrap iterations
#' ss=12345 #SEED
#' set.seed(ss)
#'
#' ##Bootstrap parametric BIAS and MSE
#' mse.pboot=mseb(pp,datar$Xk,datar$X,datar$Z,datar$n,datar$N,result,B,mod)



mseb=function (pp, Xk, X, Z, M, MM, resul, B, mod) 
{
    beta.new2 = matrix(resul[[8]][, 1], length(resul[[8]][, 1]), 
        1)
    D = length(M)
    k = ncol(resul$mean) + 1 
       if (mod == 1) {
        phi.new = matrix(resul[[9]][, 1], length(resul[[9]][, 
            1]), 1)
        phiasterisco = phi.new
        uasterisco = matrix(0, D, (k - 1))
    }
    if (mod == 2 | mod == 3) {
        d = nrow(resul$u1)
        t = D/d
        phi1.new = as.matrix(resul[[9]][1:(k - 1), 1])
        phi2.new = as.matrix(resul[[9]][k:(2 * k - 2), 1])
        phiasterisco1 = phi1.new
        phiasterisco2 = phi2.new
        uasterisco1 = matrix(0, d, (k - 1))
        uasterisco2 = matrix(0, D, (k - 1))
    }
    if (mod == 3) {
        rho.new = resul$rho
        rhoasterisco = rho.new
        ceros = matrix(rep(0, t), t, 1)
    }
    yasterisco = matrix(0, D, k)
    yasteriscor = matrix(0, D, k)
    ppp1 = list()
    for (i in 1:D) {
        ppp1[[i]] = matrix(0, (k - 1), (k - 1))
    }
    sesgo = matrix(0, D, (k - 1))
    suma = matrix(0, D, (k - 1))
    media = matrix(0, D, (k - 1))
    mse.bootstrap = matrix(0, D, (k - 1))
    mse1 = list()
    mse2 = list()
    ps = cumsum(pp + 1)
    beta.new = list()
    beta.new[[1]] = matrix(beta.new2[1:ps[1], ], (pp[1] + 1), 
        1)
    for (i in 2:(k - 1)) {
        beta.new[[i]] = matrix(beta.new2[(ps[i - 1] + 1):ps[i], 
            ], (pp[i] + 1), 1)
    }
    betaasterisco = beta.new
    it = 1
    j = 1
    valido = 0
    while (j < (B + 1)) {
    	      for (l in 1:(k - 1)) {
        	
            if (mod == 1) {
                uasterisco[, l] = rnorm(D, mean = 0, sd = phiasterisco[l]^(1/2))
            }
            if (mod == 2) {
                uasterisco1[, l] = rnorm(d, mean = 0, sd = phiasterisco1[l]^(1/2))
                uasterisco2[, l] = rnorm(D, mean = 0, sd = phiasterisco2[l]^(1/2))
            }
            if (mod == 3) {
                a = omega(t, k, rhoasterisco, phiasterisco2)
                uasterisco1[, l] = rnorm(d, mean = 0, sd = phiasterisco1[l]^(1/2))
                datos = mvrnorm(d, ceros, ((phiasterisco2[l]) * 
                  (a[[1]][[l]])))
                uasterisco2[, l] = matrix(t(datos), D, 1)
            }
        }
        if (mod == 1 ) {
            prmul = prmu(M, Xk, betaasterisco, uasterisco)
        }
        if (mod == 2 | mod == 3) {
            prmul = prmu.time(M, Xk, betaasterisco, uasterisco1, 
                uasterisco2)
        }
        theta = prmul[[3]]
        pr = prmul[[1]]
        mu = prmul[[2]]
        for (i in 1:D) {
            yasterisco[i, ] = t(rmultinom(1, M[i], pr[i, ]))
            yasteriscor[i, ] = t(rmultinom(1, MM[i] - M[i], pr[i, 
                ]))
        }
        yast = yasterisco[, 1:(k - 1)]
        yastt = yasterisco[, 1:(k - 1)] + yasteriscor[, 1:(k - 
            1)]
        if (mod == 1) {
            u.newi = matrix(rep(0.5, ((k - 1) * D)), D, (k - 
                1))
            initial = list(beta.new, as.vector(phi.new), u.newi)
            par = modelfit1(pp, Xk, X, Z, initial, yast, M, MM)
        }
        if (mod == 2) {
            u1.new = rep(0.01, ((k - 1) * d))
            dim(u1.new) = c(d, k - 1)
            u2.new = rep(0.01, ((k - 1) * D))
            dim(u2.new) = c(D, k - 1)
            initial = list(beta.new, as.vector(phi1.new), as.vector(phi2.new), 
                u1.new, u2.new)
            par = modelfit2(d, t, pp, Xk, X, Z, initial, yast, 
                M, MM)
        }
        if (mod == 3) {
            u1.new = rep(0.01, ((k - 1) * d))
            dim(u1.new) = c(d, k - 1)
            u2.new = rep(0.01, ((k - 1) * D))
            dim(u2.new) = c(D, k - 1)
            initials = list(beta.new, as.vector(phi1.new), as.vector(phi2.new), 
                u1.new, u2.new, as.vector(rho.new))
            par = modelfit3(d, t, pp, Xk, X, Z, initials, yast, 
                M, MM, 1)
        }
        valido = par$warning1
       
        if (valido == 0) {
            betaasterisco2 = matrix(par[[8]][, 1], length(par[[8]][, 
                1]), 1)
            if (mod == 1) {
                phiasterisco2 = matrix(par[[9]][, 1], length(par[[9]][, 
                  1]), 1)
            }
            if (mod == 2 | mod == 3) {
                phiasterisco12 = par[[9]][1:(k - 1), 1]
                phiasterisco22 = par[[9]][k:(2 * k - 2), 1]
            }
            if (mod == 3) {
                rhoasterisco2 = par$rho
            }
            muasterisco = par$mean
            prasterisco = par$Estimated.probabilities
            for (i in 1:D) {
                ppp1[[i]] = ppp1[[i]] + (yastt[i, ] - muasterisco[i, 
                  ]) %*% t(yastt[i, ] - muasterisco[i, ])
                sesgo[i, ] = sesgo[i, ] + (yastt[i, ] - muasterisco[i, 
                  ])
                suma[i, ] = suma[i, ] + muasterisco[i, ]
            }

            it = it + 1
        }
        if (j == B) {
            if (it < B) {
                j = j - (B - it) - 1
            }
        }
        j = j + 1
    }
    for (i in 1:D) {
        mse1[[i]] = ppp1[[i]]/it
        sesgo[i, ] = sesgo[i, ]/it
        media[i, ] = suma[i, ]/it
    }
    for (j in 1:D) {
        mse.bootstrap[j, ] = diag(mse1[[j]])
    }
    rmse.bootstrap = sqrt(mse.bootstrap)/media * 100
    colnames(mse.bootstrap) = paste("mse.pboot.", 1:(k - 1), 
        sep = "")
    colnames(sesgo) = paste("bias.pboot.", 1:(k - 1), sep = "")
    colnames(rmse.bootstrap) = paste("rmse.pboot.", 1:(k - 1), 
        sep = "")
    bootstrap = list(bias.pboot = sesgo, mse.pboot = mse.bootstrap, 
        rmse.pboot = rmse.bootstrap)
    return(bootstrap)
}


####################################################################
#' Choose between the three models
#'
#' This function chooses one of the three models. 
#'
#' @param d number of areas.
#' @param t number of time periods.
#' @param pp vector with the number of the auxiliary variables per category.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix of random effects obtained from \code{\link[mme]{data.mme}}.
#' @param initial output of the function \code{\link[mme]{initial.values}}.
#' @param y matrix with the response variable obtained from \code{\link[mme]{data.mme}}. The rows are the domains and the columns are the categories of the response variable.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @param mod a number specifying the type of models: 1=multinomial mixed model with one independent random effect in each category of the response variable
#' (Model 1), 2=multinomial mixed model with two independent random effects in each category of the response variable: one domain random effect and another independent time and domain random effect (Model 2) and
#'  3= multinomial model with two independent random effects in each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3).
#'
#' @keywords models
#' @export
#' @return  the output of the function \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}}.
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13 ,153-178.
#' @references
#' Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)  #data
#' mod=1 #Model 1
#' datar=data.mme(simdata,k,pp,mod)
#' result=model(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#' datar$n,datar$N, mod)

model <- function(d, t, pp, Xk, X, Z, initial, y, M, MM, mod){	
	
         if (mod==1) MODEL = modelfit1(pp, Xk, X, Z, initial, y, M, MM)
         if (mod==2) MODEL = modelfit2(d, t, pp, Xk, X, Z, initial, y, M, MM)
         if (mod==3) MODEL = modelfit3(d, t, pp, Xk, X, Z, initial, y, M, MM, 0)   
         return(MODEL)      
}


####################################################################
#' Print objects of class mme
#'
#' This function prints objects of class mme. 
#'
#' @param x a list with the output of \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}}. 
#' @param ... further information.
#'
#' @keywords class
#' @export
#' @seealso \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}},
#' \code{\link[mme]{modelfit3}}
#' @method print mme
#' @examples
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=1  # Type of model
#' data(simdata)
#' datar=data.mme(simdata,k,pp,mod)
#' ##Model fit
#' result=modelfit1(pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],datar$n,datar$N)
#' result


print.mme=function(x,...){
	cat("\nMultinomial mixed effects model\n\nCall:\n")

	if(x$warning1==0){
	cat("\nCoefficients\n")
	if (length(x$beta.Stddev.p.value)){
	print(zapsmall(x$beta.Stddev.p.value),digits=3)
	cat("\nRandom effects\n")
	print(zapsmall(x$phi.Stddev.p.value),digits=3)}
	if (length(x$rho.Stddev.p.value)){
		cat("\nCorrelation random effects\n")
		print(zapsmall(x$rho.Stddev.p.value),digits=3)
	cat("\n")
		
	}		
	}
	if(x$warning1==1){
		cat("\nThe model could not be fitted")
	}
		if(x$warning2==1){
		cat("\nThe value of the variance component is negative: the initial value is taken")
	}
	
	
}

####################################################################
#' Create objects of class mmedata
#'
#' This function creates objects of class mmedata. 
#'
#' @param fi input data set with (d X t) rows and 4+k+sum(pp) columns. The first four columns of the data set are, in this order: the area indicator (integer), the time indicator (integer), the sample size of each domain and the population size of each domain. The following k columns are the categories of the response variable. Then, the auxiliary variables: first, the auxiliary variables of the first category, second, the auxiliary variables of the second category, and so on. Examples of input data set are in \code{\link[mme]{simdata}}, \code{\link[mme]{simdata2}} and \code{\link[mme]{simdata3}}.
#' @param k number of categories of the response variable.
#' @param pp vector with the number of auxiliary variables per category.
#'
#' @keywords data
#' @export
#' @seealso \code{\link[mme]{modelfit1}}, \code{\link[mme]{modelfit2}},
#' \code{\link[mme]{modelfit3}}
#' @examples
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)
#' r=mmedata(simdata,k,pp)

mmedata=function(fi,k,pp){
	c=4+k+sum(pp)
	if (length(fi)!=c) stop("Dataset number of columns are incorrect")
	ppp=c-4-k
	if(sum(pp)!=ppp) stop("There is a problem with the number of the auxiliary variables")
	d=length(levels(factor(fi[,1])))
	t=length(levels(factor(fi[,2])))
	D=d*t
	if (D!=nrow(fi)) stop("The number of rows of the dataset must be the number of areas(d) times the number of time periods(t) (dXt)")
	if (!isTRUE(all(fi[,1] == floor(fi[,1])))) stop("Area indicator must only contain integer values")
	if (!isTRUE(all(fi[,2] == floor(fi[,2])))) stop("Time indicator must only contain integer values")
	if (!isTRUE(all(fi[,3] <= (fi[,4])))) stop("Population size must be higher than the sample size for all the areas")
	s=rowSums(fi[5:(4+k)])
	if(!isTRUE(all(fi[,3] ==s))) stop("The total sample size for each of the categories does not match the total sample size ")
	fi=as.matrix(fi)
	class(fi)="mmedata"
	return(fi)
	}




