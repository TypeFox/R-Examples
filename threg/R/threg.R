### The following is the definition of the threg function

"threg" <-
function (formula,data) 
{
	### Read in all arguments
	cl <- match.call()

	indx <- match(c("formula", "data"),
		  names(cl), nomatch=0) 
	if (indx[1] ==0) stop("A formula argument is required")

	mf<- cl[c(1, indx)]
	f <- Formula(formula)
	f1<-formula(f,lhs=1)
	f1<-Formula(f1)
	mf[[1]] <- as.name("model.frame")

	mf$formula <- if(missing(data)) terms(f1)
		    else              terms(f1, data=data)  
	mf$formula <- f1
	mf <- eval(mf, parent.frame())


	if (nrow(mf) ==0) stop("No (non-missing) observations")
	Terms <- terms(mf)
	Y <- model.extract(mf, "response")
	if (!inherits(Y, "Surv")) stop("Response must be a survival object")
	type <- attr(Y, "type")
	if (type!='right')
	stop(paste("threg package current doesn't support \"", type,
			  "\" survival data", sep=''))

	f2<-formula(f1, lhs = 0)
	if (length(f2[[2]])!=3)
	stop(paste("Predictors for both lny0 and mu should be specified"))

	f_lny0 <-formula(f1,lhs=0,rhs=1)
	f_mu <-formula(f1,lhs=0,rhs=2)
        x_mu<-model.matrix(f_mu,data)
        x_lny0<-model.matrix(f_lny0,data)
        



        #len_lny0=length(dimnames(x_lny0)[[2]])
        #len_mu=length(dimnames(x_mu)[[2]])
        #lny0<-function(para,len_lny0){para[1]*lkr$f.treatment2+para[2]}
               # browser()
                #para(1:len_lny0)%*%x_lny0
        lny0<-function(para_lny0){x_lny0%*%para_lny0}

                 
                
        #mu<-function(para,len_lny0,len_mu){para((len_lny0+1):(len_lny0+len_mu))%*%x_mu}
        mu<-function(para_mu){x_mu%*%para_mu}





        #d<-function(para,len_lny0,len_mu) {-mu(para,len_lny0,len_mu)/exp(lny0(para,len_lny0))}
        #v<-function(para,len_lny0,len_mu) {exp(-2*lny0(para,len_lny0))}
        d<-function(para) 
        {
          para_lny0=para[1:length(dimnames(x_lny0)[[2]])]
          para_mu=para[(length(dimnames(x_lny0)[[2]])+1):(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]]))]
          -mu(para_mu)/exp(lny0(para_lny0))
        }
        v<-function(para) 
        {
          para_lny0=para[1:length(dimnames(x_lny0)[[2]])]
          exp(-2*lny0(para_lny0))
        }

        #logf<-function(para,len_lny0,len_mu) {
#-sum(failure*(-.5*(log(2*pi*v(para,len_lny0,len_mu)*(time^3))+(d(para,len_lny0,len_mu)*time-1)^2/(v(para,len_lny0,len_mu)*time))))- sum((1-failure)*log(pnorm((1-d(para,len_lny0,len_mu)*time)/sqrt(v(para,len_lny0,len_mu)*time))-exp(2*d(para,len_lny0,len_mu)/v(para,len_lny0,len_mu))*pnorm(-(1+d(para,len_lny0,len_mu)*time)/sqrt(v(para,len_lny0,len_mu)*time))))
#}

        logf<-function(para) {
          -sum(failure*(-.5*(log(2*pi*v(para)*(time^3))+(d(para)*time-1)^2/(v(para)*time))))- sum((1-failure)*log(pnorm((1-d(para)*time)/sqrt(v(para)*time))-exp(2*d(para)/v(para))*pnorm(-(1+d(para)*time)/sqrt(v(para)*time))))
        }



        p<-rep(0,(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])))
        time<-Y[,1]
        failure<-Y[,2]

        est<-nlm(logf, p, hessian = TRUE)

        names(est$estimate) <-c(paste("lny0:",dimnames(x_lny0)[[2]]),paste("  mu:",dimnames(x_mu)[[2]]))


	loglik = (-1)*est$minimum
	fit<-list(coefficients  = est$estimate,
		    var    = solve(est$hessian),
		    loglik = loglik,
                    AIC    = (-2)*loglik+2*(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])),
		    iter   = est$iterations,
                    call   = cl,
		    mf     = mf,
                    lny0   = dimnames(x_lny0)[[2]],
                    mu     = dimnames(x_mu)[[2]])
        class(fit) <- 'threg'
        fit


}


