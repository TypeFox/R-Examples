##' Fits Clayton-Oakes or bivariate Plackett (OR) models for binary data 
##' using marginals that are on logistic form. 
##' If clusters contain more than two times, the algoritm uses a compososite likelihood
##' based on all pairwise bivariate models.
##'
##' The pairwise pairwise odds ratio model provides an alternative to the alternating logistic
##' regression (ALR).
##'
##' The reported standard errors are based on a cluster corrected score equations from the 
##' pairwise likelihoods assuming that the marginals are known. This gives correct standard errors
##' in the case of the Plackett distribution (OR model for dependence), but incorrect standard
##' errors for the Clayton-Oakes types model.
##'
##' @export
##' @aliases binomial.twostage binomial.twostage.time
##' @references
##' Two-stage binomial modelling 
##' @examples
##' data(twinstut)
##' twinstut0 <- subset(twinstut, tvparnr<2300000)
##' twinstut <- twinstut0
##' theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
##' margbin <- glm(stutter~factor(sex)+age,data=twinstut,family=binomial())
##' bin <- binomial.twostage(margbin,data=twinstut,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bin)
##' 
##' twinstut$cage <- scale(twinstut$age)
##' theta.des <- model.matrix( ~-1+factor(zyg)+cage,data=twinstut)
##' bina <- binomial.twostage(margbin,data=twinstut,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bina)
##' 
##' theta.des <- model.matrix( ~-1+factor(zyg)+factor(zyg)*cage,data=twinstut)
##' bina <- binomial.twostage(margbin,data=twinstut,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bina)
##' 
##' twinstut$binstut <- (twinstut$stutter=="yes")*1
##' ## refers to zygosity of first subject in eash pair : zyg1
##' ## could also use zyg2 (since zyg2=zyg1 within twinpair's))
##' out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut,
##'                           response="binstut",id="tvparnr",
##' 	             	      theta.formula=~-1+factor(zyg1),
##'                           score.method="fisher.scoring")
##' summary(out)
##' 
##' ## refers to zygosity of first subject in eash pair : zyg1
##' ## could also use zyg2 (since zyg2=zyg1 within twinpair's))
##' desfs<-function(x,num1="zyg1",num2="zyg2")
##'     c(x[num1]=="dz",x[num1]=="mz",x[num1]=="os")*1
##'
##' out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
##'       data=twinstut,response="binstut",id="tvparnr",
##'       score.method="fisher.scoring",theta.formula=desfs,desnames=c("mz","dz","os"))
##' summary(out3)
##'
##' @keywords binomial regression 
##' @author Thomas Scheike
##' @export
##' @param margbin Marginal binomial model 
##' @param data data frame
##' @param score.method Scoring method
##' @param Nit Number of iterations
##' @param detail Detail
##' @param clusters Cluster variable
##' @param silent Debug information
##' @param weights Weights for log-likelihood, can be used for each type of outcome in 2x2 tables. 
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.des Variance component design
##' @param var.link Link function for variance 
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size
##' @param notaylor Taylor expansion
##' @param model model
##' @param marginal.p vector of marginal probabilities 
##' @param strata strata for fitting: considers only pairs where both are from same strata 
##' @param max.clust max clusters
##' @param se.clusters clusters for iid decomposition for roubst standard errors
##' @param numDeriv uses Fisher scoring aprox of second derivative if 0, otherwise numerical derivatives 
binomial.twostage <- function(margbin,data=sys.parent(),score.method="nlminb",
                              Nit=60,detail=0,clusters=NULL,silent=1,weights=NULL,
                              control=list(),theta=NULL,theta.des=NULL,var.link=1,iid=1,
                              step=0.5,notaylor=1,model="plackett",marginal.p=NULL,strata=NULL,
                              max.clust=NULL,se.clusters=NULL,numDeriv=0)
{ ## {{{
    ## {{{ seting up design and variables
    rate.sim <- 1; sym=1; 
    if (model=="clayton.oakes") dep.model <- 1 else if (model=="plackett") dep.model <- 2 else stop("Model must by either clayton.oakes or plackett \n"); 
    antpers <- NROW(data); 

### marginal prediction and binomial response, two types of calls ## {{{
    if (class(margbin)[1]=="glm") {
        ps <- predict(margbin,type="response")
        cause <- margbin$y
    }
    else if (class(margbin)[1]=="formula") {
        margbin <- glm(margbin,data=data,family=binomial())
        ps <- predict(margbin,type="response")
        cause <- margbin$y
    }  else if (is.null(marginal.p))
        stop("without marginal model, marginal p's must be given\n"); 

    if (!is.null(marginal.p)) {
        if (length(margbin)!=antpers) 
            stop("with marginal margbin is reseponse \n")
        else cause <- margbin
        if (length(marginal.p)!=antpers) 
            stop("length same as data dimension  \n")
        else ps <- marginal.p
    }
    ## }}}

    notaylor <- 1
    if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
    if (is.null(strata)==TRUE) strata<- rep(1,antpers); 
    if (length(strata)!=antpers) stop("Strata must have length equal to number of data points \n"); 

    out.clust <- cluster.index(clusters);  
    clusters <- out.clust$clusters
    maxclust <- out.clust$maxclust 
    antclust <- out.clust$antclust
    clusterindex <- out.clust$idclust
    clustsize <- out.clust$cluster.size
    call.secluster <- se.clusters

    if (is.null(se.clusters)) { se.clusters <- clusters; 
                                antiid <- nrow(clusterindex);} else  {
                                    iids <-  unique(se.clusters); 
                                    antiid <- length(iids); 
                                    if (is.numeric(se.clusters)) se.clusters <-  fast.approx(iids,se.clusters)-1
                                    else se.clusters <- as.integer(factor(se.clusters, labels = seq(antiid)))-1
                                }
    if (length(se.clusters)!=length(clusters)) stop("Length of seclusters and clusters must be same\n"); 

    if ((!is.null(max.clust))) if (max.clust< antiid) {
        coarse.clust <- TRUE
	qq <- unique(quantile(se.clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(se.clusters, breaks = qq, include.lowest = TRUE)    
	se.clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(se.clusters))
	maxclust <- max.clust    
	antiid <- max.clusters
    }                                                         

    ratesim<-rate.sim; 
    if (is.null(theta.des)==TRUE) ptheta<-1; 
    if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) else
    theta.des<-as.matrix(theta.des); 
    ptheta<-ncol(theta.des); 
    if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");

    if (is.null(theta)==TRUE) {
        if (var.link==1) theta<- rep(-0.7,ptheta);  
        if (var.link==0) theta<- rep(exp(-0.7),ptheta);   
    }       
    if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
    theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

    if (maxclust==1) stop("No clusters, maxclust size=1\n"); 
    ## }}}

    loglike <- function(par) 
        { ## {{{
            Xtheta <- theta.des %*% matrix(c(par),ptheta,1); 
            DXtheta <- array(0,c(1,1,1));

###      dyn.load("twostage.so")
            ptrunc <- rep(1,antpers); 

            outl<-.Call("twostageloglikebin", ## {{{
                        icause=cause,ipmargsurv=ps, 
                        itheta=c(par),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
                        icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
                        ivarlink=var.link,iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
                        itrunkp=ptrunc,istrata=strata,iseclusters=se.clusters,iantiid=antiid) 
            ## }}}

            if (detail==3) print(c(par,outl$loglike))

            attr(outl,"gradient") <-outl$score 
            if (oout==0) ret <- c(-1*outl$loglike) else if (oout==1) ret <- sum(outl$score^2) else if (oout==3) ret <- outl$score else ret <- outl
            return(ret)
        } ## }}}

    if (score.method=="optimize" && ptheta!=1) {cat("optimize only works for d==1, score.mehod set to nlminb \n"); score.method <- "nlminb";}

    theta.iid <- NULL
    logl <- NULL
    p <- theta
    if (score.method=="fisher.scoring") { ## {{{
        oout <- 2;  ### output control for obj
        if (Nit>0) 
            for (i in 1:Nit)
                {
                    out <- loglike(p)
                    hess <- out$Dscore
                    if (!is.na(sum(hess))) hessi <- lava::Inverse(out$Dscore) else hessi <- hess 
                    if (detail==1) {## {{{
                        cat(paste("Fisher-Scoring ===================: it=",i,"\n")); 
                        cat("theta:");print(c(p))
                        cat("loglike:");cat(c(out$loglike),"\n"); 
                        cat("score:");cat(c(out$score),"\n"); 
                        cat("hess:\n"); cat(out$Dscore,"\n"); 
                    }## }}}
                    delta <- hessi %*% out$score *step 
                    p <- p+delta* step
                    theta <- p; 
                    if (is.nan(sum(out$score))) break; 
                    if (sum(abs(out$score))<0.00001) break; 
                    if (max(theta)>20) break; 
                }
        if (!is.nan(sum(p))) { 
            if (detail==1 && iid==1) cat("iid decomposition\n"); 
            out <- loglike(p) 
            logl <- out$loglike
            score1 <- score <- out$score
            oout <- 0; 
            hess1 <- hess <- - out$Dscore 
            if (iid==1) theta.iid <- out$theta.iid
        }
        if (numDeriv==1) {
            oout <- 3
            hess <- numDeriv::jacobian(loglike,p)
        }
        if (detail==1 & Nit==0) {## {{{
            cat(paste("Fisher-Scoring ===================: final","\n")); 
            cat("theta:");print(c(p))
            cat("loglike:");cat(c(out$loglike),"\n"); 
            cat("score:");cat(c(out$score),"\n"); 
            cat("hess:\n"); cat(out$Dscore,"\n"); 
        }## }}}
        if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- diag(nrow(hess))
        ## }}}
    } else if (score.method=="nlminb") { ## {{{ nlminb optimizer
        oout <- 0; 
        tryCatch(opt <- nlminb(theta,loglike,control=control),error=function(x) NA)
        if (detail==1) print(opt); 

        if (detail==1 && iid==1) cat("iid decomposition\n"); 
        oout <- 2
        theta <- opt$par
        out <- loglike(opt$par)
        logl <- out$loglike
        score1 <- score <- out$score
        hess1 <- hess <- - out$Dscore
        if (iid==1) theta.iid <- out$theta.iid
        if (numDeriv==1) {
            oout <- 3; 
            p <- theta
            hess <- numDeriv::jacobian(loglike,theta)
        }
        hessi <- lava::Inverse(hess); 
        ## }}}
    } else if (score.method=="optimize" && ptheta==1) { ## {{{  optimizer
        oout <- 0; 
        if (var.link==1) {mino <- -20; maxo <- 10;} else {mino <- 0.001; maxo <- 100;}
        tryCatch(opt <- optimize(loglike,c(mino,maxo)));
        if (detail==1) print(opt); 

        opt$par <- opt$minimum
        theta <- opt$par
        if (detail==1 && iid==1) cat("iid decomposition\n"); 
        oout <- 2
        out <- loglike(opt$par)
        logl <- out$loglike
        score1 <- score <- out$score
        hess1 <- hess <- - out$Dscore
        if (iid==1) theta.iid <- out$theta.iid
        if (numDeriv==1) {
            oout <- 3; 
            p <- opt$par
            hess <-  numDeriv::jacobian(loglike,p)
        }
        hessi <- lava::Inverse(hess); 
        ## }}}
    } else if (score.method=="nlm") { ## {{{ nlm optimizer
        iid <- 0; oout <- 0; 
        tryCatch(opt <- nlm(loglike,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
        iid <- 1; 
        hess <-  opt$hessian
        score <- opt$gradient
        if (detail==1) print(opt); 
        hessi <-  lava::Inverse(hess); 
        theta <- opt$estimate
        if (detail==1 && iid==1) cat("iid decomposition\n"); 
        oout <- 2
        out <- loglike(opt$estimate)
        logl <- out$loglike
        score1 <- out$score
        hess1 <- out$Dscore
        if (iid==1) theta.iid <- out$theta.iid
        ## }}}
    }  else stop("score.methods = optimize(dim=1) nlm nlminb fisher.scoring\n"); 


    ## {{{ handling output
    robvar.theta <- NULL
    var.theta <- hessi
    if (iid==1) {
        theta.iid <- out$theta.iid %*% hessi
        robvar.theta  <- (t(theta.iid) %*% theta.iid) 
    }
###  if (iid==1) var.theta <- robvar.theta else var.theta <- -hessi
    if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:ptheta,sep="")
    theta <- matrix(theta,ptheta,1)
    if (length(thetanames)==nrow(theta)) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }
    ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,model=model,robvar.theta=robvar.theta,
               theta.iid=theta.iid,thetanames=thetanames,loglike=-logl,score1=score1,Dscore=out$Dscore,margsurv=ps); 
    class(ud)<-"twostage" 
    attr(ud, "Formula") <- formula
    attr(ud, "Clusters") <- clusters
    attr(ud,"sym")<-sym; 
    attr(ud,"var.link")<-var.link; 
    attr(ud,"antpers")<-antpers; 
    attr(ud,"antclust")<-antclust; 
    attr(ud, "Type") <- model
    attr(ud, "response") <- "binomial"
    return(ud);
    ## }}}

} ## }}}

##' @export
binomial.twostage.time <- function(formula,data,id,...,silent=1,fix.censweights=1,
breaks=Inf,pairsonly=TRUE,fix.marg=NULL,cens.formula,cens.model="aalen",weights="w") {
   ## {{{ 
    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

    censtime <- model.extract(M, "response")
    status <- censtime[,2]
    time <- censtime[,1]
    timevar <- colnames(censtime)[1]
    outcome <- as.character(terms(formula)[[2]])    
    if (is.null(breaks)) breaks <-  quantile(time,c(0.25,0.5,0.75,1))

    outcome0 <- paste(outcome,"_dummy")
    res <- list()
    logor <- cif <- conc <- c()
    k <- 0
    for (tau in rev(breaks)) {
        if ((length(breaks)>1) & (silent==0)) message(tau)
        ### construct min(T_i,tau) or T_i and related censoring variable, 
        ### thus G_c(min(T_i,tau)) or G_c(T_i) as weights
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) {
        data0 <- data
        time0 <- time 
        status0 <- status 
        }
        cond0 <- (time>tau)
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) status0[cond0 & status==1] <- 3 
        if (fix.censweights==0) status0[cond0 & status==1] <- 3 
        data0[,outcome] <- data[,outcome]
        data0[cond0,outcome] <- FALSE
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) time0[cond0] <- tau
###        if (fix.censweights==0 ) time0[cond0] <- tau
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) {
		data0$S <- survival::Surv(time0,status0==1)        
	}
	if ((fix.censweights==1 & k==0) | (fix.censweights==0))
        dataw <- ipw(update(cens.formula,S~.),data=data0,cens.model=cens.model,
		     obsonly=TRUE)
        if ((fix.censweights==1)) 
		dataw[,outcome] <- (dataw[,outcome])*(dataw[,timevar]<tau)
	marg.bin <- glm(formula,data=dataw,weights=1/dataw[,weights],family="quasibinomial")
        pudz <- predict(marg.bin,newdata=dataw,type="response")
	dataw$pudz <- pudz
	datawdob <- fast.reshape(dataw,id=id)
        datawdob$minw <- pmin(datawdob$w1,datawdob$w2)
        dataw2 <- fast.reshape(datawdob)
	### removes second row of singletons 
	dataw2  <- subset(dataw2,!is.na(dataw2$minw)) 
	k <- k+1
	if (!is.null(fix.marg)) dataw2$pudz <- fix.marg[k]
        suppressWarnings( b <- binomial.twostage(dataw2[,outcome],data=dataw2,clusters=dataw2[,id],marginal.p=dataw2$pudz,weights=1/dataw2$minw,...))
        theta0 <- b$theta[1,1]
        prev <- prev0 <- exp(coef(marg.bin)[1])/(1+exp(coef(marg.bin)[1]))
	if (!is.null(fix.marg)) prev <- fix.marg[k]
        concordance <- plack.cif2(prev,prev,theta0)
	conc <- c(conc,concordance)
	cif <- c(cif,prev0)
	logor <- rbind(logor,coef(b))
###     res <- c(res,list(coef(b),concordance=concordance,cif=prev0))
    }
###    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,concordance=rev(conc),cif=rev(cif),
		time=breaks,call=m,type="time",logor=logor[k:1,])
###	coef=lapply(res,function(x) x$all),
###    class(res) <- ""
    return(res)    
} ## }}} 


##' Fits two-stage binomial for describing depdendence in binomial data
##' using marginals that are on logistic form using the binomial.twostage funcion, but
##' call is different and easier and the data manipulation is build into the function.
##' Useful in particular for family design data. 
##'
##' If clusters contain more than two times, the algoritm uses a compososite likelihood
##' based on the pairwise bivariate models.
##'
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. This gives correct standard errors
##' in the case of the plackett distribution (OR model for dependence), but incorrect for
##' the clayton-oakes types model. The OR model is often known as the ALR model, but our fitting
##' procedures gives correct standard errors and is quite a bit quicker.
##'
##' @examples
##' data(twinstut)
##' twinstut0 <- subset(twinstut, tvparnr<2300000)
##' twinstut <- twinstut0
##' theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut0)
##' margbin <- glm(stutter~factor(sex)+age,data=twinstut0,family=binomial())
##' bin <- binomial.twostage(margbin,data=twinstut0,
##' 		         clusters=twinstut0$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bin)
##' 
##' twinstut0$cage <- scale(twinstut0$age)
##' theta.des <- model.matrix( ~-1+factor(zyg)+cage,data=twinstut0)
##' bina <- binomial.twostage(margbin,data=twinstut0,
##' 		         clusters=twinstut0$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bina)
##' 
##' theta.des <- model.matrix( ~-1+factor(zyg)+factor(zyg)*cage,data=twinstut0)
##' bina <- binomial.twostage(margbin,data=twinstut0,
##' 		         clusters=twinstut0$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bina)
##' 
##' twinstut0$binstut <- (twinstut0$stutter=="yes")*1
##' out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut0,
##'                               response="binstut",id="tvparnr",
##' 			      theta.formula=~-1+factor(zyg1),
##'                               score.method="fisher.scoring")
##' summary(out)
##' 
##' ## refers to zygosity of first subject in eash pair : zyg1
##' ## could also use zyg2 (since zyg2=zyg1 within twinpair's))
##' desfs <- function(x,num1="zyg1",namesdes=c("mz","dz","os"))
##'     c(x[num1]=="mz",x[num1]=="dz",x[num1]=="os")*1
##' 
##' out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
##'                                data=twinstut0, response="binstut",id="tvparnr",
##'                                score.method="fisher.scoring",
##'                                theta.formula=desfs,desnames=c("mz","dz","os"))
##' summary(out3)
##'
##' \donttest{ ## Reduce Ex.Timings
##' n <- 10000
##' set.seed(100)
##' dd <- simBinFam(n,beta=0.3) 
##' binfam <- fast.reshape(dd,varying=c("age","x","y"))
##' ## mother, father, children  (ordered)
##' head(binfam)
##'
##' ########### ########### ########### ########### ########### ###########
##' ####  simple analyses of binomial family data 
##' ########### ########### ########### ########### ########### ###########
##' desfs <- function(x,num1="num1",num2="num2")
##' {  
##'      pp <- 1*(((x[num1]=="m")*(x[num2]=="f"))|(x[num1]=="f")*(x[num2]=="m"))
##'      pc <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
##'      cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
##'      c(pp,pc,cc)
##' } 
##'
##' ud <- easy.binomial.twostage(y~+1,data=binfam,
##'      response="y",id="id",
##'      score.method="fisher.scoring",deshelp=0,
##'      theta.formula=desfs,desnames=c("pp","pc","cc"))
##' summary(ud)
##'
##' udx <- easy.binomial.twostage(y~+x,data=binfam,
##'      response="y",id="id",
##'      score.method="fisher.scoring",
##'      theta.formula=desfs,desnames=c("pp","pc","cc"))
##' summary(udx)
##'
##' ########### ########### ########### ########### ########### ###########
##' ####  now allowing parent child POR to be different for mother and father 
##' ########### ########### ########### ########### ########### ###########
##'
##' desfsi <- function(x,num1="num1",num2="num2")
##' { 
##'     pp <- (x[num1]=="m")*(x[num2]=="f")*1
##'     mc <- (x[num1]=="m")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     fc <- (x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     c(pp,mc,fc,cc)
##' }
##'
##' udi <- easy.binomial.twostage(y~+1,data=binfam,
##'      response="y",id="id",
##'      score.method="fisher.scoring",
##'      theta.formula=desfsi,desnames=c("pp","mother-child","father-child","cc"))
##' summary(udi)
##'
##' ##now looking to see if interactions with age or age influences marginal models 
##' ##converting factors to numeric to make all involved covariates numeric
##' ##to use desfai2 rather then desfai that works on binfam 
##'
##' nbinfam <- binfam
##' nbinfam$num <- as.numeric(binfam$num)
##' head(nbinfam)
##'
##' desfsai <- function(x,num1="num1",num2="num2")
##' { 
##'     pp <- (x[num1]=="m")*(x[num2]=="f")*1
##' ### av age for pp=1 i.e parent pairs
##'     agepp <- ((as.numeric(x["age1"])+as.numeric(x["age2"]))/2-30)*pp 
##'     mc <- (x[num1]=="m")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     fc <- (x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     agecc <- ((as.numeric(x["age1"])+as.numeric(x["age2"]))/2-12)*cc 
##'     c(pp,agepp,mc,fc,cc,agecc)
##' } 
##'
##' desfsai2 <- function(x,num1="num1",num2="num2")
##' { 
##'     pp <- (x[num1]==1)*(x[num2]==2)*1
##'     agepp <- (((x["age1"]+x["age2"]))/2-30)*pp ### av age for pp=1 i.e parent pairs
##'     mc <- (x[num1]==1)*(x[num2]==3 | x[num2]==4)*1
##'     fc <- (x[num1]==2)*(x[num2]==3 | x[num2]==4)*1
##'     cc <- (x[num1]==3)*(x[num2]==3 | x[num2]==4)*1
##'     agecc <- ((x["age1"]+x["age2"])/2-12)*cc ### av age for children 
##'     c(pp,agepp,mc,fc,cc,agecc)
##' } 
##'
##' udxai2 <- easy.binomial.twostage(y~+x+age,data=binfam,
##'      response="y",id="id",
##'      score.method="fisher.scoring",deshelp=0,detail=0,
##'      theta.formula=desfsai,
##'      desnames=c("pp","pp-age","mother-child","father-child","cc","cc-age"))
##' summary(udxai2)
##' }
##' @keywords binomial regression 
##' @export
##' @param margbin Marginal binomial model 
##' @param data data frame
##' @param response name of response variable in data frame
##' @param id name of cluster variable in data frame
##' @param score.method Scoring method
##' @param Nit Number of iterations
##' @param detail Detail for more output for iterations 
##' @param silent Debug information
##' @param weights Weights for log-likelihood, can be used for each type of outcome in 2x2 tables. 
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.formula design for depedence, either formula or design function
##' @param desnames names for dependence parameters
##' @param deshelp if 1 then prints out some data sets that are used, on on which the design function operates
##' @param var.link Link function for variance 
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size
##' @param model model
##' @param marginal.p vector of marginal probabilities 
##' @param strata strata for fitting 
##' @param max.clust max clusters
##' @param se.clusters clusters for iid decomposition for roubst standard errors
easy.binomial.twostage <- function(margbin=NULL,data=sys.parent(),score.method="nlminb",
                                   response="response",id="id",
                                   Nit=60,detail=0, silent=1,weights=NULL, control=list(),
                                   theta=NULL,theta.formula=NULL,desnames=NULL,deshelp=0,var.link=1,iid=1,
                                   step=0.5,model="plackett",marginal.p=NULL,strata=NULL,max.clust=NULL,se.clusters=NULL)
{ ## {{{
    if (class(margbin)[1]=="glm") ps <- predict(margbin,type="response") 
    else if (class(margbin)=="formula") {
        margbin <- glm(margbin,data=data,family=binomial())
        ps <- predict(margbin,type="response")
    }  else if (is.null(marginal.p)) stop("without marginal model, marginal p's must be given\n"); 

    if (!is.null(marginal.p)) ps <- marginal.p

    data <- cbind(data,ps)

### make all pairs in the families,
    fam <- familycluster.index(data[,id])
    data.fam <- data[fam$familypairindex,]
    data.fam$subfam <- fam$subfamilyindex

### make dependency design using wide format for all pairs 
    data.fam.clust <- fast.reshape(data.fam,id="subfam")
    if (is.function(theta.formula)) {
        desfunction <- compiler::cmpfun(theta.formula)
	if (deshelp==1){
            cat("These names appear in wide version of pairs for dependence \n")
            cat("design function must be defined in terms of these: \n")
            cat(names(data.fam.clust)); cat("\n")
            cat("Here is head of wide version with pairs\n")
            print(head(data.fam.clust)); cat("\n")
	}
###	des.theta <- Reduce("rbind",lapply(seq(nrow(data.fam.clust)),function(i) unlist(desfunction(data.fam.clust[i,] ))))
        des.theta <- t(apply(data.fam.clust,1, function(x) desfunction(x)))
	colnames(des.theta) <- desnames
	desnames <- desnames
    } else {
        if (is.null(theta.formula)) theta.formula <- ~+1
        des.theta <- model.matrix(theta.formula,data=data.fam.clust)
        desnames <- colnames(des.theta); 
    }
    data.fam.clust <- cbind(data.fam.clust,des.theta)
    if (deshelp==1) {
        cat("These names appear in wide version of pairs for dependence \n")
        print(head(data.fam.clust))
    }

### back to long format keeping only needed variables
    data.fam <- fast.reshape(data.fam.clust,varying=c(response,id,"ps"))
    if (deshelp==1) {
	cat("Back to long format for binomial.twostage (head)\n"); 
        print(head(data.fam)); 
	cat("\n")
	cat(paste("binomial.twostage, called with reponse",response,"\n")); 
	cat(paste("cluster=",id,",  subcluster (pairs)=subfam \n")); 
	cat(paste("design variables =")); 
	cat(desnames)
	cat("\n")
    } 

    out <- binomial.twostage(data.fam[,response],data=data.fam,
                             clusters=data.fam$subfam,
                             theta.des=data.fam[,desnames],
                             detail=detail, score.method=score.method, Nit=Nit,step=step,
                             iid=iid,theta=theta, var.link=var.link,model=model, 
                             max.clust=max.clust,
                             marginal.p=data.fam[,"ps"], se.clusters=data.fam[,id])
    return(out)
} ## }}}

##' @export
simBinPlack <- function(n,beta=0.3,theta=1,...) { ## {{{ 
    x1 <- rbinom(n,1,0.5)
    x2 <- rbinom(n,1,0.5)
###
    p1 <- exp(0.5+x1*beta)
    p2 <- exp(0.5+x2*beta)
    p1 <- p1/(1+p1)
    p2 <- p2/(1+p2)
###
    p11 <- plack.cif2(p1,p2,theta)
    p10 <- p1-p11
    p01 <- p2-p11
    p00 <- 1- p10-p01-p11
###
    y1 <- rbinom(n,1,p1)
    y2 <- (y1==1)*rbinom(n,1,p11/p1)+(y1==0)*rbinom(n,1,p01/(1-p1))
    list(x1=x1,x2=x2,y1=y1,y2=y2,id=1:n)
} ## }}} 

##' @export 
simBinFam <- function(n,beta=0.0,rhopp=0.1,rhomb=0.7,rhofb=0.1,rhobb=0.7) { ## {{{ 
    xc <- runif(n)*0.5
    xm <- rbinom(n,1,0.5+xc); 
    xf <- rbinom(n,1,0.5+xc); 
    xb1 <- rbinom(n,1,0.3+xc); 
    xb2 <- rbinom(n,1,0.3+xc); 
###
    rn <- matrix(rnorm(n*4),n,4)
    corm <- matrix( c(1,rhopp,rhomb,rhomb, rhopp,1,rhofb,rhofb, rhomb,rhofb,1,rhobb, rhomb,rhofb,rhobb,1),4,4)
    rnn <- t( corm %*% t(rn))
    zm <- exp(rnn[,1]); zf <- exp(rnn[,2]); zb1 <- exp(rnn[,3]); zb2 <- exp(rnn[,4]); 
    pm <- exp(0.5+xm*beta+zm)
    pf <- exp(0.5+xf*beta+zf)
    pf <- pf/(1+pf)
    pm <- pm/(1+pm)
    pb1 <- exp(0.5+xb1*beta+zb1)
    pb1 <- pb1/(1+pb1)
    pb2 <- exp(0.5+xb2*beta+zb2)
    pb2 <- pb2/(1+pb2)
    ym <- rbinom(n,1,pm)
    yf <- rbinom(n,1,pf)
    yb1 <- rbinom(n,1,pb1)
    yb2 <- rbinom(n,1,pb2)
                                        #
    agem <- 20+runif(n)*10
    ageb1 <- 5+runif(n)*10
    data.frame(agem=agem,agef=agem+3+rnorm(n)*2,
               ageb1=ageb1,ageb2=ageb1+1+runif(n)*3,xm=xm,xf=xf,xb1=xb1,xb2=xb2,ym=ym,yf=yf,yb1=yb1,yb2=yb2,id=1:n)
} ## }}} 

##' @export
simBinFam2 <- function(n,beta=0.0,lam1=1,lam2=1,...) { ## {{{ 
    x1 <- rbinom(n,1,0.5); x2 <- rbinom(n,1,0.5); 
    x3 <- rbinom(n,1,0.5); x4 <- rbinom(n,1,0.5); 
###
    zf <- rgamma(n,shape=lam1); zb <- rgamma(n,shape=lam2); 
    pm <- exp(0.5+x1*beta+zf)
    pf <- exp(0.5+x2*beta+zf)
    pf <- pf/(1+pf)
    pm <- pm/(1+pm)
    pb1 <- exp(0.5+x1*beta+zf+zb)
    pb1 <- pb1/(1+pb1)
    ym <- rbinom(n,1,pm)
    yf <- rbinom(n,1,pf)
    yb1 <- rbinom(n,1,pb1)
    yb2 <- rbinom(n,1,pb1)
                                        #
    data.frame(x1=x1,x2=x2,ym=ym,yf=yf,yb1=yb1,yb2=yb2,id=1:n)
} ## }}} 


#####' @S3method summary twostage
###summary.twostage <-function (object,digits = 3,...) { ## {{{
###  if (!(inherits(object,"twostage"))) stop("Must be a Two-Stage object")
###  
###  var.link<-attr(object,"var.link");
###  if (object$model=="plackett") cat("Dependence parameter for Plackett model \n"); 
###  if (object$model=="clayton.oakes") cat("Dependence parameter for Clayton-Oakes model \n"); 
###
###  if (sum(abs(object$score)>0.0001) ) {
###	  cat("    Variance parameters did not converge, allow more iterations.\n"); 
###	  cat(paste("    Score:",object$score,"  \n")); 
###  }
###
###  coefs <- coef.twostage(object,...);
###
###  res <- list(estimates=coefs, type=attr(object,"Type"))
###  class(res) <- "summary.twostage"
###  res
###} ## }}}
###
#####' @S3method coef twostage
###coef.twostage <- function(object,var.link=NULL,...)
###{ ## {{{
###  theta <- object$theta
###  if (is.null(var.link))
###     if (attr(object,"var.link")==1) vlink <- 1 else vlink <- 0
###     else vlink <- var.link
###  se<-diag(object$var.theta)^0.5
###  res <- cbind(theta, se )
###  wald <- theta/se
###  waldp <- (1 - pnorm(abs(wald))) * 2
###  library(numDeriv)
###  if (object$model=="plackett") {
###  spearman <- alpha2spear(theta,link=vlink)
###  Dspear <- numDeriv::jacobian(alpha2spear,theta,link=vlink) 
###  var.spearman <- Dspear %*% object$var.theta %*%  Dspear
###  se.spearman <- diag(var.spearman)^.5
###  res <- as.matrix(cbind(res, wald, waldp,spearman,se.spearman))
###  if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Spearman Corr.","SE")
###  else colnames(res) <- c("Coef.", "SE","z", "P-val","Spearman Corr.","SE")
###  if (!is.null(object$thetanames)) rownames(res)<-object$thetanames
###  }
###  if (object$model=="clayton.oakes") {
###  kendall <- alpha2kendall(theta,link=vlink)
###  Dken <- numDeriv::jacobian(alpha2kendall,theta,link=vlink) 
###  var.kendall<- Dken %*% object$var.theta %*%  Dken
###  se.kendall <- diag(var.kendall)^.5
###  res <- as.matrix(cbind(res, wald, waldp,kendall,se.kendall))
###  if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Kendall tau","SE")
###  else colnames(res) <- c("Coef.", "SE","z", "P-val","Kendall tau","SE")
###  if (!is.null(object$thetanames)) rownames(res)<-object$thetanames
###  }
###
###  return(res)
###} ## }}}
###
#####' @export
###alpha2spear <- function(theta,link=1) {
###   if (link==1) theta <- exp(theta)
###   if (theta!=1) return( (theta+1)/(theta-1) -2* theta* log(theta)/ (theta-1)^2)
###   else return(0)
###}
###
#####' @export
###alpha2kendall <- function(theta,link=0) { 
###   if (link==1) theta <- exp(theta)
###   return(1/(1+2/theta)) 
###}
###
#####' @S3method print twostage
###print.twostage<-function(x,digits=3,...)
###{ ## {{{
###  print(x$call); 
###  cat("\n")
###  print(summary(x)); 
###} ## }}}
###
#####' @S3method plot twostage
###plot.twostage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
###		level=0.05, 
###		start.time=0,stop.time=0,add.to.plot=FALSE,mains=TRUE,
###                xlab="Time",ylab ="Cumulative regression function",...) 
###{ ## {{{
###  if (!(inherits(x, 'two.stage'))) stop("Must be a Two-Stage object")
###  object <- x; rm(x);  
### 
###  B<-object$cum; V<-object$var.cum; p<-dim(B)[[2]]; 
###  if (robust>=1) V<-object$robvar.cum; 
###
###  if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1
###  if (stop.time==0) stop.time<-max(B[,1]);
###
###  med<-B[,1]<=stop.time & B[,1]>=start.time
###  B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
###  V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
###  Vrob<-object$robvar.cum; 
###  Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 
###
###  c.alpha<- qnorm(1-level/2)
###  for (v in comp) { 
###    c.alpha<- qnorm(1-level/2)
###    est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
###    if (add.to.plot==FALSE) 
###      {
###        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
###        if (mains==TRUE) title(main=colnames(B)[v]); }
###    else lines(B[,1],est,type="s"); 
###    if (pointwise.ci>=1) {
###      lines(B[,1],ul,lty=pointwise.ci,type="s");
###      lines(B[,1],nl,lty=pointwise.ci,type="s"); }
###    if (robust>=1) {
###      lines(B[,1],ul,lty=robust,type="s"); 
###      lines(B[,1],nl,lty=robust,type="s"); }
###    abline(h=0); 
###  }
###}  ## }}}
###
#####' @S3method predict twostage
###predict.twostage <- function(object,X=NULL,Z=NULL,times=NULL,times2=NULL,theta.des=NULL,diag=TRUE,...)
###{ ## {{{
###time.coef <- data.frame(object$cum)
###if (!is.null(times)) {
###cum <- Cpred(object$cum,times);
###cum2 <- Cpred(object$cum,times);
###} else { cum <- object$cum; cum2 <- object$cum }
###if (!is.null(times2)) cum2 <- Cpred(object$cum,times2);
###
###if (is.null(X)) X <- 1;
###if (is.null(X) & (!is.null(Z))) { Z <- as.matrix(Z);  X <- matrix(1,nrow(Z),1)}
###if (is.null(Z) & (!is.null(X)))  {X <- as.matrix(X);  Z <- matrix(0,nrow(X),1); gamma <- 0}
###
###if (diag==FALSE) {
###   time.part <-  X %*% t(cum[,-1]) 
###   time.part2 <-  X %*% t(cum2[,-1]) 
###   if (!is.null(object$gamma)) { RR <- exp( Z %*% gamma ); 
###       cumhaz <- t( t(time.part) * RR ); cumhaz2 <- t( t(time.part2) * RR )}
###	    else { cumhaz <- time.part;  cumhaz2 <- time.part2;   }
###} else { 
###	time.part <-  apply(as.matrix(X*cum[,-1]),1,sum) 
###	time.part2 <-  apply(as.matrix(X*cum2[,-1]),1,sum) 
###}
###
###if (!is.null(object$gamma)) {
###	RR<- exp(Z%*%gamma); 
###	cumhaz <- t( t(time.part) * RR );  
###	cumhaz2 <- t( t(time.part2) * RR )} else {
###		cumhaz <- time.part;  cumhaz2 <- time.part2; 
###} 
###S1 <- exp(-cumhaz); S2 <- exp(-cumhaz2)
###
###if (attr(object,"var.link")==1) theta  <- exp(object$theta) else theta <- object$theta
###if (!is.null(theta.des)) theta <- c(theta.des %*% object$theta)
###
###if (diag==FALSE) St1t2<- (outer(c(S1)^{-(1/theta)},c(S2)^{-(1/theta)},FUN="+") - 1)^(-(theta)) else 
###St1t2<- ((S1^{-(1/theta)}+S2^{-(1/theta)})-1)^(-(theta))
###
###out=list(St1t2=St1t2,S1=S1,S2=S2,times=times,times2=times2,theta=theta)
###return(out)
###} ## }}}
###
#####' @export
###piecewise.twostage <- function(cut1,cut2,data=sys.parent(),timevar="time",status="status",id="id",covars=NULL,num=NULL,
###            score.method="optimize",Nit=100,detail=0,clusters=NULL,silent=1,weights=NULL,
###            control=list(),theta=NULL,theta.des=NULL,var.link=1,iid=0,step=0.5,model="plackett",data.return=0)
###{ ## {{{
###
###ud <- list()
###if (missing(cut2)) cut2 <- cut1; 
###nc1 <- length(cut1); nc2 <- length(cut2)
###names1 <- names2 <- c()
###theta.mat <- se.theta.mat <- cor.mat <- score.mat <- se.cor.mat <- matrix(0,nc1-1,nc2-1); 
###idi <- unique(data[,id]); 
###if (iid==1) theta.iid <- matrix(0,length(idi),(nc1-1)*(nc2-1)) else theta.iid <- NULL
###
###k <- 0; 
###for (i1 in 2:nc1)
###for (i2 in 2:nc2)
###{
###k <-(i1-2)*(nc2-1)+(i2-1)
###if (silent<=0) cat(paste("Data-set ",k,"out of ",(nc1-1)*(nc2-1)),"\n"); 
### datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data,timevar=timevar,
###			status=status,id=id,covars=covars,num=num,silent=silent) 
###if (silent<=-1) print(summary(datalr)); 
### boxlr <- list(left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]))
###### marg1 <- aalen(Surv(datalr$left,datalr[,timevar],datalr[,status])~+1,data=datalr,n.sim=0,max.clust=NULL,robust=0)
###datalr$tstime <- datalr[,timevar]
###datalr$tsstatus <- datalr[,status]
###datalr$tsid <- datalr[,id]
######
###f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(",num,")")))
######f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(num)")))
###marg1 <- aalen(f,data=datalr,n.sim=0,max.clust=NULL,robust=0)
###fitlr<-  twostage(marg1,data=datalr,clusters=datalr$tsid,model=model,score.method=score.method,
###              Nit=Nit,detail=detail,silent=silent,weights=weights,
###              control=control,theta=theta,theta.des=theta.des,var.link=var.link,iid=iid,step=step)
#######
###coef <- coef(fitlr)
###theta.mat[i1-1,i2-1] <- fitlr$theta
###se.theta.mat[i1-1,i2-1] <- fitlr$var.theta^.5
###cor.mat[i1-1,i2-1] <- coef[1,5]
###se.cor.mat[i1-1,i2-1] <- coef[1,6]
###score.mat[i1-1,i2-1] <- fitlr$score
###if (data.return==0) 
###ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr)
###if (data.return==1) 
###ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr,data=datalr)
###if (i2==2) names1 <- c(names1, paste(cut1[i1-1],"-",cut1[i1]))
###if (i1==2) names2 <- c(names2, paste(cut2[i2-1],"-",cut2[i2]))
###theta <- c(theta,fitlr$theta)
###
###if (iid==1) theta.iid[idi %in% unique(datalr$tsid),k] <-  fitlr$theta.iid 
###}
###
###var.thetal <- NULL
###if (iid==1)  var.thetal <- t(theta.iid) %*% theta.iid
###
###colnames(score.mat) <- colnames(cor.mat) <-  colnames(se.cor.mat)  <- colnames(se.theta.mat) <- colnames(theta.mat) <- names1; 
###rownames(score.mat) <- rownames(cor.mat) <-  rownames(se.cor.mat) <-  rownames(se.theta.mat) <- rownames(theta.mat) <- names2; 
###
###ud <- list(model.fits=ud,theta=theta.mat,var.theta=se.theta.mat^2,
###	   se.theta=se.theta.mat,thetal=theta,thetal.iid=theta.iid,var.thetal=var.thetal,model=model,
###	   cor=cor.mat,se.cor=se.cor.mat,score=score.mat); 
###class(ud)<-"pc.twostage" 
###attr(ud,"var.link")<-var.link; 
###attr(ud, "Type") <- model
###return(ud);
###} ## }}}
###
#####' @S3method summary pc.twostage
###summary.pc.twostage <- function(object,var.link=NULL,...)
###{ ## {{{
###  if (!(inherits(object,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
###  
###  res <- list(estimates=object$theta,se=object$se.theta,cor=object$cor,se.cor=object$se.cor,
###	      model=object$model,score=object$score)
###  class(res) <- "summary.pc.twostage"
###  attr(res,"var.link")<-attr(object,"var.link"); 
###  attr(res, "Type") <- object$model
###  res
###} ## }}}
###
#####' @S3method print pc.twostage
###print.pc.twostage <- function(x,var.link=NULL,...)
###{ ## {{{
###   if (!(inherits(x,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
###   print( summary(x,var.link=var.link,...))
###} ## }}}
###
#####' @S3method print summary.pc.twostage
###print.summary.pc.twostage <- function(x,var.link=NULL, digits=3,...)
###{ ## {{{
###  
###  if (is.null(var.link)) { if (attr(x,"var.link")==1) vlink <- 1 else vlink <- 0; } else vlink <- var.link
###  print(vlink)
###
###  if (x$model=="plackett") cat("Dependence parameter for Plackett model \n"); 
###  if (x$model=="clayton.oakes") cat("Dependence parameter for Clayton-Oakes model \n"); 
### 
###  if (max(x$score)>0.001) { cat("Score of log-likelihood for parameter estimates (too large?)\n"); print(x$score);cat("\n\n");}
###
###  if (vlink==1) cat("log-coefficient for dependence parameter (SE) \n")  else cat("Dependence parameter (SE) \n");
###  print(coefmat(x$estimate,x$se,digits=digits,...))
###  cat("\n") 
###
###  if (x$model=="plackett") {cat("Spearman Correlation (SE) \n");cor.type <- "Spearman Correlation"; }
###  if (x$model=="clayton.oakes") {cat("Kendall's tau (SE) \n"); cor.type <- "Kendall's tau";}
###
###  print(coefmat(x$cor,x$se.cor,digits,...))
###  cat("\n") 
###} ## }}}
###
#####' @export
###coefmat <- function(est,stderr,digits=3,...) { ## {{{
###  myest <- round(10^digits*(est))/10^digits;
###  myest <- paste(ifelse(myest<0,""," "),myest,sep="")
###  mysd <- round(10^digits*(stderr))/10^digits;  
###  res <- matrix(paste(format(myest)," (",format(mysd),")",sep=""),ncol=ncol(est))
###  dimnames(res) <- dimnames(est)
###  colnames(res) <- paste("",colnames(res))
###  noquote(res)
###} ## }}}
###
