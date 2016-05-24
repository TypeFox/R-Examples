frailtyHL <-
function(formula,data,weights,subset,na.action,RandDist="Normal",mord=0,dord=1,Maxiter=200,convergence=10^-6, varfixed=FALSE, varinit=c(0.1)){
    require(Matrix)
    require(numDeriv)
    require(survival)

    Call <- match.call()
    mc <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "na.action"), 
        names(Call), nomatch = 0)
    if (indx[1] == 0) 
        stop("A formula argument is required")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    temp$formula <- terms(subbars(formula), special)
    m <- eval(temp)
    Terms <- attr(m, "terms")
    Y <- model.extract(m, "response")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    temp$formula <- terms(formula, special)
    Terms <- temp[[2]]
    formula1<-paste(paste(Terms[[2]][[2]],Terms[[3]],sep="~")[[2]],paste(Terms[[3]])[3],sep="+")
    formula1<-formula(formula1)
    fr <- FrailtyFrames(mc, formula1, contrasts)
    namesX <- names(fr$fixef)
    namesX <- namesX[-1]
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formula1, fr, 0L, 0L)
    namesRE <- FL$namesRE
    y <- matrix(Y[,1], length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n<-nrow(x)
    p<-ncol(x)
    x1 <- x[1:n,2:p]
    x2 <- matrix(x1,n,p-1)
    x <- x2
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    del <-matrix(0,n,1)
    del[,1] <- censor<-Y[,2]
    SS <- FL$Subject
    res1<-FrailtyMakeData(y,x,del,z)
    y<-res1[1][[1]]
    x<-res1[2][[1]]
    del<-res1[3][[1]]
    z<-res1[4][[1]]
    Mi<-res1[5][[1]]
    idx2<-res1[6][[1]]
    t2<-res1[7][[1]]
    di<-res1[8][[1]]
    beta_h<-matrix(0,p,1)
    qcum <- cumsum(c(0, q))
    v_h<-matrix(0,qcum[nrand+1],1)
    alpha_h <- rep(0, nrand)
    varinit1 <- rep(0.1, nrand)
    length_init<-length(varinit)
    for (i in 1:length_init) varinit1[i]<-varinit[i]
    for (i in 1:nrand) alpha_h[i] <- varinit1[i]
    Max_iter<-Maxiter
    err<-1
    for ( i in 1:Max_iter) {
        if (err>=convergence) {
        if (RandDist=="Normal") res2<-PNfrailtyHL(x,z,y, del,Mi,idx2,t2, di, beta_h,v_h, alpha_h,mord,dord,varfixed=varfixed)
        if (RandDist=="Gamma") res2<-PGfrailtyHL(x,z,y, del,Mi,idx2,t2, di, beta_h,v_h, alpha_h,mord,dord,varfixed=varfixed)
        alpha_h<-res2[13][[1]]
        alpha_h1<-res2[14][[1]]
        beta_h<-res2[11][[1]]
        beta_h1<-res2[9][[1]]
        v_h<-res2[12][[1]]
        v_h1<-res2[10][[1]]
        temp4<-sum(abs(alpha_h-alpha_h1))+sum(abs(v_h-v_h1))+sum(abs(beta_h-beta_h1))
        err<-temp4
        alpha_h<-alpha_h1
        se_beta<-res2[20][[1]]
        print_i<-i
        print_err<-err
        }
    }
    names(print_i) <- "iteration : "
    print(print_i)
    names(print_err) <- "convergence : "
    print(print_err)
    if (err<convergence) print("converged") 
    if (err>convergence) print("did not converge")
###############################################################
############# print estimates ###########################
###############################################################
    result<-list(0)
    names(result)[1]<-"Model"
    sum_init<-sum(abs(varinit))
    if (varfixed==TRUE && sum_init<0.00001) {
         print("Results from the Cox model")
         result$Model<-"Cox model"
    } else {
    if (RandDist=="Gamma") {
         print("Results from the gamma frailty model")
         result$Model<-"gamma frailty model"
    }
    if (RandDist=="Normal") {
         print("Results from the log-normal frailty model")
         result$Model<-"log-normal frailty model"
    }
    }
    nevent<-sum(censor)
    print("Number of data : ")
    print(n)
    print("Number of event : ")
    print(nevent)
    print("Model for conditional hazard : ")
    result$formula<-formula
    print(formula)
    if (mord==0 && dord==1) {
         print("Method : HL(0,1)")   
         result$Method<-"HL(0,1)"
    }
    if (mord==0 && dord==2) {
         print("Method : HL(0,2)")     
         result$Method<-"HL(0,2)"
    }
    if (mord==1 && dord==1) {
         print("Method : HL(1,1)")   
         result$Method<-"HL(1,1)"
    }
    if (mord==1 && dord==2) {
         print("Method : HL(1,2)")
         result$Method<-"HL(1,2)"
    }
    print("Estimates from the mean model")
    z_beta<-beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff<-cbind(beta_h,se_beta,z_beta,pval)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value", "p-value")
    rownames(beta_coeff) <- namesX
    result$FixCoef<-beta_coeff
    print(beta_coeff,4)
###############################################################
############# se for lambda         ###########################
###############################################################
    if (RandDist=="Normal") res3<-PNFrailty_SE.h(res2,nrand,q,qcum,dord,varfixed=varfixed)
    if (RandDist=="Gamma") res3<-PGFrailty_SE.h(res2,nrand,q,qcum,dord,varfixed=varfixed)
    print("Estimates from the dispersion model")
    se_alpha_h<-res3[1][[1]]
    hlike<--2*res3[2][[1]]
    p1<--2*res3[3][[1]]
    p2<--2*res3[4][[1]]
    p3<--2*res3[5][[1]]
    p0<--2*res3[6][[1]]
    p4<-p3-(p1-p2)
    df1<-res3[7][[1]]
    if (varfixed==FALSE) z_lam<-alpha_h/se_alpha_h
    for (i in 1:nrand) {
       if (alpha_h[i]<0.00001) alpha_h[i]<-0
    }
    lam_coeff<-cbind(alpha_h,se_alpha_h)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
    print(lam_coeff,4)
    result$RandCoef<-lam_coeff
###############################################################
############# -2*Likelihoods         ##########################
###############################################################
    if (mord==0 && dord==1) like_value<-cbind(p0,hlike,p1)
    if (mord==0 && dord==1) colnames(like_value) <- c("-2h0","-2*hp","-2*p_b,v(hp)")
    if (mord==0 && dord==2) like_value<-cbind(p0,hlike,p1,p3)
    if (mord==0 && dord==2) colnames(like_value) <- c("-2h0","-2*hp","-2*p_b,v(hp)","-2*s_b,v(hp)")
    if (mord==1 && dord==1) like_value<-cbind(p0,hlike,p2,p1)
    if (mord==1 && dord==1) colnames(like_value) <- c("-2h0","-2*hp","-2*p_v(hp)","-2*p_b,v(hp)")
    if (mord==1 && dord==2) like_value<-cbind(p0,hlike,p2,p4,p1,p3)
    if (mord==1 && dord==2) colnames(like_value) <- c("-2h0","-2*hp","-2*p_v(hp)","-2*s_v(hp)","-2*p_b,v(hp)","-2*s_b,v(hp)")
    result$likelihood<-like_value
    result$iter<-print_i
    if (print_err<convergence) result$convergence<-"converged"
    if (print_err>convergence) result$convergence<-"did not converge"
    names(result$convergence) <- "convergence : "
    print(like_value,5)
    res4<-list(res2,res3)
    caic<-p0+2*df1
    n_lam<-nrow(lam_coeff)
    if (varfixed==TRUE) n_lam<-0
    maic<-p2+2*nrow(beta_coeff)+2*n_lam
    if (varfixed==TRUE) maic<-hlike+2*nrow(beta_coeff)+2*n_lam
    raic<-p1+2*n_lam
    if (RandDist=="Gamma" && mord==1 && dord==2) maic<-p4+2*nrow(beta_coeff)+2*n_lam
    if (RandDist=="Gamma" && mord==1 && dord==2) raic<-p3+2*n_lam
    aic<-cbind(caic,maic,raic)
    colnames(aic)<-c("cAIC","mAIC","rAIC")
    print(aic,5)
    result$aic<-aic
    return(result)   
}

