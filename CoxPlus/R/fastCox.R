#' Cox Regression (Proportional Hazards Model) with Multiple Causes and Mixed Effects
#' @docType package
#' @name CoxPlus
#' @description A high performance package estimating Proportional Hazards Model when an even can have more than one causes, including support for random and fixed effects, tied events, and time-varying variables.
#' @param head A data frame with 4~5 columns: start, stop, event, weight, strata (optional).
#' @param formula A formula specifying the independent variables
#' @param par A optional list of parameters controlling the estimation process
#' @param data The dataset, a data frame containing observations on the independent variables
#' @return A list containing the estimated parameters
#' @export
#' @importFrom Rcpp evalCpp Module initialize cpp_object_initializer
#' @importFrom stats model.matrix rnorm pchisq
#' @importFrom utils tail
#' @importFrom methods new
#' @useDynLib CoxPlus
#' @examples
#' # Simulate a dataset. lam=exp(x), suvtime depends on lam
#' x = rnorm(5000)
#' suvtime = -log(runif(length(x)))/exp(x)
#' # Censor 80% of events
#' thd = quantile(suvtime, 0.2)
#' event = as.numeric(suvtime <= thd)
#' suvtime[suvtime>thd] = thd
#'
#' # The estimates of beta should be very close to 1, the true value
#' head = cbind(start=0,stop=suvtime,event=event,weight=1)
#' est = fastCox(head,~x)
#' print(est$result)
#' @references 1. Jing Peng, Ashish Agarwal, Kartik Hosanagar, and Raghuram Iyengar. Towards Effective Information Diffusion on Social Media Platforms: A Dyadic Analysis of Network Embeddedness. Working Paper.
#' @references 2. Jing Peng, Ashish Agarwal, Kartik Hosanagar, and Raghuram Iyengar. Toward Effective Social Contagion: A Micro Level Analysis of the Impact of Dyadic Network Relationship. In Proceedings of the 2014 International Conference on Information Systems.
fastCox = function(head, formula, par=list(), data=NULL){
    # require(Rcpp)
    begin = Sys.time()
    # check on start and stop time
    if(any(head[,'start']>=head[,'stop'])){
        print(head[head[,'start']>=head[,'stop'],])
        stop('Error: start time larger than or equal to stop time!!!')
    }
    # order, should be used on all vectors of length n in par!!!!
    event_ix = head[,'event']>0
    if("strata" %in% colnames(head)) {
        ord = order(head[,'strata'], head[,'event'], head[,'stop'], head[,'start'], decreasing=T)
        nEvents = length(unique(paste(head[event_ix,'strata'], head[event_ix,'stop'])))
    }else{
        ord = order(head[,'event'], head[,'stop'], head[,'start'], decreasing=T)
        nEvents = length(unique(head[event_ix, 'stop']))
    }
    head = head[ord,]
    if(!is.null(par$sender)) par$sender=par$sender[ord]
    if(!is.null(par$receiver)) par$receiver=par$receiver[ord]

    # groups
    if(!is.null(par$egroups)) {
        if(any(is.na(par$egroups))) stop("There should be no NAs in par$egroups!")
        par$egroups = par$egroups[ord]
        egroups = factor(par$egroups)
        lvls = levels(egroups)
        ix = 1:length(lvls)
        names(ix) = lvls
        par$egroups = ix[egroups] - 1
    }
    if(!is.null(par$fgroups)) {
        par$fgroups = as.matrix(par$fgroups)
        if(any(is.na(par$fgroups))) stop("There should be no NAs in par$fgroups!")
        if(ncol(par$fgroups)<1) stop("Number of frailty groups should be >=1!")
        par$fgroups = par$fgroups[ord,]
        par$fgroups = as.matrix(par$fgroups)

        par$fglen = rep(0, ncol(par$fgroups))
        fgroups = matrix(0,nrow(par$fgroups),ncol(par$fgroups))
        lvls = list()
        accIx = 0
        for(i in 1:ncol(par$fgroups)){
            grp = as.character(par$fgroups[,i])
            # To avoid negative se, aggregate all rare groups into one group
            if(!is.null(par$fgroupThd)){
                grpCnt = table(grp)
                rareGrp = names(grpCnt)[which(grpCnt<par$fgroupThd)]
                grp[grp %in% rareGrp] = "_rareGroup_"
                writeLines(paste("Level of group",i,"reduced by",length(rareGrp)-1,"from",length(grpCnt),", total rare group obs",sum(grpCnt[rareGrp])))
            }
            grp = factor(grp)
            lvls[[i]] = levels(grp)
            ix = 1:length(lvls[[i]])
            names(ix) = lvls[[i]]
            fgroups[,i] = ix[grp] - 1 + accIx
            par$fglen[i] = length(lvls[[i]])
            accIx = accIx + par$fglen[i]
        }

        par$fgroups = fgroups
        par$nf = sum(par$fglen)
        frailNames = paste("Group", 1:par$nf, sep=":")
        # By default, assume each group has the same variance, repeat the corresponding entry respective times
        # par$thetagroups = rep(0:(length(par$fglen)-1),par$fglen)
        if(is.null(par$nFixedGroups) || par$nFixedGroups==0) {
            par$nFixedGroups=0
            par$isFixedEffect = -1
        }else if(par$nFixedGroups>=length(par$fglen)){
            par$nFixedGroups=length(par$fglen)
            par$isFixedEffect = 1
        } else par$isFixedEffect = 0
        par$thetagroups = NULL
        par$randStarts = ifelse(par$nFixedGroups>0,sum(par$fglen[1:par$nFixedGroups]), 0)
        k = 0
        par$validFrailty = rep(1, par$nf) # reduce one degree of freedom for each fixed effect term
        for(i in 1:length(par$fglen)){
            tmp = rep(k, par$fglen[i])
            if(!is.null(par$specialFrailty) && length(par$specialFrailty)>=i && length(par$specialFrailty[[i]])>=1) {
                spf = unique(par$specialFrailty[[i]])
                for(j in 1:length(spf)) tmp[lvls[[i]]==spf[j]] = k+j
            }
            par$thetagroups = c(par$thetagroups, tmp)
            k = max(par$thetagroups) + 1
            # set the last level of fixed effect to be the reference
            if(i<=par$nFixedGroups) par$validFrailty[sum(par$fglen[1:i])] = 0
        }
    }
    # offset
    if(!is.null(par$offset)) par$offset = par$offset[ord]
    # delta
    if(!is.null(par$delta)) par$delta = par$delta[ord]
    if(!is.null(par$beta_last) && !is.null(par$delta)) print("Delta will not be used when beta_last is specified")

    # strata
    if("strata" %in% colnames(head)){
        s = factor(head[,'strata'])
        nStrata = length(unique(s))
        strats = sort(unique(s), decreasing=T)
        counts = table(s)
        offset = 0;
        strataBounds = matrix(rep(0, 2*nStrata), ncol=2)
        for(i in 1:nStrata){
            strataBounds[i,]=c(offset,offset+counts[strats[i]]-1)
            offset = offset + counts[strats[i]]
        }
        par$strataBounds = strataBounds
    }


    # model.matrix
    if(is.null(data)) {
        X = model.matrix(formula)
    }else X = model.matrix(formula, data)
    if(nrow(X) != length(ord)) stop("nrow(X)!=length(ord), there might be missing values in the regressors")
    X = X[ord,]
    X = X[, -1, drop=F] # drop intercept
    # deleted not so useful demean because it may affect interaction terms

    # figure out the appropriate names for beta
    fnames = colnames(X)
    realnames = fnames
    if(!is.null(par$dropCols)){
        par$dropFlag = as.integer(fnames %in% c(par$dropCols, par$fixedCol))
        realnames = fnames[!fnames %in% c(par$dropCols, par$fixedCol)]
    }
    if(!is.null(par$fgroups)) realnames = c(realnames, frailNames)
    # time varying cols
    ix = 1:ncol(X)
    names(ix) = fnames
    if(!is.null(par$fixedCol)) par$fixedCoefIndex = ix[par$fixedCol] - 1
    if(!is.null(par$timeVarCols) && length(intersect(par$timeVarCols,fnames))>0){
        par$timeVarCols = intersect(par$timeVarCols,fnames)
        par$timeVarIndex = ix[par$timeVarCols] -1
        # print(par$timeVarIndex)
        if("age" %in% par$timeVarCols) par$ageIndex=which(par$timeVarCols=="age")-1
        if("logDiggNum" %in% par$timeVarCols) par$diggNumIndex=which(par$timeVarCols=="logDiggNum")-1
        # construction a map for interaction terms
        subnames = strsplit(fnames, split=":")
        nInteractions = 0
        for(i in 1:ncol(X)){
            snames = unlist(subnames[i])
            if(length(snames)>1 && any(snames %in% par$timeVarCols)){
                nInteractions = nInteractions+1
            }
        }
        # does not support three way interactions at this moment
        interMap = matrix(rep(-1, 3*nInteractions), ncol=3)
        j = 1
        for(i in 1:ncol(X)){
            snames = unlist(subnames[i])
            if(length(snames)>1 && any(snames %in% par$timeVarCols)){
                interMap[j,] = c(i, ix[snames]) - 1
                j = j+1
            }
        }
        par$interMap = interMap
    }

    # verbose
    if(is.null(par$verbose)) par$verbose = 1
    if(nrow(X)>100000) par$verbose = max(2, par$verbose)
    # show ties or not distinguishable stop times
    stoptime = head[,'stop'][head[,'event']>0]
    time_diff = abs(diff(stoptime))
    isDistinguishable = time_diff < 1e-9 * pmin(stoptime[1:(length(stoptime)-1)], stoptime[2:length(stoptime)])
    if(!is.null(par$showties) && any(isDistinguishable)) {
        print(paste('Warning: The difference between ', sum(isDistinguishable), ' stop times are not distinguishable!!!'))
        false_ix = which(isDistinguishable)
        row1 = stoptime[false_ix]
        row2 = stoptime[1+false_ix]
        print(cbind(false_ix, row1, row2, abs(row2-row1)/pmin(row1,row2)))
    }

    # estimate
    if(is.null(par$method)) par$method = 1
    mod=Rcpp::Module("cox_module", PACKAGE='CoxPlus')
    X = t(X)
    if(!is.null(par$isRandBeta) && par$isRandBeta) {
        writeLines("Random initial values are provided for beta")
        par$beta=par$randScale*rnorm(nrow(X) + ifelse(!is.null(par$nf), par$nf, 0))
        if(!is.null(par$isNBetaOnly) && par$isNBetaOnly && !is.null(par$nf))
            par$beta[(nrow(X)+1):length(par$beta)] = 0
        print(head(par$beta))
        print(tail(par$beta))
    }
    head = as.matrix(head[,c("start","stop","event","weight")])

    if(is.null(par$recursive)){
        if(!is.null(par$savePath)) save(head,X,par,file=par$savePath)
        inst = new(mod$CoxReg,head,X,par)
        fit = inst$estimate()
    } else if (par$recursive %in% c("Twostage","TwostageMin","TwostageMax")){
        if(!is.null(par$delta) || !is.null(par$deltaType) ) writeLines("No need to specifiy delta and deltaType for two stage model.")
        if(par$method!=3) writeLines("Forcing to use method 3 for two stage model")
        par$delta = NULL
        par$deltaType=ifelse(par$recursive=="TwostageMin",-1,1)
        par$method=3
        if(!is.null(par$savePath)) save(head,X,par,file=par$savePath)
        inst = new(mod$CoxReg,head,X,par)
        fit = inst$estimate()
        fit = inst$estimate()
        print("Done with two-stage estimation!")
    }else if (par$recursive=="Multistage"){
        if(!is.null(par$delta) || !is.null(par$deltaType) ) writeLines("No need to specifiy delta and deltaType for multi-stage model.")
        if(par$method!=3) writeLines("Forcing to use method 3 for multi-stage model")
        par$delta = NULL
        par$method=3
        if(!is.null(par$savePath)) save(head,X,par,file=par$savePath)
        inst = new(mod$CoxReg,head,X,par)
        if(is.null(par$stages)) par$stages = 5
        if(is.null(par$lik_tol)) par$lik_tol = 1e-6
        if(is.null(par$beta_tol)) par$beta_tol = 1e-3
        beta_last = 0
        Ln = -Inf
        for(iter in 1:par$stages) {
            print(paste("Working on the",iter,"/",par$stages,"iteration..."))
            # use two stage weights to speed up may not be a good idea because w and beta are not compatible in the model. The large likelihood is meaningless because it's not the likelihood of the model of interest
            # The two-stage speed up method may have a reversed interpretation where the likelihood decreases until compatible with the model
            # starting with equal weight is a good idea because it is compatible with beta=0
            if(is.null(par$isEnableSpeedup)) par$isEnableSpeedup=0
            if(par$isEnableSpeedup==1){ inst$deltaType = ifelse(iter==2,1,0) }
            else{ inst$deltaType = 0 }
            fit = inst$estimate()
            Lc = fit$likelihood
            # The first two rounds for initializing beta if speedup
            if(iter>=2+par$isEnableSpeedup*2 && ( all( abs(fit$beta-beta_last)<abs(fit$beta)*par$beta_tol )||(Lc-Ln)<abs(Ln)*par$lik_tol ) ){
                print(paste("Early stop at iteration",iter,", L=",Lc,", Ln=",Ln,sep=""))
                break
            }
            if(iter==par$stages) writeLines(paste("---@@@ The EM algorithm has not fully converged yet, relative changes in L=",-(Lc-Ln)/Ln,"\nRelative changes in beta:",paste(format((fit$beta-beta_last)/beta_last,digit=3),collapse=" "),sep=""))
            beta_last = fit$beta
            Ln = Lc
        }

        #         print("Done with multi-stage estimation!")
    }

    # print results
    beta = fit$beta
    se = fit$se
    rse = fit$rse
    if(!is.null(par$fgroups) || length(rse)==0) {
        print(paste("length beta=",length(beta),", length se=",length(se)))
        if(length(se)==0) se = NA
        z = beta/se
        res = cbind(beta, se, exp(beta), 1 - pchisq(z^2, 1))
        colnames(res) = c("coef", "se(coef)", "exp(beta)", "Pr(>|z|)")
    }else{
        if(length(rse)==0) rse = NA
        z = beta/rse
        res = cbind(beta, se, rse, 1 - pchisq(z^2, 1))
        colnames(res) = c("coef", "se(coef)", "robust se", "Pr(>|z|)")
    }

    #     writeLines(paste("res length:",nrow(res),", name length: ", length(realnames)))
    rownames(res) = realnames

    out = list()
    out$par = par
    out$formula = formula
    out$result = res
    out$likelihood = fit$likelihood
	out$converged = fit$converged
    out$tests = cbind(LR=fit$LR, wald=fit$wald, score=fit$score)
    if(!is.null(par$detail) && par$detail==T){
        out$R = fit$R
        out$U = fit$U
        out$meta = fit$meta
    }
    if(!is.null(par$fgroups)) out$theta_se = fit$theta_se
    out$I = fit$I
    out$V = fit$V

    # AIC and BIC
    out$df = fit$p
	if(!is.null(par$fgroups)){
		# if #special frailty terms = #subjects, needs to subtract one (to be fixed in the future)
	    out$df = out$df + ncol(par$fgroups) - par$nFixedGroups # Number of variance terms for RE
		if(length(par$specialFrailty)>par$nFixedGroups) # Additional variance terms for RE
		    out$df = out$df + length(unlist(par$specialFrailty[(par$nFixedGroups+1):length(par$specialFrailty)]))
		if(par$nFixedGroups>=1) out$df = out$df + sum(par$fglen[1:par$nFixedGroups]-1) # Levels for each FE
    }
	out$AIC = 2 * out$df - 2 * out$likelihood
	if(!is.null(par$nEvents) && par$nEvents!=nEvents){
	    print('Provide nEvents and calculated nEvents do not match, will use provided')
	    nEvents = par$nEvents
	}
    out$BIC = log(nEvents) * out$df - 2 * out$likelihood
    out$nEvents = nEvents


    if(!is.null(fit$Vr)) out$Vr = fit$Vr
    if(!is.null(fit$theta)) out$theta = fit$theta
    print(Sys.time() - begin)
    rm(X, head, par, fit, inst, mod)
    out
}
