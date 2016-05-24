spfeml<-function(formula, data=list(), index=NULL, listw, listw2 = NULL, na.action, model = c("lag","error", "sarar"), effects = c('spfe','tpfe','sptpfe'), method="eigen", quiet = TRUE, zero.policy = NULL, interval1 = NULL, interval2 = NULL, trs1 = NULL, trs2 = NULL, tol.solve = 1e-10, control = list(), legacy = FALSE, llprof = NULL, cl = NULL, Hess = FALSE, LeeYu = FALSE, ...){

	  
        # timings <- list()
       # .ptime_start <- proc.time()

model<-match.arg(model)
effects <- match.arg(effects)


if (model == "sarar") con <- list(LAPACK = FALSE,  Imult = 2L, cheb_q = 5L, MC_p = 16L, MC_m=30L, super=FALSE, opt_method = "nlminb", opt_control = list(), pars = NULL, npars = 4L, pre_eig1 = NULL, pre_eig2 = NULL)

else     con <- list(tol.opt = .Machine$double.eps^0.5,  Imult = 2, cheb_q = 5, MC_p = 16, MC_m = 30, super = NULL, spamPivot = "MMD", in_coef = 0.1, type = "MC", correct = TRUE, trunc = TRUE, SE_method = "LU", nrho = 200, interpn = 2000, SElndet = NULL, LU_order = FALSE, pre_eig = NULL)
	


nmsC <- names(con)
con[(namc <- names(control))] <- control
    
    if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))

##    if (is.null(quiet)) # now this has a default in spml(), hence it never is
##	quiet <- !get("verbose", envir = spdep:::.spdepOptions)
##    stopifnot(is.logical(quiet))

	if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
	  
	  
	  ## reorder data if needed
  if(!is.null(index)) {
    #require(plm)
    data <- plm.data(data, index)
    }

  index <- data[,1]
  tindex <- data[,2]

	  ## record call
        ## now passed on from spml() but to be sure:
        if(is.null(cl)) cl <- match.call()

#check the model
# model<-match.arg(model)


#check the effects
effects<-match.arg(effects)


  ## check
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")

  ## reduce X,y to model matrix values (no NAs)
  x<-model.matrix(formula,data=data)
  clnames <- colnames(x)
  rwnames <- rownames(x)

  y<-model.response(model.frame(formula,data=data))
  ## reduce index accordingly
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]

  ## reorder data by cross-sections, then time
  oo<-order(tind,ind)
  x<-x[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]


  #make sure that the model has no intercept if effects !=pooled
  if (attr(attributes(model.frame(formula,data=data))$terms, "intercept") == 1) {
  	x <- as.matrix(x[,-1])
    colnames(x)<-clnames[-1]
    dimnames(x)[[1]]<-rwnames
    clnames <- clnames[-1]
  }
	x <- as.matrix(x)
    k <- dim(x)[2]
      
    

  ## det. number of groups and df
  N<-length(unique(ind))
  n<-N
  ## det. max. group numerosity
  T<-max(tapply(tind,ind,length))

  ## det. total number of obs. (robust vs. unbalanced panels)
  NT<-length(ind)


  mt<-terms(formula,data=data)
  mf<-lm(formula,data,method="model.frame")#,na.action=na.fail

  na.act<-attr(mf,'na.action')


##checks on listw
  if(is.matrix(listw)) {
    if(dim(listw)[[1]] !=N ) stop("Non conformable spatial weights")
    #require(spdep)
    listw <- mat2listw(listw)
   }
  if (!inherits(listw, "listw"))
        stop("No neighbourhood list")
     
 		can.sim <- FALSE
    if (listw$style %in% c("W", "S")) 
        can.sim <- can.be.simmed(listw)
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
}
     
###specific checks for the SARAR        
if(model == "sarar"){
	
	 if (is.null(listw2)) 
        listw2 <- listw
     if (!is.null(con$pre_eig1) && is.null(con$pre_eig2)) 
        con$pre_eig2 <- con$pre_eig1
    else if (!inherits(listw2, "listw")) 
        stop("No 2nd neighbourhood list")

# if (is.null(con$fdHess)) con$fdHess <- method != "eigen"
        # stopifnot(is.logical(con$fdHess))

     if (!is.null(con$pars)) {
        stopifnot(is.numeric(con$pars))
    }
    stopifnot(is.integer(con$npars))
    # stopifnot(is.logical(con$fdHess))
    stopifnot(is.logical(con$LAPACK))
    stopifnot(is.logical(con$super))

    can.sim2 <- FALSE
    if (listw2$style %in% c("W", "S")) 
        can.sim2 <- can.be.simmed(listw2)
    if (!is.null(na.act)) {
        subset <- !(1:length(listw2$neighbours) %in% na.act)
        listw2 <- subset(listw2, subset, zero.policy = zero.policy)
    }

	}


switch(model, lag = if (!quiet) cat("\n Spatial Lag Fixed Effects Model \n"),
	    error = if (!quiet) cat("\n Spatial Error Fixed Effects Model\n"),
	    sarar = if (!quiet) cat("\n Spatial SARAR Fixed Effects Model\n"),
	    stop("\nUnknown model type\n"))


  ## check whether the panel is balanced
  balanced<-N*T==NT
  if(!balanced) stop("Estimation method unavailable for unbalanced panels")


	indic<-seq(1,T)
	inde<-as.numeric(rep(indic,each=N)) ####takes the first n observations
	indic1<-seq(1,N)
	inde1<-rep(indic1,T) ####takes observations 1,  n+1, 2n+1...
  ### generates Wy if model=='lag'
  
  
env <- new.env(parent=globalenv())
assign("y",y, envir=env)
assign("x",x, envir=env)
assign("listw",listw, envir=env)
assign("NT",NT, envir=env)
assign("T",T, envir=env)
assign("k",k, envir=env)
assign("n",n, envir=env)


wy<-unlist(tapply(y,inde, function(u) lag.listw(listw,u, zero.policy = zero.policy), simplify=TRUE))
	

#demeaning of the y and x variables depending both on model and effects

if (effects=="tpfe" | effects=="sptpfe"){
	ytms<-tapply(y,inde,mean) ####for each time period takes the mean for the cross section observations
	tpms<-function(q) tapply(q,inde,mean)
	xtms<-apply(x,2,tpms)   ###same thing for the X variable
	ytm<-rep(ytms,each=N) ##generate the NT variables
	xtm<-matrix(,NT,k)
	for (i in 1:k) xtm[,i]<-rep(xtms[,i],each=N)

if (model %in% c("lag", "sarar")) {
		wytms<-tapply(wy,inde,mean) ###same thing for Wy
		wytm<-rep(wytms,each=N)
assign("wytms",wytms, envir=env)
				
}

assign("ytms",ytms, envir=env)
assign("xtms",xtms, envir=env)
	}


if (effects=="spfe" | effects=="sptpfe"){
	ysms<-tapply(y,inde1,mean) ###for each cross-sectional unit takes the mean over the time periods
	spms<-function(q) tapply(q,inde1,mean)
	xsms<-apply(x,2,spms)
	ysm<-rep(ysms,T)
	xsm<-matrix(,NT,k)
	for (i in 1:k) xsm[,i]<-rep(xsms[,i],T)

if (model %in% c("lag", "sarar")){
			wysms<-tapply(wy,inde1,mean)
			wysm<-rep(wysms,T)
			assign("wysms",wysms, envir=env)
			}
			
assign("ysms",ysms, envir=env)
assign("xsms",xsms, envir=env)		
	}
	
	
# if (effects=='pooled'){
	# yt<-y  	###keep the variables with no transformation
	# xt<-x
	# }


if (effects=="tpfe"){ ####generate the demeaned variables for tpfe
	yt<-y-ytm
	xt<-x-xtm
						}


if(effects=="spfe"){ ####generate the demeaned variables for spfe
	yt<-y-ysm
	xt<-x-xsm

	 					}

if (effects=="sptpfe"){ ####generate the demeaned variables for both types of FE
	yt<-y - ysm - ytm + rep(mean(y),NT)
	xmm<-matrix(,NT,(k))
	for (i in 1:(k)) xmm[,i]<-rep(mean(x[,i]),NT)
	xt<-x - xsm - xtm + xmm
								}
								

	wyt<-unlist(tapply(yt,inde, function(u) lag.listw(listw,u), simplify=TRUE))

if(model=="sarar")	{
	w2yt<-unlist(tapply(yt,inde, function(u) lag.listw(listw2,u), simplify=TRUE))
	w2wyt<-unlist(tapply(wyt,inde, function(u) lag.listw(listw2,u), simplify=TRUE))
	
	}
								
								
if 	(model == "error"){
	dm<-function(A) trash<-unlist(tapply(A,inde,function(TT) lag.listw(listw,TT), simplify=TRUE))
   wxt<-apply(xt,2,dm)
   # colnames(wxt)<-paste('Lag.',colnames(x), sep="")
   wx<-apply(x,2,dm)
   # colnames(wx)<-paste('lag.',colnames(x), sep="")
	}

if 	(model == "sarar"){
	dm<-function(A) trash<-unlist(tapply(A,inde,function(TT) lag.listw(listw2,TT), simplify=TRUE))
   wxt<-apply(xt,2,dm)
   # colnames(wxt)<-paste('Lag.',colnames(x), sep="")
   wx<-apply(x,2,dm)
   # colnames(wx)<-paste('lag.',colnames(x), sep="")
	}

# print(clnames)
colnames(xt)<- clnames

	

assign("yt",yt, envir=env)
assign("xt",xt, envir=env)
assign("wyt",wyt, envir=env)
assign("wy",wy, envir=env)


if (model %in% c("error", "sarar")){
	assign("wx",wx, envir=env)
	assign("wxt",wxt, envir=env)
if (model == "sarar")	{
	assign("w2yt",w2yt, envir=env)
	assign("w2wyt",w2wyt, envir=env)
	assign("listw2", listw2, envir = env)
	assign("can.sim2", can.sim2, envir=env)
	assign("similar2", FALSE, envir = env)
	}
	} 

if(model %in% c("lag", "error") ){
	
	assign("compiled_sse", con$compiled_sse, envir=env)	
	assign("f_calls", 0L, envir = env)
    assign("hf_calls", 0L, envir = env)

}



assign("verbose", !quiet, envir = env)
# assign("first_time", TRUE, envir = env)
assign("LAPACK", con$LAPACK, envir = env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir = env)
assign("family", "SAR", envir = env)
assign("inde",inde, envir=env)
assign("con", con, envir=env)



    if (!quiet) 
        cat(paste("\nSpatial fixed effects model\n", "Jacobian calculated using "))

if(model == "lag"){
    interval1 <- jacobianSetup(method, env, con, pre_eig = con$pre_eig, trs = trs1, interval = interval1)
    assign("interval1", interval1, envir = env)


    RES<- splaglm(env = env, zero.policy = zero.policy, interval = interval1, con = con, llprof = llprof, tol.solve= tol.solve, Hess = Hess, method = method, LeeYu = LeeYu, effects = effects)
    
    res.eff<-felag(env = env, beta = RES$coeff, sige = RES$s2, effects = effects, method = method, lambda = RES$lambda, legacy = legacy, zero.policy = zero.policy)    

	}

if(model == "sarar"){
	
    interval1 <- jacobianSetup(method, env, con, pre_eig = con$pre_eig1, trs = trs1, interval = interval1, which = 1)
    assign("interval1", interval1, envir = env)
    interval2 <- jacobianSetup(method, env, con, pre_eig = con$pre_eig2, trs = trs2, interval = interval2, which = 2)
    assign("interval2", interval2, envir = env)
    # nm <- paste(method, "set_up", sep = "_")
    # timings[[nm]] <- proc.time() - .ptime_start
    # .ptime_start <- proc.time()
    
      RES<- spsararlm(env = env, zero.policy = zero.policy, con = con, llprof = llprof, tol.solve = tol.solve, Hess = Hess, LeeYu = LeeYu, effects = effects)
  
  
res.eff<-felag(env = env, beta=RES$coeff, sige=RES$s2, effects = effects ,method = method, lambda = RES$lambda, legacy = legacy, zero.policy = zero.policy)    	


		}


if (model=='error'){

    interval1 <- jacobianSetup(method, env, con, pre_eig = con$pre_eig, trs = trs1, interval = interval1)
    assign("interval1", interval1, envir = env)
    # nm <- paste(method, "set_up", sep = "_")
    # timings[[nm]] <- proc.time() - .ptime_start
    # .ptime_start <- proc.time()	

  RES <- sperrorlm(env = env, zero.policy = zero.policy, interval = interval1, Hess = Hess, LeeYu = LeeYu, effects = effects)	
    	res.eff<-feerror(env = env, beta=RES$coeff, sige=RES$s2, effects = effects ,method =method, rho=RES$rho, legacy = legacy)
    	
    }
    
	
	

    ##calculate the R-squared
    yme <- y-mean(y)
    rsqr2 <- crossprod(yme)
    rsqr1 <- crossprod(res.eff[[1]]$res.e)
    res.R2<- 1- rsqr1/rsqr2

	#generate fixed values (from fixed_effect)
	y.hat <- res.eff[[1]]$xhat
	res <- as.numeric(res.eff[[1]]$res.e)
	N.vars<-res.eff$N.vars



	nam.rows <- dimnames(x)[[1]]
	
	
   names(y.hat) <- nam.rows
   names(res) <- nam.rows


	## make model data
   model.data <- data.frame(cbind(y,x))
   dimnames(model.data)[[1]] <- nam.rows




if (model == "lag")   spat.coef<-RES$lambda
if (model == "error") spat.coef<-RES$rho
if (model == "sarar") spat.coef <- c(RES$lambda, RES$rho)


Coeff<-c(spat.coef, RES$coeff)


type <- paste("fixed effects", model)

if (Hess){

	if(model == "lag" ){
   		var<-matrix(0,(ncol(RES$asyvar1)+1),(ncol(RES$asyvar1)+1))
		var[1,1]<-	RES$lambda.se
		var[(2:ncol(var)),(2:ncol(var))]<-RES$asyvar1
	}
	
	if(model == "error" ){
	 	var<-matrix(0,(ncol(RES$asyvar1)+1),(ncol(RES$asyvar1)+1))
    	var[1,1]<-	RES$rho.se
    	var[(2:ncol(var)),(2:ncol(var))]<-RES$asyvar1
	}
	
	if(model == "sarar"){
		var <- matrix(0,(ncol(RES$asyvar1)+2),(ncol(RES$asyvar1)+2))
	    var[1,1] <-	RES$lambda.se
	    var[2,2] <-	RES$rho.se
	    var[(3:ncol(var)),(3:ncol(var))] <- RES$asyvar1
	}
	
} 

else{

if(model == "lag" ){
   var<-matrix(0,(ncol(RES$asyvar1)+1),(ncol(RES$asyvar1)+1))
   var[1,1]<-	RES$lambda.se
   var[(2:ncol(var)),(2:ncol(var))]<-RES$asyvar1
	}

if(model == "error" ){
	var<-matrix(0,(ncol(RES$asyvar1)+1),(ncol(RES$asyvar1)+1))
   var[1,1]<-	RES$rho.se
   var[(2:ncol(var)),(2:ncol(var))]<-RES$asyvar1
	}
	
if(model == "sarar"){
   var<-matrix(0,(ncol(RES$asyvar1)+2),(ncol(RES$asyvar1)+2))
   var[1,1]<-	RES$lambda.se
   var[2,2]<-	RES$rho.se
   var[(3:ncol(var)),(3:ncol(var))]<-RES$asyvar1

	}

}


spmod <- list(coefficients=Coeff, errcomp=NULL,
                vcov = var ,spat.coef=spat.coef,
                vcov.errcomp=NULL,
                residuals=res, fitted.values=y.hat,
                sigma2=RES$s2, type=type, model = model.data,
                call=cl, logLik=RES$ll, method = method, effects=effects, 
                res.eff=res.eff)
                
if (!is.null(na.act)) 
        spmod$na.action <- na.act
                
  class(spmod) <- "splm"
  return(spmod)

	}
       
       
       
       
