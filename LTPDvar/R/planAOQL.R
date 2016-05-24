planAOQL<- function (N, pbar, pL, method = c("exact", "napprox","ewmaSK","ewma2"), cm = 1,
	intdif = 20,lam=1)
{
	Gderivacex = function(x_, n_, k_) {
    	A_ = ((1/n_) + (k_^2/(2 * n_ - 2)))^0.5
    	xg_ = pnorm((x_ - k_)/A_) - (pnorm(-x_)/A_) * exp(-1 *
        	((1 - A_^2) * x_^2 - 2 * k_ * x_ + k_^2)/(2 * A_^2))
    	return(xg_)
	}
	Mnk = function(n_, k_) {
    	A_ = ((1/n_) + (k_^2/(2 * n_ - 2)))^0.5
    	xm_ = uniroot(function(x_) Gderivacex(x_, n_, k_), c(k_/(1 +
        	A_), (k_ + A_ * sqrt(k_^2 - 2 * (1 - A_^2) * log(A_)))/(1 -
        	A_^2)))$root
    	return(xm_)
	}
	k0AOQL = function(n_, pl_, nbig_) {
    	A_ = function(k_) ((1/n_) + (k_^2/(2 * n_ - 2)))^0.5
    	delta = function(k_) {
        	kMnk_ = Mnk(n_, k_)
        	-(pnorm(-kMnk_) * pnorm((kMnk_ - k_)/(((1/n_) + k_^2/(2 *
            	n_ - 2))^0.5)) - pl_/(1 - n_/nbig_))/((-pnorm(-kMnk_)/((A_(k_))^3 *
            	sqrt(2 * pi))) * (1/n_ + k_ * kMnk_/(2 * (n_ -
            	1)) * exp(-(kMnk_ - k_)^2/(2 * (A_(k_))^2))))
    	}
    	fra2 = function(i_) {
        	K_ = 1.6
        	for (i__ in (1:i_)) {
            	kK_ = K_
            	K_ = kK_ + delta(kK_)
        	}
        	return(K_)
    	}
    	fra2(25)
	}
	alpha0 = function(n_, pl_, nbig_, pbar_) pnorm((k0AOQL(n_,
    	pl_, nbig_) - qnorm(1 - pbar_))/((1/n_) + (k0AOQL(n_,
    	pl_, nbig_)^2/(2 * n_ - 2)))^0.5)
	ImsAOQL0 = function(n_, cm_, nbig_, pbar_, pl_) {
    	n_ * cm_ + (nbig_ - n_) * alpha0(n_, pl_, nbig_, pbar_)
	}
	fMinSearch0 = function(nl_, nu_, cm_, nbig_, pbar_, pl_) {
    	nl_init_ = nl_
    	nu_init_ = nu_
    	fMS = function(nl_, nu_, cm_, nbig_, pbar_, pl_) {
        	ifelse(nl_ == nu_, nl_, ifelse(ImsAOQL0(nl_ + floor((nu_ -
            	nl_)/2), cm_, nbig_, pbar_, pl_) <= ImsAOQL0(nl_ +
            	floor((nu_ - nl_)/2) + 1, cm_, nbig_, pbar_,
            	pl_), fMS(floor(nl_), floor(nl_) + floor((nu_ -
            	nl_)/2), cm_, nbig_, pbar_, pl_), fMS(floor(nl_) +
            	floor((nu_ - nl_)/2) + 1, ceiling(nu_), cm_,
            	nbig_, pbar_, pl_)))
    	}
    	out_fMS0_ = fMS(nl_, nu_, cm_, nbig_, pbar_, pl_)
    	if (out_fMS0_ == nu_init_)
        	print("\n out_fMS0_: upper search interval limit reached")
    	if (out_fMS0_ == nl_init_)
        	print("\n out_fMS0_: lower search interval limit reached")
    	return(out_fMS0_)
	}
	init_ = fMinSearch0(7, N/2, cm, N, pbar, pL)
	method = match.arg(method)
	if (method == "napprox")
    	return(new("ACSPlan", n = init_, k = k0AOQL(init_, pL,
        	N)))
	if (method %in% c("ewmaSK","ewma2","exact")) {
    	fMSmodq2 = function(pl0_, pu0_, n_, k_, nbig_,type, lam) {
        	fMSmodqOpt = function(p_) AOQ(p_, n_, k_, nbig_,type, lam)
        	p_centre_index_init_ = which.max(sapply(seq(pl0_,
            	pu0_, length = 500), function(p_) AOQ(p_, n_,
            	k_, nbig_,type, lam)))
        	pl_init_ = seq(pl0_, pu0_, length = 500)[max(1, p_centre_index_init_ -
            	1)]
        	pu_init_ = seq(pl0_, pu0_, length = 500)[min(500,
            	p_centre_index_init_ + 1)]
        	outp = optimize(f = fMSmodqOpt, interval = c(pl_init_,
            	pu_init_), maximum = T)
        	outpx = outp$maximum
        	if ((pu_init_ - outpx) < 1e-06)
            	print("in fMS: upper search interval limit reached")
        	if ((outpx - pl_init_) < 1e-06)
            	print("in fMS: lower search interval limit reached")
        	if (min(AOQ(pl_init_, n_, k_, nbig_), AOQ(pu_init_,
            	n_, k_, nbig_,type, lam)) == AOQ(seq(pl0_, pu0_, length = 500)[p_centre_index_init_],
            	n_, k_, nbig_,type, lam))
            	print("in fMS: constant objective; unsuitable interval?")
        	return(outp$objective)
    	}
    	kAOQL = function(n_, pl_, nbig_,type, lam) {
        	k1_ = uniroot(function(k_) fMSmodq2(1e-06, 0.3, n_,
            	k_, nbig_,type, lam) - pl_, c(max(k0AOQL(n_, pl_, nbig_) -
            	0.06, 1.2), 3.2), tol = .Machine$double.eps,
            	maxiter = 5000)$root
        	return(k1_)
    	}
    	fMinSearch = function(nl_, nu_, cm_, nbig_, pbar_, pl_,type, lam) {
        	ImsAOQL = function(n_, cm_, nbig_, pbar_, pl_,type, lam) {
            	kk_ = kAOQL(n_, pl_, nbig_,type, lam)
            	n_ * cm_ + (nbig_ - n_) * (1 - OC(pbar_, n_,
              	kk_, type,lam))
        	}
        	nl_init_ = nl_
        	nu_init_ = nu_
        	fMS = function(nl_, nu_, cm_, nbig_, pbar_, pl_,type, lam) {
            	ifelse(nl_ == nu_, nl_, ifelse(ImsAOQL(nl_ +
              	floor((nu_ - nl_)/2), cm_, nbig_, pbar_, pl_,type, lam) <=
              	ImsAOQL(nl_ + floor((nu_ - nl_)/2) + 1, cm_,
                	nbig_, pbar_, pl_,type, lam), fMS(floor(nl_), floor(nl_) +
              	floor((nu_ - nl_)/2), cm_, nbig_, pbar_, pl_,type, lam),
              	fMS(floor(nl_) + floor((nu_ - nl_)/2) + 1,
                	ceiling(nu_), cm_, nbig_, pbar_, pl_,type, lam)))
        	}
        	out_fMS_ = fMS(nl_, nu_, cm_, nbig_, pbar_, pl_,type, lam=lam)
        	if (out_fMS_ == nu_init_)
            	print("upper search interval limit reached")
        	if (out_fMS_ == nl_init_)
            	print("lower search interval limit reached")
        	return(out_fMS_)
    	}
    	kfMS = fMinSearch(max(7, init_ - intdif), min(N, init_ +
        	intdif), cm, N, pbar, pL,type=method, lam=lam)
    	return(new("ACSPlan", n = kfMS, k = kAOQL(kfMS, pL, N,type=method,lam=lam)))
	}
}

