lca.model <-
function(ped,probs,param,optim.param,fit=TRUE,optim.probs.indic=c(TRUE,TRUE,TRUE,TRUE),tol=0.001,x=NULL,var.list=NULL,famdep=TRUE,modify.init=NULL)
{
    fam <- ped[,1]
    id <- ped[,2]
    dad <- ped[,3]
    mom <- ped[,4]
    sex <- ped[,5]
    status <- ped[,6]

	y <- as.matrix(ped[,-(1:6)])

    K <- length(probs$p)

	if(famdep==TRUE&(is.null(probs$p)|is.null(probs$p0)|is.null(probs$p.trans))) stop("one component of probs is missing, case familial dependence")
	
	if(famdep==TRUE&is.null(probs$p.found))
	{
		#if all founds are affected
		if(all(status[dad==0]==2)) probs$p.found <- 1
		else
		{
			#if all founds are unaffected
			if(all(status[dad==0]==1)) probs$p.found <- 0
			else
			#if at least two different statuses for founders
			stop("component p.found of probs is needed")
		}
	}
	
	if(famdep==TRUE&is.null(probs$p.child))
	{
		#if all children are affected
		if(all(status[dad>0]==2)) probs$p.child <- 1
		else
		{
			#if all children are unaffected
			if(all(status[dad>0]==1)) probs$p.child <- 0
			else
			#if at least two different statuses for children
			stop("component p.child of probs is needed")
		}
	}

	if(famdep==FALSE&(is.null(probs$p)|is.null(probs$p0)|is.null(probs$p.aff))) stop("one component of probs is missing, case indpendent individuals")

	children <- which(dad>0)
	connects <- intersect(children,union(which(id%in%dad),which(id%in%mom)))
	if(length(connects)>0&famdep==TRUE&is.null(probs$p0connect)) stop("p0connect must be given when analysing extended pedigrees")

	optim.param <- attrib.dens(optim.param)
	if(is.null(attr(optim.param,"type"))) stop("function optim.param must have an attribute of name \'type\'")
	else
	{
		if(attr(optim.param,"type")=="ordi") dens <- dens.prod.ordi
		if(attr(optim.param,"type")=="norm") dens <- dens.norm
	}

    if(any(dad!=0&mom==0)|any(dad==0&mom!=0)) stop("there is individuals with a single parent\n")
    if(all(dad==0)&all(mom==0)) stop("it needs at least one child for computing p.trans\n")
    if(length(id)<ncol(y)) warning("not enough individuals to have invertible variance-covariance matrix\n")

	if(any(status[!id%in%c(dad,mom)]!=2)) stop("All children without descendants must have a symptom status = 2")

    peel <- list()
    for(f in unique(fam))
    {
        id.fam <- id[fam==f]
        dad.fam <- dad[fam==f]
        mom.fam <- mom[fam==f]
        sex.fam <- sex[fam==f]
        status.fam <- status[fam==f]

        couple <- cbind(unique(dad.fam[dad.fam>0]),unique(mom.fam[mom.fam>0]))
        couple <- rbind(couple,cbind(couple[,2],couple[,1]))
        # Changement par rapport aux versions 1.1 et antérieur pour s'adapter au module kinship2
        depth <- kindepth(id.fam,dad.fam,mom.fam,align=TRUE)
        generation <- max(depth)
        out <- NULL
        for(generat in 1:generation) out[generat] <- sum(depth==generat-1)/2
        peel.connect <- matrix(0,nrow=generation,ncol=max(out))
        peel.connect[generation,1] <- id.fam[depth==0][1]
        if(generation>1) for(generat in 1:(generation-1))
        {
            connect.generat <- intersect(id.fam[depth==generat][(id.fam[depth==generat])&(dad.fam[depth==generat]>0)],union(dad.fam,mom.fam))
            peel.connect[generation-generat,1:length(connect.generat)] <- connect.generat
        }
        peel[[which(unique(fam)==f)]] <- list("generation"=generation,"peel.connect"=peel.connect,"couple"=couple)
    }
    cat("CHECK OF ALL PEDIGREES DONE\n")

    if(any(probs$p<0)|abs(sum(probs$p)-1)>.Machine$double.eps) stop("there is a problem with the initail value of founder calss probability p.\n")
    if(famdep) if(any(probs$p.trans<0)|any(abs(apply(array(probs$p.trans[,,1:K],dim=rep(K,times=3)),2:3,sum)-1)>.Machine$double.eps)|
       any(abs(apply(matrix(probs$p.trans[,1:K,K+1],nrow=K,ncol=K),2,sum)-1)>.Machine$double.eps))
    stop("there is a problem with the initail value of transition probability p.trans.\n")

	y.x.aff <- as.matrix(y[status==2,])
    if(!is.null(x)) y.x.aff <- cbind(y.x.aff,as.matrix(x[status==2,]))

    fyc <- matrix(1,nrow=length(id),ncol=K+1)
    fyc[status==2,1:K] <- t(apply(y.x.aff,1,dens,param,var.list))
    if(K>1) prob0 <- apply(fyc[status==2,1:K],2,function(vec) all(vec<.Machine$double.eps))
    else prob0 <- fyc[status==2,1]<.Machine$double.eps
    if(fit&any(prob0))
    {
        if(!is.null(modify.init)) while(any(prob0))
        {
            param <- modify.init(prob0,param)
            fyc[status==2,1:K] <- t(apply(y.x.aff,1,dens,param,var.list))
            if(K>1) prob0 <- apply(fyc[status==2,1:K],2,function(vec) all(vec<.Machine$double.eps))
            else prob0 <- fyc[status==2,1]<.Machine$double.eps
        }
        else stop("initial parameter values produce zero probabilities for at least one subject.\n")
    }

    iter <- n.warning <- 0
    ll <- -Inf
    continue <- 1
    while(continue)
    {
        iter <- iter+1
        cat("%%%%%%%%%%%%%%%%%%%%%%%%%%% iteration",iter,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
        res.weight <- e.step(ped,probs,param,dens,peel,x,var.list,famdep)
        ll.new <- res.weight$ll
        if(!fit) continue <- 0
        else
        {
            probs <- optim.probs(ped,probs,optim.probs.indic,res.weight,famdep)

		y.aff <- as.matrix(y[status==2,])
            weight <- as.matrix(res.weight$w[,1,1:K])

            param <- optim.param(y.aff,status,weight,param,x=x,var.list=var.list)

            if(ll.new<ll)
            {
                cat("log-likelihood is decreasing!!!!!\n")
                n.warning <- n.warning+1
            }
            else if(ll.new-ll<tol|n.warning>=20) continue <- 0
        }
        ll <- ll.new
    }
    res <- list("probs"=probs,"param"=param,"weight"=res.weight$w,"ll"=ll)
    res
}

