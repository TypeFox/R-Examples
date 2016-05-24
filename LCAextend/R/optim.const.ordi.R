optim.const.ordi <-
function(y,status,weight,param,x=NULL,var.list=NULL)
{
	#if K=1 use optim.noconst instead
	if(ncol(weight)==1) param <- optim.noconst.ordi(y,status,weight,param,x,var.list)
	else
	{
		if(is.null(var.list)) var.list <- rep(list(NULL),times=ncol(y))
		else if(is.null(x)) stop("covariates data frame is missing in function optim.const.ordi\n")
		S <- apply(y,2,max,na.rm=TRUE)
		alpha <- list(NULL)
		for(j in 1:ncol(y))
		{
			val.y.j <- as.numeric(levels(factor(y[,j])))
			miss.val <- which(!1:S[j]%in%val.y.j)
			if(length(miss.val)>0) S[j] <- S[j]-length(miss.val)

			z <- factor(rep(0:(ncol(weight)-1),rep(sum(status==2)+S[j]*sum(status==0),ncol(weight))))
			sympt <- rep(c(y[,j],rep(val.y.j,rep(sum(status==0),S[j]))),ncol(weight))

			if(!is.null(var.list[[j]]))
			{
				xrep1 <- as.matrix(x[status==2,])
				for(s in 1:S[j]) xrep1 <- rbind(xrep1,as.matrix(x[status==0,]))
				xrep <- NULL
				for(i in 1:ncol(weight)) xrep <- rbind(xrep,xrep1)
			}

			if(is.null(var.list[[j]])) data.lrm <- data.frame(sympt=sympt,z=z)
			else data.lrm <- data.frame(xrep,sympt=sympt,z=z)

			w <- matrix(NA,nrow=sum(status==2),ncol=ncol(weight))

			w[1:sum(status==2),] <- weight[status==2,]
            # Version 1.1: section modifiée pour tenir compte des covariables.
			if(sum(status==0)>0)
			{
			  if(length(miss.val)>0) alpha.j <- as.matrix(param$alpha[[j]][,-miss.val])
			  else alpha.j <- param$alpha[[j]]

			  if(!is.null(var.list[[j]]))
			    {  
			        S.cov <- length(var.list[[j]])
                    S.alp <- ncol(param$alpha[[j]])-S.cov+1
					
                    covar.x <- ifelse(S.cov==0,0,sum(param$alpha[[j]][k,S.alp:(S.alp+S.cov-1)]*x[status=0,var.list[[j]]]))
					# Boucle sur les sujets avec données manquantes
					for (m in 1:sum(status==0))
					{
					f.j.s.k <- t(apply(alpha.j,1,p.compute,decal=covar.x[m]))
					# On va chercher les poids du m^e sujet avec données manquantes
					for(s in 1:S[j]) w <- rbind(w,as.matrix(weight[status==0,][m,])*f.j.s.k[,s])
					}
			    }
			  else
				{
				f.j.s.k <- t(apply(alpha.j,1,p.compute))
				for(s in 1:S[j]) w <- rbind(w,as.matrix(weight[status==0,])*f.j.s.k[,s])
				}
			}
			if(is.null(var.list[[j]])) formula.lrm <- sympt~z
			else formula.lrm <- eval(parse(text=paste(c("sympt~z",names(data.lrm)[var.list[[j]]]),collapse="+")))

			lrm.coef <- coef(lrm(formula.lrm,weights=as.vector(w),data=data.lrm))

			alpha.vector <- -c(lrm.coef[1]+c(0,lrm.coef[S[j]:(S[j]+ncol(weight)-2)]),diff(lrm.coef[1:(S[j]-1)]))
			if(!is.null(var.list[[j]])) alpha.vector <- c(alpha.vector,-lrm.coef[(S[j]+ncol(weight)-1):length(lrm.coef)])  

			alpha[[j]] <- matrix(NA,nrow=ncol(weight),ncol=S[j]-1+length(var.list[[j]]))
			alpha[[j]][,1] <- alpha.vector[1:ncol(weight)]

			if(S[j]>=3) alpha[[j]][,2:(S[j]-1)] <- matrix(alpha.vector[(ncol(weight)+1):(ncol(weight)+S[j]-2)],nrow=ncol(weight),ncol=S[j]-2,byrow=TRUE)
			if(!is.null(var.list[[j]])) alpha[[j]][,S[j]:(S[j]+length(var.list[[j]])-1)] <- matrix(alpha.vector[(ncol(weight)+S[j]-1):length(alpha.vector)],
                                                                                                nrow=ncol(weight),ncol=length(var.list[[j]]),byrow=TRUE)
			#f there are jumps in y values, remplace alpha with 0
			if(length(miss.val)>0)
			{
				alpha.j <- matrix(NA,nrow=ncol(weight),ncol=max(y[,j])-1)

				alpha.j[,-miss.val] <- alpha[[j]]
				alpha.j[,miss.val] <- matrix(0,nrow=ncol(weight),ncol=length(miss.val))
				alpha[[j]] <- alpha.j
			}
		}
		param <- list("alpha"=alpha)
	}
    param
}

