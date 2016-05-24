###Plots####

####DENSITY PLOTS######
###function to draw density plots and color the 95% equal tailed credible interval###
#if(FALSE){
#plot.posterior.density=function(chain, title, color){
#    plot(density(chain), main=title, ylab="", xlab="")
#     Qs=quantile(chain, probs=c(0.025, 0.975, 0.5))
#        x1 <- min(which(density(chain)$x >= Qs[1]))  
#x2 <- max(which(density(chain)$x <  Qs[2]))
#with(density(chain), polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), density=NA, col=color))
#abline(v=Qs[3], col="gray25") 
#}
#}


######"BOX" PLOTS########
###chain.list is a list MCMC output, either a vector for one parameter or matrix for several parameters###
plotHLSM.random.fit=function(fitted.model, parameter, burnin = 0, thin = 1){
	colors=c("navy", "orange", "darkmagenta")
#	niter = fitted.model$call$niter
#	if(is.null(burnin)){burnin = round(0.1*niter)}
#	if(is.null(thin)&) {thin = round(0.9*niter/200)}

	if(parameter=="Beta"){
		chain=getBeta(fitted.model, burnin=burnin, thin=thin)
		if(length(dim(chain)) == 3){
                p=dim(chain)[2]
                n.nets=dim(chain)[3] 
               }else{
			p = 1
			n.nets = dim(chain)[2] }
                x=c(1:n.nets)
		for(pie in 1:p){
			int.quantiles=matrix(NA, n.nets,5)
			for(kappa in 1:n.nets){
				if(p > 1){
 			    int.quantiles[kappa,]=quantile(chain[,pie,kappa],
						 probs=c(0.025, 0.25, 0.5, 0.75, 0.975)) 
				}else{ int.quantiles[kappa,]=quantile(chain[,kappa],
						 probs=c(0.025, 0.25, 0.5, 0.75, 0.975)) }
			}
			plot(x, int.quantiles[,3], pch=".",
				 ylim=c(min(int.quantiles), max(int.quantiles)),
				 col="white", xlab="Network",
				 ylab=expression(beta[pie]),
				 main=paste("Beta", pie))
			rect((x-0.05), int.quantiles[,1], (x+0.05),
				int.quantiles[,5], col=colors[3], border=NA)
			rect((x-0.15), int.quantiles[,2], (x+0.15), 
				int.quantiles[,4], col=colors[2], border=NA)
			segments((x-0.15), int.quantiles[,3], (x+0.15),
				int.quantiles[,3], lwd=2, col=colors[1])
    		}
	}
	if(parameter=="Intercept"){
		chain=getIntercept(fitted.model, burnin=burnin, thin=thin)
                p=1
                n.nets=dim(chain)[2]
                x=c(1:n.nets)
		int.quantiles=matrix(NA, n.nets,5)
		for(kappa in 1:n.nets){
			int.quantiles[kappa,]=quantile(chain[,kappa],
					 probs=c(0.025, 0.25, 0.5, 0.75, 0.975))}
			plot(x, int.quantiles[,3], pch=".", 
				ylim=c(min(int.quantiles), max(int.quantiles)),col="white",
				xlab="Network", ylab=expression(beta[0]), main=expression(beta[0]))
			rect((x-0.05), int.quantiles[,1], (x+0.05), int.quantiles[,5], col=colors[3], border=NA)
			rect((x-0.15), int.quantiles[,2], (x+0.15), int.quantiles[,4], col=colors[2], border=NA)
			segments((x-0.15), int.quantiles[,3], (x+0.15), int.quantiles[,3], lwd=2, col=colors[1])
    }
	
	if(parameter=="Alpha"){
		chain=getAlpha(fitted.model, burnin=burnin,thin=thin)
		x=1
		int.quantiles=matrix(NA, 1,5)
	        int.quantiles=quantile(chain, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
		plot(x, int.quantiles[3], pch=".",  ylim=c(min(int.quantiles),
			max(int.quantiles)), col="white", xlab="Treatment Effect",
			 ylab=expression(alpha), main=expression(alpha))
		rect((x-0.05), int.quantiles[1], (x+0.05), int.quantiles[5], col=colors[3], border=NA)
		rect((x-0.15), int.quantiles[2], (x+0.15), int.quantiles[4], col=colors[2], border=NA)
		segments((x-0.15), int.quantiles[3], (x+0.15), int.quantiles[3], lwd=2, col=colors[1])
    	}
}

###############

plotHLSM.fixed.fit=function(fitted.model, parameter, burnin =0, thin = 1){
	colors=c("navy", "orange", "darkmagenta")
#	if(is.null(burnin)){burnin = 0.1*fitted.model$call$niter}
#	if(is.null(thin)) {thin = round(0.9*fitted.model$call$niter/200)}
	if(parameter=="Beta"){
		chain=getBeta(fitted.model,burnin=burnin,thin=thin)
                if(is.null(dim(chain))==TRUE){p=1
		}else {p=dim(chain)[2]}
                x=c(1:p)
		int.quantiles=matrix(NA, p,5)
		if(p>1){
			for(kappa in 1:p){
				int.quantiles[kappa,]=quantile(chain[,kappa], 
					probs=c(0.025, 0.25, 0.5, 0.75, 0.975))}
		}else{
			kappa=1
			int.quantiles[kappa,]=quantile(chain, 
				probs=c(0.025, 0.25, 0.5, 0.75, 0.975))}
                      
		plot(x, int.quantiles[,3], pch=".", lab=c(p, 5, 7), 
			xlim=c(0.75, p+0.25), ylim=c(min(int.quantiles), max(int.quantiles)),
			col="white", xlab="Network", ylab=expression(beta), main=paste("Betas"))
		rect((x-0.05), int.quantiles[,1], (x+0.05), int.quantiles[,5], col=colors[3], border=NA)
		rect((x-0.15), int.quantiles[,2], (x+0.15), int.quantiles[,4], col=colors[2], border=NA)
		segments((x-0.15), int.quantiles[,3], (x+0.15), int.quantiles[,3], lwd=2, col=colors[1])
    	}

	if(parameter=="Intercept"){
		chain=getIntercept(fitted.model,burnin=burnin, thin=thin)
                x=1
		int.quantiles=quantile(chain, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
		plot(x, int.quantiles[3], pch=".", lab=c(1,5,7), 
			ylim=c(min(int.quantiles), max(int.quantiles)),
			col="white", xlab="Network", ylab=expression(beta[0]),main=expression(beta[0]))
		rect((x-0.05), int.quantiles[1], (x+0.05), int.quantiles[5], col=colors[3], border=NA)
		rect((x-0.15), int.quantiles[2], (x+0.15), int.quantiles[4], col=colors[2], border=NA)
		segments((x-0.15), int.quantiles[3], (x+0.15), int.quantiles[3], lwd=2, col=colors[1])
    }

	if(parameter=="Alpha"){
		chain=getAlpha(fitted.model,burnin=burnin,thin=thin)
		x=1
		int.quantiles=matrix(NA, 1,5)
		int.quantiles=quantile(chain, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
		plot(x, int.quantiles[3], pch=".", lab=c(1,5,7), 
			ylim=c(min(int.quantiles), max(int.quantiles)), col="white",
			xlab="Treatment Effect", ylab=expression(alpha), main=expression(alpha))
		rect((x-0.05), int.quantiles[1], (x+0.05), int.quantiles[5], col=colors[3], border=NA)
		rect((x-0.15), int.quantiles[2], (x+0.15), int.quantiles[4], col=colors[2], border=NA)
		segments((x-0.15), int.quantiles[3], (x+0.15), int.quantiles[3], lwd=2, col=colors[1])
    }
}



######LATENT SPACE POSTIONS########


plotLS = function(LS,xx,fitted.model, node.name = FALSE,nodenames = NULL){
	xcor = apply(LS[[xx]][,1,],1,mean)
	ycor = apply(LS[[xx]][,2,],1,mean)
	plot(xcor,ycor,main = paste('Network',xx),cex = 0.5,pch = 19,cex.lab = 1.5,col = 'red')
	if(node.name == TRUE){
		if(!is.null(nodenames)){
			text(xcor,ycor,lab = nodenames) }
	}
}

plotHLSM.LS = function(fitted.model,pdfname = NULL,burnin=0,thin=1,...){
	if(class(fitted.model)!='HLSM'){
		stop('fitted.model must be of class HLSM')}
#	niter = length(fitted.model$draws$Alpha)
#	if(is.null(burnin)){burnin=0.1*niter}
#	if(is.null(thin)){thin = round(0.9*niter/200)}
	LS = getLS(fitted.model,burnin=burnin, thin=thin)
	if(!is.null(pdfname)){pdf(file=paste(pdfname,'.pdf',sep=''))}else(dev.new(height = 10,width = 10))
	if(length(LS) > 5){
		mat = matrix(1:length(LS),ceiling(length(LS)/5),5,byrow = TRUE)
		layout(mat,widths = rep.int(2.5,ncol(mat)), heights = rep.int(1,nrow(mat)))
	}else{
		x1 = ceiling(length(LS)/2)
		par(mfrow = c(2,x1))
	}
	lapply(1:length(LS),function(x) plotLS(LS,x,fitted.model,...))
	if(!is.null(pdfname))dev.off()
}







        




    
