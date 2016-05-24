



For.MC3.REG<- function(i, g, Ys, Xs, PI, K, nu, lambda, phi, outs.list)
{

    if (g$flag == 1) {
        if (sum(g$M0.var) != 0) 
            g$M0.1 <- sum(2^((0:(length(g$M0.var) - 1))[g$M0.var])) + 1
        else g$M0.1 <- 1
        if (sum(g$M0.out) != 0) 
            g$M0.2 <- sum(2^((0:(length(g$M0.out) - 1))[g$M0.out])) + 1
        else g$M0.2 <- 1
    }
    
    M1 <- MC3.REG.choose(g$M0.var, g$M0.out)


    if (sum(M1$var) != 0) 
        M1.1 <- sum( 2^(  (0:(length(g$M0.var) - 1)) [M1$var]) ) + 1
    else M1.1 <- 1
    if (sum(M1$out) != 0) 
        M1.2<- sum(2^((0:(length(g$M0.out) - 1))[M1$out])) + 1
    else M1.2 <- 1


    if (sum(g$big.list[, 1] == M1.1 & g$big.list[, 2] == M1.2) == 
        0) 
    {


        if (M1.1 == 1)			#null model
        {
            if (g$outcnt != 0) 
                a <- (dim(Ys)[1] - sum(M1$out)) * log(1 - PI) + sum(M1$out) * log(PI) + MC3.REG.logpost(Ys, Xs, 0, 0, outs.list[M1$out], K, nu, lambda, 
                  phi)
            else a <- MC3.REG.logpost(Ys, Xs, 0, 0, outs.list[M1$out], 
                K, nu, lambda, phi)
        }
        else {

            if (g$outcnt != 0) 
                a <- (dim(Ys)[1] - sum(M1$out)) * log(1 - PI) + 
                  sum(M1$out) * log(PI) + MC3.REG.logpost(Ys, 
                  Xs, M1$var, sum(M1$var), outs.list[M1$out], 
                  K, nu, lambda, phi)
            else a <- MC3.REG.logpost(Ys, Xs, M1$var, sum(M1$var), 
                outs.list[M1$out], K, nu, lambda, phi)
        }



        g$big.list<- rbind(g$big.list, c(M1.1, M1.2, a, 0))
    }

    BF <- exp(g$big.list[g$big.list[, 1] == M1.1 & g$big.list[, 2] == 
        M1.2, 3] - g$big.list[g$big.list[, 1] == g$M0.1 & g$big.list[, 
        2] == g$M0.2, 3])

#print("")
#print("g = ")
#print(g)

    if (BF >= 1) 
        g$flag <- 1
    else g$flag <- rbinom(1, 1, BF)
    if (g$flag == 1) {
        g$M0.var <- M1$var
        g$M0.out <- M1$out
        g$M0.1 <- M1.1
        g$M0.2 <- M1.2
    }
    g$big.list[g$big.list[, 1] == g$M0.1 & g$big.list[, 2] == g$M0.2, 4]<- g$big.list[g$big.list[, 
        1] == g$M0.1 & g$big.list[, 2] == g$M0.2, 4] + 1
    return(g)
}







MC3.REG<- function(all.y, all.x, num.its, M0.var = NULL, M0.out = NULL, outs.list = NULL, outliers = TRUE,
		   PI=.1*(length(all.y) <50) + .02*(length(all.y) >= 50), 
		   K=7, nu= NULL, lambda= NULL, phi= NULL)
{


    cl <- match.call()

    all.x<- data.frame(all.x)

    if (is.null(M0.var))
	M0.var<- rep(TRUE, ncol(all.x))
    if ((sum(M0.var) == 0)) 
        stop("\nInput error:  M0.var cannot be null model")

    if (outliers)
    {
    	if (is.null(outs.list))
    	{
    		outs.list<- out.ltsreg(all.x, all.y, 2)
	}

	if (length(outs.list) == 0)
		outliers<- FALSE
    }


    if (outliers)
    {
    	if (is.null(M0.out))
    	{
    		M0.out<- rep(TRUE, length(outs.list))
	}

    }


    if ((length(M0.out) != length(outs.list)) || (length(M0.var) != 
        dim(all.x)[2])) 
        stop("\nInput error: M0.*** is not the right length")



    # calculate R for full model to determine defaults for nu, lambda and phi

    if (is.null(nu) | is.null(lambda) | is.null(phi) )
    {
     	r2<- summary(lm(all.y ~ ., data = data.frame(all.y, all.x)))$r.squared
	if (r2 < 0.9)
	{
		new.nu<- 2.58
		new.lambda<- 0.28
		new.phi<- 2.85
	}
	else
	{
		new.nu<- 0.2
		new.lambda<- 0.1684
		new.phi<- 9.20
	}
	if (is.null(nu)) nu<- new.nu
	if (is.null(lambda)) lambda<- new.lambda
	if (is.null(phi)) phi<- new.phi
    }	

    var.names<- colnames(all.x)
    var.numbers<- 1:ncol(all.x)
    outlier.numbers<- outs.list
    

    Ys <- scale(all.y)
    Xs <- scale(all.x)

    #global variables
    g<- list()

    g$flag<- 1
    g$M0.var<- M0.var
    g$M0.out<- M0.out

    if (is.null(outs.list)) g$outcnt<- 0
    else g$outcnt<- sum(outs.list)

    g$big.list<- matrix(0, 1, 4)
    g$big.list[1, 1]<- sum(2^((0:(length(g$M0.var) - 1))[g$M0.var])) + 1
    if (sum(g$M0.out) != 0) 
        g$big.list[1, 2]<- sum(2^((0:(length(g$M0.out) - 1))[g$M0.out])) + 1
    else g$big.list[1, 2]<- 1

    if (g$outcnt != 0) 
        g$big.list[1, 3] <- (dim(Ys)[1] - sum(g$M0.out)) * log(1 - 
            PI) + sum(g$M0.out) * log(PI) + MC3.REG.logpost(Ys, 
            Xs, g$M0.var, sum(g$M0.var), outs.list[g$M0.out], K, nu, 
            lambda, phi)
    else g$big.list[1, 3] <- MC3.REG.logpost(Ys, Xs, g$M0.var, sum(g$M0.var), 
        outs.list[g$M0.out], K, nu, lambda, phi)


    for(i in 1:num.its) g<- For.MC3.REG(i, g, Ys, Xs, PI, K, nu, lambda, phi, outs.list)




    var.matrix<- matrix(as.logical(rep(g$big.list[, 1] - 1, rep(length(g$M0.var), 
        length(g$big.list[, 1])))%/%2^(0:(length(g$M0.var) - 1))%%2), 
        ncol = length(g$M0.var), byrow = TRUE)
    n.var <- length(g$M0.var)
    ndx <- 1:n.var

    Xn <- rep("X", n.var)

    labs <- paste(Xn, ndx, sep = "")




    colnames(var.matrix) <- var.names



    postprob<- exp(g$big.list[, 3])/(sum(exp(g$big.list[, 3])))
    visits <- g$big.list[, 4]




    if (length(outs.list) != 0) 
    {
        out.matrix <- matrix(as.logical(rep(g$big.list[, 2] - 1, 
            rep(length(outs.list), length(g$big.list[, 2])))%/%2^(0:(length(outs.list) - 
            1))%%2), ncol = length(outs.list), byrow = TRUE)
        colnames(out.matrix) <- outs.list
        
    }
    else out.matrix<- NULL 



    ordr<- order(-postprob)

    result<- list(post.prob = postprob[ordr], 
		  variables = var.matrix[ordr, ,drop=FALSE], 
		  outliers = out.matrix[ordr, ,drop=FALSE], 
		  visit.count = visits[ordr],
		  var.numbers = var.numbers, 
		  outlier.numbers = outlier.numbers,
		  var.names = var.names,
		  n.models = length(postprob[ordr]),
		  PI = PI, K=K, nu=nu, lambda=lambda, phi=phi,
		  call = cl)


    class(result)<- "mc3"
    return(result)
}






MC3.REG.choose<-function(M0.var,M0.out) 
{
    var <- M0.var
    in.or.out <- sample(c(1:length(M0.var), rep(0, length(M0.out))), 
        1)
    if (in.or.out == 0) {
        out <- M0.out
        in.or.out2 <- sample(1:length(M0.out), 1)
        out[in.or.out2] <- !M0.out[in.or.out2]
    }
    else {
        var[in.or.out] <- !M0.var[in.or.out]
        out <- M0.out
    }
    return(list(var=var, out=out))
}




MC3.REG.logpost<- function(Y, X, model.vect, p, i, K, nu, lambda, phi)
{
#print(list(Y=Y,X=X,model.vect=model.vect, p=p, i=i, K=K, nu=nu, lambda=lambda,phi=phi))

    n <- dim(Y)[1]
    ones <- rep(1, n)
    A <- cbind(ones, X[, model.vect])
    V <- diag(c(1, rep(phi^2, p)))
    ones[i] <- K^2
    det <- diag(ones) + A %*% V %*% t(A)
    divs <- prod(eigen(det, TRUE, TRUE)$values)^0.5
    denom <- (t(Y) %*% solve(det, Y)) + (nu * lambda)
    lgamma((n + nu)/2) + log(nu * lambda) * (nu/2) - log(pi) * 
        (n/2) - lgamma(nu/2) - log(divs) - ((nu + n)/2) * log(denom)
}




out.ltsreg <- function(x,y,delta)  
{
#   require(rrcov)
   abc<-ltsReg(x,y)$residuals
   (1:length(y))[abs(as.vector(abc)/mad(abc))>=delta]
}


print.mc3<- function(x, digits = max(3, getOption("digits") - 3), n.models = nrow(x$variables), ...)
{
	mc3.out<- x
	n.best<- n.models
 	cat("\nCall:\n", deparse(mc3.out$call), "\n\n", sep = "")
	cat(paste("Model parameters: PI = ", mc3.out$PI, 
		  " K = ", mc3.out$K, 
		  " nu = ",mc3.out$nu, 
		  " lambda = ",mc3.out$lambda,
		  " phi = ", mc3.out$phi, "\n",sep=""))
    	cat("\nModels visited: \n")
	
        n.var<- ncol(mc3.out$variables)
	n.out<- ncol(mc3.out$outliers)
	n.rows<- n.best
	pretty.var<- apply(mc3.out$variables+0,1,paste,sep=" ",collapse=" ")

	if (!is.null(mc3.out$outliers))
	{
		pretty.out<- apply(mc3.out$outliers+0, 1,paste,sep=" ",collapse=" ")
	
		pretty<- format(data.frame(posterior=mc3.out$post.prob, 
			    n.visits= mc3.out$visit.count, 
			    variables = pretty.var, outliers = pretty.out, row.names=NULL),
		        digits=digits, row.names=NULL)
	}
	else

		pretty<- format(data.frame(posterior=mc3.out$post.prob, 
			    n.visits= mc3.out$visit.count, 
			    variables = pretty.var, row.names=NULL),
		        digits=digits, row.names=NULL)


	print.default(pretty[1:n.rows,], ...)
}


summary.mc3<- function(object,  n.models = 5, digits = max(3, getOption("digits") - 3), ...)
{

	mc3.out<- object

 	cat("\nCall:\n", deparse(mc3.out$call), "\n\n", sep = "")
	cat(paste("Model parameters: PI = ", mc3.out$PI, 
		  " K = ", mc3.out$K, 
		  " nu = ",mc3.out$nu, 
		  " lambda = ",mc3.out$lambda,
		  " phi = ", mc3.out$phi, "\n",sep=""))
   	n.models <- min(n.models, mc3.out$n.models)
    	sel <- 1:n.models
    	cat("\n ", mc3.out$n.models, " models were selected")
    	cat("\n Best ", n.models, " models (cumulative posterior probability = ", 
        round(sum(mc3.out$post.prob[sel]), digits), "): \n\n")

	# calculate marginal posteriors
	var.probs<- mc3.out$post.prob %*% mc3.out$variables



	if (!is.null(mc3.out$outliers))
	{

	out.probs<- mc3.out$post.prob %*% mc3.out$outliers

	marginals<- rbind(format(t(var.probs),digits=digits),"",format(t(out.probs),digits=digits))

	probs<- format(mc3.out$post.prob[sel], digits=digits)
#	vars<- format(t(mc3.out$variables[sel, ,drop=FALSE]) + 0,digits=1)
#	outs<- format(t(mc3.out$outliers[sel,]) + 0,digits=1)

	vars<- matrix(".", ncol=ncol(mc3.out$variables[sel, ,drop=FALSE]), nrow=nrow(mc3.out$variables[sel, ,drop=FALSE]))
	vars[mc3.out$variables[sel, ,drop=FALSE]]<- "x"
	vars<- t(vars)
	outs<- matrix(".", ncol=ncol(mc3.out$outliers[sel, ,drop=FALSE]), nrow=nrow(mc3.out$outliers[sel, ,drop=FALSE]))
	outs[mc3.out$outliers[sel, ,drop=FALSE]]<- "x"
	outs<- t(outs)


    	decpos <- nchar(unlist(strsplit(probs[1], "\\."))[1])
    	offset2 <- paste(rep(" ", times = decpos ), sep = "", collapse = "")

	# now loop through vars and outs pasting offset, since R 'paste' does not do vectors :-(
	for (i in 1:ncol(vars))
		vars[,i]<- paste(offset2,vars[,i],sep="")
	for (i in 1:ncol(outs))
		outs[,i]<- paste(offset2,outs[,i],sep="")

	seprow<- rep("",times=ncol(vars))

	rght<- rbind(seprow,vars,seprow,outs,seprow,probs)
	lft<- rbind("",marginals,"","")

	all<- cbind(lft,rght)
	colnames(all)<- c("prob", paste("model ",1:n.models,sep=""))
	rownames(all)<- c("variables",
			  paste(" ",mc3.out$var.names),
			  "outliers", 
			  paste(" ", mc3.out$outlier.numbers),
		          "","post prob")

	}
	else
	{

	marginals<- rbind(format(t(var.probs),digits=digits))

	probs<- format(mc3.out$post.prob[sel], digits=digits)
#	vars<- format(t(mc3.out$variables[sel,]) + 0,digits=1)

	vars<- matrix(".", ncol=ncol(mc3.out$variables[sel,]), nrow=nrow(mc3.out$variables[sel,]))
	vars[mc3.out$variables[sel,]]<- "x"
	vars<- t(vars)

    	decpos <- nchar(unlist(strsplit(probs[1], "\\."))[1])
    	offset2 <- paste(rep(" ", times = decpos ), sep = "", collapse = "")

	# now loop through vars and outs pasting offset, since R 'paste' does not do vectors :-(
	for (i in 1:ncol(vars))
		vars[,i]<- paste(offset2,vars[,i],sep="")

	seprow<- rep("",times=ncol(vars))

	rght<- rbind(seprow,vars,seprow,probs)
	lft<- rbind("",marginals,"","")

	all<- cbind(lft,rght)
	colnames(all)<- c("prob", paste("model ",1:n.models,sep=""))
	rownames(all)<- c("variables",
			  paste(" ",mc3.out$var.names),
		          "","post prob")
	}

        print.default(all, print.gap = 2, quote = FALSE, ...)

}







"[.mc3"<- function(x,  ... )
{
	as.data.frame.mc3(x)[...]
}


as.data.frame.mc3<- function(x, ...)
	{

	y<- list()
	y$post.prob<- x$post.prob
	y$variables<- x$variables + 0
	y$outliers<- x$outliers + 0
	y$visit.count<- x$visit.count

	colnames(y$variables)<- x$var.names
	colnames(y$outliers)<- x$outlier.numbers

	outliers<- !is.null(dim(y$outliers))

	if (outliers)
		yy<- data.frame(y$post.prob, 
		       y$visit.count,
		       y$variables,
		       y$outliers, ...)
	else
		yy<- data.frame(y$post.prob, 
		       y$visit.count,
		       y$variables, ...)

	nms<- c("post.prob","visit.count",x$var.names)
	if (outliers)
		nms<- c(nms, x$outlier.numbers)
	names(yy)<- nms
	
	return(yy)
}






