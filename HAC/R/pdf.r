# pdf.r ##################################################################################################################
# FUNCTION:               	DESCRIPTION:
#  .dAC						Computes the values of the bivariate copula density. (Internal function)
#  .gumb.12.density			Bivariate density of the Gumbel copula. (Internal function)
#  .clay.12.density			Bivariate density of the Clayton copula. (Internal function)
#  .frank.12.density		    Bivariate density of the Frank copula. (Internal function)
#  .joe.12.density			Bivariate density of the Joe copula. (Internal function)
#  .amh.12.density			Bivariate density of the Ali-Mikhail-Haq copula. (Internal function)
#  .d.multi.AC				Computes the values of the multivariate Archimedean copula density. (Internal function)
#  dHAC						Returns the values of the an arbitrary HAC density.
#  .cop.pdf					Derives a function for the copula density or evalutes the derived function instantaneously. (Internal function)
#  .d.dell                  Derives the copula expression given by .constr.expr with respect to the arguments of the copula, which are defined on [0,1]. (Internal function)
#  .constr.expr             Returns an expression of the HAC for a given copula type. (Internal function)
#  to.logLik                Returns the log-Likelihood function or evalutes the log-likelihood instantaneously.
#  .tree.without.params     Tranforms a tree of a 'hac' object with numeric values as parameters into a tree with symbolic parameters. (Internal function)
#  .curviture				Fisher-Scoring contributions. (Internal function)
#  .score					Score contributions of a bivariate log-likelihood. (Internal function)
#  .score.amh				Score contributions of a bivariate log-likelihood for AMH. (Internal function)
#  .score.joe				Score contributions of a bivariate log-likelihood for Joe. (Internal function)
#  .score.frank				Score contributions of a bivariate log-likelihood for Frank. (Internal function)
#  .score.clayton			Score contributions of a bivariate log-likelihood for Clayton. (Internal function)
#  .score.gumbel			Score contributions of a bivariate log-likelihood for Gumbel. (Internal function)
##########################################################################################################################

.dAC = function(x, y, theta, type){	
	if((type == 1) | (type == 2)){
		.gumb.12.density(x, y, theta)
	}else 
	if((type == 3) | (type == 4)){
		.clay.12.density(x, y, theta)       
	}else 
	if((type == 5) | (type == 6)){
		.frank.12.density(x, y, theta)
	}else 
	if((type == 7) | (type == 8)){
		.joe.12.density(x, y, theta)
	}else 
	if((type == 9) | (type == 10)){
		.amh.12.density(x, y, theta)
	}
}

#-------------------------------------------------------------------------------------------------------------------------------

.gumb.12.density = function(x, y, theta){
	lu1 = -log(x)
	lu2 = -log(y)
	(lu1^(-1 + theta) * (-1 + theta + (lu1^theta + lu2^theta)^(1/theta)) * (lu1^theta + lu2^theta)^(-2 + 1/theta) * lu2^(-1 + theta))/(exp((lu1^theta + lu2^theta)^(1/theta)) *x*y)
}
	
#-------------------------------------------------------------------------------------------------------------------------------
	
.clay.12.density = function(x, y, theta){
	u1pt = x^(-theta)
	u2pt = y^(-theta)
	(1 + theta) * u1pt * u2pt * ((u1pt + u2pt - 1)^( -1/theta - 2))/(x * y)
}

#-------------------------------------------------------------------------------------------------------------------------------
	
.frank.12.density = function(x, y, theta){
	u1pt = exp(theta)
	(u1pt^(1+x+y)*(u1pt-1)*theta)/(u1pt-u1pt^(1+x)+u1pt^(x+y)-u1pt^(1+y))^2
}

#-------------------------------------------------------------------------------------------------------------------------------
	
.joe.12.density = function(x, y, theta){
	x1 = (1 - x)^theta
	y1 = (1 - y)^theta
	(1 - x)^(-1 + theta)*(theta - (-1 + x1)*(-1 + y1))*(x1 - (-1 + x1)*y1)^(-2 + 1/theta)*(1 - y)^(-1 + theta)
}

#-------------------------------------------------------------------------------------------------------------------------------
	
.amh.12.density = function(x, y, theta){
	(theta^2*(x+y-x*y-1)-theta*(x+y+x*y-2)-1)/(theta*(x-1)*(y-1)-1)^3
}

#-------------------------------------------------------------------------------------------------------------------------------

.d.multi.AC = function(X, theta, type, log = TRUE){	
            colnames(X) = c();
            if(type == 2){
                copGumbel@dacopula(X, theta, log = log)
            }else
            if(type == 4){
            	copClayton@dacopula(X, theta, log = log)
        		}else
            if(type == 6){
            	copFrank@dacopula(X, theta, log = log)
        		}else
            if(type == 8){
            	copJoe@dacopula(X, theta, log = log)
        		}else
            if(type == 10){
            	copAMH@dacopula(X, theta, log = log)
        		}
}
	
#-------------------------------------------------------------------------------------------------------------------------------

dHAC = function(X, hac, eval = TRUE, margins = NULL, na.rm = FALSE, ...){ 

    X = .one.ob(X, margins)
    names = colnames(X)    
    if(any(!(names %in% .get.leaves(hac$tree)))){stop("The colnames of X have to coincide with the specifications of the copula model hac.")}
			
	if(na.rm){X = na.omit(X, ...)}
    
    type = hac$type; d = NCOL(X);
        if((d >= 3) & ((type == 1) | (type == 3) | (type == 5) | (type == 7) | (type == 9))){
           return(.cop.pdf(tree = hac$tree, sample = X, type = type, d = d, names = names, eval = eval))
        }else{
        	if((type == 1) | (type == 3) | (type == 5) | (type == 7) | (type == 9)){type = type + 1}
			.d.multi.AC(X, theta = hac$tree[[length(hac$tree)]], type = type, log = FALSE)[-1]
        }
}

#-------------------------------------------------------------------------------------------------------------------------------

.cop.pdf = function(tree, sample, type, d, names, eval){
	if(is.null(colnames(sample))){stop("Specify colnames for X.")}
	string.expr = .constr.expr(tree, type)
    Dd = .d.dell(parse(text = string.expr), names, d)
    
    if(eval){
        for(i in 1:d){formals(Dd)[[i]]=sample[-1 ,i]}
        c(attr(Dd(), "gradient"))
    }else{
        Dd
    }
}

#-------------------------------------------------------------------------------------------------------------------------------

.d.dell = function(expr, names, order){
   if(order==1){
        deriv(expr, names[order], function.arg = names)
   }else{
        .d.dell(D(expr, names[order]), names, order-1)}
}

#---------------------------------------------------------------------------------------------------

.constr.expr = function(tree, type){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
 
     if(any(s)){
         if(any(!s)){
           if(type==1){
                 paste("exp(-(", paste("(-log(", unlist(tree[which(s)]),"))^", tree[[n]], collapse="+", sep = ""),"+", paste("(-log(",sapply(tree[which(!s)], .constr.expr, type=type),"))^", tree[[n]], collapse="+", sep = ""),")^(1/", tree[[n]],"))", sep="")
             }else
           if(type==3){
                 paste("(", paste("(", unlist(tree[which(s)]),"^(-", tree[[n]],")-1)", collapse="+", sep = ""),"+", paste("((", sapply(tree[which(!s)], .constr.expr, type=type),")^(-", tree[[n]],")-1)", collapse="+", sep = ""), "+1)^(-1/", tree[[n]], ")", sep="")
             }else
           if(type==5){
                 paste("-log(1-(1-exp(-", tree[[n]],"))*exp(", paste("log((exp(-(",unlist(tree[which(s)]),")*", tree[[n]],")-1)/(exp(-", tree[[n]],")-1))", collapse="+", sep = ""), "+", paste("log((exp(-(", sapply(tree[which(!s)], .constr.expr, type=type),")*", tree[[n]],")-1)/(exp(-", tree[[n]],")-1))", collapse="+", sep = ""),"))/",tree[[n]], sep="")
             }else
           if(type==7){
                 paste("1-(1-exp(", paste("log(1-(1-", unlist(tree[which(s)]),")^(", tree[[n]],"))", collapse="+", sep = ""), "+", paste("log(1-(1-(", sapply(tree[which(!s)], .constr.expr, type=type),"))^(", tree[[n]],"))", collapse="+", sep = ""), "))^(1/", tree[[n]], ")", sep="")
             }else
           if(type==9){
                  paste("(1-", tree[[n]],")/(", paste("((1-", tree[[n]],")/(", unlist(tree[which(s)]),")+", tree[[n]],")", collapse="*", sep = ""),"*", paste("((1-", tree[[n]],")/(", sapply(tree[which(!s)], .constr.expr, type=type), ")+", tree[[n]],")", collapse="*", sep = ""), "-", tree[[n]],")", sep="")
             }
 }else{
             if(type==1){
                 paste("exp(-(", paste("(-log(", unlist(tree[-n]),"))^", tree[[n]], collapse="+", sep = ""),")^(1/", tree[[n]],"))", sep="")
             }else
             if(type==3){
                 paste("(", paste("(",unlist(tree[-n]),"^(-", tree[[n]],")-1)", collapse="+", sep = ""), "+1)^(-1/", tree[[n]], ")", sep="")
             }else
             if(type==5){
                 paste("-log(1-(1-exp(-", tree[[n]],"))*exp(", paste("log((exp(-(",unlist(tree[-n]),")*", tree[[n]],")-1)/(exp(-", tree[[n]],")-1))", collapse="+", sep = ""), "))/",tree[[n]], sep="")
             }else
             if(type==7){
                 paste("1-(1-exp(", paste("log(1-(1-", unlist(tree[-n]),")^(", tree[[n]],"))", collapse="+", sep = ""), "))^(1/", tree[[n]], ")", sep="")
             }else
             if(type==9){
                 paste("(1-", tree[[n]],")/(", paste("((1-", tree[[n]],")/(", unlist(tree[-n]),")+", tree[[n]],")", collapse="*", sep = ""), "-", tree[[n]],")", sep="")
             }
 }}else{
             if(type==1){
                 paste("exp(-(", paste("(-log(", sapply(tree[-n], .constr.expr, type=type),"))^", tree[[n]], collapse="+", sep = ""),")^(1/", tree[[n]],"))", sep="")           
             }else
             if(type==3){
                 paste("(", paste("((", sapply(tree[-n], .constr.expr, type=type),")^(-", tree[[n]],")-1)", collapse="+", sep = ""), "+1)^(-1/", tree[[n]], ")", sep="")
             }else
           if(type==5){
                 paste("-log(1-(1-exp(-", tree[[n]],"))*exp(", paste("log((exp(-(", sapply(tree[-n], .constr.expr, type=type),")*", tree[[n]],")-1)/(exp(-", tree[[n]],")-1))", collapse="+", sep = ""),"))/",tree[[n]], sep="")
             }else
           if(type==7){
                 paste("1-(1-exp(",paste("log(1-(1-(", sapply(tree[-n], .constr.expr, type=type),"))^(", tree[[n]],"))", collapse="+", sep = ""), "))^(1/", tree[[n]], ")", sep="")
             }else
           if(type==9){
                  paste("(1-", tree[[n]],")/(", paste("((1-", tree[[n]],")/(", sapply(tree[-n], .constr.expr, type=type),")+", tree[[n]],")", collapse="*", sep = ""), "-", tree[[n]],")", sep="")
}}
}

#---------------------------------------------------------------------------------------------------

to.logLik = function(X, hac, eval = FALSE, margins = NULL, sum.log = TRUE, na.rm = FALSE, ...){
	X = .margins(X, margins)
			
	if(na.rm){X = na.omit(X, ...)}
    
    tree = .tree.without.params(hac$tree)
    thetas = .read.params(tree); values = get.params(hac); d = NCOL(X)
    expr = .constr.expr(tree, hac$type)
    f = .d.dell(parse(text=expr), c(colnames(X), thetas[order(values)]), order=d)
    for(i in 1:d){formals(f)[[i]]=X[,i]}
    
    g = function(theta, density=f){
            n.par = length(theta)
            for(i in 1:n.par){formals(density)[[length(formals(density))-n.par+i]]=theta[i]}
            if(sum.log){
	            sum(log(c(attr(density(), "gradient"))))
	        }else{
	        	log(c(attr(density(), "gradient")))
	        }
    }
        
    if(!eval){g}else{g(values[order(values)])}
}
 
#---------------------------------------------------------------------------------------------------
 
.tree.without.params = function(tree, k=1, l=1){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
     tree[[n]] = paste("theta",k,".",l, sep="")
     
     if(any(s)){
         if(any(!s)){
            for(i in which(!s)){
                tree[[i]]=.tree.without.params(tree[[i]], k=k+1,l=i)
            }       
         }else{
            tree = tree
         }}else{
         for(i in 1:(n-1)){
                tree[[i]]=.tree.without.params(tree[[i]], k=k+1,l=i)
         }}
    return(tree)        
}

#-------------------------------------------------------------------------------------------------------------------------------
#
#.d.dell.score = function(expr, names, order){
#   if(order==1){
#        deriv(parse(text = paste("log(",paste(deparse(expr), sep = "", collapse = ""), ")", sep = "", collapse = "")), names[order], function.arg = names)
#   }else{
#        .d.dell.score(D(expr, names[order]), names, order-1)}
#}

#-----------------------------------------------------------------------------------------------------------------

#.score = function(X, tree, type){
#	hac.names = .get.leaves(tree)
#	values = .read.params(tree);
	  
#	tree = .tree.without.params(tree)
#    thetas = .read.params(tree); 
#    expr = .constr.expr(tree, type)
#    f = .d.dell.score(parse(text=expr), c(thetas, hac.names), order = length(hac.names) + 1)
 	
# 	values = cbind(values, X)
#    for(i in 1:NCOL(values)){formals(f)[[i]]=values[,i]}
#    c(attr(f(), "gradient"))
#}

#-------------------------------------------------------------------------------------------------------------------------------

.curviture = function(X, theta, type){
   d = NCOL(X)
  .names = colnames(X)
  if(d < 3){
     #.AC = vector("list", d+1)
     #for(j in 1:d){.AC[[j]] = .names[[j]]}; .AC[[d+1]] = theta;
     #.score(X = X, tree = .AC, type = type)^2
     .score(X[,1], X[,2], theta, type + 1)^2
  }else{
      jacobian(function(r){.d.multi.AC(X = X, theta = r, type = type + 1)}, x = theta, method = "simple", method.args = list(eps = 1e-8))^2
  }
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.score = function(x, y, theta, type){
	        if(type == 2){
                .score.gumbel(x, y, theta)
            }else
            if(type == 4){
            	.score.clayton(x, y, theta)
       		}else
            if(type == 6){
            	.score.frank(x, y, theta)
        		}else
            if(type == 8){
            	.score.joe(x, y, theta)
        		}else
            if(type == 10){
            	.score.amh(x, y, theta)
        		}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.score.amh = function(x, y, theta){
	m.x = -1 + x
	m.y = -1 + y
	(-1 + x*(2 - 4*y) - theta^2*(-1 + x)^2*m.y^2 + 2*y - 2*theta*m.x*m.y*(m.x + y + x*y))/((-1 + theta*m.x*m.y)*(1 + theta*(m.x + m.y + theta*m.x*m.y + x*y)))
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.score.joe = function(x, y, theta){
	x.1 = (1-x)^theta
	y.1 = (1-y)^theta
	-((theta*((-(-1 + theta)^2 + x.1)*x.1 + ((-1 + theta)*theta + (2 + (-1 + theta)*theta)*x.1 - 2*x.1^2)*y.1 + (-1 + x.1)*(-theta + x.1)*y.1^2)*log1p(-x) - (-theta + (-1 + x.1)*(-1 + y.1))*(-x.1 + (-1 + x.1)*y.1)*log(x.1 - (-1 + x.1)*y.1) + theta*((-1 + x.1)^2*y.1^2*log(1 - y) + theta*x.1*(1 + (-1 + theta + x.1)*log1p(-y)) - (-1 + x.1)*y.1*(theta + (-(-1 + theta)^2 + (1 + theta)*x.1)*log1p(-y))))/(theta^2*(-theta + (-1 + x.1)*(-1 + y.1))*(x.1 - (-1 + x.1)*y.1)))
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.score.frank = function(x, y, theta){
	dif.xy = x-y
	sum.xy = x+y
(exp(theta*(2 + y))*(1 + theta*dif.xy) + exp(theta + theta*x)*(-1 + theta*(1 + dif.xy)) + exp(theta*(2 + x))*(1 - theta*dif.xy) + exp(theta + theta*y)*(-1 - theta*(-1 + dif.xy)) + exp(theta*(1 + sum.xy))*(-1 + theta*(-2 + sum.xy)) + exp(theta*sum.xy)*(1 - theta*(-1 + sum.xy)) + exp(theta)*(1 + theta*(-1 + sum.xy)) - exp(2*theta)*(1 + theta*sum.xy))/((-1 + exp(theta))*(-exp(theta*sum.xy) + exp(theta)*(-1 + exp(theta*x) + exp(theta*y)))*theta)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.score.clayton = function(x, y, theta){
	x.theta = x^theta 
	y.theta = y^theta
	theta.1 = 1 + theta
	(theta*theta.1*(-(theta*x.theta) + theta.1*y.theta)*log(x) + theta*theta.1*(theta.1*x.theta - theta*y.theta)*log(y) + (x.theta + y.theta)*(theta^2 + theta.1*log(x^(-theta) + y^(-theta))))/(theta^2*theta.1*(x.theta + y.theta))
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
.score.gumbel = function(x, y, theta){
	lx = (-log1p(-x))^theta
	ly = (-log1p(-y))^theta
	xy = lx + ly
	xy.t = xy^(1/theta)
(theta*(-(lx*((-1+theta)^2+(-3+2*theta)*xy.t+(lx+ly)^(2/theta)))+theta*(-1+theta+xy.t)*ly)*log(-log1p(-x))+(1-theta+(-3+theta)*xy.t+xy.t^2)*xy*log(xy)+theta*(theta*lx*(1+(-1 +theta+xy.t)*log(-log1p(-y)))+ly*(theta-((-1+theta)^2+(-3+2*theta)*xy.t+xy^(2/theta))*log(-log1p(-y)))))/(theta^2*(-1+theta+xy.t)*xy)
}