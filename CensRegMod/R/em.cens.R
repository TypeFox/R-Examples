em.cens <- function(cc,x,y,nu="NULL",dist="Normal",diagnostic="FALSE",typediag="NULL"){

  type=dist
  
	if(ncol(as.matrix(y)) > 1) stop("Only univariate linear regression supported!")
	if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
	if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
	if( (type != "T") && (type != "Normal") ) stop("Distribution family not supported. Check documentation!") 
	if( type == "T" )
  {
  	if(length(nu) > 1) stop("Nu parameters must have only one parameter")
  	if(length(nu) == 0) stop("Nu parameters must be provided.")
  	if(nu <= 0) stop("Nu parameters must be positive.")
  } 
	if( (diagnostic=="TRUE") && (length(typediag)==0) ) stop("typediag parameter must be provided.")
	if( (diagnostic=="TRUE") && ( (typediag!=1) && (typediag!=2) && (typediag!=3) )) stop("typediag must be 1, 2 or 3")

	p <- ncol(x)
	n <- nrow(x)
	reg <- lm(y ~ x[,2:p])

	ERRO<-1e-6
	TOLERANCIA<-1e-6
	MAX_NU<-150
	MIN_NU <- 1.01
	
	beta<- as.vector(coefficients(reg),mode="numeric")
	sigma2 <- sum((y-x%*%beta)^2)/(n-p)
	
	if(type=="T")
  {
		teta_velho <- matrix(c(beta,sigma2,nu),ncol=1)
 		cont <- 0
	
		repeat
    {
      mu<-x%*%beta
			d <- ((y-mu)^2)/sigma2
		
			aux1t<-CalMom(mu,sigma2,nu,y,type)
      		
			u0<- (nu+1)/(nu+d)
  		u1<- y*(nu+1)/(nu+d) 
      u2<- y^2*(nu+1)/(nu+d)
		
			u1[cc==1]<- aux1t$Ey1[cc==1]
      u2[cc==1]<- aux1t$Ey2[cc==1]
      u0[cc==1]<- aux1t$Ey0[cc==1]           									 		
				           			      
  		suma1<-matrix(0,p,p)
      suma2<-matrix(0,p,1)
      		
  		for(i in 1:n)
      {
  		  suma1<-suma1+u0[i]*((x[i,])%*%t(x[i,]))
      	suma2<-suma2+(x[i,])*u1[i]
     	}

  		beta<-solve(suma1)%*%suma2
  		sigma2<-sum(u2-2*u1*mu+mu^2*u0)/n
		 
			auxf<-(y-x%*%beta)/sqrt(sigma2)

			if(sum(cc)>0)
			{
			  ft <- function(nu){sum(log(dt(auxf[cc==0],df=nu)/sqrt(sigma2)))+ sum(log(pt(-auxf[cc==1],df=nu))) }
			}else{
			  ft <- function(nu){sum(log(dt(auxf[cc==0],df=nu)/sqrt(sigma2)))}
			}
			nu <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum 
      
    	teta_novo<-matrix(c(beta,sigma2,nu),ncol=1)

			if (sqrt(sum((teta_velho-teta_novo)^2)) < ERRO ){break}
		
   		teta_velho <- teta_novo
		}
  
		suma3 <- IM.CensGrad(cc=cc,x=x,y=y,BETA=beta,SIGMA2=sigma2,nu=nu,type=type)
		logver<-ft(nu)
		AIC <- (-2)*logver + 2*(p+2)
		BIC <- (-2)*logver + (p+2)*log(n)
		EDC <- (-2)*logver + (p+2)*0.2*sqrt(n)
	}
	
	
	if(type=="Normal")
  {
	  teta_velho <- matrix(c(beta,sigma2),ncol=1)
 		cont <- 0
	
		repeat
    {
      mu<-x%*%beta
      u1<-y 
      u2<-y^2  
      u0<-rep(1,n)       
      aux1<-CalMom(mu=mu,sigma2=sigma2,a=y,type="Normal")

      u1[cc==1]<- aux1$Ey1[cc==1]
      u2[cc==1]<- aux1$Ey2[cc==1]
      u0[cc==1]<- aux1$Ey0[cc==1]

			suma1<-matrix(0,p,p)
      suma2<-matrix(0,p,1)

  		for(i in 1:n)
      {
  		  suma1<-suma1+u0[i]*((x[i,])%*%t(x[i,]))
        suma2<-suma2+(x[i,])*u1[i]
  		}

			beta<-solve(suma1)%*%suma2
  		sigma2<-sum(u2-2*u1*mu+mu^2*u0)/n

  		teta_novo<-matrix(c(beta,sigma2),ncol=1)

			if (sqrt(sum((teta_velho-teta_novo)^2)) < ERRO ){break}
		
			teta_velho <- teta_novo	
		} 
    
		suma3 <- IM.CensGrad(cc=cc,x=x,y=y,BETA=beta,SIGMA2=sigma2,nu=nu,type=type) 
		auxpdf<-dnorm(y,x%*%beta,sqrt(sigma2))
		auxcdf<-pnorm(-(y-x%*%beta)/sqrt(sigma2))
		logver<-sum(log(auxpdf[cc==0]))+sum(log(auxcdf[cc==1]))
		AIC <- (-2)*logver + 2*(p+1)  
		BIC <- (-2)*logver + (p+1)*log(n)
		EDC <- (-2)*logver + (p+1)*0.2*sqrt(n)
	}
  
  
  namesx <- ('x1     ')
  if(ncol(as.matrix(x))>1)
  {
    for(i in 2:ncol(as.matrix(x))){namesx <- cbind(namesx, paste("x",i,"     ",sep=""))}
  }
	SE <- t(suma3)
	SE <- round(t(SE),digits=5)
	param <- round(cbind(rbind(beta,sigma2),SE),digits=5)
	namespar <- colnames(x) 
	colx <- ncol(as.matrix(x))
	if(length(namespar)==0)namespar <- namesx[1:colx]
	dimnames(param) <- list(c(namespar,expression(sigma^2)),c("Estimates", "SE"))
	if( (type=="T"))
	{
	  nu1 <- matrix(round(nu,digits=5),ncol=1,nrow=1)
	  row.names(nu1) <- "nu"
	  colnames(nu1) <- " "
	}
	cat('\n') 
	cat('-------------------------------------------\n')
	cat('         EM estimates and SE \n')
	cat('-------------------------------------------\n')
	print(param)
	if(type=="T")
	{
	  print(nu1)
	}
	cat('------------------------------------------\n')
	cat('\r \n')
	critFin <- c(logver, AIC, BIC, EDC)
	critFin <- round(t(as.matrix(critFin)),digits=3)
	dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","EDC"))
	cat('\n') 
	cat('Model selection criteria\n')
	cat('-------------------------------------------\n')
	print(critFin)
	cat('-------------------------------------------\n')
	cat('\r \n')
  
	if(diagnostic=="TRUE")
  {	
		obs<-c(rep(0,n))
		
		Matriz<-HessianaQ(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2)
 		MatrizQ<-Matriz$MatrizQ
 		MatrizQbeta<-Matriz$MatrizQbeta
 		MatrizQsigma2<-Matriz$MatrizQsigma2

		if(typediag==1)
    {

			CDm<-CD(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2,MatrizQ=MatrizQ,MatrizQbeta=MatrizQbeta,MatrizQsigma2=MatrizQsigma2)
			GD<-CDm$GD
			GDbeta<-CDm$GDbeta
			GDsigma2<-CDm$GDsigma2

			medida <- data.frame(GD,GDbeta,GDsigma2)

			for(i in 1:n)
      {
				if(GD[i] > (2*(p+1)/n)){obs[i]<-i}
			}
			influentGD<-obs[obs!=0]

			obs<-c(rep(0,n))
			for(i in 1:n)
      {
				if(GDbeta[i]> (2*p/n) ){obs[i]<-i}
			}
			influentGDbeta<-obs[obs!=0]

			obs<-c(rep(0,n))
			for(i in 1:n)
      {
				if(GDsigma2[i]>(2/n)){obs[i]<-i}
			}
			influentGDsigma2<-obs[obs!=0]

			if(length(influentGD)==0){influentGD<-NULL}
			if(length(influentGDbeta)==0){influentGDbeta<-NULL}
			if(length(influentGDsigma2)==0){influentGDsigma2<-NULL}

			influent <- list(GD=influentGD,GDbeta=influentGDbeta,GDsigma2=influentGDsigma2)
		}

		if(typediag==2)
    {
			If1<-If(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2,MatrizQ=MatrizQ,Perturb=1)
			
			medida<-diag(If1)/sum(diag(If1))
			benchmark<-mean(medida)+3.5*sd(medida)
			
			for(i in 1:n)
      {
				if(medida[i]>benchmark){obs[i]<-i}
			}
			influent<-obs[obs!=0]
		}

		if(typediag==3)  
    {

			If2<-If(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2,MatrizQ=MatrizQ,Perturb=2)
			medida<-diag(If2)/sum(diag(If2))
			benchmark<-mean(medida)+3.5*sd(medida)
			
			for(i in 1:n)
      {
				if(medida[i]>benchmark){obs[i]<-i}
			}
			influent<-obs[obs!=0]
		}

    if(length(influent)==0){influent<-NULL}

    
		if(type=="T")
    {
			return(list(betas=beta,sigma2=sigma2,nu=nu,logver=logver,SE=suma3,influents=influent,measure=medida,AIC=AIC,BIC=BIC,EDC=EDC))
		}

		if(type=="Normal")
    {
			return(list(betas=beta,sigma2=sigma2,logver=logver,SE=suma3,influents=influent,measure=medida,AIC=AIC,BIC=BIC,EDC=EDC))
		}

	}else{
		if(type=="T")
    {
			return(list(betas=beta,sigma2=sigma2,nu=nu,logver=logver,SE=suma3,AIC=AIC,BIC=BIC,EDC=EDC  ))
		}
		if(type=="Normal")
    {
			return(list(betas=beta,sigma2=sigma2,logver=logver,SE=suma3,AIC=AIC,BIC=BIC,EDC=EDC  ))
		}
	}
  
  
  
}













