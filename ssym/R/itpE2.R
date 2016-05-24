itpE2 <-
function(vP,objeto){
	epsilon <- objeto$epsilon
	maxiter <- objeto$maxiter
	response <- objeto$y
	n <- length(response)
	orig <- objeto$orig
	pspm <- objeto$pspm
    penm <- objeto$penm
    pspp <- objeto$pspp
    penp <- objeto$penp
	family <- objeto$family
	p <- objeto$p
	qm <- objeto$qm
	q <- objeto$q
	l <- objeto$l
	xi <- objeto$xi	
	theta_new <- vP
	tol <- 1
	cont <- 0
	l1.mu <- objeto$l1.mu
	l1.phi <- objeto$l1.phi
	l.mu.i <- objeto$l.mu.i
	l.phi.i <- objeto$l.phi.i
	u <- function(z) (xi[2] + 1)/(xi[2] + 4*sinh(z)*sinh(z)/(xi[1]^2))
	v <- function(z,u){
		  4*sinh(z)*cosh(z)*u/(xi[1]^2*z) - tanh(z)/z
	}
	dg <- function(xis){
		  2 + 4/(xis^2) - (sqrt(2*pi)/xis)*(1-2*(pnorm(sqrt(2)/xis,mean=0,sd=sqrt(2)/2)-0.5))*exp(2/(xis^2))
	}
	fg <- function(xis){
		  (2.071 + xis*(1.936)      + xis^2*(-6.003e-02) + xis^3*(1.337e-03)  + xis^4*(-1.748e-05) + xis^5*(1.202e-07) + xis^6*(-3.343e-10))*(xis>2) + 
		  (3.001 + xis*(-1.030e-02) + xis^2*(1.050)      + xis^3*(-8.828e-02) + xis^4*(-1.859e-01) + xis^5*(8.899e-02) + xis^6*(-1.327e-02))*(xis<=2)
	}
	while(tol > epsilon){
		 theta <- theta_new
		 if(orig=="nonlinear"){
		   mu_work <-  objeto$mu(theta)
		 }else{mu_work <-  l.mu.i(pspm%*%theta[1:(p+sum(qm))])}
		 phi_work <- l.phi.i(pspp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
		 z_work <- (response-mu_work)/sqrt(phi_work)
		 tol_EM <- 1
		 theta_EM_new <- theta
		 cont_EM <- 0
		 dgs <- dg(xi[1]/sqrt(u(z_work)))
		 fgs <- fg(xi[1]/sqrt(u(z_work)))		 
		 while(tol_EM > epsilon){
		     theta_EM <- theta_EM_new
			 if(orig=="nonlinear"){
		       mu_es <-  objeto$mu(theta_EM)
		       pspm <- objeto$GradD(theta_EM)
		     }else{mu_es <-  l.mu.i(pspm%*%theta_EM[1:(p+sum(qm))])}
		     phi_es <- l.phi.i(pspp%*%theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
		     z_es <- (response-mu_es)/sqrt(phi_es)
			 pspmw <- pspm*matrix(dgs*l1.mu(mu_es)^2/phi_es,n,p+sum(qm))
			 v_es_work <- v(z_es,u(z_work))
			 thetam <- theta_EM[1:(p+sum(qm))] + solve(crossprod(pspm,pspmw) + penm)%*%(crossprod(pspm,v_es_work*(z_es)*l1.mu(mu_es)/sqrt(phi_es))-penm%*%theta_EM[1:(p+sum(qm))])
 			 psppw <- pspp*matrix(l1.phi(phi_es)^2*(fgs-1)/4,n,l+sum(q))
			 thetap <- theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)] + solve(crossprod(pspp,psppw) + penp)%*%(crossprod(pspp,l1.phi(phi_es)*(v_es_work*z_es^2 - 1)/2) - penp%*%theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
			 theta_EM_new <- c(thetam,thetap)
		 	 theta_pd <- ifelse(abs(theta_EM)<=epsilon,1,theta_EM)
		     tol_EM <- max(abs(theta_EM_new-theta_EM)/abs(theta_pd))
		     cont_EM <- cont_EM + 1
		     if(cont_EM > maxiter) stop("no convergence was obtained!!",call.=FALSE)
		 }
		 theta_new <- theta_EM_new
		 theta_pd <- ifelse(abs(theta)<=epsilon,1,theta)
	     tol <- max(abs(theta_new-theta)/abs(theta_pd))
		 cont <- cont + 1
		 if(cont > maxiter) stop("no convergence was obtained!!",call.=FALSE)
    }
	theta_new
}
