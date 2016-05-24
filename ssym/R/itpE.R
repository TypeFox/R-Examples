itpE <-
function(vP,objeto){
	epsilon <- objeto$epsilon
	maxiter <- objeto$maxiter
	response <- objeto$y
	orig <- objeto$orig
	n <- length(response)
	pspm <- objeto$pspm
    penm <- objeto$penm
    pspp <- objeto$pspp
    penp <- objeto$penp
	p <- objeto$p
	qm <- objeto$qm
	q <- objeto$q
	l <- objeto$l
	l.mu.i <- objeto$l.mu.i
	l.phi.i <- objeto$l.phi.i
	l1.mu <- objeto$l1.mu
	l1.phi <- objeto$l1.phi
	theta_new <- vP
	v <- objeto$v
	tol <- 1
	cont <- 0
	if(attr(l1.phi,"link")=="logarithmic") efip <- solve(crossprod(pspp) + 2*penp)
	while(tol > epsilon){
		 theta <- theta_new
		 if(orig=="nonlinear"){
		   mu_work <-  objeto$mu(theta)
		 }else{mu_work <-  l.mu.i(pspm%*%theta[1:(p+sum(qm))])}
		 phi_work <- l.phi.i(pspp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
		 z_work <- (response-mu_work)/sqrt(phi_work)
		 v_work <- v(z_work)
		 tol_EM <- 1
		 theta_EM_new <- theta
		 cont_EM <- 0
		 while(tol_EM > epsilon){
		     theta_EM <- theta_EM_new
			 if(orig=="nonlinear"){
				mu_es <-  objeto$mu(theta_EM)
			 	pspm <- objeto$GradD(theta_EM)
			 }else{mu_es <-  l.mu.i(pspm%*%theta_EM[1:(p+sum(qm))])}
		     phi_es <- l.phi.i(pspp%*%theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
			 z_es <- (response - mu_es)/sqrt(phi_es)
			 pspmw <- pspm*matrix(v_work*(l1.mu(mu_es)^2)/phi_es,n,ncol(pspm))
			 thetam <- theta_EM[1:(p+sum(qm))] + solve(crossprod(pspmw,pspm) + penm)%*%(crossprod(pspm,v_work*l1.mu(mu_es)*z_es/sqrt(phi_es))-penm%*%theta_EM[1:(p+sum(qm))])
			 if(attr(l1.phi,"link")=="logarithmic") 
			 	thetap <- theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)] + efip%*%(crossprod(pspp,v_work*(z_es)^2 - 1) - 2*penp%*%theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
			 else
			 	thetap <- theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)] + solve(crossprod(pspp*matrix(l1.phi(phi_es)^2,n,ncol(pspp))) + 2*penp)%*%(crossprod(pspp,l1.phi(phi_es)*(v_work*(z_es)^2 - 1)) - 2*penp%*%theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
			 theta_EM_new <- c(thetam,thetap)
		 	 theta_pd <- ifelse(abs(theta_EM_new) <= epsilon,1,theta_EM_new)
		     tol_EM <- max(abs(theta_EM_new-theta_EM)/abs(theta_pd))
		     cont_EM <- cont_EM + 1
		     if(cont_EM > maxiter) stop("no convergence was obtained!!",call.=FALSE)
		 }
		 theta_new <- theta_EM_new
		 theta_pd <- ifelse(abs(theta) <= epsilon,1,theta)
	     tol <- max(abs(theta_new-theta)/abs(theta_pd))
		 tol <- ifelse(objeto$family=="Normal",epsilon/2,tol)
		 cont <- cont + 1
		 if(cont > maxiter) stop("no convergence was obtained!!",call.=FALSE)
    }
	theta_new
}
