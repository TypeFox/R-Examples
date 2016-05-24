itpE3 <-
function(vP,objeto){
	epsilon <- objeto$epsilon
	maxiter <- objeto$maxiter
	fg <- objeto$fg
	dg <- objeto$dg	
	response <- objeto$y
	n <- length(response)
	orig <- objeto$orig
	pspm <- objeto$pspm
    penm <- objeto$penm
    pspp <- objeto$pspp
    penp <- objeto$penp
	p <- objeto$p
	l1.mu <- objeto$l1.mu
	l1.phi <- objeto$l1.phi
	l.mu.i <- objeto$l.mu.i
	l.phi.i <- objeto$l.phi.i
	qm <- objeto$qm
	q <- objeto$q
	l <- objeto$l
	theta_new <- vP
	v <- objeto$v
	tol <- 1
	cont <- 0
	if(attr(l1.phi,"link")=="logarithmic") efip <- solve(((fg-1)/4)*crossprod(pspp) + penp)
	while(tol > epsilon){
		     theta <- theta_new
			 if(orig=="nonlinear"){
				mu_es <-  objeto$mu(theta)
			 	pspm <- objeto$GradD(theta)
			 }else{mu_es <-  l.mu.i(pspm%*%theta[1:(p+sum(qm))])}
		     phi_es <- l.phi.i(pspp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
			 z_es <- (response - mu_es)/sqrt(phi_es)
			 v_es <- v(z_es)
			 pspmw <- pspm*matrix(dg*l1.mu(mu_es)^2/phi_es,n,p+sum(qm))
			 thetam <- theta[1:(p+sum(qm))] + (1/2)*solve(crossprod(pspm,pspmw) + penm)%*%(crossprod(pspm,l1.mu(mu_es)*v_es*(z_es)/sqrt(phi_es))-penm%*%theta[1:(p+sum(qm))])
		 	 if(attr(l1.phi,"link")=="logarithmic")
				 thetap <- theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)] + (1/2)*efip%*%(crossprod(pspp,(v_es*(z_es)^2 - 1)/2) - penp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
			 else
				 thetap <- theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)] + (1/2)*solve(((fg-1)/4)*crossprod(pspp*matrix(l1.phi(phi_es)^2,n,ncol(pspp))) + penp)%*%(crossprod(pspp,l1.phi(phi_es)*(v_es*(z_es)^2 - 1)/2) - penp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
			 theta_new <- c(thetam,thetap)
		 	 theta_pd <- ifelse(abs(theta)<=epsilon,1,theta)
		     tol <- max(abs(theta_new-theta)/abs(theta_pd))
			 cont <- cont + 1
		     if(cont > maxiter) stop("no convergence was obtained!!",call.=FALSE)
	}
	theta_new
}
