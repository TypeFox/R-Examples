itpEC2 <-
function(vP,objeto){
	epsilon <- objeto$epsilon
	maxiter <- objeto$maxiter
	y <- objeto$y
	event <- objeto$event
	n <- length(y)
	pspm <- objeto$pspm
    penm <- objeto$penm
    pspp <- objeto$pspp
    penp <- objeto$penp
	p <- objeto$p
	qm <- objeto$qm
	q <- objeto$q
	l <- objeto$l
	theta_new <- vP
	v <- objeto$v
	vp <- objeto$vp	
	h <- objeto$h
	l.mu.i <- objeto$l.mu.i
	l.phi.i <- objeto$l.phi.i
	l1.mu <- objeto$l1.mu
	l1.phi <- objeto$l1.phi
	l2.mu <- objeto$l2.mu
	l2.phi <- objeto$l2.phi	
	tol <- 1
	cont <- 0
	Ltt <- matrix(0,ncol(pspm)+ncol(pspp),ncol(pspm)+ncol(pspp))
	Ut <- matrix(0,ncol(pspm)+ncol(pspp),1)
	while(tol > epsilon){
		theta <- theta_new
  		mu_es <-  l.mu.i(pspm%*%theta[1:(p + sum(qm))])
		phi_es <- l.phi.i(pspp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
		z_es <- (y - mu_es)/sqrt(phi_es)
		h_es <- ifelse(event==1,h(z_es),0)
		v_es <- v(z_es)
		vp_es <- ifelse(event==1,0,vp(z_es))
		Dc   <- ifelse(event==1,h_es*(h_es-v_es*z_es),vp_es*z_es + v_es)
		Dcdot   <- l1.mu(mu_es)*l2.mu(mu_es)*ifelse(event==1,h_es,v_es*z_es)
		Dcdot2   <- (1 + phi_es*phi_es*l1.phi(phi_es)*l2.phi(phi_es))*ifelse(event==1,h_es*z_es,v_es*z_es^2 - 1)/2
		Dcar <- ifelse(event==1,h_es + z_es*h_es*(h_es-v_es*z_es),vp_es*z_es^2 + 2*v_es*z_es)/2
		Dcab <- Dcar*z_es/2
		Ltt[1:ncol(pspm),1:ncol(pspm)] <- -crossprod(pspm,pspm*matrix(l1.mu(mu_es)^2*(Dc/phi_es + Dcdot/sqrt(phi_es)),n,ncol(pspm))) - penm
		Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)] <- -crossprod(pspp,pspm*matrix(l1.mu(mu_es)*l1.phi(phi_es)*(Dcar/sqrt(phi_es)),n,ncol(pspm)))
		Ltt[1:ncol(pspm),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- t(Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)])
		Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- -crossprod(pspp,pspp*matrix(l1.phi(phi_es)^2*(Dcab + Dcdot2),n,ncol(pspp))) - penp
		vt_es <- ifelse(event==1,h_es/z_es,v_es)
		s_es <- ifelse(event==1,h_es*z_es+1,v_es*z_es^2)		
		Ut[1:ncol(pspm)] <- crossprod(pspm,l1.mu(mu_es)*vt_es*z_es/sqrt(phi_es)) - penm%*%theta[1:(p+sum(qm))]
		Ut[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- crossprod(pspp,l1.phi(phi_es)*(s_es - 1)/2) - penp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)]
		theta_new <- theta + crossprod(solve(-Ltt),Ut)
		theta_pd <- ifelse(abs(theta_new) <= epsilon,1,theta_new)
		tol <- max(abs(theta_new-theta)/abs(theta_pd))
		cont <- cont + 1
		if(cont > maxiter) stop("no convergence was obtained!!",call.=FALSE)
    }
	theta_new
}
