itpEC <-
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
	l.mu.i <- objeto$l.mu.i
	l.phi.i <- objeto$l.phi.i
	l1.mu <- objeto$l1.mu
	l1.phi <- objeto$l1.phi
	l2.mu <- objeto$l2.mu
	l2.phi <- objeto$l2.phi	
	p <- objeto$p
	qm <- objeto$qm
	q <- objeto$q
	l <- objeto$l
	theta_new <- vP
	v <- objeto$v
	vb <- objeto$vb
	h <- objeto$h	
	tol <- 1
	cont <- 0
	while(tol > epsilon){
		 theta <- theta_new
  		 mu_work <-  pspm%*%theta[1:(p+sum(qm))]
		 phi_work <- exp(pspp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
		 z_work <- (y-mu_work)/sqrt(phi_work)
		 v_work <- ifelse(event==1,vb(z_work),v(z_work))
		 y_work <- ifelse(event==1,mu_work + h(z_work)*sqrt(phi_work)/v_work,y)
		 mt_work <- ifelse(event==1,mu_work^2 + sqrt(phi_work)*(h(z_work)*(y + mu_work) + sqrt(phi_work))/v_work,0)		 
		 theta_EM_new <- theta
		 tol_EM <- 1
		 cont_EM <- 0
		 Ltt <- matrix(0,ncol(pspm)+ncol(pspp),ncol(pspm)+ncol(pspp))
	     Ut <- matrix(0,ncol(pspm)+ncol(pspp),1)
		 while(tol_EM > epsilon){  
		 	 theta_EM <- theta_EM_new
  		     mu_es <-  l.mu.i(pspm%*%theta_EM[1:(p + sum(qm))])
		     phi_es <- l.phi.i(pspp%*%theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
		     z_es <- (y_work-mu_es)/sqrt(phi_es)
			 lres <- v_work*ifelse(event==1,(y_work-mu_es)^2 + mt_work - y_work^2,(y-mu_es)^2)/phi_es
			 Ltt[1:ncol(pspm),1:ncol(pspm)] <- -crossprod(pspm,pspm*matrix(l1.mu(mu_es)^2*(v_work/phi_es)*(1 + sqrt(phi_es)*z_es*l1.mu(mu_es)*l2.mu(mu_es)),n,ncol(pspm))) - penm
			 Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)] <- -crossprod(pspp,pspm*matrix(l1.mu(mu_es)*l1.phi(phi_es)*v_work*z_es/sqrt(phi_es),n,ncol(pspm)))
			 Ltt[1:ncol(pspm),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- t(Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)])
			 Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- -crossprod(pspp,pspp*matrix(l1.phi(phi_es)^2*(lres/2 + (1 + phi_es^2*l1.phi(phi_es)*l2.phi(phi_es))*(lres - 1)/2),n,ncol(pspp))) - penp
			 Ut[1:ncol(pspm)] <- crossprod(pspm,l1.mu(mu_es)*v_work*z_es/sqrt(phi_es)) - penm%*%theta_EM[1:(p+sum(qm))]
			 Ut[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- crossprod(pspp,l1.phi(phi_es)*(lres - 1)/2) - penp%*%theta_EM[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)]
			 theta_EM_new <- theta_EM + crossprod(solve(-Ltt),Ut)
		 	 theta_pd <- ifelse(abs(theta_EM) <= epsilon,1,theta_EM)
		     tol_EM <- max(abs(theta_EM_new-theta_EM)/abs(theta_pd))
		     cont_EM <- cont_EM + 1
		     if(cont_EM > maxiter) stop("no convergence was obtained!!",call.=FALSE)
		 }
		 theta_new <- theta_EM_new
		 theta_pd <- ifelse(abs(theta) <= epsilon,1,theta)
	     tol <- max(abs(theta_new-theta)/abs(theta_pd))
		 cont <- cont + 1
		 if(cont > maxiter) stop("no convergence was obtained!!",call.=FALSE)
    }
	theta_new
}
