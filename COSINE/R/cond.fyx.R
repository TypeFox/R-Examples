cond.fyx <-
function(data.y,data.x,type)
{
	cl <- as.integer(as.factor(type))
	if(nlevels(as.factor(cl))<3)
	{
		n1 <- length(cl[cl==1])
		n2 <- length(cl[cl==2])
		fx <- c(data.x[cl==1],data.x[cl==2])
		fy <- c(data.y[cl==1],data.y[cl==2])
		mux1 <- mean(fx[1:n1])
		mux2 <- mean(fx[(n1+1):(n1+n2)])
		muy1 <- mean(fy[1:n1])
		muy2 <- mean(fy[(n1+1):(n1+n2)])
		varx1 <- var(fx[1:n1])
		varx2 <- var(fx[(n1+1):(n1+n2)])
		vary1 <- var(fy[1:n1])
		vary2 <- var(fy[(n1+1):(n1+n2)])
		pho1 <- cor(fx[1:n1],fy[1:n1])
		pho2 <- cor(fx[(n1+1):(n1+n2)],fy[(n1+1):(n1+n2)])
		x <- ( (muy1-muy2)-(mux1-mux2)*pho2*sqrt(vary2)/sqrt(varx2) )^2
		x <- x + varx1*(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))^2
		x <- x*n1/(n1+n2)
		y <- ( (muy1-muy2)-(mux1-mux2)*pho1*sqrt(vary1)/sqrt(varx1) )^2
		y <- y + varx2*(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))^2
		y <- y*n2/(n1+n2)
		z <- (x+y)/((n1+n2)*(vary1*(1-pho1*pho1)/n2+vary2*(1-pho2*pho2)/n1))
		if(is.na(z))
			z <- 0
		return(z)
	}

	if(nlevels(as.factor(cl))>2)
	{
                n1 <- length(cl[cl==1])
                n2 <- length(cl[cl==2])
                n3 <- length(cl[cl==3])
                p1 <- n1/(n1+n2+n3)
                p2 <- n2/(n1+n2+n3)
                p3 <- n3/(n1+n2+n3)
		     fx <- c(data.x[cl==1],data.x[cl==2],data.x[cl==3])
		     fy <- c(data.y[cl==1],data.y[cl==2],data.y[cl==3])
		     mux1 <- mean(fx[1:n1])
		     mux2 <- mean(fx[(n1+1):(n1+n2)])
                     mux3 <- mean(fx[(n1+n2+1):(n1+n2+n3)])
		     muy1 <- mean(fy[1:n1])
		     muy2 <- mean(fy[(n1+1):(n1+n2)])
                     muy3 <- mean(fy[(n1+n2+1):(n1+n2+n3)])
		     varx1 <- var(fx[1:n1])
		     varx2 <- var(fx[(n1+1):(n1+n2)])
                     varx3 <- var(fx[(n1+n2+1):(n1+n2+n3)])
		     vary1 <- var(fy[1:n1])
		     vary2 <- var(fy[(n1+1):(n1+n2)])
                     vary3 <- var(fy[(n1+n2+1):(n1+n2+n3)])
		     pho1 <- cor(fx[1:n1],fy[1:n1])
		     pho2 <- cor(fx[(n1+1):(n1+n2)],fy[(n1+1):(n1+n2)])
                     pho3 <- cor(fx[(n1+n2+1):(n1+n2+n3)],fy[(n1+n2+1):(n1+n2+n3)])

                ###k=1,i=1,j=2   
                item.1 <- p1*p1*p2*(((muy1-muy2)-(mux1*pho1*sqrt(vary1)/sqrt(varx1)-mux2*pho2*sqrt(vary2)/sqrt(varx2))+(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))*mux1)^2+(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))^2*varx1)

                ###k=1,i=1,j=3
                item.2 <- p1*p3*p1*(((muy1-muy3)-(mux1*pho1*sqrt(vary1)/sqrt(varx1)-mux3*pho3*sqrt(vary3)/sqrt(varx3))+(pho1*sqrt(vary1)/sqrt(varx1)-pho3*sqrt(vary3)/sqrt(varx3))*mux1)^2+(pho1*sqrt(vary1)/sqrt(varx1)-pho3*sqrt(vary3)/sqrt(varx3))^2*varx1)

                ###k=1,i=2,j=3
                item.3 <- p2*p3*p1*(((muy2-muy3)-(mux2*pho2*sqrt(vary2)/sqrt(varx2)-mux3*pho3*sqrt(vary3)/sqrt(varx3))+(pho2*sqrt(vary2)/sqrt(varx2)-pho3*sqrt(vary3)/sqrt(varx3))*mux1)^2+(pho2*sqrt(vary2)/sqrt(varx2)-pho3*sqrt(vary3)/sqrt(varx3))^2*varx1)

                ###k=2,i=1,j=2
                item.4 <- p1*p2*p2*(((muy1-muy2)-(mux1*pho1*sqrt(vary1)/sqrt(varx1)-mux2*pho2*sqrt(vary2)/sqrt(varx2))+(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))*mux2)^2+(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))^2*varx2)

                ###k=2,i=1,j=3
                item.5 <- p1*p3*p2*(((muy1-muy3)-(mux1*pho1*sqrt(vary1)/sqrt(varx1)-mux3*pho3*sqrt(vary3)/sqrt(varx3))+(pho1*sqrt(vary1)/sqrt(varx1)-pho3*sqrt(vary3)/sqrt(varx3))*mux2)^2+(pho1*sqrt(vary1)/sqrt(varx1)-pho3*sqrt(vary3)/sqrt(varx3))^2*varx2) 

                ###k=2,i=2,j=3
                item.6 <- p2*p3*p2*(((muy2-muy3)-(mux2*pho2*sqrt(vary2)/sqrt(varx2)-mux3*pho3*sqrt(vary3)/sqrt(varx3))+(pho2*sqrt(vary2)/sqrt(varx2)-pho3*sqrt(vary3)/sqrt(varx3))*mux2)^2+(pho2*sqrt(vary2)/sqrt(varx2)-pho3*sqrt(vary3)/sqrt(varx3))^2*varx2)

                ###k=3,i=1,j=2
                item.7 <- p1*p2*p3*(((muy1-muy2)-(mux1*pho1*sqrt(vary1)/sqrt(varx1)-mux2*pho2*sqrt(vary2)/sqrt(varx2))+(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))*mux3)^2+(pho1*sqrt(vary1)/sqrt(varx1)-pho2*sqrt(vary2)/sqrt(varx2))^2*varx3)

                ###k=3,i=1,j=3
                item.8 <- p1*p3*p3*(((muy1-muy3)-(mux1*pho1*sqrt(vary1)/sqrt(varx1)-mux3*pho3*sqrt(vary3)/sqrt(varx3))+(pho1*sqrt(vary1)/sqrt(varx1)-pho3*sqrt(vary3)/sqrt(varx3))*mux3)^2+(pho1*sqrt(vary1)/sqrt(varx1)-pho3*sqrt(vary3)/sqrt(varx3))^2*varx3)

                ###k=3,i=2,j=3
                item.9 <- p2*p3*p3*(((muy2-muy3)-(mux2*pho2*sqrt(vary2)/sqrt(varx2)-mux3*pho3*sqrt(vary3)/sqrt(varx3))+(pho2*sqrt(vary2)/sqrt(varx2)-pho3*sqrt(vary3)/sqrt(varx3))*mux3)^2+(pho2*sqrt(vary2)/sqrt(varx2)-pho3*sqrt(vary3)/sqrt(varx3))^2*varx3)

                z <- (item.1+item.2+item.3+item.4+item.5+item.6+item.7+item.8+item.9)/(p1*vary1*(1-pho1^2)+p2*vary2*(1-pho2^2)+p3*vary3*(1-pho3^2))

		if(is.na(z))
			z <- 0
		return(z)   
		
	}	
}

