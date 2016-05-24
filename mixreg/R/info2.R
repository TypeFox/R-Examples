info2 <- function(x,y,theta,gma) {
#
# Note: The matrix x includes the column of 1's corresponding
# to the intercept term if an intercept is being fitted.
#

K    <- length(theta)
lK   <- theta[[K]]$lambda
gKK  <- aux1(gma,K,K)
pp1  <- ncol(x) + 1
pp2  <- pp1 + 1
nd   <- K*pp2-1
rslt <- matrix(NA,nd,nd)
for(j in 1:K) {
	tj <- theta[[j]]
	bj <- tj$beta
	sj <- tj$sigsq
	lj <- tj$lambda
	rj <- drop(y - x%*%bj)
	mj <- rj*x/sj
	vj <- 0.5*(rj**2/sj-1)/sj
	gjK  <- aux1(gma,j,K)
	for(k in j:K) {
		tk <- theta[[k]]
		bk <- tk$beta
		sk <- tk$sigsq
		lk <- tk$lambda
		rk <- drop(y - x%*%bk)
		mk <- rk*x/sk
		vk <- 0.5*(rk**2/sk-1)/sk
		gjk  <- aux1(gma,j,k)
		gkK  <- aux1(gma,k,K)
# beta-beta:
		t11 <- t(mj)%*%gjk%*%mk
# sigsq-beta:
		t21 <- drop(vj%*%gjk%*%mk)
# beta-sigsq:
		t12 <- drop(vk%*%t(gjk)%*%mj)
# sigsq-sigsq:
		t22 <- sum(gjk*vj%o%vk)
		if(j < K) {
# lambda-beta:
			t31 <- apply(gjk%*%mk,2,sum)/lj -
					apply(t(gkK)%*%mk,2,sum)/lK
# lambda-sigsq:
			t32 <- sum(gjk%*%vk)/lj - sum(vk%*%gkK)/lK
		}
		else t31 <- t32 <- NULL
		if(k < K) {
# beta-lambda:
			t13 <- apply(t(gjk)%*%mj,2,sum)/lk -
					apply(t(gjK)%*%mj,2,sum)/lK
# sigsq-lambda:
			t23 <- sum(t(gjk)%*%vj)/lk - sum(vj%*%gjK)/lK
		}
		else t13 <- t23 <- NULL
# lambda-lambda:
		t33 <- if(k < K && j < K)
			sum(gjk)/(lj*lk) - sum(gjK)/(lj*lK) -
				sum(gkK)/(lk*lK) + sum(gKK)/(lK*lK)
			else
				NULL

		tmp <- aux2(t11,t12,t13,t21,t22,t23,t31,t32,t33)
		indj <- if(j < K) (j-1)*pp2 + (1:pp2) else (j-1)*pp2 + (1:pp1)
		indk <- if(k < K) (k-1)*pp2 + (1:pp2) else (k-1)*pp2 + (1:pp1)
		rslt[indj,indk] <- tmp
		if(j!=k) rslt[indk,indj] <- t(tmp)
	}
}
rslt
}
