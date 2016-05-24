
	###########################
	#	third expansion	  #
	###########################

	#bi:r.size*division	#aMat:k.size*r.size*division
	#Sigma:k.size*k.size

	#env:d.size,r.size,k.size,division,delta,state,pars,aMat,mu,Sigma,
	#    invSigma,invY,block,my.range,Diff

	#first:k*k*k*t, second:k*t

	I_123 <- function(b1,b2,b3,env){

		r.size <- env$r.size
		k.size <- env$k.size
		delta <- env$delta
		aMat <- env$aMat
		invSigma <- env$invSigma
		block <- env$block
		my.range <- env$my.range
		Diff <- env$Diff

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n5 <- length(b3)
		n6 <- sum(b3 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) * (n5 - n6 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,k.size,block))
			second <- matrix(0,k.size,block)

			return(list(first=first,second=second))
		}

		b1a <- matrix(0,k.size,block)
		b2a <- matrix(0,k.size,block)
		b3a <- matrix(0,k.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
			b2 <- matrix(b2,1,block)
			b3 <- matrix(b3,1,block)
		}

		for(t in 1:block){
			b1a[,t] <- matrix(b1[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
			b2a[,t] <- matrix(b2[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
			b3a[,t] <- matrix(b3[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
		}

		first <- array(0,dim=c(k.size,k.size,k.size,block))

		tmp1 <- matrix(0,k.size,block)
		tmp2 <- array(0,dim=c(k.size,k.size,block))
		tmp3 <- array(0,dim=c(k.size,k.size,k.size,block))

		for(t in 2:block){
			tmp1[,t] <- b3a %*% Diff[t,] * delta
		}

		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(t in 2:block){
			tmp2[k1,k2,t] <- (b2a[k1,] * tmp1[k2,]) %*%
					     Diff[t,] * delta
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			for(t in 2:block){
		        tmp3[k1,k2,k3,t] <- (b1a[k1,] * tmp2[k2,k3,]) %*% 
						    Diff[t,] * delta
			}
		    }
		  }
		}


		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(t in 1:block){
			first[k1,k2,,t] <- matrix(tmp3[k1,k2,,t],1,k.size) %*%
						 invSigma
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(k3 in 1:k.size){
		    for(t in 1:block){
			first[k1,,k3,t] <- matrix(first[k1,,k3,t],1,k.size) %*%
						 invSigma
		    }
		  }
		}

		for(k2 in 1:k.size){
		  for(k3 in 1:k.size){
		    for(t in 1:block){
			first[,k2,k3,t] <- matrix(first[,k2,k3,t],1,k.size) %*%
						 invSigma
		    }
		  }
		}


		second1 <- matrix(0,k.size,block)

		for(t in 1:block){
		  for(k1 in 1:k.size){
		    for(k2 in 1:k.size){
			second1[,t] <- second1[,t] + tmp3[k1,k2,,t] * invSigma[k1,k2]
		    }
		  }
		  second1[,t] <- matrix(second1[,t],1,k.size) %*% invSigma
		}

		second2 <- matrix(0,k.size,block)

		for(t in 1:block){
		  for(k1 in 1:k.size){
		    for(k3 in 1:k.size){
			second2[,t] <- second2[,t] + tmp3[k1,,k3,t] * invSigma[k1,k3]
		    }
		  }
		  second2[,t] <- matrix(second2[,t],1,k.size) %*% invSigma
		}

		second3 <- matrix(0,k.size,block)

		for(t in 1:block){
		  for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			second3[,t] <- second3[,t] + tmp3[,k2,k3,t] * invSigma[k2,k3]
		    }
		  }
		  second3[,t] <- matrix(second3[,t],1,k.size) %*% invSigma
		}

		second <- - second1 - second2 - second3
		return(list(first=first,second=second))
	}


	I_123_x <- function(x,get_I_123,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(get_I_123$first[,,,block],dim=c(k.size,k.size,k.size))

		for(k1 in 1:k.size){
		for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			result1 <- result1 + first[k1,k2,k3] *
				     x[k1] * x[k2] * x[k3]
		}
		}
		}

		second <- get_I_123$second[,block]

		result2 <- matrix(x,1,k.size) %*%
			     matrix(second,k.size,1)

		result <- result1 + result2

		return(result)
	}


	#first:k*k*k*k*t, second:k*k*t, third:t

	I_1234 <- function(b1,b2,b3,b4,env){

		r.size <- env$r.size
		k.size <- env$k.size
		delta <- env$delta
		aMat <- env$aMat
		invSigma <- env$invSigma
		block <- env$block
		my.range <- env$my.range
		Diff <- env$Diff

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n5 <- length(b3)
		n6 <- sum(b3 == 0)

		n7 <- length(b4)
		n8 <- sum(b4 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) *
		     (n5 - n6 != 0) * (n7 - n8 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
			second <- array(0,dim=c(k.size,k.size,block))
			third <- double(block)

			return(list(first=first,second=second,third=third))
		}

		b1a <- matrix(0,k.size,block)
		b2a <- matrix(0,k.size,block)
		b3a <- matrix(0,k.size,block)
		b4a <- matrix(0,k.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
			b2 <- matrix(b2,1,block)
			b3 <- matrix(b3,1,block)
			b4 <- matrix(b4,1,block)
		}

		for(t in 1:block){
			b1a[,t] <- matrix(b1[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
			b2a[,t] <- matrix(b2[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
			b3a[,t] <- matrix(b3[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
			b4a[,t] <- matrix(b4[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
		}

		first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))

		tmp1 <- matrix(0,k.size,block)
		tmp2 <- array(0,dim=c(k.size,k.size,block))
		tmp3 <- array(0,dim=c(k.size,k.size,k.size,block))
		tmp4 <- array(0,dim=c(k.size,k.size,k.size,k.size,block))

		for(t in 2:block){
			tmp1[,t] <- b4a %*% Diff[t,] * delta
		}

		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(t in 2:block){
			tmp2[k1,k2,t] <- (b3a[k1,] * tmp1[k2,]) %*%
					     Diff[t,] * delta
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			for(t in 2:block){
			  tmp3[k1,k2,k3,t] <- (b2a[k1,] * tmp2[k2,k3,]) %*%
						    Diff[t,]  * delta
			}
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			for(k4 in 1:k.size){
			  for(t in 2:block){
			    tmp4[k1,k2,k3,k4,t] <- (b1a[k1,] * tmp3[k2,k3,k4,]) %*%
							   Diff[t,]* delta
			  }
			}
		    }
		  }
		}


		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			for(t in 1:block){
			  first[k1,k2,k3,,t] <- matrix(tmp4[k1,k2,k3,,t],1,k.size) %*%
							invSigma
			}
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(k4 in 1:k.size){
			for(t in 1:block){
			  first[k1,k2,,k4,t] <- matrix(first[k1,k2,,k4,t],1,k.size) %*%
							invSigma
			}
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(k3 in 1:k.size){
		    for(k4 in 1:k.size){
			for(t in 1:block){
			  first[k1,,k3,k4,t] <- matrix(first[k1,,k3,k4,t],1,k.size) %*%
							invSigma
			}
		    }
		  }
		}

		for(k2 in 1:k.size){
		  for(k3 in 1:k.size){
		    for(k4 in 1:k.size){
			for(t in 1:block){
			  first[,k2,k3,k4,t] <- matrix(first[,k2,k3,k4,t],1,k.size) %*%
							invSigma
			}
		    }
		  }
		}


		second1 <- array(0,dim=c(k.size,k.size,block))

		tmp5 <- array(0,dim=c(k.size,k.size,block))

		for(t in 1:block){
		  for(k1 in 1:k.size){
		    for(k2 in 1:k.size){
			tmp5[,,t] <- tmp5[,,t] + tmp4[k1,k2,,,t] *
					    invSigma[k1,k2]
		    }
		  }
		}

		for(k3 in 1:k.size){
		  for(t in 1:block){
		    second1[k3,,t] <- matrix(tmp5[k3,,t],1,k.size) %*%
					    invSigma
		  }
		}

		for(k4 in 1:k.size){
		  for(t in 1:block){
		    second1[,k4,t] <- matrix(second1[,k4,t],1,k.size) %*%
					    invSigma
		  }
		}


		second2 <- array(0,dim=c(k.size,k.size,block))

		tmp6 <- array(0,dim=c(k.size,k.size,block))

		for(t in 1:block){
		  for(k1 in 1:k.size){
		    for(k3 in 1:k.size){
			tmp6[,,t] <- tmp6[,,t] + tmp4[k1,,k3,,t] *
					    invSigma[k1,k3]
		    }
		  }
		}

		for(k2 in 1:k.size){
		  for(t in 1:block){
		    second2[k2,,t] <- matrix(tmp6[k2,,t],1,k.size) %*%
					    invSigma
		  }
		}

		for(k4 in 1:k.size){
		  for(t in 1:block){
		    second2[,k4,t] <- matrix(second2[,k4,t],1,k.size) %*%
					    invSigma
		  }
		}


		second3 <- array(0,dim=c(k.size,k.size,block))

		tmp7 <- array(0,dim=c(k.size,k.size,block))

		for(t in 1:block){
		  for(k1 in 1:k.size){
		    for(k4 in 1:k.size){
			tmp7[,,t] <- tmp7[,,t] + tmp4[k1,,,k4,t] *
					    invSigma[k1,k4]
		    }
		  }
		}

		for(k2 in 1:k.size){
		  for(t in 1:block){
		    second3[k2,,t] <- matrix(tmp7[k2,,t],1,k.size) %*%
					    invSigma
		  }
		}

		for(k3 in 1:k.size){
		  for(t in 1:block){
		    second3[,k3,t] <- matrix(second3[,k3,t],1,k.size) %*%
					    invSigma
		  }
		}


		second4 <- array(0,dim=c(k.size,k.size,block))

		for(t in 1:block){
		  for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			second4[,,t] <- second4[,,t] + tmp4[,k2,k3,,t] *
					    invSigma[k2,k3]
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(t in 1:block){
		    second4[k1,,t] <- matrix(second4[k1,,t],1,k.size) %*%
					    invSigma
		  }
		}

		for(k4 in 1:k.size){
		  for(t in 1:block){
		    second4[,k4,t] <- matrix(second4[,k4,t],1,k.size) %*%
					    invSigma
		  }
		}


		second5 <- array(0,dim=c(k.size,k.size,block))

		for(t in 1:block){
		  for(k2 in 1:k.size){
		    for(k4 in 1:k.size){
			second5[,,t] <- second5[,,t] + tmp4[,k2,,k4,t] *
					    invSigma[k2,k4]
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(t in 1:block){
		    second5[k1,,t] <- matrix(second5[k1,,t],1,k.size) %*%
					    invSigma
		  }
		}

		for(k3 in 1:k.size){
		  for(t in 1:block){
		    second5[,k3,t] <- matrix(second5[,k3,t],1,k.size) %*%
					    invSigma
		  }
		}


		second6 <- array(0,dim=c(k.size,k.size,block))

		for(t in 1:block){
		  for(k3 in 1:k.size){
		    for(k4 in 1:k.size){
			second6[,,t] <- second6[,,t] + tmp4[,,k3,k4,t] *
					    invSigma[k3,k4]
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(t in 1:block){
		    second6[k1,,t] <- matrix(second6[k1,,t],1,k.size) %*%
					    invSigma
		  }
		}

		for(k2 in 1:k.size){
		  for(t in 1:block){
		    second6[,k2,t] <- matrix(second6[,k2,t],1,k.size) %*%
					    invSigma
		  }
		}

		second <- - second1 - second2 - second3 - second4 - second5 - second6


		third1 <- double(block)

		for(t in 1:block){
		  for(k3 in 1:k.size){
		    for(k4 in 1:k.size){
			third1[t] <- third1[t] + tmp5[k3,k4,t] *
					 invSigma[k3,k4]
		    }
		  }
		}


		third2 <- double(block)

		for(t in 1:block){
		  for(k2 in 1:k.size){
		    for(k4 in 1:k.size){
			third2[t] <- third2[t] + tmp6[k2,k4,t] *
					 invSigma[k2,k4]
		    }
		  }
		}

		third3 <- double(block)


		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		for(t in 1:block){
		  for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			third3[t] <- third3[t] + tmp6[k2,k3,t] *
					 invSigma[k2,k3]
		    }
		  }
		}

		third <- third1 + third2 + third3

		return(list(first=first,second=second,third=third))
	}


	I_1234_x <- function(x,get_I_1234,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(get_I_1234$first[,,,,block],
				   dim=c(k.size,k.size,k.size,k.size))

		for(k1 in 1:k.size){
		for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			for(k4 in 1:k.size){
			  result1 <- result1 + first[k1,k2,k3,k4] *
					 x[k1] * x[k2] * x[k3] * x[k4]
			}
		}
		}
		}

		second <- matrix(get_I_1234$second[,,block],k.size,k.size)

		result2 <- t(matrix(x,k.size,1)) %*% second %*% 
			     matrix(x,k.size,1)

		third <- get_I_1234$third[block]

		result <- result1 + result2 + third

		return(result)
	}


	b1b2 <- function(b1,b2,env){

		r.size <- env$r.size
		block <- env$block

		result <- double(block)

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		if(n == 0){
			return(result)
		}

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
			b2 <- matrix(b2,1,block)
		}

		for(t in 1:block){
			result[t] <- matrix(b1[,t],1,r.size) %*%
					 matrix(b2[,t],r.size,1)
		}

		return(result)
	}


	b1_b2 <- function(b1,b2,env){

		delta <- env$delta
		block <- env$block
		Diff <- env$Diff

		tmp <- b1b2(b1,b2,env)

		result <- double(block)

		for(t in 1:block){
			result[t] <- tmp %*% Diff[t,] * delta
		}

		return(result)
	}
	

	I0 <- function(f,env){

		delta <- env$delta
		block <- env$block
		Diff <- env$Diff

		result <- double(block)

		n1 <- length(f)
		n2 <- sum(f == 0)

		n <- n1 - n2

		if(n == 0){
			return(result)
		}

		for(t in 1:block){
			result[t] <- f %*% Diff[t,] * delta
		}

		return(result)
	}


	#k*t

	I_1 <- function(b1,env){

		r.size <- env$r.size
		k.size <- env$k.size
		delta <- env$delta
		aMat <- env$aMat
		invSigma <- env$invSigma
		block <- env$block
		my.range <- env$my.range
		Diff <- env$Diff

		result <- matrix(0,k.size,block)

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		if(n1 == n2){
			return(result)
		}

		b1a <- matrix(0,k.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
		}

		for(t in 1:block){
			b1a[,t] <- matrix(b1[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
		}

		for(t in 2:block){
			result[,t] <- t(invSigma) %*% b1a %*% Diff[t,] * delta
		}

		return(result)
	}


	I_1_x <- function(x,get_I_1,env){

		k.size <- env$k.size
		block <- env$block

            tmp <- get_I_1[,block]

		result <- matrix(tmp,1,k.size) %*% matrix(x,k.size,1)

		return(result)
	}


	#first:k*k*t, second:t

	I_12 <- function(b1,b2,env){

		r.size <- env$r.size
		k.size <- env$k.size
		delta <- env$delta
		aMat <- env$aMat
		Sigma <- env$Sigma
		invSigma <- env$invSigma
		block <- env$block
		my.range <- env$my.range
		Diff <- env$Diff

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,block))
			second <- double(block)

			return(list(first=first,second=second))
		}

		b1a <- matrix(0,k.size,block)
		b2a <- matrix(0,k.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
			b2 <- matrix(b2,1,block)
		}

		for(t in 1:block){
			b1a[,t] <- matrix(b1[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
			b2a[,t] <- matrix(b2[,t],1,r.size) %*%
				     t(matrix(aMat[,,my.range[t]],k.size,r.size))
		}

		first <- array(0,dim=c(k.size,k.size,block))

		tmp1 <- matrix(0,k.size,block)
		tmp2 <- array(0,dim=c(k.size,k.size,block))

		for(t in 2:block){
			tmp1[,t] <- b2a %*% Diff[t,] * delta
		}

		for(k1 in 1:k.size){
		  for(k2 in 1:k.size){
		    for(t in 2:block){
			tmp2[k1,k2,t] <- (b1a[k1,] * tmp1[k2,]) %*%
					     Diff[t,] * delta
		    }
		  }
		}

		for(k1 in 1:k.size){
		  for(t in 1:block){
		    first[k1,,t] <- matrix(tmp2[k1,,t],1,k.size) %*%
					  invSigma
		  }
		}

		for(k2 in 1:k.size){
		  for(t in 1:block){
		    first[,k2,t] <- matrix(first[,k2,t],1,k.size) %*%
					  invSigma
		  }
		}


		second <- double(block)

		for(t in 1:block){
		  for(k1 in 1:k.size){
		    for(k2 in 1:k.size){
			second[t] <- second[t] + tmp2[k1,k2,t] * invSigma[k1,k2]
		    }
		  }
		}

		second <- - second

		return(list(first=first,second=second))
	}
		

	I_12_x <- function(x,get_I_12,env){

		k.size <- env$k.size
		block <- env$block

            first <- matrix(get_I_12$first[,,block],k.size,k.size)
		second <- get_I_12$second[block]

		result <- matrix(x,1,k.size) %*% first %*%
			    matrix(x,k.size,1) + second

		return(result)
	}


	#first:k*k*t, second:t

	I_1_2 <- function(b1,b2,env){

		k.size <- env$k.size
		block <- env$block

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,block))
			second <- double(block)

			return(list(first=first,second=second))
		}

		first.tmp <- I_12(b1,b2,env)
		second.tmp <- I_12(b2,b1,env)
		third.tmp <- b1_b2(b1,b2,env)

		first <- first.tmp$first + second.tmp$first
		second <- first.tmp$second + second.tmp$second + third.tmp

		return(list(first=first,second=second))
	}


	I_1_2_x <- function(x,get_I_1_2,env){

		result <- I_12_x(x,get_I_1_2,env)

		return(result)
	}


	#first:k*k*t, second:t

	I_1_23 <- function(b1,b2,b3,env){

		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block
		my.range <- env$my.range
		Diff <- env$Diff

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n5 <- length(b3)
		n6 <- sum(b3 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) * (n5 - n6 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,k.size,block))
			second <- matrix(0,k.size,block)

			return(list(first=first,second=second))
		}

		first.tmp <- I_123(b1,b2,b3,env)
		second.tmp <- I_123(b2,b1,b3,env)
		third.tmp <- I_123(b2,b3,b1,env)

		tmp1 <- b1_b2(b1,b3,env)
		tmp2 <- matrix(0,r.size,block)

		if(r.size == 1){
			b2 <- matrix(b2,1,block)
		}

		for(t in 1:block){
			tmp2[,t] <- tmp1[t] * b2[,t]
		}

		fourth.tmp <- I_1(tmp2,env)

		tmp3 <- b1b2(b1,b2,env)
		tmp4 <- I_1(b3,env)
		tmp5 <- matrix(0,k.size,block)

		for(t in 1:block){
			tmp5[,t] <- tmp3[t] * tmp4[,t]
		}

		fifth.tmp <- matrix(0,k.size,block)

		for(k in 1:k.size){
			fifth.tmp[k,] <- I0(tmp5[k,],env)
		}

		first <- first.tmp$first + second.tmp$first + third.tmp$first
		second <- first.tmp$second + second.tmp$second + third.tmp$second +
			    fourth.tmp + fifth.tmp

		return(list(first=first,second=second))
	}


	I_1_23_x <- function(x,get_I_1_23,env){

		result <- I_123_x(x,get_I_1_23,env)

		return(result)
	}


	#first:k*k*t, second:t

	I_12_3 <- function(b1,b2,b3,env){

		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n5 <- length(b3)
		n6 <- sum(b3 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) * (n5 - n6 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,k.size,block))
			second <- matrix(0,k.size,block)

			return(list(first=first,second=second))
		}

		first.tmp <- I_123(b1,b2,b3,env)
		second.tmp <- I_123(b1,b3,b2,env)

		tmp1 <- b1_b2(b2,b3,env)
		tmp2 <- matrix(0,r.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
		}

		for(t in 1:block){
			tmp2[,t] <- tmp1[t] * b1[,t]
		}

		third.tmp <- I_1(tmp2,env)

		first <- first.tmp$first + second.tmp$first
		second <- first.tmp$second + second.tmp$second + third.tmp

		return(list(first=first,second=second))
	}


	I_12_3_x <- function(x,get_I_12_3,env){

		result <- I_123_x(x,get_I_12_3,env)

		return(result)
	}


	#first:k*k*t, second:t
	
	I_1_2_3 <- function(b1,b2,b3,env){

		k.size <- env$k.size
		block <- env$block

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n5 <- length(b3)
		n6 <- sum(b3 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) * (n5 - n6 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,k.size,block))
			second <- matrix(0,k.size,block)

			return(list(first=first,second=second))
		}

		first.tmp <- I_123(b1,b2,b3,env)
		second.tmp <- I_123(b2,b1,b3,env)
		third.tmp <- I_123(b2,b3,b1,env)
		fourth.tmp <- I_123(b1,b3,b2,env)
		fifth.tmp <- I_123(b3,b1,b2,env)
		sixth.tmp <- I_123(b3,b2,b1,env)

		tmp1 <- b1_b2(b1,b3,env)
		tmp2 <- I_1(b2,env)
		seventh.tmp <- matrix(0,k.size,block)

		for(t in 1:block){
			seventh.tmp[,t] <- tmp1[t] * tmp2[,t]
		}

		tmp3 <- b1_b2(b1,b2,env)
		tmp4 <- I_1(b3,env)
		eighth.tmp <- matrix(0,k.size,block)

		for(t in 1:block){
			eighth.tmp[,t] <- tmp3[t] * tmp4[,t]
		}

		tmp5 <- b1_b2(b2,b3,env)
		tmp6 <- I_1(b1,env)
		ninth.tmp <- matrix(0,k.size,block)

		for(t in 1:block){
			ninth.tmp[,t] <- tmp5[t] * tmp6[,t]
		}

		first <- first.tmp$first + second.tmp$first + third.tmp$first +
			   fourth.tmp$first + fifth.tmp$first + sixth.tmp$first

		second <- first.tmp$second + second.tmp$second + third.tmp$second +
			    fourth.tmp$second + fifth.tmp$second + sixth.tmp$second +
			    seventh.tmp + eighth.tmp + ninth.tmp

		return(list(first=first,second=second))
	}


	I_1_2_3_x <- function(x,get_I_1_2_3,env){

		result <- I_123_x(x,get_I_1_2_3,env)

		return(result)
	}


	#first:k*k*k*k*t, second:k*k*t, third:t

	I_12_34 <- function(b1,b2,b3,b4,env){

		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		n1 <- length(b1)
		n2 <- sum(b1 == 0)

		n3 <- length(b2)
		n4 <- sum(b2 == 0)

		n5 <- length(b3)
		n6 <- sum(b3 == 0)

		n7 <- length(b4)
		n8 <- sum(b4 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) *
		     (n5 - n6 != 0) * (n7 - n8 != 0)

		if(n == 0){
			first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
			second <- array(0,dim=c(k.size,k.size,block))
			third <- double(block)

			return(list(first=first,second=second,third=third))
		}

		first.tmp <- I_1234(b1,b2,b3,b4,env)
		second.tmp <- I_1234(b1,b3,b4,b2,env)
		third.tmp <- I_1234(b1,b3,b2,b4,env)

		tmp1 <- b1_b2(b2,b4,env)
		tmp2 <- matrix(0,r.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
			b2 <- matrix(b2,1,block)
			b3 <- matrix(b3,1,block)
			b4 <- matrix(b4,1,block)
		}

		for(t in 1:block){
			tmp2[,t] <- tmp1[t] * b3[,t]
		}

		fourth.tmp <- I_12(b1,tmp2,env)

		tmp3 <- b1_b2(b2,b3,env)
		tmp4 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp4[,t] <- tmp3[t] * b1[,t]
		}

		fifth.tmp <- I_12(tmp4,b4,env)

		tmp5 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp5[,t] <- tmp3[t] * b4[,t]
		}

		sixth.tmp <- I_12(b1,tmp5,env)

		seventh.tmp <- I_1234(b3,b4,b1,b2,env)
		eighth.tmp <- I_1234(b3,b1,b2,b4,env)
		ninth.tmp <- I_1234(b3,b1,b4,b2,env)

		tmp6 <- b1_b2(b4,b2,env)
		tmp7 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp7[,t] <- tmp6[t] * b1[,t]
		}

		tenth.tmp <- I_12(b3,tmp7,env)

		tmp8 <- b1_b2(b4,b1,env)
		tmp9 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp9[,t] <- tmp8[t] * b3[,t]
		}

		eleventh.tmp <- I_12(tmp9,b2,env)

		tmp10 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp10[,t] <- tmp8[t] * b2[,t]
		}

		twelfth.tmp <- I_12(b3,tmp10,env)

		tmp11 <- b1_b2(b1,b3,env)
		tmp12 <- I_12(b2,b4,env)
		thirteenth.tmp <- tmp12

		for(t in 1:block){
			thirteenth.tmp$first[,,t] <- tmp11[t] * tmp12$first[,,t]
			thirteenth.tmp$second[t] <- tmp11[t] * tmp12$second[t]
		}

		tmp13 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp13[,t] <- tmp11[t] * b2[,t]
		}

		fourteenth.tmp <- I_12(tmp13,b4,env)

		tmp14 <- I_12(b4,b2,env)
		fifteenth.tmp <- tmp14

		for(t in 1:block){
			fifteenth.tmp$first[,,t] <- tmp11[t] * tmp14$first[,,t]
			fifteenth.tmp$second[t] <- tmp11[t] * tmp14$second[t]
		}

		tmp15 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp15[,t] <- tmp11[t] * b4[,t]
		}

		sixteenth.tmp <- I_12(tmp15,b2,env)

		tmp16 <- b1_b2(b4,b2,env)
		tmp17 <- matrix(0,r.size,block)

		for(t in 1:block){
			tmp17[,t] <- tmp16[t] * b3[,t]
		}

		tmp18 <- b1b2(b1,tmp17,env)

		seventeenth.tmp <- I0(tmp18,env)

		first <- first.tmp$first + second.tmp$first + third.tmp$first +
			   seventh.tmp$first + eighth.tmp$first + ninth.tmp$first

		second <- first.tmp$second + second.tmp$second + third.tmp$second +
			    fourth.tmp$first + fifth.tmp$first - sixth.tmp$first +
			    seventh.tmp$second + eighth.tmp$second + ninth.tmp$second +
			    tenth.tmp$first + eleventh.tmp$first - twelfth.tmp$first +
			    thirteenth.tmp$first - fourteenth.tmp$first + 
			    fifteenth.tmp$first - sixteenth.tmp$first

		third <- first.tmp$third + second.tmp$third + third.tmp$third +
			   fourth.tmp$second + fifth.tmp$second - sixth.tmp$second +
			   seventh.tmp$third + eighth.tmp$third + ninth.tmp$third +
			   tenth.tmp$second + eleventh.tmp$second - twelfth.tmp$second +
			   thirteenth.tmp$second - fourteenth.tmp$second + 
			   fifteenth.tmp$second - sixteenth.tmp$second +
			   seventeenth.tmp

		return(list(first=first,second=second,third=third))
	}


	I_12_34_x <- function(x,get_I_12_34,env){

		result <- I_1234_x(x,get_I_12_34,env)

		return(result)
	}


	#first:k*k*k*k*t, second:k*k*t, third:t

	I_1_2_34 <- function(b1,b2,b3,b4,env){

		division <- env$division

		first.tmp <- I_12_34(b1,b2,b3,b4,env)
		second.tmp <- I_12_34(b2,b1,b3,b4,env)

		tmp1 <- b1_b2(b1,b2,env)
		tmp2 <- I_12(b3,b4,env)
		third.tmp <- tmp2

		for(t in 1:division){
			third.tmp$first[,,t] <- tmp1[t] * tmp2$first[,,t]
			third.tmp$second[t] <- tmp1[t] * tmp2$second[t]
		}

		first <- first.tmp$first + second.tmp$first
		second <- first.tmp$second + second.tmp$second + third.tmp$first
		third <- first.tmp$third + second.tmp$third + third.tmp$second

		return(list(first=first,second=second,third=third))
	}


	I_1_2_34_x <- function(x,get_I_1_2_34,env){

		result <- I_1234_x(x,get_I_1_2_34,env)

		return(result)
	}


	I_1_2_3_4 <- function(b1,b2,b3,b4){


	}


##p16 Lemma3

	I0_a <- function(b1,c1,env){

		r.size <- env$r.size
		block <- env$block

		first.coef <- I0(c1,env)
		first <- b1

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
		}

		second <- matrix(0,r.size,block) 

		for(r in 1:r.size){
			second[r,] <- first.coef * b1[r,]
		}

		return(list(first.coef=first.coef,first=first,second=second))
	}


	I0_b <- function(b1,b2,c1,env){

		r.size <- env$r.size
		block <- env$block

		first.coef <- I0(c1,env)
		first <- list()

		first[[1]] <- b1
		first[[2]] <- b2

		second <- list() 

		second[[1]] <- matrix(0,r.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
		}

		for(r in 1:r.size){
			second[[1]][r,] <- first.coef * b1[r,]
		}

		second[[2]] <- b2

		return(list(first.coef=first.coef,first=first,second=second))
	}


	I0_c <- function(b1,b2,c1,c2,env){

		r.size <- env$r.size
		block <- env$block

		tmp1 <- I0(c2,env)

		tmp2 <- c1 * tmp1

		first.coef <- I0(tmp2,env)
		first <- list()

		first[[1]] <- b1
		first[[2]] <- b2

		second <- list() 

		second[[1]] <- matrix(0,r.size,block)

		if(r.size == 1){
			b1 <- matrix(b1,1,block)
		}

		for(r in 1:r.size){
			second[[1]][r,] <- first.coef * b1[r,]
		}

		second[[2]] <- b2

		third.coef <- I0(c1,env)
		third <- list()

		third[[1]] <- matrix(0,r.size,block)

		for(r in 1:r.size){
			third[[1]][r,] <- tmp1 * b1[r,]
		}

		third[[2]] <- b2

		fourth <- list()

		fourth[[1]] <- matrix(0,r.size,block)

		for(r in 1:r.size){
			fourth[[1]][r,] <- third.coef * tmp1 * b1[r,]
		}

		fourth[[2]] <- b2

		return(list(first.coef=first.coef,first=first,second=second,
				third.coef=third.coef,third=third,fourth=fourth))
	}


	#d*t

	Y_e_V0 <- function(X.t0,de.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,d.size,block)

		tmp <- c()

		assign(pars[1],0)

		for(t in 1:block){
			for(d in 1:d.size){
				assign(state[d],X.t0[my.range[t],d])
			}
			for(d in 1:d.size){
				tmp[d] <- eval(de.drift[[d]])
			}
			result[,t] <- invY[,,t] %*% tmp
		}

		return(result)
	}


	#d*r*t

	Y_e_V <- function(X.t0,de.diffusion,env){

		d.size <- env$d.size
		r.size <- env$r.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,r.size,block))

		tmp <- matrix(0,d.size,r.size)

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i in 1:d.size){
		    for(r in 1:r.size){
			tmp[i,r] <- eval(de.diffusion[[i]][[r]])
		    }
		  }
		  result[,,t] <- invY[,,t] %*% tmp
		}

		return(result)
	}


	#d*t

	Y_D <- function(X.t0,tmpY,get_Y_e_V0,env){

		d.size <- env$d.size
		delta <- env$delta
		block <- env$block
		Diff <- env$Diff

		result <- matrix(0,d.size,block)

		for(i in 1:d.size){
		  for(t in 2:block){
		    result[i,] <- get_Y_e_V0[i,] %*% Diff[t,] * delta
		  }
		}

		for(t in 2:block){
			result[,t] <- tmpY[,,t] %*% result[,t]
		}

		return(result)
	}


	#d*d*d*t

	Y_x1_x2_V0 <- function(X.t0,dxdx.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,d.size,block))

		tmp <- c()

		assign(pars[1],0)


		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			tmp.drift <- dxdx.drift[[i1]][[i2]] 

			for(d in 1:d.size){
			  tmp[d] <- eval(tmp.drift[[d]])
			}

			result[i1,i2,,t] <- invY[,,t] %*% tmp
		    }
		  }
		}

		return(result)
	}


	#d*d*d*r*t

	Y_x1_x2_V <- function(X.t0,dxdx.diffusion,env){

		d.size <- env$d.size
		r.size <- env$r.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,d.size,r.size,block))

		tmp <- matrix(0,d.size,r.size)

		assign(pars[1],0)

		for(t in 1:block){
		for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		  for(i2 in 1:d.size){
			tmp.diffusion <- dxdx.diffusion[[i1]][[i2]] 

			for(i in 1:d.size){
			for(r in 1:r.size){
			tmp[i,r] <- eval(tmp.diffusion[[i]][[r]])
			}
			}
			result[i1,i2,,,t] <- invY[,,t] %*% tmp
		}
		}
		}

		return(result)
	}


	#d*d*t

	Y_x_e_V0 <- function(X.t0,dxde.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,block))

		tmp <- c()

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    tmp.drift <- dxde.drift[[i1]]

		    for(d in 1:d.size){
			tmp[d] <- eval(tmp.drift[[d]])
		    }

		    result[i1,,t] <- invY[,,t] %*% tmp
		  }
		}

		return(result)
	}


	#d*d*r*t

	Y_x_e_V <- function(X.t0,dxde.diffusion,env){

		d.size <- env$d.size
		r.size <- env$r.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,r.size,block))

		tmp <- matrix(0,d.size,r.size)

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    tmp.diffusion <- dxde.diffusion[[i1]] 

		    for(i in 1:d.size){
			for(r in 1:r.size){
			  tmp[i,r] <- eval(tmp.diffusion[[i]][[r]])
			}
		    }
		    result[i1,,,t] <- invY[,,t] %*% tmp
		  }
		}

		return(result)
	}


	#d*t

	Y_e_e_V0 <- function(X.t0,dede.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,d.size,block)

		tmp <- c()

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(d in 1:d.size){
		    tmp[d] <- eval(dede.drift[[d]])
		  }

		  result[,t] <- invY[,,t] %*% tmp
		}

		return(result)
	}


	#d*r*t

	Y_e_e_V <- function(X.t0,dede.diffusion,env){

		d.size <- env$d.size
		r.size <- env$r.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,r.size,block))

		tmp <- matrix(0,d.size,r.size)

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i in 1:d.size){
		    for(r in 1:r.size){
			tmp[i,r] <- eval(dede.diffusion[[i]][[r]])
		    }
		  }
		  result[,,t] <- invY[,,t] %*% tmp
		}

		return(result)
	}


	#d*d*d*d*t

	Y_x1_x2_x3_V0 <- function(X.t0,dxdxdx.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,d.size,d.size,block))

		tmp <- c()

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			for(i3 in 1:d.size){
			  tmp.drift <- dxdxdx.drift[[i1]][[i2]][[i3]]

			  for(d in 1:d.size){
			    tmp[d] <- eval(tmp.drift[[d]])
			  }

			  result[i1,i2,i3,,t] <- invY[,,t] %*% tmp
			}
		    }
		  }
		}

		return(result)
	}


	#d*d*d*t

	Y_x1_x2_e_V0 <- function(X.t0,dxdxde.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,d.size,block))

		tmp <- c()

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			tmp.drift <- dxdxde.drift[[i1]][[i2]] 

			for(d in 1:d.size){
			  tmp[d] <- eval(tmp.drift[[d]])
			}

			result[i1,i2,,t] <- invY[,,t] %*% tmp
		    }
		  }
		}

		return(result)
	}


	#d*d*d*r*t

	Y_x1_x2_e_V <- function(X.t0,dxdxde.diffusion,env){

		d.size <- env$d.size
		r.size <- env$r.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,d.size,r.size,block))

		tmp <- matrix(0,d.size,r.size)

		assign(pars[1],0)

		for(t in 1:block){
		for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		 }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			tmp.diffusion <- dxdxde.diffusion[[i1]][[i2]]

			for(i in 1:d.size){
			for(r in 1:r.size){
			tmp[i,r] <- eval(tmp.diffusion[[i]][[r]])
			}
			}
			result[i1,i2,,,t] <- invY[,,t] %*% tmp
		}
		}
		}

		return(result)
	}


	#d*d*t

	Y_x_e_e_V0 <- function(X.t0,dxdede.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,block))

		tmp <- c()

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    tmp.drift <- dxdede.drift[[i1]]

		    for(d in 1:d.size){
			tmp[d] <- eval(tmp.drift[[d]])
		    }

		    result[i1,,t] <- invY[,,t] %*% tmp
		  }
		}

		return(result)
	}


	#d*d*r*t

	Y_x_e_e_V <- function(X.t0,dxdede.diffusion,env){

		d.size <- env$d.size
		r.size <- env$r.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,r.size,block))

		tmp <- matrix(0,d.size,r.size)

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    tmp.diffusion <- dxdede.diffusion[[i1]] 

		    for(i in 1:d.size){
			for(r in 1:r.size){
			  tmp[i,r] <- eval(tmp.diffusion[[i]][[r]])
			}
		    }
		    result[i1,,,t] <- invY[,,t] %*% tmp
		  }
		}

		return(result)
	}


	#d*t

	Y_e_e_e_V0 <- function(X.t0,dedede.drift,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,d.size,block)

		tmp <- c()

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(d in 1:d.size){
		    tmp[d] <- eval(dedede.drift[[d]])
		  }

		  result[,t] <- invY[,,t] %*% tmp
		}

		return(result)
	}


	#d*r*t

	Y_e_e_e_V <- function(X.t0,dedede.diffusion,env){

		d.size <- env$d.size
		r.size <- env$r.size
		state <- env$state
		pars <- env$pars
		invY <- env$invY
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,r.size,block))

		tmp <- matrix(0,d.size,r.size)

		assign(pars[1],0)

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i in 1:d.size){
		    for(r in 1:r.size){
			tmp[i,r] <- eval(dedede.diffusion[[i]][[r]])
		    }
		  }
		  result[,,t] <- invY[,,t] %*% tmp
		}

		return(result)
	}


	D0_t <- function(tmpY, get_Y_D, get_Y_e_V){

		first <- get_Y_D

		second.coef <- tmpY
		second <- get_Y_e_V

		return(list(first=first,second.coef=second.coef,second=second))
	}


	#i*t

	e_t <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x1_x2_V0, get_Y_x_e_V0, get_Y_e_e_V0, env){

		d.size <- env$d.size
		delta <- env$delta
		block <- env$block
		Diff <- env$Diff

		result <- matrix(0,d.size,block)

		first <- matrix(0,d.size,block)
		second <- matrix(0,d.size,block)
		third <- matrix(0,d.size,block)
		fourth <- matrix(0,d.size,block)

		first.tmp <- array(0,dim=c(d.size,d.size,block))
		second.tmp <- array(0,dim=c(d.size,d.size,block))
		third.tmp <- array(0,dim=c(d.size,d.size,block))

		for(i in 1:d.size){
		  for(j in 1:d.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(j1 in 1:d.size){
			    for(j2 in 1:d.size){
				tmp1 <- b1_b2(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)
				tmp2 <- double(block)

				for(t in 2:block){
				  tmp2[t] <- (get_Y_x1_x2_V0[i1,i2,j,] * 
						  tmpY[i1,j1,] * tmpY[i2,j2,] * 
						  tmp1) %*% Diff[t,] * delta
				}

				first.tmp[i,j,] <- first.tmp[i,j,] + tmp2
			    }
			  }

			  tmp3 <- double(block)

			  for(t in 2:block){
			    tmp3[t] <- (get_Y_x1_x2_V0[i1,i2,j,] * get_Y_D[i1,] *
					    get_Y_D[i2,]) %*% Diff[t,] * delta
			  }

			  second.tmp[i,j,] <- second.tmp[i,j,] + tmp3
			}

			tmp4 <- double(block)

			for(t in 2:block){
			  tmp4[t] <- (get_Y_x_e_V0[i1,j,] * 
					  get_Y_D[i1,]) %*% Diff[t,] * delta
			}

			third.tmp[i,j,] <- third.tmp[i,j,] + tmp4
		    }

		    tmp5 <- double(block)

		    for(t in 2:block){
			tmp5[t] <- get_Y_e_e_V0[j,] %*% Diff[t,] * delta
		    }

		    first[i,] <- first[i,] + first.tmp[i,j,] * tmpY[i,j,]
		    second[i,] <- second[i,] + second.tmp[i,j,] * tmpY[i,j,]
		    third[i,] <- third[i,] + 2 * third.tmp[i,j,] * tmpY[i,j,]
		    fourth[i,] <- fourth[i,] +  tmp5 * tmpY[i,j,]
		  }
		}

		result <- first + second + third + fourth

		return(result)
	}


	#j1*j*s

	U_t <- function(tmpY, get_Y_D, get_Y_x1_x2_V0, get_Y_x_e_V0, env){

		d.size <- env$d.size
		block <- env$block

		result <- array(0,dim=c(d.size,d.size,block))

		first <- array(0,dim=c(d.size,d.size,block))
		second <- array(0,dim=c(d.size,d.size,block))

		for(j1 in 1:d.size){
		  for(j in 1:d.size){
		    for(s in 1:block){
			for(i1 in 1:d.size){
			  tmp.Y <- as.matrix(tmpY[,,s])

			  for(i2 in 1:d.size){
			    first[j1,j,s] <- first[j1,j,s] + get_Y_x1_x2_V0[i1,i2,j,s] *
						   get_Y_D[i1,s] * tmp.Y[i2,j1]
			  }

			  second[j1,j,s] <- second[j1,j,s] + get_Y_x_e_V0[i1,j,s] *
						  tmp.Y[i1,j1]
			}
		    }
		  }
		}

		result <- first + second

		return(result)
	}


	#j1*j*r*t

	U_hat_t <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x1_x2_V0, get_Y_x_e_V, env){

		d.size <- env$d.size
		r.size <- env$r.size
		block <- env$block

		result <- array(0,dim=c(d.size,d.size,r.size,block))

		first <- array(0,dim=c(d.size,d.size,r.size,block))
		second <- array(0,dim=c(d.size,d.size,r.size,block))

		for(j1 in 1:d.size){
		  for(j in 1:d.size){
		    for(s in 1:block){
			tmp.Y <- as.matrix(tmpY[,,s])

			for(i1 in 1:d.size){
			  first[j1,j,,s] <- first[j1,j,,s] + get_Y_x_e_V[i1,j,,s] *
						  tmp.Y[i1,j1]
			}
		    }
		  }
		}

		for(j1 in 1:d.size){
		  for(j in 1:d.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(j2 in 1:d.size){
			    tmp1 <- get_Y_x1_x2_V0[i1,i2,j,] *
					tmpY[i1,j1,] * tmpY[i2,j2,]

			    tmp2 <- I0(tmp1,env)
			    tmp3 <- matrix(0,r.size,block)

			    for(r in 1:r.size){
				tmp3[r,] <- tmp2 * get_Y_e_V[j2,r,]
			    }

			    second[j1,j,,] <- second[j1,j,,] + tmp3
			  }
			}
		    }
		  }
		}

		result <- first - second

		return(result)
	}


	E0_t <- function(tmpY, get_Y_e_V, get_Y_x1_x2_V0, get_Y_e_e_V, get_e_t, get_U_t, get_U_hat_t, env){

		d.size <- env$d.size
		r.size <- env$r.size
		block <- env$block

		first <- get_e_t

		second.coef <- array(0,dim=c(d.size,d.size,d.size,block))	# i, j1, j2, t

		for(i in 1:d.size){
		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){
			for(j in 1:d.size){

			  tmp2 <- double(block)

			  for(i1 in 1:d.size){
			    for(i2 in 1:d.size){

				tmp1 <- get_Y_x1_x2_V0[i1,i2,j,] *
					  tmpY[i1,j1,] * tmpY[i2,j2,]

				tmp2 <- tmp2 + tmp1
			    }
			  }

			  second.coef[i,j1,j2,] <- second.coef[i,j1,j2,] +
							   2 * tmpY[i,j,] * I0(tmp2,env)
			}
		    }
		  }
		}

		second <- list()
		second[[1]] <- get_Y_e_V
		second[[2]] <- get_Y_e_V

		third.coef <- array(0,dim=c(d.size,d.size,block))	# i, j1, t

		for(j1 in 1:d.size){
		  for(j in 1:d.size){
		    tmp3 <- I0(get_U_t[j1,j,],env)

		    for(i in 1:d.size){
			third.coef[i,j1,] <- third.coef[i,j1,] +
						   2 * tmpY[i,j,] * tmp3
		    }
		  }
		}

		third <- get_Y_e_V

		fourth.coef <- - 2 * tmpY
		fourth <- array(0,dim=c(d.size,r.size,block))

		for(j in 1:d.size){
		  for(j1 in 1:d.size){
		    tmp4 <- I0(get_U_t[j1,j,],env)

		    for(r in 1:r.size){
			fourth[j,r,] <- fourth[j,r,] + tmp4 * get_Y_e_V[j1,r,]
		    }
		  }
		}

		fifth.coef <- 2 * tmpY
		fifth <- list()
		fifth[[1]] <- get_U_hat_t		
		fifth[[2]] <- get_Y_e_V

		sixth.coef <- tmpY
		sixth <- get_Y_e_e_V

		return(list(first=first,second.coef=second.coef,second=second,
				third.coef=third.coef,third=third,
				fourth.coef=fourth.coef,fourth=fourth,
				fifth.coef=fifth.coef,fifth=fifth,
				sixth.coef=sixth.coef,sixth=sixth))
	}


	#l*t

	e_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,k.size,block)

		assign(pars[1],0)

		de.f0 <- list()

		for(k in 1:k.size){
		  de.f0[[k]] <- parse(text=deparse(D(tmp.f[[1]][k],pars[1])))
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    result[k,t] <- eval(de.f0[[k]])
		  }
		}

		return(result)
	}


	#l*r*t

	e_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,r.size,block))

		assign(pars[1],0)

		de.f <- list()

		for(k in 1:k.size){
		  de.f[[k]] <- list()

		  for(r in 1:r.size){
		    de.f[[k]][r] <- parse(text=deparse(D(tmp.f[[r+1]][k],pars[1])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(r in 1:r.size){
			result[k,r,t] <- eval(de.f[[k]][[r]])
		    }
		  }
		}

		return(result)
	}


	#l*i*t

	x_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,block))

		assign(pars[1],0)

		dx.f0 <- list()

		for(k in 1:k.size){
		  dx.f0[[k]] <- list()

		  for(i in 1:d.size){
		    dx.f0[[k]][i] <- parse(text=deparse(D(tmp.f[[1]][k],state[i])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i in 1:d.size){
			result[k,i,t] <- eval(dx.f0[[k]][[i]])
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*t

	x1_x2_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,block))

		assign(pars[1],0)

		dxdx.f0 <- list()

		for(k in 1:k.size){
		  dxdx.f0[[k]] <- list()

		  for(i1 in 1:d.size){
		    dxdx.f0[[k]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdx.f0[[k]][[i1]][i2] <- parse(text=deparse(D(D(tmp.f[[1]][k],state[i2]),state[i1])))
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  result[k,i1,i2,t] <- eval(dxdx.f0[[k]][[i1]][[i2]])
			}
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*r*t

	x1_x2_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,r.size,block))

		assign(pars[1],0)

		dxdx.f <- list()

		for(k in 1:k.size){
		  dxdx.f[[k]] <- list()

		  for(i1 in 1:d.size){
		    dxdx.f[[k]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdx.f[[k]][[i1]][[i2]] <- list()

			for(r in 1:r.size){
			  dxdx.f[[k]][[i1]][[i2]][r] <- parse(text=deparse(D(D(tmp.f[[r+1]][k],state[i2]),state[i1])))
			}
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(r in 1:r.size){
			    result[k,i1,i2,r,t] <- eval(dxdx.f[[k]][[i1]][[i2]][[r]])
			  }
			}
		    }
		  }
		}

		return(result)
	}


	#l*i*t

	x_e_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,block))

		assign(pars[1],0)

		dxde.f0 <- list()

		for(k in 1:k.size){
		  dxde.f0[[k]] <- list()

		  for(i in 1:d.size){
		    dxde.f0[[k]][i] <- parse(text=deparse(D(D(tmp.f[[1]][k],pars[1]),state[i])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i in 1:d.size){
			result[k,i,t] <- eval(dxde.f0[[k]][[i]])
		    }
		  }
		}

		return(result)
	}


	#l*i*r*t

	x_e_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,r.size,block))

		assign(pars[1],0)

		dxde.f <- list()

		for(k in 1:k.size){
		  dxde.f[[k]] <- list()

		  for(i in 1:d.size){
		    dxde.f[[k]][[i]] <- list()

		    for(r in 1:r.size){
			dxde.f[[k]][[i]][r] <- parse(text=deparse(D(D(tmp.f[[r+1]][k],pars[1]),state[i])))
		    }
		  }
		}

		 for(t in 1:block){
		   for(d in 1:d.size){
		     assign(state[d],X.t0[my.range[t],d])
		   }

		  for(k in 1:k.size){
		    for(i in 1:d.size){
			for(r in 1:r.size){
			  result[k,i,r,t] <- eval(dxde.f[[k]][[i]][[r]])
			}
		    }
		  }
		}

		return(result)
	}


	#l*t

	e_e_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,k.size,block)

		assign(pars[1],0)

		dede.f0 <- list()

		for(k in 1:k.size){
		  dede.f0[[k]] <- parse(text=deparse(D(D(tmp.f[[1]][k],pars[1]),pars[1])))
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    result[k,t] <- eval(dede.f0[[k]])
		  }
		}

		return(result)
	}


	#l*r*t

	e_e_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,r.size,block))

		assign(pars[1],0)

		dede.f <- list()

		for(k in 1:k.size){
		  dede.f[[k]] <- list()

		  for(r in 1:r.size){
		    dede.f[[k]][r] <- parse(text=deparse(D(D(tmp.f[[r+1]][k],pars[1]),pars[1])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(r in 1:r.size){
			result[k,r,t] <- eval(dede.f[[k]][[r]])
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*i3*t

	x1_x2_x3_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,d.size,block))

		assign(pars[1],0)

		dxdxdx.f0 <- list()

		for(k in 1:k.size){
		  dxdxdx.f0[[k]] <- list()

		  for(i1 in 1:d.size){
		    dxdxdx.f0[[k]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdxdx.f0[[k]][[i1]][[i2]] <- list()

			for(i3 in 1:d.size){
			  dxdxdx.f0[[k]][[i1]][[i2]][i3] <- parse(text=deparse(D(D(D(tmp.f[[1]][k],state[i3]),state[i2]),state[i1])))
			}
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(i3 in 1:d.size){
			    result[k,i1,i2,i3,t] <- eval(dxdxdx.f0[[k]][[i1]][[i2]][[i3]])
			  }
			}
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*i3*r*t

	x1_x2_x3_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,d.size,r.size,block))

		assign(pars[1],0)

		dxdxdx.f <- list()

		for(k in 1:k.size){
		  dxdxdx.f[[k]] <- list()

		  for(i1 in 1:d.size){
		    dxdxdx.f[[k]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdxdx.f[[k]][[i1]][[i2]] <- list()

			for(i3 in 1:d.size){
			  dxdxdx.f[[k]][[i1]][[i2]][[i3]] <- list()

			  for(r in 1:r.size){
			    dxdxdx.f[[k]][[i1]][[i2]][[i3]][r] <- parse(text=deparse(D(D(D(tmp.f[[r+1]][k],state[i3]),state[i2]),state[i1])))
			  }
			}
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(i3 in 1:d.size){
			    for(r in 1:r.size){
				result[k,i1,i2,i3,r,t] <- eval(dxdxdx.f[[k]][[i1]][[i2]][[i3]][[r]])
			    }
			  }
			}
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*t

	x1_x2_e_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,block))

		assign(pars[1],0)

		dxdxde.f0 <- list()

		for(k in 1:k.size){
		  dxdxde.f0[[k]] <- list()

		  for(i1 in 1:d.size){
		    dxdxde.f0[[k]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdxde.f0[[k]][[i1]][i2] <- parse(text=deparse(D(D(D(tmp.f[[1]][k],pars[1]),state[i2]),state[i1])))
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  result[k,i1,i2,t] <- eval(dxdxde.f0[[k]][[i1]][[i2]])
			}
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*r*t

	x1_x2_e_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,r.size,block))

		assign(pars[1],0)

		dxdxde.f <- list()

		for(k in 1:k.size){
		  dxdxde.f[[k]] <- list()

		  for(i1 in 1:d.size){
		    dxdxde.f[[k]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdxde.f[[k]][[i1]][[i2]] <- list()

			for(r in 1:r.size){
			  dxdxde.f[[k]][[i1]][[i2]][r] <- parse(text=deparse(D(D(D(tmp.f[[r+1]][k],pars[1]),state[i2]),state[i1])))
			}
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(r in 1:r.size){
			    result[k,i1,i2,r,t] <- eval(dxdxde.f[[k]][[i1]][[i2]][[r]])
			  }
			}
		    }
		  }
		}

		return(result)
	}


	#l*i*t

	x_e_e_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,block))

		assign(pars[1],0)

		dxdede.f0 <- list()

		for(k in 1:k.size){
		  dxdede.f0[[k]] <- list()

		  for(i in 1:d.size){
		    dxdede.f0[[k]][i] <- parse(text=deparse(D(D(D(tmp.f[[1]][k],pars[1]),pars[1]),state[i])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i in 1:d.size){
			result[k,i,t] <- eval(dxdede.f0[[k]][[i]])
		    }
		  }
		}

		return(result)
	}


	#l*i*r*t

	x_e_e_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,r.size,block))

		assign(pars[1],0)

		dxdede.f <- list()

		for(k in 1:k.size){
		  dxdede.f[[k]] <- list()

		  for(i in 1:d.size){
		    dxdede.f[[k]][[i]] <- list()

		    for(r in 1:r.size){
			dxdede.f[[k]][[i]][r] <- parse(text=deparse(D(D(D(tmp.f[[r+1]][k],pars[1]),pars[1]),state[i])))
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(i in 1:d.size){
			for(r in 1:r.size){
			  result[k,i,r,t] <- eval(dxdede.f[[k]][[i]][[r]])
			}
		    }
		  }
		}

		return(result)
	}


	#l*t

	e_e_e_f0 <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,k.size,block)

		assign(pars[1],0)

		dedede.f0 <- list()

		for(k in 1:k.size){
		  dedede.f0[[k]] <- parse(text=deparse(D(D(D(tmp.f[[1]][k],pars[1]),pars[1]),pars[1])))
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    result[k,t] <- eval(dedede.f0[[k]])
		  }
		}

		return(result)
	}


	#l*r*t

	e_e_e_f <- function(X.t0,tmp.f,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,r.size,block))

		dedede.f <- list()

		for(k in 1:k.size){
		  dedede.f[[k]] <- list()

		  for(r in 1:r.size){
		    dedede.f[[k]][r] <- parse(text=deparse(D(D(D(tmp.f[[r+1]][k],pars[1]),pars[1]),pars[1])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(k in 1:k.size){
		    for(r in 1:r.size){
			result[k,r,t] <- eval(dedede.f[[k]][[r]])
		    }
		  }
		}

		return(result)
	}


##p.19

	#l*t

	F_t <- function(tmpY, get_Y_e_V, get_Y_D, get_x1_x2_f0, get_x_e_f0, get_e_e_f0, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		delta <- env$delta
		block <- env$block
		Diff <- env$Diff

		result <- matrix(0,k.size,block)

		first <- matrix(0,k.size,block)
		second <- matrix(0,k.size,block)
		third <- matrix(0,k.size,block)
		fourth <- matrix(0,k.size,block)

		for(l in 1:k.size){
		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			for(j1 in 1:d.size){
			  for(j2 in 1:d.size){
			    tmp1 <- b1_b2(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)
			    tmp2 <- double(block)

			    for(t in 2:block){
				tmp2[t] <- (get_x1_x2_f0[l,i1,i2,] * 
						tmpY[i1,j1,] * tmpY[i2,j2,] *
						tmp1) %*% Diff[t,] * delta
			    }

			    first[l,] <- first[l,] + tmp2
			  }
			}

			tmp3 <- double(block)

			for(t in 2:block){
			  tmp3[t] <- (get_x1_x2_f0[l,i1,i2,] * get_Y_D[i1,] *
					  get_Y_D[i2,]) %*% Diff[t,] * delta
			}

			second[l,] <- second[l,] + tmp3
		    }

		    tmp4 <- double(block)

		    for(t in 2:block){
			tmp4[t] <- (get_x_e_f0[l,i1,] * 
					get_Y_D[i1,]) %*% Diff[t,] * delta
		    }

		    third[l,] <- third[l,] + 2 * tmp4
		  }

		  tmp5 <- double(block)

		  for(t in 2:block){
		    tmp5[t] <- get_e_e_f0[l,] %*% Diff[t,] * delta
		  }

		  fourth[l,] <- fourth[l,] + tmp5
		}

		result <- first + second + third + fourth

		return(result)
	}


	#l*j1*t

	W_t <- function(tmpY, get_Y_D, get_x1_x2_f0, get_x_e_f0, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		result <- array(0,dim=c(k.size,d.size,block))

		first <- array(0,dim=c(k.size,d.size,block))
		second <- array(0,dim=c(k.size,d.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  first[l,j1,] <- first[l,j1,] + get_x1_x2_f0[l,i1,i2,] *
						get_Y_D[i1,] * tmpY[i2,j1,]
			}

			second[l,j1,] <- second[l,j1,] + get_x_e_f0[l,i1,] *
					     tmpY[i1,j1,]
		    }
		  }
		}

		result <- first + second

		return(result)
	}


	#l*j1*r*t

	W_hat_t <- function(tmpY, get_Y_e_V, get_x1_x2_f0, get_x_e_f, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		result <- array(0,dim=c(k.size,d.size,r.size,block))

		first <- array(0,dim=c(k.size,d.size,r.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(r in 1:r.size){
			for(i1 in 1:d.size){
			  first[l,j1,r,] <- first[l,j1,r,] + get_x_e_f[l,i1,r,] *
						  tmpY[i1,j1,]
			}
		    }
		  }
		}

		second <- array(0,dim=c(k.size,d.size,r.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			tmp2 <- double(block)
			tmp3 <- matrix(0,r.size,block)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){

			    tmp1 <- get_x1_x2_f0[l,i1,i2,] *
					tmpY[i1,j1,] * tmpY[i2,j2,]

			    tmp2 <- tmp2 + tmp1
			  }
			}

			for(r in 1:r.size){
			  tmp3[r,] <- I0(tmp2,env) * get_Y_e_V[j2,r,]
			}

			second[l,j1,,] <- second[l,j1,,] + tmp3
		    }
		  }
		}

		result <- first - second

		return(result)
	}


	#first:l, second:l*j2*r*t&j2*r*t, third:l*r*t

	F_tilde1_1 <- function(tmpY, get_Y_e_V, get_x1_x2_f0, get_e_e_f, get_F_t, get_W_t, get_W_hat_t, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first <- get_F_t[block]

		second <- list()
		second[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))
		second[[2]] <- get_Y_e_V

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			tmp2 <- double(block)
			tmp3 <- double(1)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){
			    tmp1 <- get_x1_x2_f0[l,i1,i2,] *
					tmpY[i1,j1,] * tmpY[i2,j2,]

			    tmp2 <- tmp2 + tmp1
			  }
			}

			tmp3 <- tmp3 + 2 * I0(tmp2,env)[block]

			for(r in 1:r.size){
			  second[[1]][l,j2,r,] <- second[[1]][l,j2,r,] +
							  tmp3 * get_Y_e_V[j1,r,]
			}
		    }
		  }
		}

		third.coef <- array(0,dim=c(k.size,d.size,block))	#l, j1, t
		third <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){

		    third.coef[l,j1,] <- 2 * I0(get_W_t[l,j1,],env)

		    for(r in 1:r.size){
			third[l,r,] <- third[l,r,] +
					   third.coef[l,j1,block] * get_Y_e_V[j1,r,]
		    }
		  }
		}

		fourth <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(r in 1:r.size){
			fourth[l,r,] <- fourth[l,r,] -
					    third.coef[l,j1,] * get_Y_e_V[j1,r,]
		    }
		  }
		}

		fifth <- list()
		fifth[[1]] <- 2 * get_W_hat_t
		fifth[[2]] <- get_Y_e_V

		sixth <- get_e_e_f

		second[[1]] <- second[[1]] + fifth[[1]]
		third <- third + fourth + sixth

		return(list(first=first,second=second,third=third))
	}


	#first:d*k*k*t, second:d*k*t, third:d*t

	E0_bar <- function(get_E0_t,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- get_E0_t$first

		second.coef <- get_E0_t$second.coef
		second.wtmp <- get_E0_t$second

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,d.size,block)

		n1 <- length(second.coef)
		n2 <- sum(second.coef == 0)

		n3 <- length(second.wtmp[[1]])
		n4 <- sum(second.wtmp[[1]] == 0)

		n5 <- length(second.wtmp[[2]])
		n6 <- sum(second.wtmp[[2]] == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) * (n5 - n6 != 0)

		for(j1 in 1:1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j2 in 1:d.size){
		    tmp1 <- I_12(second.wtmp[[1]][j1,,],second.wtmp[[2]][j2,,],env)

		    for(i in 1:d.size){
			for(t in 1:block){
			  second.tmp$first[i,,,t] <- second.tmp$first[i,,,t] + 
							     second.coef[i,j1,j2,t] * tmp1$first[,,t]
			  second.tmp$second[i,t] <- second.tmp$second[i,t] + 
							    second.coef[i,j1,j2,t] * tmp1$second[t]
			}
		    }
		  }
		}

		third.coef <- get_E0_t$third.coef
		third.wtmp <- get_E0_t$third

		third.tmp <- array(0,dim=c(d.size,k.size,block))

		n7 <- length(third.coef)
		n8 <- sum(third.coef == 0)

		n9 <- length(third.wtmp)
		n10 <- sum(third.wtmp == 0)

		n <- (n7 - n8 != 0) * (n9 - n10 != 0)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){
		    tmp2 <- I_1(third.wtmp[j1,,],env)

		    for(t in 1:block){
			third.tmp[i,,t] <- third.tmp[i,,t] + third.coef[i,j1,t] * tmp2[,t]
		    }
		  }
		}

		fourth.coef <- get_E0_t$fourth.coef
		fourth.wtmp <- get_E0_t$fourth

		fourth.tmp <- array(0,dim=c(d.size,k.size,block))

		fifth.coef <- get_E0_t$fifth.coef
		fifth.wtmp <- get_E0_t$fifth

		fifth.tmp <- list()
		fifth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		fifth.tmp$second <- matrix(0,d.size,block)

		sixth.coef <- get_E0_t$sixth.coef
		sixth.wtmp <- get_E0_t$sixth

		sixth.tmp <- array(0,dim=c(d.size,k.size,block))

		n11 <- length(sixth.wtmp)
		n12 <- sum(sixth.wtmp == 0)

		n <- n11 - n12

		for(i in 1:d.size){
		  for(j in 1:d.size){

		    tmp3 <- I_1(fourth.wtmp[j,,],env)

		    for(k in 1:k.size){
			fourth.tmp[i,k,] <- fourth.tmp[i,k,] + 
						  fourth.coef[i,j,] * tmp3[k,]
		    }

		    for(j1 in 1:d.size){

			tmp4 <- I_12(fifth.wtmp[[1]][j1,j,,],fifth.wtmp[[2]][j1,,],env)

			for(t in 1:block){
			  fifth.tmp$first[i,,,t] <- fifth.tmp$first[i,,,t] + 
							    fifth.coef[i,j,t] * tmp4$first[,,t]
			  fifth.tmp$second[i,t] <- fifth.tmp$second[i,t] + 
							   fifth.coef[i,j,t] * tmp4$second[t]
			}
		    }
		  }

		  if(n == 0){
		    next
		  }

		  tmp5 <- I_1(sixth.wtmp[j,,],env)

		  for(i in 1:d.size){
		    for(k in 1:k.size){
			sixth.tmp[i,k,] <- sixth.tmp[i,k,] + 
						 sixth.coef[i,j,] * tmp5[k,]
		    }
		  }
		}

		first <- second.tmp$first + fifth.tmp$first
		second <- third.tmp + fourth.tmp + sixth.tmp
		third <- first.tmp + second.tmp$second + fifth.tmp$second

		return(list(first=first,second=second,third=third))
	}


	E0_bar_x <- function(x,get_E0_bar,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size

		result <- double(d.size)

		for(i in 1:d.size){
		  first.tmp <- list()
		  first.tmp$first <- get_E0_bar$first[i,,,]
		  first.tmp$second <- get_E0_bar$third[i,]
		  second.tmp <- get_E0_bar$second[i,,]

		  result1 <- I_12_x(x,first.tmp,env)
		  result2 <- I_1_x(x,second.tmp,env)

		  result[i] <- result1 + result2
		}

		return(result)
	}


	#first:l, second:l*j2*r*t&j2*r*t, third:l*r*t

	F_tilde1_2 <- function(get_E0_t,get_x_f0,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		result1 <- get_E0_t$first

		first <- double(k.size)

		for(l in 1:k.size){
		  for(i in 1:d.size){
		    tmp1 <- get_x_f0[l,i,] * result1[i,]
		    first[l] <- first[l] + I0(tmp1,env)[block]
		  }
		}

		result2.coef <- get_E0_t$second.coef
		result2 <- get_E0_t$second[[1]]

		second <- list()
		second[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))	#l,j2,r,t
		second[[2]] <- result2

		third <- list()
		third[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))
		third[[2]] <- result2

		for(l in 1:k.size){
		  for(j2 in 1:d.size){
		    for(j1 in 1:d.size){

			tmp3 <- double(block)

			for(i in 1:d.size){
			  tmp2 <- get_x_f0[l,i,] * result2.coef[i,j1,j2,]

			  tmp3 <- tmp3 + tmp2
			}

			tmp4 <- I0(tmp3,env)

			for(r in 1:r.size){
			  second[[1]][l,j2,r,] <- second[[1]][l,j2,r,] +
							  tmp4[block] * result2[j1,r,]

			  third[[1]][l,j2,r,] <- third[[1]][l,j2,r,] -
							 tmp4 * result2[j1,r,]
			}
		    }
		  }
		}

		result4.coef <- get_E0_t$third.coef

		fourth <- array(0,dim=c(k.size,r.size,block))

		fifth <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){

		    tmp6 <- double(block)

		    for(i in 1:d.size){
			tmp5 <- get_x_f0[l,i,] * result4.coef[i,j1,]

			tmp6 <- tmp6 + tmp5
		    }

		    tmp7 <- I0(tmp6,env)

		    for(r in 1:r.size){
			fourth[l,r,] <- fourth[l,r,] + tmp7[block] * result2[j1,r,]

			fifth[l,r,] <- fifth[l,r,] - tmp7 * result2[j1,r,]
		    }
		  }
		}

		result6.coef <- get_E0_t$fourth.coef
		result6 <- get_E0_t$fourth

		result8 <- get_E0_t$fifth[[1]]

		result10 <- get_E0_t$sixth

		sixth <- array(0,dim=c(k.size,r.size,block))

		seventh <- array(0,dim=c(k.size,r.size,block))

		eighth <- list()
		eighth[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))
		eighth[[2]] <- result2

		ninth <- list()
		ninth[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))
		ninth[[2]] <- result2

		tenth <- array(0,dim=c(k.size,r.size,block))

		eleventh <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j in 1:d.size){
		    tmp9 <- double(block)

		    for(i in 1:d.size){
			tmp8 <- get_x_f0[l,i,] * result6.coef[i,j,]

			tmp9 <- tmp9 + tmp8
		    }

		    tmp10 <- I0(tmp9,env)

		    for(r in 1:r.size){
			sixth[l,r,] <- sixth[l,r,] + tmp10[block] * result6[j,r,]

			seventh[l,r,] <- seventh[l,r,] - tmp10 * result6[j,r,]

			for(j1 in 1:d.size){
			  eighth[[1]][l,j1,r,] <- eighth[[1]][l,j1,r,] -
							  tmp10[block] * result8[j1,j,r,]

			  ninth[[1]][l,j1,r,] <- ninth[[1]][l,j1,r,] +
							 tmp10 * result8[j1,j,r,]
			}

			tenth[l,r,] <- tenth[l,r,] -
					   tmp10[block]/2 * result10[j,r,]

			eleventh[l,r,] <- eleventh[l,r,] +
						tmp10/2 * result10[j,r,]
		    }
		  }
		}

		second[[1]] <- second[[1]] + third[[1]] + eighth[[1]] +
				   ninth[[1]]

		third <- fourth + fifth + sixth + seventh + tenth + eleventh

		return(list(first=first,second=second,third=third))
	}


	#first:l*t, second.coef:l*j1*j2*t, second:j1*r*t&j2*r*t,
	#third.coef:l*j1*t, third:j1*r*t, fourth.coef:l*j*t, fourth:j*r*t,
	#fifth.coef:l*j*t, fifth:j1*j*r*t&j1*r*t, sixth.coef:l*j*t,
	#sixth:j*r*t


##p.22

	#l*i*t

	x_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,block))

		assign(pars[1],0)

		dx.F <- list()

		for(l in 1:k.size){
		  dx.F[[l]] <- list()

		  for(i in 1:d.size){
		    dx.F[[l]][i] <- parse(text=deparse(D(F[l],state[i])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    for(i in 1:d.size){
			result[l,i,t] <- eval(dx.F[[l]][[i]])
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*t

	x1_x2_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,block))

		assign(pars[1],0)

		dxdx.F <- list()

		for(l in 1:k.size){
		  dxdx.F[[l]] <- list()

		  for(i1 in 1:d.size){
		    dxdx.F[[l]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdx.F[[l]][[i1]][i2] <- parse(text=deparse(D(D(F[l],state[i2]),state[i1])))
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  result[l,i1,i2,t] <- eval(dxdx.F[[l]][[i1]][[i2]])
			}
		    }
		  }
		}

		return(result)
	}


	#l*i*t

	x_e_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,block))

		assign(pars[1],0)

		dxde.F <- list()

		for(l in 1:k.size){
		  dxde.F[[l]] <- list()

		  for(i in 1:d.size){
		    dxde.F[[l]][i] <- parse(text=deparse(D(D(F[l],pars[1]),state[i])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    for(i in 1:d.size){
			result[l,i,t] <- eval(dxde.F[[l]][[i]])
		    }
		  }
		}

		return(result)
	}


	#l*t

	e_e_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,k.size,block)

		assign(pars[1],0)

		dede.F <- list()

		for(l in 1:k.size){
		  dede.F[[l]] <- parse(text=deparse(D(D(F[l],pars[1]),pars[1])))
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    result[l,t] <- eval(dede.F[[l]])
		  }
		}

		return(result)
	}


	#l*i1*i2*i3*t

	x1_x2_x3_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,d.size,block))

		assign(pars[1],0)

		dxdxdx.F <- list()

		for(l in 1:k.size){
		  dxdxdx.F[[l]] <- list()

		  for(i1 in 1:d.size){
		    dxdxdx.F[[l]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdxdx.F[[l]][[i1]][[i2]] <- list()

			for(i3 in 1:d.size){
			  dxdxdx.F[[l]][[i1]][[i2]][i3] <- parse(text=deparse(D(D(D(F[l],state[i3]),state[i2]),state[i1])))
			}
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(i3 in 1:d.size){
			    result[l,i1,i2,i3,t] <- eval(dxdxdx.F[[l]][[i1]][[i2]][[i3]])
			  }
			}
		    }
		  }
		}

		return(result)
	}


	#l*i1*i2*t

	x1_x2_e_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,d.size,block))

		assign(pars[1],0)

		dxdxde.F <- list()

		for(l in 1:k.size){
		  dxdxde.F[[l]] <- list()

		  for(i1 in 1:d.size){
		    dxdxde.F[[l]][[i1]] <- list()

		    for(i2 in 1:d.size){
			dxdxde.F[[l]][[i1]][i2] <- parse(text=deparse(D(D(D(F[l],pars[1]),state[i2]),state[i1])))
		    }
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  result[l,i1,i2,t] <- eval(dxdxde.F[[l]][[i1]][[i2]])
			}
		    }
		  }
		}

		return(result)
	}


	#l*i*t

	x_e_e_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(k.size,d.size,block))

		assign(pars[1],0)

		dxdede.F <- list()

		for(l in 1:k.size){
		  dxdede.F[[l]] <- list()

		  for(i in 1:d.size){
		    dxdede.F[[l]][i] <- parse(text=deparse(D(D(D(F[l],pars[1]),pars[1]),state[i])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    for(i in 1:d.size){
			result[l,i,t] <- eval(dxdede.F[[l]][[i]])
		    }
		  }
		}

		return(result)
	}


	#l*t

	e_e_e_F <- function(X.t0,F,env){

		d.size <- env$d.size
		k.size <- env$k.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- matrix(0,k.size,block)

		assign(pars[1],0)

		dedede.F <- list()

		for(l in 1:k.size){
		  dedede.F[[l]] <- parse(text=deparse(D(D(D(F[l],pars[1]),pars[1]),pars[1])))
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(l in 1:k.size){
		    result[l,t] <- eval(dedede.F[[l]])
		  }
		}

		return(result)
	}


	#first:l,second:l*j2*r*t&j2*r*t,third:l*r*t

	F_tilde1_3 <- function(get_E0_t,get_x_F,env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		result1 <- get_E0_t$first

		first <- double(k.size)

		for(l in 1:k.size){
		  for(i in 1:d.size){
		    tmp1 <- get_x_F[l,i,block] * result1[i,block]
		    first[l] <- first[l] + tmp1
		  }
		}

		result2.coef <- get_E0_t$second.coef
		result2 <- get_E0_t$second[[1]]

		second <- list()
		second[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))
		second[[2]] <- result2

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			tmp3 <- double(1)

			for(i in 1:d.size){
			  tmp2 <- get_x_F[l,i,block] * result2.coef[i,j1,j2,block]

			  tmp3 <- tmp3 + tmp2
			}

			for(r in 1:r.size){
			  second[[1]][l,j2,r,] <- second[[1]][l,j2,r,] +
							  tmp3 * result2[j1,r,]
			}
		    }
		  }
		}

		result3.coef <- get_E0_t$third.coef

		third <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){

		    tmp5 <- double(1)

		    for(i in 1:d.size){
			tmp4 <- get_x_F[l,i,block] * result3.coef[i,j1,block]

			tmp5 <- tmp5 + tmp4
		    }

		    for(r in 1:r.size){
			third[l,r,] <- third[l,r,] + tmp5 * result2[j1,r,]
		    }
		  }
		}

		result4.coef <- get_E0_t$fourth.coef
		result4 <- get_E0_t$fourth

		result5 <- get_E0_t$fifth[[1]]

		result6 <- get_E0_t$sixth

		fourth <- array(0,dim=c(k.size,r.size,block))

		fifth <- list()
		fifth[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))
		fifth[[2]] <- result2

		sixth <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j in 1:d.size){

		    tmp7 <- double(1)

		    for(i in 1:d.size){
			tmp6 <- get_x_F[l,i,block] * result4.coef[i,j,block]

			tmp7 <- tmp7 + tmp6
		    }

		    for(r in 1:r.size){
			fourth[l,r,] <- fourth[l,r,] + tmp7 * result4[j,r,]

			for(j1 in 1:d.size){
			  fifth[[1]][l,j1,r,] <- fifth[[1]][l,j1,r,] -
							 tmp7 * result5[j1,j,r,]
			}

			sixth[l,r,] <- sixth[l,r,] - tmp7/2 * result6[j1,r,]
		    }
		  }
		}

		second[[1]] <- second[[1]] + fifth[[1]]

		third <- third + fourth + sixth

		return(list(first=first,second=second,third=third))
	}


	#first:l,second:l*j2*r*t&j2*r*t,third:l*r*t

	F_tilde1_4 <- function(tmpY, get_Y_e_V, get_Y_D, get_x_F, get_x1_x2_F, get_x_e_F, get_e_e_F, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first <- double(k.size)

		for(l in 1:k.size){
		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			first[l] <- first[l] + get_x1_x2_F[l,i1,i2,block] *
					get_Y_D[i1,block] * get_Y_D[i2,block]
		    }
		  }
		}

		second <- double(k.size)

		for(l in 1:k.size){
		  for(i in 1:d.size){
		    second[l] <- second[l] + 2 * get_x_e_F[l,i,block] *
				     get_Y_D[i1,block]
		  }
		}

		third <- get_e_e_F[,block]

		fourth <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j1 in 1:d.size){

		    tmp1 <- double(1)

		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  tmp1 <- tmp1 + 2 * get_x1_x2_F[l,i1,i2,block] *
				    get_Y_D[i1,block] * tmpY[i2,j1,block]
			}
		    }

		    for(r in 1:r.size){
			fourth[l,r,] <- fourth[l,r,] + tmp1 * get_Y_e_V[j1,r,]
		    }
		  }
		}

		fifth <- array(0,dim=c(k.size,r.size,block))

		for(l in 1:k.size){
		  for(j in 1:d.size){

		    tmp2 <- double(1)

		    for(i in 1:d.size){
			tmp2 <- tmp2 + 2 * get_x_e_F[l,i,block] * tmpY[i,j,block]
		    }

		    for(r in 1:r.size){
			fifth[l,r,] <- fifth[l,r,] + tmp2 * get_Y_e_V[j,r,]
		    }
		  }
		}

		sixth <- list()
		sixth[[1]] <- array(0,dim=c(k.size,d.size,r.size,block))
		sixth[[2]] <- get_Y_e_V

		seventh <- double(k.size)

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			tmp3 <- double(1)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){
			    tmp3 <- tmp3 + 2 * get_x1_x2_F[l,i1,i2,block] *
					tmpY[i1,j1,block] * tmpY[i2,j2,block]
			  }
			}

			for(r in 1:r.size){
			  sixth[[1]][l,j2,r,] <- sixth[[1]][l,j2,r,] +
							 tmp3 * get_Y_e_V[j1,r,]
			}

			tmp4 <- b1_b2(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)[block]

			seventh[l] <- seventh[l] + tmp3/2 * tmp4
		    }
		  }
		}

		first <- first + second + third + seventh

		second <- list()
		second[[1]] <- sixth[[1]]
		second[[2]] <- sixth[[2]]

		third <- fourth + fifth

		return(list(first=first,second=second,third=third))
	}


	#result1:l, result2:l*j2*r*t/j2*r*t, result3:l*r*t

	F_tilde1 <- function(get_F_tilde1_1, get_F_tilde1_2, get_F_tilde1_3, get_F_tilde1_4){

		result1 <- get_F_tilde1_1$first + get_F_tilde1_2$first +
			     get_F_tilde1_3$first + get_F_tilde1_4$first

		result1 <- result1/2

		result2 <- list()
		result2[[1]] <- get_F_tilde1_1$second[[1]] +
				    get_F_tilde1_2$second[[1]] +
				    get_F_tilde1_3$second[[1]] +
				    get_F_tilde1_4$second[[1]]

		result2[[1]] <- result2[[1]]/2

		result2[[2]] <- get_F_tilde1_1$second[[2]]

		result3 <- get_F_tilde1_1$third + get_F_tilde1_2$third +
			     get_F_tilde1_3$third + get_F_tilde1_4$third

		result3 <- result3/2

		return(list(result1=result1,result2=result2,result3=result3))
	}


##p.31(I)

	#a:l, b:l*j*r*t/j*r*t, c:l*r*t

	#l1*l2

	ft_1 <- function(a1,a2,env){

		k.size <- env$k.size

		result <- matrix(0,k.size,k.size)

		for(l1 in 1:k.size){
			for(l2 in 1:k.size){
				result[l1,l2] <- a1[l1] * a2[l2]
			}
		}

		return(result)
	}


	#first:l1*l2*k*k, second:l1*l2

	ft_2 <- function(a1,b1,env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size))
		second <- matrix(0,k.size,k.size)

		for(l2 in 1:k.size){
		  for(j in 1:d.size){
		    tmp <- I_12(b1[[1]][l2,j,,],b1[[2]][j,,],env)

		    for(l1 in 1:k.size){
			first[l1,l2,,] <- first[l1,l2,,] + a1[l1] * tmp$first[,,block]

			second[l1,l2] <- second[l1,l2] + a1[l1] * tmp$second[block]
		    }
		  }
		}

		return(list(first=first,second=second))
	}


	#l1*l2*k

	ft_3 <- function(a1,c1,env){

		k.size <- env$k.size
		block <- env$block

		result <- array(0,dim=c(k.size,k.size,k.size))

		for(l2 in 1:k.size){
		  tmp <- I_1(c1[l2,,],env)

		  for(l1 in 1:k.size){
		    result[l1,l2,] <- a1[l1] * tmp[,block]
		  }
		}

		return(result)
	}


	#first:l1*l2*k*k*k*k, second:l1*l2*k*k*, third:l1*l2

	ft_4 <- function(b1,b2,env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size,k.size,k.size))
		second <- array(0,dim=c(k.size,k.size,k.size,k.size))
		third <- matrix(0,k.size,k.size)

		for(l1 in 1:k.size){
		  for(l2 in 1:k.size){
		    for(j1 in 1:d.size){
			for(j2 in 1:d.size){

			  tmp <- I_12_34(b1[[1]][l1,j1,,],b1[[2]][j1,,],
					     b2[[1]][l2,j2,,],b2[[2]][j2,,],env)

			  first[l1,l2,,,,] <- first[l1,l2,,,,] + tmp$first[,,,,block]

			  second[l1,l2,,] <- second[l1,l2,,] + tmp$second[,,block]

			  third[l1,l2] <- third[l1,l2] + tmp$third[block]
			}
		    }
		  }
		}

		return(list(first=first,second=second,third=third))
	}


	#first:l1*l2*k*k*k, second:l1*l2*k

	ft_5 <- function(b1,c1,env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size,k.size))
		second <- array(0,dim=c(k.size,k.size,k.size))

		for(l1 in 1:k.size){
		  for(l2 in 1:k.size){
		    for(j in 1:d.size){
			tmp <- I_1_23(c1[l2,,],b1[[1]][l1,j,,],b1[[2]][j,,],env)

			first[l1,l2,,,] <- first[l1,l2,,,] + tmp$first[,,,block]

			second[l1,l2,] <- second[l1,l2,] + tmp$second[,block]
		    }
		  }
		}

		return(list(first=first,second=second))
	}


	#first:l1*l2*k*k, second:l1*l2

	ft_6 <- function(c1,c2,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size))
		second <- matrix(0,k.size,k.size)

		for(l1 in 1:k.size){
		  for(l2 in 1:k.size){
		    tmp <- I_1_2(c1[l1,,],c2[l2,,],env)

		    first[l1,l2,,] <- tmp$first[,,block]

		    second[l1,l2] <- tmp$second[block]
		  }
		}

		return(list(first=first,second=second))
	}


	#first:l1*l2*k*k*k*k, second:l1*l2*k*k*k,
	#third:l1*l2*k*k, fourth:l1*l2*k, fifth:l1*l2

	F_tilde1__2 <- function(get_F_tilde1,env){

		k.size <- env$k.size

		temp <- list()
		temp$first <- array(0,dim=c(k.size,k.size,k.size,k.size,k.size,k.size))
		temp$second <- array(0,dim=c(k.size,k.size,k.size,k.size,k.size))
		temp$third <- array(0,dim=c(k.size,k.size,k.size,k.size))
		temp$fourth <- array(0,dim=c(k.size,k.size,k.size))
		temp$fifth <- matrix(0,k.size,k.size)

		calc.range <- c(1:3)

		tmp1 <- get_F_tilde1$result1

		n1 <- length(tmp1)
		n2 <- sum(tmp1 == 0)

		if(n1 == n2){
			calc.range <- calc.range[calc.range != 1]
		}

		tmp2 <- get_F_tilde1$result2[[1]]

		n3 <- length(tmp2)
		n4 <- sum(tmp2 == 0)

		tmp3 <- get_F_tilde1$result2[[2]]

		n5 <- length(tmp3)
		n6 <- sum(tmp3 == 0)

		n <- (n3 - n4 != 0) * (n5 - n6 != 0)

		if(n == 0){
			calc.range <- calc.range[calc.range != 2]
		}

		tmp3 <- get_F_tilde1$result3

		n7 <- length(tmp3)
		n8 <- sum(tmp3 == 0)

		if(n7 == n8){
			calc.range <- calc.range[calc.range != 3]
		}

		for(i1 in calc.range){
		  for(i2 in calc.range){

		    tmp1 <- switch(i1,"a","b","c")

		    tmp2 <- switch(i2,"a","b","c")

		    tmp <- paste(tmp1,tmp2,sep="")

		    result <- switch(tmp,"aa"=ft_1(get_F_tilde1[[i1]],get_F_tilde1[[i2]],env),
						 "ab"=ft_2(get_F_tilde1[[i1]],get_F_tilde1[[i2]],env),
						 "ac"=ft_3(get_F_tilde1[[i1]],get_F_tilde1[[i2]],env),
						 "ba"=ft_2(get_F_tilde1[[i2]],get_F_tilde1[[i1]],env),
						 "bb"=ft_4(get_F_tilde1[[i1]],get_F_tilde1[[i2]],env),
						 "bc"=ft_5(get_F_tilde1[[i1]],get_F_tilde1[[i2]],env),
						 "ca"=ft_3(get_F_tilde1[[i2]],get_F_tilde1[[i1]],env),
						 "cb"=ft_5(get_F_tilde1[[i2]],get_F_tilde1[[i1]],env),
						 "cc"=ft_6(get_F_tilde1[[i1]],get_F_tilde1[[i2]],env))

		    nlist <- length(result)

		    if(nlist != 2 && nlist != 3){

			tmp3 <- 7 - length(dim(result))	#4 or 5

			temp[[tmp3]] <- temp[[tmp3]] + result

		    }else if(nlist == 2){

			tmp4 <- 7 - length(dim(result[[1]]))	#2 or 3
			temp[[tmp4]] <- temp[[tmp4]] + result[[1]]

			tmp5 <- 7 - length(dim(result[[2]]))	#4 or 5
			temp[[tmp5]] <- temp[[tmp5]] + result[[2]]

		    }else{

			temp[[1]] <- temp[[1]] + result[[1]]

			temp[[3]] <- temp[[3]] + result[[2]]

			temp[[5]] <- temp[[5]] + result[[3]]
		    }
		  }
		}

		first <- temp$first
		second <- temp$second
		third <- temp$third
		fourth <- temp$fourth
		fifth <- temp$fifth

		return(list(first=first,second=second,third=third,
				fourth=fourth,fifth=fifth))
	}


	F_tilde1__2_x <- function(l1,l2,x,get_F_tilde1__2,env){

		k.size <- env$k.size

		first <- array(get_F_tilde1__2$first[l1,l2,,,,],
				   dim=c(k.size,k.size,k.size,k.size))

		result1 <- 0

		for(k1 in 1:k.size){
		for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			for(k4 in 1:k.size){
			  result1 <- result1 + first[k1,k2,k3,k4] *
					 x[k1] * x[k2] * x[k3] * x[k4]
			}
		}
		}
		}


		second <- array(get_F_tilde1__2$second[l1,l2,,,],
				    dim=c(k.size,k.size,k.size))

		result2 <- 0

		for(k1 in 1:k.size){
		for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			result2 <- result2 + second[k1,k2,k3] *
				     x[k1] * x[k2] * x[k3]
		    }
		  }
		}


		third <- matrix(get_F_tilde1__2$third[l1,l2,,],
				    k.size,k.size)

		result3 <- 0

		for(k1 in 1:k.size){
		for(k2 in 1:k.size){
		    result3 <- result3 + third[k1,k2] *
				   x[k1] * x[k2]
		  }
		}


		fourth <- get_F_tilde1__2$fourth[l1,l2,]

		result4 <- 0

		for(k1 in 1:k.size){
		  result4 <- result4 + fourth[k1] * x[k1]
		}


		fifth <- get_F_tilde1__2$fifth[l1,l2]

		result <- result1 + result2 + result3 + result4 + fifth

		return(result)
	}


##p.32(II)-1

#b:r*t, c:1*t


	#k*t

	C_b <- function(b,c1,c2,env){

		k.size <- env$k.size
		block <- env$block

		result <- matrix(0,k.size,block)

		tmp1 <- I_1(b,env)

		for(k in 1:k.size){
			tmp2 <- c2 * tmp1[k,]

			tmp3 <- I0(tmp2,env)

			result[k,] <- c1 * tmp3
		}

		return(result)
	}


	#first:k*k*t, second:t

	C_c <- function(b1,b2,c1,c2,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,block))

		tmp1 <- I_12(b1,b2,env)

		for(k1 in 1:k.size){
			for(k2 in 1:k.size){
				tmp2 <- c2 * tmp1$first[k1,k2,]

				tmp3 <- I0(tmp2,env)

				first[k1,k2,] <- c1 * tmp3
			}
		}

		second <- double(block)

		tmp4 <- c2 * tmp1$second

		tmp5 <- I0(tmp4,env)

		second <- c1 * tmp5

		return(list(first=first,second=second))
	}


	#first:k*k*t, second:t

	C_d <- function(b1,b2,c1,c2,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,block))

		tmp1 <- I_1_2(b1,b2,env)

		for(k1 in 1:k.size){
			for(k2 in 1:k.size){
				tmp2 <- c2 * tmp1$first[k1,k2,]

				tmp3 <- I0(tmp2,env)

				first[k1,k2,] <- c1 * tmp3
			}
		}

		tmp4 <- c2 * tmp1$second

		tmp5 <- I0(tmp4,env)

		second <- c1 * tmp5

		return(list(first=first,second=second))
	}


	#first:k*k*k*t, second:k*t

	C_e <- function(b1,b2,b3,c1,c2,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,block))

		tmp1 <- I_1_23(b1,b2,b3,env)

		for(k1 in 1:k.size){
			for(k2 in 1:k.size){
				for(k3 in 1:k.size){
					tmp2 <- c2 * tmp1$first[k1,k2,k3,]

					tmp3 <- I0(tmp2,env)

					first[k1,k2,k3,] <- c1 * tmp3
				}
			}
		}

		second <- matrix(0,k.size,block)

		for(k in 1:k.size){
			tmp4 <- c2 * tmp1$second[k,]

			tmp5 <- I0(tmp4,env)

			second[k,] <- c1 * tmp5
		}

		return(list(first=first,second=second))
	}


	#k*t

	C_f <- function(b1,c1,env){

		k.size <- env$k.size
		block <- env$block

		result <- matrix(0,k.size,block)

		tmp1 <- I_1(b1,env)

		for(k in 1:k.size){
			result[k,] <- c1 * tmp1[k,]
		}

		return(result)
	}


	#first:k*k*t, second:t

	C_g <- function(b1,b2,c1,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,block))

		tmp1 <- I_12(b1,b2,env)

		for(k1 in 1:k.size){
			for(k2 in 1:k.size){
				first[k1,k2,] <- c1 * tmp1$first[k1,k2,]
			}
		}

		second <- c1 * tmp1$second

		return(list(first=first,second=second))
	}


	#first:k*k*k*t, second:k*t

	C_h <- function(b1,b2,b3,c1,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,block))

		tmp1 <- I_123(b1,b2,b3,env)

		for(k1 in 1:k.size){
			for(k2 in 1:k.size){
				for(k3 in 1:k.size){
					first[k1,k2,k3,] <- c1 * tmp1$first[k1,k2,k3,]
				}
			}
		}

		second <- matrix(0,k.size,block)

		for(k in 1:k.size){
			second[k,] <- c1 * tmp1$second[k,]
		}

		return(list(first=first,second=second))
	}


	#first:k*k*k*t, second:k*t

	C_i <- function(b1,b2,b3,c1,c2,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,block))

		tmp1 <- I_1_2_3(b1,b2,b3,env)

		for(k1 in 1:k.size){
			for(k2 in 1:k.size){
				for(k3 in 1:k.size){
					tmp2 <- c2 * tmp1$first[k1,k2,k3,]

					tmp3 <- I0(tmp2,env)

					first[k1,k2,k3,] <- c1 * tmp3
				}
			}
		}

		second <- matrix(0,k.size,block)

		for(k in 1:k.size){
			tmp4 <- c2 * tmp1$second[k,]

			tmp5 <- I0(tmp4,env)

			second[k,] <- c1 * tmp5
		}

		return(list(first=first,second=second))
	}


##p.24

	#first:i*k*k*k*t, second:i*k*k*t, third:i*k*t, fourth:i*t

	C_hat1 <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x1_x2_V0, get_Y_e_e_V, get_e_t, get_U_t, get_U_hat_t, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- matrix(0,d.size,block)

		n1 <- length(get_Y_x1_x2_V0)
		n2 <- sum(get_Y_x1_x2_V0 == 0)

		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n5 <- length(get_e_t)
		n6 <- sum(get_e_t == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0) * (n5 - n6 != 0)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- double(block)

		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  tmp1 <- 3 * get_Y_x1_x2_V0[i3,i4,j,] * get_Y_D[i3,] *
				    get_e_t[i4,]

			  tmp2 <- tmp2 + I0(tmp1,env)
			}
		    }

		    first.tmp[i,] <- first.tmp[i,] + tmpY[i,j,] * tmp2
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,d.size,block)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  for(j1 in 1:d.size){
			    for(j2 in 1:d.size){
				for(j4 in 1:d.size){

				  tmp4 <- double(block)

				  for(i1 in 1:d.size){
				    for(i2 in 1:d.size){
					tmp3 <- get_Y_x1_x2_V0[i1,i2,j4,] * tmpY[i1,j1,] *
						  tmpY[i2,j2,]

					tmp4 <- tmp4 + I0(tmp3,env)
				    }
				  }

				  tmp5 <- 6 * get_Y_x1_x2_V0[i3,i4,j,] * get_Y_D[i3,] *
					    tmpY[i4,j4,] * tmp4

				  tmp6 <- C_c(get_Y_e_V[j1,,],get_Y_e_V[j2,,],tmpY[i,j,],tmp5,env)

				  second.tmp$first[i,,,] <- second.tmp$first[i,,,] +
								    tmp6$first
				  second.tmp$second[i,] <- second.tmp$second[i,] +
								     tmp6$second
				}
			    }
			  }
			}
		    }
		  }
		}

		third.tmp <- array(0,dim=c(d.size,k.size,block))
		fourth.tmp <- array(0,dim=c(d.size,k.size,block))
		fifth.tmp <- list()
		fifth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		fifth.tmp$second <- matrix(0,d.size,block)
		sixth.tmp <- array(0,dim=c(d.size,k.size,block))

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		tmp11 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  for(j1 in 1:d.size){
			    for(j4 in 1:d.size){

				tmp7 <- I0(get_U_t[j1,j4,],env)

				tmp8 <- 6 * get_Y_x1_x2_V0[i3,i4,j,] * get_Y_D[i3,] *
					  tmpY[i4,j4,] * tmp7

				tmp9 <- C_b(get_Y_e_V[j1,,],tmpY[i,j,],tmp8,env)

				third.tmp[i,,] <- third.tmp[i,,] + tmp9

				tmp10 <- 6 * get_Y_x1_x2_V0[i3,i4,j,] * get_Y_D[i3,] *
					   tmpY[i4,j4,]

				for(t in 1:block){
				  tmp11[,t] <- tmp7[t] * get_Y_e_V[j1,,t]
				}

				tmp12 <- C_b(tmp11,tmpY[i,j,],tmp10,env)

				fourth.tmp[i,,] <- fourth.tmp[i,,] - tmp12

				tmp13 <- C_c(get_U_hat_t[j1,j4,,],get_Y_e_V[j1,,],tmpY[i,j,],tmp10,env)

				fifth.tmp$first[i,,,] <- fifth.tmp$first[i,,,] +
								 tmp13$first
				fifth.tmp$second[i,] <- fifth.tmp$second[i,] +
								tmp13$second

				tmp14 <- C_b(get_Y_e_e_V[j4,,],tmpY[i,j,],tmp10/2,env)

				sixth.tmp[i,,] <- sixth.tmp[i,,] + tmp14
			    }
			  }
			}
		    }
		  }
		}

		seventh.tmp <- array(0,dim=c(d.size,k.size,block))

		n <- (n1 - n2 != 0) * (n5 - n6 != 0)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  for(j3 in 1:d.size){

			    tmp15 <- 3 * get_Y_x1_x2_V0[i3,i4,j,] * tmpY[i3,j3,] *
					 get_e_t[i4,]

			    tmp16 <- C_b(get_Y_e_V[j3,,],tmpY[i,j,],tmp15,env)

			    seventh.tmp[i,,] <- seventh.tmp[i,,] + tmp16
			  }
			}
		    }
		  }
		}

		eighth.tmp <- list()
		eighth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,k.size,block))
		eighth.tmp$second <- array(0,dim=c(d.size,k.size,block))

		n <- n1 - n2

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  for(j1 in 1:d.size){
			    for(j2 in 1:d.size){
				for(j4 in 1:d.size){

				  tmp18 <- double(block)

				  for(i1 in 1:d.size){
				    for(i2 in 1:d.size){

					tmp17 <- get_Y_x1_x2_V0[i1,i2,j4,] *
						   tmpY[i1,j1,] * tmpY[i2,j2,]

					tmp18 <- tmp18 + I0(tmp17,env)
				    }
				  }

				  tmp19 <- 6 * get_Y_x1_x2_V0[i3,i4,j,] *
					     tmpY[i3,j3,] * tmpY[i4,j4,] * tmp18

				  tmp20 <- C_e(get_Y_e_V[j3,,],get_Y_e_V[j1,,],get_Y_e_V[j2,,],tmpY[i,j,],tmp19,env)

				  eighth.tmp$first[i,,,,] <- eighth.tmp$first[i,,,,] +
								     tmp20$first
				  eighth.tmp$second[i,,] <- eighth.tmp$second[i,,] +
								    tmp20$second
				}
			    }
			  }
			}
		    }
		  }
		}

		ninth.tmp <- list()
		ninth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		ninth.tmp$second <- matrix(0,d.size,block)
		tenth.tmp <- list()
		tenth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		tenth.tmp$second <- matrix(0,d.size,block)
		eleventh.tmp <- list()
		eleventh.tmp$first <- array(0,dim=c(d.size,k.size,k.size,k.size,block))
		eleventh.tmp$second <- array(0,dim=c(d.size,k.size,block))
		twelfth.tmp <- list()
		twelfth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		twelfth.tmp$second <- matrix(0,d.size,block)

		n <- n1 - n2

		tmp25 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  for(j1 in 1:d.size){
			    for(j4 in 1:d.size){

				tmp21 <- I0(get_U_t[j1,j4,],env)

				tmp22 <- 6 * get_Y_x1_x2_V0[i3,i4,j,] * tmpY[i3,j3,] *
					   tmpY[i4,j4,] * tmp21

				tmp23 <- C_d(get_Y_e_V[j3,,],get_Y_e_V[j1,,],tmpY[i,j,],tmp22,env)

				ninth.tmp$first[i,,,] <- ninth.tmp$first[i,,,] +
								 tmp23$first
				ninth.tmp$second[i,] <- ninth.tmp$second[i,] +
								tmp23$second

				tmp24 <- 6 * get_Y_x1_x2_V0[i3,i4,j,] * tmpY[i3,j3,] *
					   tmpY[i4,j4,]

				for(t in 1:block){
				  tmp25[,t] <- tmp21[t] * get_Y_e_V[j1,,t]
				}

				tmp26 <- C_d(get_Y_e_V[j3,,],tmp25,tmpY[i,j,],tmp24,env)

				tenth.tmp$first[i,,,] <- tenth.tmp$first[i,,,] -
								 tmp26$first
				tenth.tmp$second[i,] <- tenth.tmp$second[i,] -
								tmp26$second

				tmp27 <- C_e(get_Y_e_V[j3,,],get_U_hat_t[j1,j4,,],get_Y_e_V[j1,,],tmpY[i,j,],tmp24,env)

				eleventh.tmp$first[i,,,,] <- eleventh.tmp$first[i,,,,] +
								     tmp27$first
				eleventh.tmp$second[i,,] <- eleventh.tmp$second[i,,] +
								    tmp27$second

				tmp28 <- C_d(get_Y_e_V[j3,,],get_Y_e_e_V[j4,,],tmpY[i,j,],tmp24/2,env)

				twelfth.tmp$first[i,,,] <- twelfth.tmp$first[i,,,] +
								    tmp28$first
				twelfth.tmp$second[i,] <- twelfth.tmp$second[i,] +
								   tmp28$second
			    }
			  }
			}
		    }
		  }
		}

		first <- eighth.tmp$first + eleventh.tmp$first
		second <- second.tmp$first + fifth.tmp$first + ninth.tmp$first +
			    tenth.tmp$first + twelfth.tmp$first
		third <- third.tmp + fourth.tmp + sixth.tmp + seventh.tmp +
			   eighth.tmp$second + eleventh.tmp$second
		fourth <- first.tmp + second.tmp$second + fifth.tmp$second +
			    ninth.tmp$second + tenth.tmp$second +
			    twelfth.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:i*k*k*t, second:i*k*t, third:i*t

	C_hat2 <- function(tmpY, get_Y_e_V, get_Y_x1_x2_V0, get_Y_x_e_V0, get_Y_e_e_V, get_e_t, get_U_t, get_U_hat_t, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- matrix(0,d.size,block)

		n1 <- length(get_Y_x1_x2_V0)
		n2 <- sum(get_Y_x1_x2_V0 == 0)

		n3 <- length(get_e_t)
		n4 <- sum(get_e_t == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- double(block)

		    for(i3 in 1:d.size){

			tmp1 <- 3 * get_Y_x_e_V0[i3,j,] * get_e_t[i3,]

			tmp2 <- tmp2 + I0(tmp1,env)
		    }

		    first.tmp[i,] <- first.tmp[i,] + tmpY[i,j,] * tmp2
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,d.size,block)

		n5 <- length(get_Y_x_e_V0)
		n6 <- sum(get_Y_x_e_V0 == 0)

		n <- (n1 - n2 != 0) * (n5 - n6 != 0)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(j1 in 1:d.size){
			  for(j2 in 1:d.size){
			    for(j3 in 1:d.size){

				tmp4 <- double(block)

				for(i1 in 1:d.size){
				  for(i2 in 1:d.size){

				    tmp3 <- get_Y_x1_x2_V0[i1,i2,j3,] *
						tmpY[i1,j1,] * tmpY[i2,j2,]

				    tmp4 <- tmp4 + I0(tmp3,env)
				  }
				}

				tmp5 <- 6 * get_Y_x_e_V0[i3,j,] * 
					  tmpY[i3,j3,] * tmp4

				tmp6 <- C_c(get_Y_e_V[j1,,],get_Y_e_V[j2,,],tmpY[i,j,],tmp5,env)

				second.tmp$first[i,,,] <- second.tmp$first[i,,,] +
								  tmp6$first
				second.tmp$second[i,] <- second.tmp$second[i,] +
								 tmp6$second
			    }
			  }
			}
		    }
		  }
		}

		third.tmp <- array(0,dim=c(d.size,k.size,block))
		fourth.tmp <- array(0,dim=c(d.size,k.size,block))
		fifth.tmp <- list()
		fifth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		fifth.tmp$second <- matrix(0,d.size,block)
		sixth.tmp <- array(0,dim=c(d.size,k.size,block))

		n <- n1 - n2

		tmp11 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(j1 in 1:d.size){
			  for(j3 in 1:d.size){

			    tmp7 <- I0(get_U_t[j1,j3,],env)

			    tmp8 <- 6 * get_Y_x_e_V0[i3,j,] *
					tmpY[i3,j3,] * tmp7

			    tmp9 <- C_b(get_Y_e_V[j1,,],tmpY[i,j,],tmp8,env)

			    third.tmp[i,,] <- third.tmp[i,,] + tmp9

			    tmp10 <- 6 * get_Y_x_e_V0[i3,j,] * tmpY[i3,j3,]

			    for(r in 1:r.size){
				tmp11[r,] <- tmp7 * get_Y_e_V[j1,r,]
			    }

			    tmp12 <- C_b(tmp11,tmpY[i,j,],tmp10,env)

			    fourth.tmp[i,,] <- fourth.tmp[i,,] - tmp12

			    tmp13 <- C_c(get_U_hat_t[j1,j3,,],get_Y_e_V[j1,,],tmpY[i,j,],tmp10,env)

			    fifth.tmp$first[i,,,] <- fifth.tmp$first[i,,,] +
							     tmp13$first
			    fifth.tmp$second[i,] <- fifth.tmp$second[i,] +
							    tmp13$second

			    tmp14 <- C_b(get_Y_e_e_V[j3,,],tmpY[i,j,],tmp10/2,env)

			    sixth.tmp[i,,] <- sixth.tmp[i,,] + tmp14
			  }
			}
		    }
		  }
		}

		first <- second.tmp$first + fifth.tmp$first
		second <- third.tmp + fourth.tmp + sixth.tmp
		third <- first.tmp + second.tmp$second + fifth.tmp$second

		return(list(first=first,second=second,third=third))
	}


	#first:i*k*k*k*t, second:i*k*k*t, third:i*k*t, fourth:i*t

	C_hat3 <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x1_x2_x3_V0, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- matrix(0,d.size,block)

		n1 <- length(get_Y_x1_x2_x3_V0)
		n2 <- sum(get_Y_x1_x2_x3_V0 == 0)

		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- double(block)

		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(i3 in 1:d.size){

			    tmp1 <- get_Y_x1_x2_x3_V0[i1,i2,i3,j,] *
					get_Y_D[i1,] * get_Y_D[i2,] *
					get_Y_D[i3,]

			    tmp2 <- tmp2 + I0(tmp1,env)
			  }
			}
		    }

		    first.tmp[i,] <- first.tmp[i,] + tmpY[i,j,] * tmp2
		  }
		}

		second.tmp <- array(0,dim=c(d.size,k.size,block))

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j3 in 1:d.size){

			tmp4 <- double(block)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){
			    for(i3 in 1:d.size){

				tmp3 <- get_Y_x1_x2_x3_V0[i1,i2,i3,j3,] *
					  get_Y_D[i1,] * get_Y_D[i2,] *
					  tmpY[i3,j3,]

				tmp4 <- tmp4 + tmp3
			    }
			  }
			}

			tmp5 <- C_b(get_Y_e_V[j3,,],tmpY[i,j,],tmp4,env)

			second.tmp[i,,] <- second.tmp[i,,] + 3 * tmp5
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		third.tmp$second <- matrix(0,d.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j2 in 1:d.size){
			for(j3 in 1:d.size){

			  tmp7 <- double(block)

			  for(i1 in 1:d.size){
			    for(i2 in 1:d.size){
				for(i3 in 1:d.size){

				  tmp6 <- get_Y_x1_x2_x3_V0[i1,i2,i3,j,] *
					    get_Y_D[i1,] * tmpY[i2,j2,] *
					    tmpY[i3,j3,]

				  tmp7 <- tmp7 + tmp6
				}
			    }
			  }

			  tmp8 <- C_d(get_Y_e_V[j2,,],get_Y_e_V[j3,,],tmpY[i,j,],tmp7,env)

			  third.tmp$first[i,,,] <- third.tmp$first[i,,,] +
							   3 * tmp8$first
			  third.tmp$second[i,] <- third.tmp$second[i,] +
							  3 * tmp8$second
			}
		    }
		  }
		}

		fourth.tmp <- list()
		fourth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,k.size,block))
		fourth.tmp$second <- array(0,dim=c(d.size,k.size,block))

		n <- n1 - n2

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){
			for(j2 in 1:d.size){
			  for(j3 in 1:d.size){

			    tmp10 <- double(block)

			    for(i1 in 1:d.size){
				for(i2 in 1:d.size){
				  for(i3 in 1:d.size){

				    tmp9 <- get_Y_x1_x2_x3_V0[i1,i2,i3,j,] *
						tmpY[i1,j1,] * tmpY[i2,j2,] *
						tmpY[i3,j3,]

				    tmp10 <- tmp10 + tmp9
				  }
				}
			    }

			    tmp11 <- C_i(get_Y_e_V[j1,,],get_Y_e_V[j2,,],get_Y_e_V[j3,,],tmpY[i,j,],tmp10,env)

			    fourth.tmp$first[i,,,,] <- fourth.tmp$first[i,,,,] +
							       tmp11$first
			    fourth.tmp$second[i,,] <- fourth.tmp$second[i,,] +
							      tmp11$second
			  }
			}
		    }
		  }
		}

		first <- fourth.tmp$first
		second <- third.tmp$first
		third <- second.tmp + fourth.tmp$second
		fourth <- first.tmp + third.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:i*k*k*t, second:i*k*t, third:i*t

	C_hat4 <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x1_x2_e_V0, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		n1 <- length(get_Y_x1_x2_e_V0)
		n2 <- sum(get_Y_x1_x2_e_V0 == 0)

		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		first.tmp <- matrix(0,d.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- double(block)

		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){

			  tmp1 <- 3 * get_Y_x1_x2_e_V0[i1,i2,j,] * get_Y_D[i1,] *
				    get_Y_D[i2,]

			  tmp2 <- tmp2 + I0(tmp1,env)
			}
		    }

		    first.tmp[i,] <- first.tmp[i,] + tmpY[i,j,] * tmp2
		  }
		}

		second.tmp <- array(0,dim=c(d.size,k.size,block))

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j2 in 1:d.size){

			tmp4 <- double(block)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){

			    tmp3 <- 6 * get_Y_x1_x2_e_V0[i1,i2,j,] *
					get_Y_D[i1,] * tmpY[i2,j2,]

			    tmp4 <- tmp4 + tmp3
			  }
			}

			tmp5 <- C_b(get_Y_e_V[j2,,],tmpY[i,j,],tmp4,env)

			second.tmp[i,,] <- second.tmp[i,,] + tmp5
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		third.tmp$second <- matrix(0,d.size,block)

		n <- n1 - n2

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){
			for(j2 in 1:d.size){

			  tmp7 <- double(block)

			  for(i1 in 1:d.size){
			    for(i2 in 1:d.size){

				tmp6 <- 3 * get_Y_x1_x2_e_V0[i1,i2,j,] *
					  tmpY[i1,j1,] * tmpY[i2,j2,]

				tmp7 <- tmp7 + tmp6
			    }
			  }

			  tmp8 <- C_d(get_Y_e_V[j1,,],get_Y_e_V[j2,,],tmpY[i,j,],tmp7,env)

			  third.tmp$first[i,,,] <- third.tmp$first[i,,,] +
							   tmp8$first
			  third.tmp$second[i,] <- third.tmp$second[i,] +
							  tmp8$second
			}
		    }
		  }
		}

		first <- third.tmp$first
		second <- second.tmp
		third <- first.tmp + third.tmp$second

		return(list(first=first,second=second,third=third))
	}


	#first:i*k*t, second:i*t

	C_hat5 <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x_e_e_V0, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		n1 <- length(get_Y_x_e_e_V0)
		n2 <- sum(get_Y_x_e_e_V0 == 0)

		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		first.tmp <- matrix(0,d.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- double(block)

		    for(i1 in 1:d.size){

			tmp1 <- 3 * get_Y_x_e_e_V0[i1,j,] * get_Y_D[i1,]

			tmp2 <- tmp2 + I0(tmp1,env)
		    }

		    first.tmp[i,] <- first.tmp[i,] + tmpY[i,j,] * tmp2
		  }
		}

		second.tmp <- array(0,dim=c(d.size,k.size,block))

		n <- n1 - n2

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){

			tmp4 <- double(block)

			for(i1 in 1:d.size){

			  tmp3 <- 3 * get_Y_x_e_e_V0[i1,j,] * tmpY[i1,j1,]

			  tmp4 <- tmp4 + tmp3
			}

			tmp5 <- C_b(get_Y_e_V[j1,,],tmpY[i,j,],tmp4,env)

			second.tmp[i,,] <- second.tmp[i,,] + tmp5
		    }
		  }
		}

		first <- second.tmp
		second <- first.tmp

		return(list(first=first,second=second))
	}


	#i*t

	C_hat6 <- function(tmpY, get_Y_e_e_e_V0, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		result <- matrix(0,d.size,block)

		n1 <- length(get_Y_e_e_e_V0)
		n2 <- sum(get_Y_e_e_e_V0 == 0)

		n <- n1 - n2

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    tmp1 <- I0(get_Y_e_e_e_V0[j,],env)

		    result[i,] <- result[i,] + tmpY[i,j,] * tmp1
		  }
		}

		return(result)
	}


	#first:i*k*k*k*t, second:i*k*k*t, third:i*k*t, fourth:i*t

	C_hat7 <- function(tmpY, get_Y_e_V, get_Y_x1_x2_V0, get_Y_x_e_V, get_Y_e_e_V, get_e_t, get_U_t, get_U_hat_t, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(d.size,k.size,block))

		n1 <- length(get_Y_x_e_V)
		n2 <- sum(get_Y_x_e_V == 0)

		n3 <- length(get_e_t)
		n4 <- sum(get_e_t == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		tmp1 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- matrix(0,r.size,block)

		    for(i3 in 1:d.size){

			for(r in 1:r.size){
			  tmp1[r,] <- 3 * get_Y_x_e_V[i3,j,r,] * get_e_t[i3,]
			}

			tmp2 <- tmp2 + tmp1
		    }

		    tmp3 <- C_f(tmp2,tmpY[i,j,],env)

		    first.tmp[i,,] <- first.tmp[i,,] + tmp3
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(d.size,k.size,k.size,k.size,block))
		second.tmp$second <- array(0,dim=c(d.size,k.size,block))

		n5 <- length(get_Y_x1_x2_V0)
		n6 <- sum(get_Y_x1_x2_V0 == 0)

		n <- (n1 - n2 != 0) * (n5 - n6 != 0)

		tmp7 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(j1 in 1:d.size){
			  for(j2 in 1:d.size){
			    for(j3 in 1:d.size){

				tmp5 <- double(block)

				for(i1 in 1:d.size){
				  for(i2 in 1:d.size){

				    tmp4 <- get_Y_x1_x2_V0[i1,i2,j3,] * tmpY[i1,j1,] *
						tmpY[i2,j2,]

				    tmp5 <- tmp5 + tmp4
				  }
				}

				tmp6 <- I0(tmp5,env)

				for(r in 1:r.size){
				  tmp7[r,] <- 6 * get_Y_x_e_V[i3,j,r,] * 
						  tmpY[i3,j3,] * tmp6
				}

				tmp8 <- C_h(tmp7,get_Y_e_V[j1,,],get_Y_e_V[j2,,],tmpY[i,j,],env)

				second.tmp$first[i,,,,] <- second.tmp$first[i,,,,] +
								   tmp8$first
				second.tmp$second[i,,] <- second.tmp$second[i,,] +
								  tmp8$second
			    }
			  }
			}
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		third.tmp$second <- matrix(0,d.size,block)
		fourth.tmp <- list()
		fourth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		fourth.tmp$second <- matrix(0,d.size,block)
		fifth.tmp <- list()
		fifth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,k.size,block))
		fifth.tmp$second <- array(0,dim=c(d.size,k.size,block))
		sixth.tmp <- list()
		sixth.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		sixth.tmp$second <- matrix(0,d.size,block)

		tmp10 <- matrix(0,r.size,block)

		tmp12 <- matrix(0,r.size,block)
		tmp13 <- matrix(0,r.size,block)

		for(i in 1:d.size){
		  for(j in 1:d.size){
		    for(i3 in 1:d.size){
			for(j1 in 1:d.size){
			  for(j3 in 1:d.size){

			    tmp9 <- I0(get_U_t[j1,j3,],env)

			    for(r in 1:r.size){
				tmp10[r,] <- 6 * get_Y_x_e_V[i3,j,r,] *
					      tmpY[i3,j3,] * tmp9
			    }

			    tmp11 <- C_g(tmp9,get_Y_e_V[j1,,],tmpY[i,j,],env)

			    third.tmp$first[i,,,] <- third.tmp$first[i,,,] +
							     tmp11$first
			    third.tmp$second[i,] <- third.tmp$second[i,] +
							    tmp11$second

			    for(r in 1:r.size){
				tmp12[r,] <- 6 * get_Y_x_e_V[i3,j,r,] * tmpY[i3,j3,]
				tmp13[r,] <- tmp9 * get_Y_e_V[j1,r,]
			    }

			    tmp14 <- C_g(tmp12,tmp13,tmpY[i,j,],env)

			    fourth.tmp$first[i,,,] <- fourth.tmp$first[i,,,] -
							      tmp14$first
			    fourth.tmp$second[i,] <- fourth.tmp$second[i,] -
							     tmp14$second

			    tmp15 <- C_h(tmp12,get_U_hat_t[j1,j3,,],get_Y_e_V[j1,,],tmpY[i,j,],env)

			    fifth.tmp$first[i,,,,] <- fifth.tmp$first[i,,,,] +
							      tmp15$first
			    fifth.tmp$second[i,,] <- fifth.tmp$second[i,,] +
							     tmp15$second

			    tmp16 <- C_g(tmp12/2,get_Y_e_e_V[j3,,],tmpY[i,j,],env)

			    sixth.tmp$first[i,,,] <- sixth.tmp$first[i,,,] +
							     tmp16$first
			    sixth.tmp$second[i,] <- sixth.tmp$second[i,] +
							    tmp16$second
			  }
			}
		    }
		  }
		}

		first <- second.tmp$first + fifth.tmp$first
		second <- third.tmp$first + fourth.tmp$first + sixth.tmp$first
		third <- first.tmp + second.tmp$second + fifth.tmp$second
		fourth <- third.tmp$second + fourth.tmp$second + sixth.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:i*k*k*k*t, second:i*k*k*t, third:i*k*t, fourth:i*t

	C_hat8 <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x1_x2_e_V, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(d.size,k.size,block))

		n1 <- length(get_Y_x1_x2_e_V)
		n2 <- sum(get_Y_x1_x2_e_V == 0)

		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		tmp1 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- matrix(0,r.size,block)

		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){

			  for(r in 1:r.size){
			    tmp1[r,] <- 3 * get_Y_x1_x2_e_V[i1,i2,j,r,] *
					    get_Y_D[i1,] * get_Y_D[i2,]
			  }

			  tmp2 <- tmp2 + tmp1
			}
		    }

		    tmp3 <- C_f(tmp2,tmpY[i,j,],env)

		    first.tmp[i,,] <- first.tmp[i,,] + tmp3
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,d.size,block)

		tmp4 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j2 in 1:d.size){

			tmp5 <- matrix(0,r.size,block)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){

			    for(r in 1:r.size){
				tmp4[r,] <- 6 * get_Y_x1_x2_e_V[i1,i2,j,r,] *
						get_Y_D[i1,] * tmpY[i2,j2,]
			    }

			    tmp5 <- tmp5 + tmp4
			  }
			}

			tmp6 <- C_g(tmp5,get_Y_e_V[j2,,],tmpY[i,j,],env)

			second.tmp$first[i,,,] <- second.tmp$first[i,,,] +
							  tmp6$first
			second.tmp$second[i,] <- second.tmp$second[i,] +
							 tmp6$second
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(d.size,k.size,k.size,k.size,block))
		third.tmp$second <- array(0,dim=c(d.size,k.size,block))

		n <- n1 - n2

		tmp7 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){
			for(j2 in 1:d.size){

			  tmp8 <- matrix(0,r.size,block)

			  for(i1 in 1:d.size){
			    for(i2 in 1:d.size){

				for(r in 1:r.size){
				  tmp7[r,] <- 6 * get_Y_x1_x2_e_V[i1,i2,j,r,] *
						  tmpY[i1,j1,] * tmpY[i2,j2,]
				}

				tmp8 <- tmp8 + tmp7
			    }
			  }

			  tmp9 <- C_h(tmp8,get_Y_e_V[j1,,],get_Y_e_V[j2,,],tmpY[i,j,],env)

			  third.tmp$first[i,,,,] <- third.tmp$first[i,,,,] +
							    tmp9$first
			  third.tmp$second[i,,] <- third.tmp$second[i,,] +
							   tmp9$second
			}
		    }
		  }
		}

		fourth.tmp <- array(0,dim=c(d.size,k.size,block))

		tmp11 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){
			for(j2 in 1:d.size){

			  tmp10 <- b1_b2(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)
			  tmp12 <- matrix(0,r.size,block)

			  for(i1 in 1:d.size){
			    for(i2 in 1:d.size){

				for(r in 1:r.size){
				  tmp11[r,] <- 3 * get_Y_x1_x2_e_V[i1,i2,j,r,] *
						   tmpY[i1,j1,] * tmpY[i2,j2,] *
						   tmp10
				}

				tmp12 <- tmp12 + tmp11
			    }
			  }

			  tmp13 <- C_f(tmp12,tmpY[i,j,],env)

			  fourth.tmp[i,,] <- fourth.tmp[i,,] + tmp13
			}
		    }
		  }
		}

		first <- third.tmp$first
		second <- second.tmp$first
		third <- first.tmp + third.tmp$second + fourth.tmp
		fourth <- second.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:i*k*k*t, second:i*k*t, third:i*t

	C_hat9 <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x_e_e_V, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		n1 <- length(get_Y_x_e_e_V)
		n2 <- sum(get_Y_x_e_e_V == 0)

		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		first.tmp <- array(0,dim=c(d.size,k.size,block))

		tmp1 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp2 <- matrix(0,r.size,block)

		    for(i1 in 1:d.size){

			for(r in 1:r.size){
			  tmp1[r,] <- get_Y_x_e_e_V[i1,j,r,] * get_Y_D[i1,]
			}

			tmp2 <- tmp2 + tmp1
		    }

		    tmp3 <- C_f(tmp2,tmpY[i,j,],env)

		    first.tmp[i,,] <- first.tmp[i,,] + tmp3
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(d.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,d.size,block)

		n <- n1 - n2

		tmp4 <- matrix(0,r.size,block)

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){

			tmp5 <- matrix(0,r.size,block)

			for(i1 in 1:d.size){

			  for(r in 1:r.size){
			    tmp4[r,] <- get_Y_x_e_e_V[i1,j,r,] * tmpY[i1,j1,]
			  }

			  tmp5 <- tmp5 + tmp4
			}

			tmp6 <- C_g(tmp5,get_Y_e_V[j1,,],tmpY[i,j,],env)

			second.tmp$first[i,,,] <- second.tmp$first[i,,,] +
							  tmp6$first
			second.tmp$second[i,] <- second.tmp$second[i,] +
							 tmp6$second
		    }
		  }
		}

		first <- 3 * second.tmp$first
		second <- 3 * first.tmp
		third <- 3 * second.tmp$second

		return(list(first=first,second=second,third=third))
	}


	#first:i*k*t

	C_hat10 <- function(tmpY, get_Y_e_e_e_V, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		result <- array(0,dim=c(d.size,k.size,block))

		n1 <- length(get_Y_e_e_e_V)
		n2 <- sum(get_Y_e_e_e_V == 0)

		n <- n1 - n2

		for(i in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    result[i,,] <- result[i,,] + C_f(get_Y_e_e_e_V[j,,],tmpY[i,j,],env)
		  }
		}

		return(result)
	}


	#first:i*k*k*k*t, second:i*k*k*t, third:i*k*t, fourth:i*t

	C0_t <- function(get_C_hat1, get_C_hat2, get_C_hat3, get_C_hat4,
			     get_C_hat5, get_C_hat6, get_C_hat7, get_C_hat8,
			     get_C_hat9, get_C_hat10){

		result1 <- get_C_hat1
		result2 <- get_C_hat2
		result3 <- get_C_hat3
		result4 <- get_C_hat4
		result5 <- get_C_hat5
		result6 <- get_C_hat6
		result7 <- get_C_hat7
		result8 <- get_C_hat8
		result9 <- get_C_hat9
		result10 <- get_C_hat10

		first <- result1$first + result3$first + result7$first +
			   result8$first

		second <- result1$second + result2$first + result3$second +
			    result4$first + result7$second + result8$second +
			    result9$first

		third <- result1$third + result2$second + result3$third +
			   result4$second + result5$first + result7$third +
			   result8$third + result9$second + result10

		fourth <- result1$fourth + result2$third + result3$fourth +
			    result4$third + result5$second + result6 +
			    result7$fourth + result8$fourth + result9$third

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


##p.32(II)-1

	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_1_1 <- function(get_x_f0,get_C0_t,env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		second <- array(0,dim=c(k.size,k.size,k.size,block))
		third <- array(0,dim=c(k.size,k.size,block))
		fourth <- matrix(0,k.size,block)

		n1 <- length(get_x_f0)
		n2 <- sum(get_x_f0 == 0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i in 1:d.size){
		    for(k1 in 1:k.size){
			for(k2 in 1: k.size){
			  for(k3 in 1:k.size){

			    tmp1 <- get_x_f0[l,i,] *
					get_C0_t$first[i,k1,k2,k3,]

			    first[l,k1,k2,k3,] <- first[l,k1,k2,k3,] + I0(tmp1,env)/6
			  }

			  tmp2 <- get_x_f0[l,i,] *
				    get_C0_t$second[i,k1,k2,]

			  second[l,k1,k2,] <- second[l,k1,k2,] + I0(tmp2,env)/6
			}

			tmp3 <- get_x_f0[l,i,] *
				  get_C0_t$third[i,k1,]

			third[l,k1,] <- third[l,k1,] + I0(tmp3,env)/6
		    }

		    tmp4 <- get_x_f0[l,i,] *
				get_C0_t$fourth[i,]

		    fourth[l,] <- fourth[l,] + I0(tmp4,env)/6
		  }
		}

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_1_2 <- function(get_x_F,get_C0_t,env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		second <- array(0,dim=c(k.size,k.size,k.size,block))
		third <- array(0,dim=c(k.size,k.size,block))
		fourth <- matrix(0,k.size,block)

		n1 <- length(get_x_F)
		n2 <- sum(get_x_F == 0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i in 1:d.size){
		    for(k1 in 1:k.size){
			for(k2 in 1: k.size){
			  for(k3 in 1:k.size){

			    tmp1 <- get_x_F[l,i,] *
					get_C0_t$first[i,k1,k2,k3,]

			    first[l,k1,k2,k3,] <- first[l,k1,k2,k3,] + tmp1/6
			  }

			  tmp2 <- get_x_F[l,i,] *
				    get_C0_t$second[i,k1,k2,]

			  second[l,k1,k2,] <- second[l,k1,k2,] + tmp2/6
			}

			tmp3 <- get_x_F[l,i,] *
				  get_C0_t$third[i,k1,]

			third[l,k1,] <- third[l,k1,] + tmp3/6
		    }

		    tmp4 <- get_x_F[l,i,] *
				get_C0_t$fourth[i,]

		    fourth[l,] <- fourth[l,] + tmp4/6
		  }
		}

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


##p.33(II)-2

	#first:i3*i4*k*k*k*t, second:i3*i4*k*k*t, third:i3*i4*k*t, fourth:i3*i4*t

	D_E <- function(tmpY, get_Y_e_V, get_Y_D, get_Y_x1_x2_V0, get_Y_e_e_V, get_e_t, get_U_t, get_U_hat_t, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(d.size,d.size,block))

		n1 <- length(get_Y_D)
		n2 <- sum(get_Y_D == 0)

		n <- n1 - n2

		for(i3 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(i4 in 1:d.size){
		    first.tmp[i3,i4,] <- get_Y_D[i3,] * get_e_t[i4,]
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,block))
		second.tmp$second <- array(0,dim=c(d.size,d.size,block))

		n3 <- length(get_Y_x1_x2_V0)
		n4 <- sum(get_Y_x1_x2_V0 == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		I0_V0_Y_Y <- array(0,dim=c(d.size,d.size,d.size,block))	#j1,j2,j4,t

		for(j1 in 1:d.size){

		  if(n3 == n4){
		    break
		  }

		  for(j2 in 1:d.size){
		    for(j4 in 1:d.size){

			tmp2 <- double(block)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){
			    tmp1 <- get_Y_x1_x2_V0[i1,i2,j4,] * tmpY[i1,j1,] *
					tmpY[i2,j2,]

			    tmp2 <- tmp2 + tmp1
			  }
			}

			I0_V0_Y_Y[j1,j2,j4,] <- I0_V0_Y_Y[j1,j2,j4,] + I0(tmp2,env)
		    }
		  }
		}

		for(j1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j2 in 1:d.size){

		    tmp3 <- I_12(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)

		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){

			  tmp4 <- double(block)

			  for(j4 in 1:d.size){
			    tmp4 <- tmp4 + 2 * get_Y_D[i3,] * tmpY[i4,j4,] *
					I0_V0_Y_Y[j1,j2,j4,]
			  }

			  for(t in 1:block){
			    second.tmp$first[i3,i4,,,t] <- second.tmp$first[i3,i4,,,t] +
								     tmp4[t] * tmp3$first[,,t]

			    second.tmp$second[i3,i4,t] <- second.tmp$second[i3,i4,t] +
								    tmp4[t] * tmp3$second[t]
			  }
			}
		    }
		  }
		}

		third.tmp  <- array(0,dim=c(d.size,d.size,k.size,block))

		n5 <- length(get_U_t)
		n6 <- sum(get_U_t == 0)

		n <- (n1 - n2 != 0) * (n5 - n6 != 0)

		I0_U <- array(0,dim=c(d.size,d.size,block))	#j1,j4,t

		for(j1 in 1:d.size){

		  if(n5 == n6){
		    break
		  }

		  for(j4 in 1:d.size){
		    I0_U[j1,j4,] <- I0(get_U_t[j1,j4,],env)
		  }
		}

		for(j1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  tmp5 <- I_1(get_Y_e_V[j1,,],env)

		  for(i3 in 1:d.size){
		    for(i4 in 1:d.size){

			tmp6 <- double(block)

			for(j4 in 1:d.size){
			  tmp6 <- tmp6 + 2 * get_Y_D[i3,] * tmpY[i4,j4,] *
				    I0_U[j1,j4,]
			}

			for(k in 1:k.size){
			  third.tmp[i3,i4,k,] <- third.tmp[i3,i4,k,] +
							 tmp6 * tmp5[k,]
			}
		    }
		  }
		}

		fourth.tmp  <- array(0,dim=c(d.size,d.size,k.size,block))

		tmp7 <- matrix(0,r.size,block)

		for(j4 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  tmp8 <- matrix(0,r.size,block)

		  for(j1 in 1:d.size){

		    for(r in 1:r.size){
			tmp7[r,] <- I0_U[j1,j4,] * get_Y_e_V[j1,r,]
		    }

		    tmp8 <- tmp8 + tmp7
		  }

		  tmp9 <- I_1(tmp8,env)


		  for(i3 in 1:d.size){
		    for(i4 in 1:d.size){
			for(k in 1:k.size){
			  fourth.tmp[i3,i4,k,] <- fourth.tmp[i3,i4,k,] - 2 *
							  get_Y_D[i3,] * tmpY[i4,j4,] *
							  tmp9[k,]
			}
		    }
		  }
		}

		fifth.tmp <- list()
		fifth.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,block))
		fifth.tmp$second <- array(0,dim=c(d.size,d.size,block))

		n <- n1 - n2

		for(j1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j4 in 1:d.size){

		    tmp10 <- I_12(get_U_hat_t[j1,j4,,],get_Y_e_V[j1,,],env)

		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  for(t in 1:block){
			    fifth.tmp$first[i3,i4,,,t] <- fifth.tmp$first[i3,i4,,,t] + 2 *
							 	    get_Y_D[i3,t] * tmpY[i4,j4,t] *
								    tmp10$first[,,t]

			    fifth.tmp$second[i3,i4,t] <- fifth.tmp$second[i3,i4,t] + 2 *
								   get_Y_D[i3,t] * tmpY[i4,j4,t] *
								   tmp10$second[t]
			  }
			}
		    }
		  }
		}

		sixth.tmp  <- array(0,dim=c(d.size,d.size,k.size,block))

		n7 <- length(get_Y_e_e_V)
		n8 <- sum(get_Y_e_e_V == 0)

		n <- (n1 - n2 != 0) * (n7 - n8 != 0)

		for(j4 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  tmp11 <- I_1(get_Y_e_e_V[j4,,],env)

		  for(i3 in 1:d.size){
		    for(i4 in 1:d.size){
			for(k in 1:k.size){
			  sixth.tmp[i3,i4,k,] <- sixth.tmp[i3,i4,k,] +
							 get_Y_D[i3,] * tmpY[i4,j4,] *
							 tmp11[k,]
			}
		    }
		  }
		}

		seventh.tmp  <- array(0,dim=c(d.size,d.size,k.size,block))

		n9 <- length(get_e_t)
		n10 <- sum(get_e_t == 0)

		n <- n9 - n10

		for(j3 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  tmp12 <- I_1(get_Y_e_V[j3,,],env)

		  for(i3 in 1:d.size){
		    for(i4 in 1:d.size){
			for(k in 1:k.size){
			  seventh.tmp[i3,i4,k,] <- seventh.tmp[i3,i4,k,] +
							   tmpY[i3,j3,] * get_e_t[i4,] *
							   tmp12[k,]
			}
		    }
		  }
		}

		eighth.tmp <- list()
		eighth.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,k.size,block))
		eighth.tmp$second <- array(0,dim=c(d.size,d.size,k.size,block))

		n <- n3 - n4

		for(j1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j2 in 1:d.size){
		    for(j3 in 1:d.size){

			tmp13 <- I_1_23(get_Y_e_V[j3,,],get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)

			for(i3 in 1:d.size){
			  for(i4 in 1:d.size){

			    tmp14 <- double(block)

			    for(j4 in 1:d.size){

				tmp14 <- tmp14 + 2 * tmpY[i3,j3,] * tmpY[i4,j4,] *
					   I0_V0_Y_Y[j1,j2,j4,]
			    }

			    for(t in 1:block){
				eighth.tmp$first[i3,i4,,,,t] <- eighth.tmp$first[i3,i4,,,,t] +
									  tmp14[t] * tmp13$first[,,,t]

				eighth.tmp$second[i3,i4,,t] <- eighth.tmp$second[i3,i4,,t] +
									 tmp14[t] * tmp13$second[,t]
			    }
			  }
			}
		    }
		  }
		}

		ninth.tmp <- list()
		ninth.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,block))
		ninth.tmp$second <- array(0,dim=c(d.size,d.size,block))

		n <- n5 - n6

		for(j1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j3 in 1:d.size){

		    tmp15 <- I_1_2(get_Y_e_V[j3,,],get_Y_e_V[j1,,],env)

		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){

			  tmp16 <- double(block)

			  for(j4 in 1:d.size){

			    tmp16 <- tmp16 + 2 * tmpY[i3,j3,] * tmpY[i4,j4,] *
				 	 I0_U[j1,j4,]
			  }

			  for(t in 1:block){
			    ninth.tmp$first[i3,i4,,,t] <- ninth.tmp$first[i3,i4,,,t] +
								    tmp16[t] * tmp15$first[,,t]

			    ninth.tmp$second[i3,i4,t] <- ninth.tmp$second[i3,i4,t] +
								   tmp16[t] * tmp15$second[t]
			  }
			}
		    }
		  }
		}

		tenth.tmp <- list()
		tenth.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,block))
		tenth.tmp$second <- array(0,dim=c(d.size,d.size,block))

		for(j3 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j4 in 1:d.size){

		    tmp17 <- matrix(0,r.size,block)

		    for(j1 in 1:d.size){
			for(r in 1:r.size){
			  tmp17[r,] <- tmp17[r,] + I0_U[j1,j4,] *
					   get_Y_e_V[j1,r,]						 
			}
		    }

		    tmp18 <- I_1_2(get_Y_e_V[j3,,],tmp17,env)

		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){
			  tmp19 <- 2 * tmpY[i3,j3,] * tmpY[i4,j4,]

			  for(t in 1:block){
			    tenth.tmp$first[i3,i4,,,t] <- tenth.tmp$first[i3,i4,,,t] -
								    tmp19[t] * tmp18$first[,,t]

			    tenth.tmp$second[i3,i4,t] <- tenth.tmp$second[i3,i4,t] -
								   tmp19[t] * tmp18$second[t]
			  }
			}
		    }
		  }
		}

		eleventh.tmp <- list()
		eleventh.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,k.size,block))
		eleventh.tmp$second <- array(0,dim=c(d.size,d.size,k.size,block))

		for(j1 in 1:d.size){
		  for(j3 in 1:d.size){
		    for(j4 in 1:d.size){

			tmp20 <- I_1_23(get_Y_e_V[j3,,],get_U_hat_t[j1,j4,,],get_Y_e_V[j1,,],env)

			for(i3 in 1:d.size){
			  for(i4 in 1:d.size){

			    tmp21 <- 2 * tmpY[i3,j3,] * tmpY[i4,j4,]

			    for(t in 1:block){
			 	eleventh.tmp$first[i3,i4,,,,t] <- eleventh.tmp$first[i3,i4,,,,t] +
								 	    tmp21[t] * tmp20$first[,,,t]

				eleventh.tmp$second[i3,i4,,t] <- eleventh.tmp$second[i3,i4,,t] +
									   tmp21[t] * tmp20$second[,t]
			    }
			  }
			}
		    }
		  }
		}

		twelfth.tmp <- list()
		twelfth.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,block))
		twelfth.tmp$second <- array(0,dim=c(d.size,d.size,block))

		n <- n7 - n8

		for(j3 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j4 in 1:d.size){

		    tmp22 <- I_1_2(get_Y_e_V[j3,,],get_Y_e_e_V[j4,,],env)

		    for(i3 in 1:d.size){
			for(i4 in 1:d.size){

			  tmp23 <- tmpY[i3,j3,] * tmpY[i4,j4,]

			  for(t in 1:block){
			    twelfth.tmp$first[i3,i4,,,t] <- twelfth.tmp$first[i3,i4,,,t] +
									tmp23[t] * tmp22$first[,,t]

			    twelfth.tmp$second[i3,i4,t] <- twelfth.tmp$second[i3,i4,t] +
								     tmp23[t] * tmp22$second[t]
			  }
			}
		    }
		  }
		}

		first <- eighth.tmp$first + eleventh.tmp$first
		second <- second.tmp$first + fifth.tmp$first + ninth.tmp$first +
			    tenth.tmp$first + twelfth.tmp$first
		third <- third.tmp + fourth.tmp + sixth.tmp + seventh.tmp +
			   eighth.tmp$second + eleventh.tmp$second
		fourth <- first.tmp + second.tmp$second + fifth.tmp$second +
			    ninth.tmp$second + tenth.tmp$second +
			    twelfth.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_2_1 <- function(get_x1_x2_f0,get_D_E,env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		second <- array(0,dim=c(k.size,k.size,k.size,block))
		third <- array(0,dim=c(k.size,k.size,block))
		fourth <- matrix(0,k.size,block)

		n1 <- length(get_x1_x2_f0)
		n2 <- sum(get_x1_x2_f0 == 0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i3 in 1:d.size){
		    for(i4 in 1:d.size){
			for(k1 in 1:k.size){
			  for(k2 in 1:k.size){
			    for(k3 in 1:k.size){

				tmp1 <- get_x1_x2_f0[l,i3,i4,] *
					  get_D_E$first[i3,i4,k1,k2,k3,]

				first[l,k1,k2,k3,] <- first[l,k1,k2,k3,] + I0(tmp1,env)/2
			    }

			    tmp2 <- get_x1_x2_f0[l,i3,i4,] *
					get_D_E$second[i3,i4,k1,k2,]

			    second[l,k1,k2,] <- second[l,k1,k2,] + I0(tmp2,env)/2
			  }

			  tmp3 <- get_x1_x2_f0[l,i3,i4,] *
				    get_D_E$third[i3,i4,k1,]

			  third[l,k1,] <- third[l,k1,] + I0(tmp3,env)/2
			}

			tmp4 <- get_x1_x2_f0[l,i3,i4,] *
				  get_D_E$fourth[i3,i4,]

			fourth[l,] <- fourth[l,] + I0(tmp4,env)/2
		    }
		  }
		}

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_2_2 <- function(get_x1_x2_F,get_D_E,env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		second <- array(0,dim=c(k.size,k.size,k.size,block))
		third <- array(0,dim=c(k.size,k.size,block))
		fourth <- matrix(0,k.size,block)

		n1 <- length(get_x1_x2_F)
		n2 <- sum(get_x1_x2_F == 0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i3 in 1:d.size){
		    for(i4 in 1:d.size){
			for(k1 in 1: k.size){
			  for(k2 in 1:k.size){
			    for(k3 in 1:k.size){

				tmp1 <- get_x1_x2_F[l,i3,i4,] *
					  get_D_E$first[i4,i3,k1,k2,k3,]

				first[l,k1,k2,k3,] <- first[l,k1,k2,k3,] + tmp1/2
			    }

			    tmp2 <- get_x1_x2_F[l,i3,i4,] *
					get_D_E$second[i4,i3,k1,k2,]

			    second[l,k1,k2,] <- second[l,k1,k2,] + tmp2/2
			  }

			  tmp3 <- get_x1_x2_F[l,i3,i4,] *
				    get_D_E$third[i4,i3,k1,]

			  third[l,k1,] <- third[l,k1,] + tmp3/2
			}

			tmp4 <- get_x1_x2_F[l,i3,i4,] *
				  get_D_E$fourth[i4,i3,]

			fourth[l,] <- fourth[l,] + tmp4/2
		    }
		  }
		}

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


##p.34(II)-3

	#first:l*k*k*t, second:l*k*t, third:l*t

	F_tilde2_3_1 <- function(tmpY, get_Y_e_V, get_Y_x1_x2_V0, get_Y_e_e_V, get_U_t, get_U_hat_t, get_e_t, get_x_e_f0, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		h.tmp <- get_x_e_f0/2

		first <- array(0,dim=c(k.size,k.size,k.size,block))
		second <- array(0,dim=c(k.size,k.size,block))
		third <- matrix(0,k.size,block)

		first.tmp <- matrix(0,k.size,block)

		n1 <- length(h.tmp)
		n2 <- sum(h.tmp == 0)
		n3 <- length(get_e_t)
		n4 <- sum(get_e_t == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  tmp1 <- double(block)

		  for(i in 1:d.size){
		    tmp1 <- tmp1 + h.tmp[l,i,] * get_e_t[i,]
		  }

		  first.tmp[l,] <- I0(tmp1,env)
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,k.size,block)

		n5 <- length(get_Y_x1_x2_V0)
		n6 <- sum(get_Y_x1_x2_V0 == 0)

		n <- (n1 - n2 != 0) * (n5 - n6 != 0)

		I0.tmp1 <- array(0,dim=c(k.size,d.size,d.size,block))	#l,j1,j2,t

		for(l in 1:k.size){
		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){
			for(i in 1:d.size){
			  for(j in 1:d.size){
			    tmp3 <- double(block)

			    for(i1 in 1:d.size){
				for(i2 in 1:d.size){
				  tmp2 <- get_Y_x1_x2_V0[i1,i2,j,] * tmpY[i1,j1,] *
					    tmpY[i2,j2,]

				  tmp3 <- tmp3 + tmp2
				}
			    }

			    tmp4 <- 2 * h.tmp[l,i,] * tmpY[i,j,] * I0(tmp3,env)
			    tmp5 <- I0(tmp4,env)

			    I0.tmp1[l,j1,j2,] <- I0.tmp1[l,j1,j2,] + tmp5
			  }
			}
		    }
		  }
		}

		
		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){
			tmp6 <- I_12(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)

			for(t in 1:block){
			  second.tmp$first[l,,,t] <- second.tmp$first[l,,,t] +
							     I0.tmp1[l,j1,j2,t] * tmp6$first[,,t]

			  second.tmp$second[l,t] <- second.tmp$second[l,t] +
							    I0.tmp1[l,j1,j2,t] * tmp6$second[t]
			}
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		third.tmp$second <- matrix(0,k.size,block)

		tmp7 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			for(r in 1:r.size){
			  tmp7[r,] <- I0.tmp1[l,j1,j2,] * get_Y_e_V[j1,r,]
			}

			tmp8 <- I_12(tmp7,get_Y_e_V[j2,,],env)

			for(t in 1:block){
			  third.tmp$first[l,,,t] <- third.tmp$first[l,,,t] -
							    tmp8$first[,,t]

			  third.tmp$second[l,t] <- third.tmp$second[l,t] -
							   tmp8$second[t]
			}
		    }
		  }
		}

		fourth.tmp  <- array(0,dim=c(k.size,k.size,block))
		fifth.tmp  <- array(0,dim=c(k.size,k.size,block))

		n7 <- length(get_U_t)
		n8 <- sum(get_U_t == 0)

		n <- (n1 - n2 != 0) * (n7 - n8 != 0)

		tmp12 <- matrix(0,r.size,block)

		I0_U <- array(0,dim=c(d.size,d.size,block))	#j1,j,t

		for(j1 in 1:d.size){
		  for(j in 1:d.size){
		    I0_U[j1,j,] <- I0_U[j1,j,] + I0(get_U_t[j1,j,],env)
		  }
		}

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){

		    tmp9 <- double(block)

		    for(i in 1:d.size){
			for(j in 1:d.size){
			  tmp9 <- tmp9 + 2 * h.tmp[l,i,] * tmpY[i,j,] *
				    I0_U[j1,j,]
			}
		    }

		    tmp10 <- I0(tmp9,env)

		    tmp11 <- I_1(get_Y_e_V[j1,,],env)

		    for(k in 1:k.size){
			fourth.tmp[l,k,] <- fourth.tmp[l,k,] +
						  tmp10 * tmp11[k,]
		    }

		    for(r in 1:r.size){
			tmp12[r,] <- tmp10 * get_Y_e_V[j1,r,]
		    }

		    tmp13 <- I_1(tmp12,env)

		    for(k in 1:k.size){
			  fifth.tmp[l,k,] <- fifth.tmp[l,k,] - tmp13[k,]
		    }
		  }
		}

		sixth.tmp  <- array(0,dim=c(k.size,k.size,block))

		tmp14 <- matrix(0,r.size,block)
		I0.tmp2 <- array(0,dim=c(k.size,d.size,block))	#l,j,t

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    for(i in 1:d.size){
			I0.tmp2[l,j,] <- I0.tmp2[l,j,] + 2 * h.tmp[l,i,] * tmpY[i,j,]
		    }

		    I0.tmp2[l,j,] <- I0(I0.tmp2[l,j,],env)

		    tmp15 <- matrix(0,k.size,block)

		    for(j1 in 1:d.size){

			for(r in 1:r.size){
			  tmp14[r,] <- I0_U[j1,j,] * get_Y_e_V[j1,r,]
			}

			tmp15 <- tmp15 + I_1(tmp14,env)
		    }

		    for(k in 1:k.size){
			sixth.tmp[l,k,] <- sixth.tmp[l,k,] -
						 I0.tmp2[l,j,] * tmp15[k,]
		    }
		  }
		}

		seventh.tmp  <- array(0,dim=c(k.size,k.size,block))

		tmp16 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  tmp17 <- matrix(0,r.size,block)

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){

			for(r in 1:r.size){
			  tmp16[r,] <- I0.tmp2[l,j,] * I0_U[j1,j,] *
					   get_Y_e_V[j1,r,]
			}

			tmp17 <- tmp17 + tmp16
		    }
		  }

		  seventh.tmp[l,,] <- I_1(tmp17,env)
		}

		eighth.tmp <- list()
		eighth.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		eighth.tmp$second <- matrix(0,k.size,block)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){

			tmp18 <- I_12(get_U_hat_t[j1,j,,],get_Y_e_V[j1,,],env)

			for(t in 1:block){
			  eighth.tmp$first[l,,,t] <- eighth.tmp$first[l,,,t] +
							     I0.tmp2[l,j,t] * tmp18$first[,,t]

			  eighth.tmp$second[l,t] <- eighth.tmp$second[l,t] +
							    I0.tmp2[l,j,t] * tmp18$second[t]

			}
		    }
		  }
		}

		ninth.tmp <- list()
		ninth.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		ninth.tmp$second <- matrix(0,k.size,block)

		tmp19 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){
		    for(j1 in 1:d.size){

			for(r in 1:r.size){
			  tmp19[r,] <- I0.tmp2[l,j,] * get_U_hat_t[j1,j,r,]
			}

			tmp20 <- I_12(tmp19,get_Y_e_V[j1,,],env)

			ninth.tmp$first[l,,,] <- ninth.tmp$first[l,,,] -
							 tmp20$first

			ninth.tmp$second[l,] <- ninth.tmp$second[l,] -
							tmp20$second
		    }
		  }
		}

		tenth.tmp <- array(0,dim=c(k.size,k.size,block))

		n9 <- length(get_Y_e_e_V)
		n10 <- sum(get_Y_e_e_V == 0)

		n <- (n1 - n2 != 0) * (n9 - n10 != 0)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j in 1:d.size){

		    tmp21 <- I0.tmp2[l,j,]/2

		    tmp22 <- I_1(get_Y_e_e_V[j,,],env)

		    for(k in 1:k.size){
			 tenth.tmp[l,k,] <- tenth.tmp[l,k,] +
						  tmp21 * tmp22[k,]
		    }
		  }
		}

		eleventh.tmp <- array(0,dim=c(k.size,k.size,block))

		tmp24 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  tmp25 <- matrix(0,r.size,block)

		  for(j in 1:d.size){

		    tmp23 <- I0.tmp2[l,j,]/2

		    for(r in 1:r.size){
			 tmp24[r,] <- tmp23 * get_Y_e_e_V[j,r,]
		    }

		    tmp25 <- tmp25 + tmp24
		  }

		  eleventh.tmp[l,,] <- - I_1(tmp25,env)
		}

		first <- second.tmp$first + third.tmp$first + eighth.tmp$first +
			   ninth.tmp$first
		second <- fourth.tmp + fifth.tmp + sixth.tmp + seventh.tmp +
			    tenth.tmp + eleventh.tmp
		third <- first.tmp + second.tmp$second + third.tmp$second +
			   eighth.tmp$second + ninth.tmp$second

		return(list(first=first,second=second,third=third))
	}
		

	#first:l*k*k*t, second:l*k*t, third:l*t

	F_tilde2_3_2 <- function(get_E0_bar, get_x_e_F, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first <- array(0,dim=c(k.size,k.size,k.size,block))
		second <- array(0,dim=c(k.size,k.size,block))
		third <- matrix(0,k.size,block)

		n1 <- length(get_x_e_F)
		n2 <- sum(get_x_e_F == 0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i in 1:d.size){

		    tmp1 <- get_x_e_F[l,i,]/2

		    for(k1 in 1:k.size){
			for(k2 in 1:k.size){

			  first[l,k1,k2,] <- first[l,k1,k2,] +
						   tmp1 * get_E0_bar$first[i,k1,k2,]
			}

			second[l,k1,] <- second[l,k1,] +
					     tmp1 * get_E0_bar$second[i,k1,]
		    }

		    third[l,] <- third[l,] +
				     tmp1 * get_E0_bar$third[i,]
		  }
		}

		return(list(first=first,second=second,third=third))
	}


##p.34(II)-4

	#first:i1*i2*i3*k*k*k*t, second:i1*i2*i3*k*k*t, third:i1*i2*i3*k*t, fourth:i1*i2*i3*t

	D_a <- function(tmpY, get_Y_e_V, get_Y_D, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(d.size,d.size,d.size,block))

		n1 <- length(get_Y_D)
		n2 <- sum(get_Y_D == 0)

		n <- n1 - n2

		for(i1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(i2 in 1:d.size){
		    for(i3 in 1:d.size){
			first.tmp[i1,i2,i3,] <- get_Y_D[i1,] * get_Y_D[i2,] *
							get_Y_D[i3,]
		    }
		  }
		}

		second.tmp <- array(0,dim=c(d.size,d.size,d.size,k.size,block))

		for(j3 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  tmp1 <- I_1(get_Y_e_V[j3,,],env)

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			for(i3 in 1:d.size){

			  tmp2 <- 3 * get_Y_D[i1,] * get_Y_D[i2,] * tmpY[i3,j3,]

			  for(k in 1:k.size){
			    second.tmp[i1,i2,i3,k,] <- second.tmp[i1,i2,i3,k,] +
								 tmp2 * tmp1[k,]
			  }
			}
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(d.size,d.size,d.size,k.size,k.size,block))
		third.tmp$second <- array(0,dim=c(d.size,d.size,d.size,block))

		for(j2 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(j3 in 1:d.size){

		    tmp3 <- I_1_2(get_Y_e_V[j2,,],get_Y_e_V[j3,,],env)

		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){
			  for(i3 in 1:d.size){

			    tmp4 <- 3 * get_Y_D[i1,] * tmpY[i2,j2,] * tmpY[i3,j3,]

			    for(t in 1:block){
				third.tmp$first[i1,i2,i3,,,t] <- third.tmp$first[i1,i2,i3,,,t] +
									   tmp4[t] * tmp3$first[,,t]
				third.tmp$second[i1,i2,i3,t] <- third.tmp$second[i1,i2,i3,t] +
									  tmp4[t] * tmp3$second[t]
			    }
			  }
			}
		    }
		  }
		}

		fourth.tmp <- list()
		fourth.tmp$first <- array(0,dim=c(d.size,d.size,d.size,k.size,k.size,k.size,block))
		fourth.tmp$second <- array(0,dim=c(d.size,d.size,d.size,k.size,block))

		for(j1 in 1:d.size){
		  for(j2 in 1:d.size){
		    for(j3 in 1:d.size){

			tmp5 <- I_1_2_3(get_Y_e_V[j1,,],get_Y_e_V[j2,,],get_Y_e_V[j3,,],env)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){
			    for(i3 in 1:d.size){

				tmp6 <- tmpY[i1,j1,] * tmpY[i2,j2,] * tmpY[i3,j3,]

				for(t in 1:block){
				  fourth.tmp$first[i1,i2,i3,,,,t] <- fourth.tmp$first[i1,i2,i3,,,,t] +
										 tmp6[t] * tmp5$first[,,,t]
				  fourth.tmp$second[i1,i2,i3,,t] <- fourth.tmp$second[i1,i2,i3,,t] +
										tmp6[t] * tmp5$second[,t]
				}
			    }
			  }
			}
		    }
		  }
		}

		first <- fourth.tmp$first
		second <- third.tmp$first
		third <- second.tmp + fourth.tmp$second
		fourth <- first.tmp + third.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:i1*i2*k*k*t, second:i1*i2*k*t, third:i1*i2*t

	D_b <- function(tmpY, get_Y_e_V, get_Y_D, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(d.size,d.size,block))

		n1 <- length(get_Y_D)
		n2 <- sum(get_Y_D == 0)

		n <- n1 - n2

		for(i1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  for(i2 in 1:d.size){
		    first.tmp[i1,i2,] <- get_Y_D[i1,] * get_Y_D[i2,]
		  }
		}

		second.tmp <- array(0,dim=c(d.size,d.size,k.size,block))

		for(j1 in 1:d.size){

		  if(n == 0){
		    break
		  }

		  tmp1 <- I_1(get_Y_e_V[j2,,],env)

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){

			tmp2 <- 2 * get_Y_D[i1,] * tmpY[i2,j2,]

			for(k in 1:k.size){
			  second.tmp[i1,i2,k,] <- second.tmp[i1,i2,k,] +
							  tmp2 * tmp1[k,]
			}
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(d.size,d.size,k.size,k.size,block))
		third.tmp$second <- array(0,dim=c(d.size,d.size,block))

		for(j1 in 1:d.size){
		  for(j2 in 1:d.size){

		    tmp3 <- I_1_2(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)

		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){

			  tmp4 <- tmpY[i1,j1,] * tmpY[i2,j2,]

			  for(t in 1:block){
			    third.tmp$first[i1,i2,,,t] <- third.tmp$first[i1,i2,,,t] +
								    tmp4[t] * tmp3$first[,,t]
			    third.tmp$second[i1,i2,t] <- third.tmp$second[i1,i2,t] +
								   tmp4[t] * tmp3$second[t]
			  }
			}
		    }
		  }
		}

		first <- third.tmp$first
		second <- second.tmp
		third <- first.tmp + third.tmp$second

		return(list(first=first,second=second,third=third))
	}


	#first:i1*k*t, second:i1*t

	D_c <- function(tmpY, get_Y_e_V, get_Y_D, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- get_Y_D

		second.tmp <- array(0,dim=c(d.size,k.size,block))

		for(j1 in 1:d.size){

		  tmp1 <- I_1(get_Y_e_V[j1,,],env)

		  for(i1 in 1:d.size){

		    tmp2 <- tmpY[i1,j1,]

		    for(k in 1:k.size){
			second.tmp[i1,k,] <- second.tmp[i1,k,] +
						   tmp2 * tmp1[k,]
		    }
		  }
		}

		first <- second.tmp
		second <- first.tmp

		return(list(first=first,second=second))
	}


	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_4_1 <- function(get_x1_x2_x3_f0, get_x1_x2_e_f0, get_x_e_e_f0, get_e_e_e_f0, get_D_a, get_D_b, get_D_c, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- list()
		first.tmp$first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		first.tmp$second <- array(0,dim=c(k.size,k.size,k.size,block))
		first.tmp$third <- array(0,dim=c(k.size,k.size,block))
		first.tmp$fourth <- matrix(0,k.size,block)

		n1 <- length(get_x1_x2_x3_f0)
		n2 <- sum(get_x1_x2_x3_f0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			for(i3 in 1:d.size){

			  tmp1 <- get_x1_x2_x3_f0[l,i1,i2,i3,]

			  for(k1 in 1:k.size){
			    for(k2 in 1:k.size){
				for(k3 in 1:k.size){

				  tmp2 <- tmp1 * get_D_a$first[i1,i2,i3,k1,k2,k3,]

				  first.tmp$first[l,k1,k2,k3,] <- first.tmp$first[l,k1,k2,k3,] +
									    I0(tmp2,env)/6
				}

				tmp3 <- tmp1 * get_D_a$second[i1,i2,i3,k1,k2,]

				first.tmp$second[l,k1,k2,] <- first.tmp$second[l,k1,k2,] +
									I0(tmp3,env)/6
			    }

			    tmp4 <- tmp1 * get_D_a$third[i1,i2,i3,k1,]

			    first.tmp$third[l,k1,] <- first.tmp$third[l,k1,] +
								I0(tmp4,env)/6
			  }

			  tmp5 <- tmp1 * get_D_a$fourth[i1,i2,i3,]

			  first.tmp$fourth[l,] <- first.tmp$fourth[l,] +
							  I0(tmp5,env)/6
			}
		    }
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		second.tmp$second <- array(0,dim=c(k.size,k.size,block))
		second.tmp$third <- matrix(0,k.size,block)

		n3 <- length(get_x1_x2_e_f0)
		n4 <- sum(get_x1_x2_e_f0 == 0)

		n <- n3 - n4

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){

			tmp6 <- get_x1_x2_e_f0[l,i1,i2,]

			for(k1 in 1:k.size){
			  for(k2 in 1:k.size){

			    tmp7 <- tmp6 * get_D_b$first[i1,i2,k1,k2,]

			    second.tmp$first[l,k1,k2,] <- second.tmp$first[l,k1,k2,] +
								    I0(tmp7,env)/2
			  }

			  tmp8 <- tmp6 * get_D_b$second[i1,i2,k1,]

			  second.tmp$second[l,k1,] <- second.tmp$second[l,k1,] +
							      I0(tmp8,env)/2
			}

			tmp9 <- tmp6 * get_D_b$third[i1,i2,]

			second.tmp$third[l,] <- second.tmp$third[l,] +
							   I0(tmp9,env)/2
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(k.size,k.size,block))
		third.tmp$second <- matrix(0,k.size,block)

		n5 <- length(get_x_e_e_f0)
		n6 <- sum(get_x_e_e_f0 == 0)

		n <- n5 - n6

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i1 in 1:d.size){

		    tmp10 <- get_x_e_e_f0[l,i1,]

		    for(k1 in 1:k.size){

			tmp11 <- tmp10 * get_D_c$first[i1,k1,]

			third.tmp$first[l,k1,] <- third.tmp$first[l,k1,] +
							  I0(tmp11,env)/2
		    }

		    tmp12 <- tmp10 * get_D_c$second[i1,]

		    third.tmp$second[l,] <- third.tmp$second[l,] +
						    I0(tmp12,env)/2
		  }
		}

		fourth.tmp <- matrix(0,k.size,block)

		for(l in 1:k.size){
			fourth.tmp[l,] <- I0(get_e_e_e_f0[l,],env)/6
		}

		first <- first.tmp$first
		second <- first.tmp$second + second.tmp$first
		third <- first.tmp$third + second.tmp$second + third.tmp$first
		fourth <- first.tmp$fourth + second.tmp$third + third.tmp$second +
			    fourth.tmp

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_4_2 <- function(get_x1_x2_x3_F, get_x1_x2_e_F, get_x_e_e_F, get_e_e_e_F, get_D_a, get_D_b, get_D_c, env){

		d.size <- env$d.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- list()
		first.tmp$first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		first.tmp$second <- array(0,dim=c(k.size,k.size,k.size,block))
		first.tmp$third <- array(0,dim=c(k.size,k.size,block))
		first.tmp$fourth <- matrix(0,k.size,block)

		n1 <- length(get_x1_x2_x3_F)
		n2 <- sum(get_x1_x2_x3_F == 0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			for(i3 in 1:d.size){

			  tmp1 <- get_x1_x2_x3_F[l,i1,i2,i3,]/6

			  for(k1 in 1:k.size){
			    for(k2 in 1:k.size){
				for(k3 in 1:k.size){

				  tmp2 <- tmp1 * get_D_a$first[i1,i2,i3,k1,k2,k3,]

				  first.tmp$first[l,k1,k2,k3,] <- first.tmp$first[l,k1,k2,k3,] +
									    tmp2
				}

				tmp3 <- tmp1 * get_D_a$second[i1,i2,i3,k1,k2,]

				first.tmp$second[l,k1,k2,] <- first.tmp$second[l,k1,k2,] +
									tmp3
			    }

			    tmp4 <- tmp1 * get_D_a$third[i1,i2,i3,k1,]

			    first.tmp$third[l,k1,] <- first.tmp$third[l,k1,] + tmp4
			  }

			  tmp5 <- tmp1 * get_D_a$fourth[i1,i2,i3,]

			  first.tmp$fourth[l,] <- first.tmp$fourth[l,] + tmp5
			}
		    }
		  }
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		second.tmp$second <- array(0,dim=c(k.size,k.size,block))
		second.tmp$third <- matrix(0,k.size,block)

		n3 <- length(get_x1_x2_e_F)
		n4 <- sum(get_x1_x2_e_F == 0)

		n <- n3 - n4

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){

			tmp6 <- get_x1_x2_e_F[l,i1,i2,]/2

			for(k1 in 1:k.size){
			  for(k2 in 1:k.size){

			    tmp7 <- tmp6 * get_D_b$first[i1,i2,k1,k2,]

			    second.tmp$first[l,k1,k2,] <- second.tmp$first[l,k1,k2,] +
								    tmp7
			  }

			  tmp8 <- tmp6 * get_D_b$second[i1,i2,k1,]

			  second.tmp$second[l,k1,] <- second.tmp$second[l,k1,] +
							      tmp8
			}

			tmp9 <- tmp6 * get_D_b$third[i1,i2,]

			second.tmp$third[l,] <- second.tmp$third[l,] +
							   tmp9
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(k.size,k.size,block))
		third.tmp$second <- matrix(0,k.size,block)

		n5 <- length(get_x_e_e_F)
		n6 <- sum(get_x_e_e_F == 0)

		n <- n5 - n6

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(i1 in 1:d.size){

		    tmp10 <- get_x_e_e_F[l,i1,]/2

		    for(k1 in 1:k.size){

			tmp11 <- tmp10 * get_D_c$first[i1,k1,]

			third.tmp$first[l,k1,] <- third.tmp$first[l,k1,] +
							  tmp11
		    }

		    tmp12 <- tmp10 * get_D_c$second[i1,]

		    third.tmp$second[l,] <- third.tmp$second[l,] +
						      tmp12
		  }
		}

		fourth.tmp <- matrix(0,k.size,block)

		for(l in 1:k.size){
			fourth.tmp[l,] <- I0(get_e_e_e_F[l,],env)/6
		}

		first <- first.tmp$first
		second <- first.tmp$second + second.tmp$first
		third <- first.tmp$third + second.tmp$second + third.tmp$first
		fourth <- first.tmp$fourth + second.tmp$third + third.tmp$second +
			    fourth.tmp

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


##p.36(II)-5

	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_5 <- function(tmpY, get_Y_e_V, get_Y_x1_x2_V0, get_Y_e_e_V, get_e_t, get_U_t, get_U_hat_t, get_x_e_f, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(k.size,k.size,block))

		n1 <- length(get_x_e_f)
		n2 <- sum(get_x_e_f == 0)
		n3 <- length(get_e_t)
		n4 <- sum(get_e_t == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		tmp1 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  tmp2 <- matrix(0,r.size,block)

		  for(i in 1:d.size){

		    for(r in 1:r.size){
			tmp1[r,] <- get_x_e_f[l,i,r,] * get_e_t[i,]
		    }

		    tmp2 <- tmp2 + tmp1
		  }

		  first.tmp[l,,] <- first.tmp[l,,] + I_1(tmp2,env)/2
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		second.tmp$second <- array(0,dim=c(k.size,k.size,block))

		n5 <- length(get_Y_x1_x2_V0)
		n6 <- sum(get_Y_x1_x2_V0 == 0)

		n <- (n1 - n2 != 0) * (n5 - n6 != 0)

		tmp5 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			tmp6 <- matrix(0,r.size,block)

			for(j3 in 1:d.size){
			  for(i3 in 1:d.size){

			    tmp4 <- double(block)

			    for(i1 in 1:d.size){
				for(i2 in 1:d.size){

				  tmp3 <- get_Y_x1_x2_V0[i1,i2,j3,] * tmpY[i1,j1,] *
					    tmpY[i2,j2,]

				  tmp4 <- tmp4 + I0(tmp3,env)
				}
			    }

			    for(r in 1:r.size){
				tmp5[r,] <- get_x_e_f[l,i3,r,] * 
						tmpY[i3,j3,] * tmp4
			    }

			    tmp6 <- tmp6 + tmp5
			  }
			}

			tmp7 <- I_123(tmp6,get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)

			second.tmp$first[l,,,,] <- second.tmp$first[l,,,,] +
							   tmp7$first
			second.tmp$second[l,,] <- second.tmp$second[l,,] +
							  tmp7$second
		    }
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		third.tmp$second <- matrix(0,k.size,block)
		fourth.tmp <- list()
		fourth.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		fourth.tmp$second <- matrix(0,k.size,block)
		fifth.tmp <- list()
		fifth.tmp$first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		fifth.tmp$second <- array(0,dim=c(k.size,k.size,block))
		sixth.tmp <- list()
		sixth.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		sixth.tmp$second <- matrix(0,k.size,block)

		n5 <- length(get_U_t)
		n6 <- sum(get_U_t == 0)

		n <- n5 - n6

		f_Y <- array(0,dim=c(k.size,d.size,r.size,block))	#l,j3,r,t

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j3 in 1:d.size){
		    for(i3 in 1:d.size){
			for(r in 1:r.size){
			  f_Y[l,j3,r,] <- f_Y[l,j3,r,] +
						get_x_e_f[l,i3,r,] * tmpY[i3,j3,]
			}
		    }
		  }
		}

		tmp9 <- matrix(0,r.size,block)
		tmp13 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){

		    tmp10 <- matrix(0,r.size,block)

		    for(j3 in 1:d.size){

			tmp8 <- I0(get_U_t[j1,j3,],env)

			for(r in 1:r.size){
			  tmp9[r,] <- f_Y[l,j3,r,] * tmp8
			}

			tmp10 <- tmp10 + tmp9
		    }

		    tmp11 <- I_12(tmp10,get_Y_e_V[j1,,],env)

		    third.tmp$first[l,,,] <- third.tmp$first[l,,,] +
						     tmp11$first
		    third.tmp$second[l,] <- third.tmp$second[l,] +
						    tmp11$second
		  }

		  for(j3 in 1:d.size){

		    tmp14 <- matrix(0,r.size,block)

		    for(j1 in 1:d.size){

			tmp12 <- I0(get_U_t[j1,j3,],env)

			for(r in 1:r.size){
			  tmp13[r,] <- tmp12 * get_Y_e_V[j1,r,]
			}

			tmp14 <- tmp14 + tmp13
		    }

		    tmp15 <- I_12(f_Y[l,j3,,],tmp14,env)

		    fourth.tmp$first[l,,,] <- fourth.tmp$first[l,,,] -
							tmp15$first
		    fourth.tmp$second[l,] <- fourth.tmp$second[l,] -
						     tmp15$second
		  }

		  for(j3 in 1:d.size){
		    for(j1 in 1:d.size){

			tmp16 <- I_123(f_Y[l,j3,,],get_U_hat_t[j1,j3,,],get_Y_e_V[j1,,],env)

			fifth.tmp$first[l,,,,] <- fifth.tmp$first[l,,,,] +
							  tmp16$first
			fifth.tmp$second[l,,] <- fifth.tmp$second[l,,] +
							 tmp16$second
		    }

		    tmp17 <- I_12(f_Y[l,j3,,]/2,get_Y_e_e_V[j3,,],env)

		    sixth.tmp$first[l,,,] <- sixth.tmp$first[l,,,] +
						     tmp17$first
		    sixth.tmp$second[l,] <- sixth.tmp$second[l,] +
						    tmp17$second
		  }
		}

		first <- second.tmp$first + fifth.tmp$first
		second <- third.tmp$first + fourth.tmp$first + sixth.tmp$first
		third <- first.tmp + second.tmp$second + fifth.tmp$second
		fourth <- third.tmp$second + fourth.tmp$second + sixth.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


##p.36(II)-6

	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2_6 <- function(tmpY, get_Y_e_V, get_Y_D, get_x1_x2_e_f, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(k.size,k.size,block))

		n1 <- length(get_x1_x2_e_f)
		n2 <- sum(get_x1_x2_e_f == 0)
		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		tmp1 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  tmp2 <- matrix(0,r.size,block)

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){

			for(r in 1:r.size){
			  tmp1[r,] <- get_x1_x2_e_f[l,i1,i2,r,] *
					  get_Y_D[i1,] * get_Y_D[i2,]
			}

			tmp2 <- tmp2 + tmp1
		    }
		  }

		  tmp3 <- I_1(tmp2,env)/2

		  first.tmp[l,,] <- first.tmp[l,,] + tmp3
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,k.size,block)

		tmp4 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j2 in 1:d.size){

		    tmp5 <- matrix(0,r.size,block)

		    for(i1 in 1:d.size){
			for(i2 in 1:d.size){

			  for(r in 1:r.size){
			    tmp4[r,] <- get_x1_x2_e_f[l,i1,i2,r,] *
					    get_Y_D[i1,] * tmpY[i2,j2,]
			  }

			  tmp5 <- tmp5 + tmp4
			}
		    }

		    tmp6 <- I_12(tmp5,get_Y_e_V[j2,,],env)

		    second.tmp$first[l,,,] <- second.tmp$first[l,,,] +
							tmp6$first
		    second.tmp$second[l,] <- second.tmp$second[l,] +
						     tmp6$second
		  }
		}

		third.tmp <- list()
		third.tmp$first <- array(0,dim=c(k.size,k.size,k.size,k.size,block))
		third.tmp$second <- array(0,dim=c(k.size,k.size,block))

		n <- n1 - n2

		tmp7 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			tmp8 <- matrix(0,r.size,block)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){

			    for(r in 1:r.size){
				tmp7[r,] <- get_x1_x2_e_f[l,i1,i2,r,] *
						tmpY[i1,j1,] * tmpY[i2,j2,]
			    }

			    tmp8 <- tmp8 + tmp7
			  }
			}

			tmp9 <- I_123(tmp8,get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)

			third.tmp$first[l,,,,] <- third.tmp$first[l,,,,] +
							  tmp9$first
			third.tmp$second[l,,] <- third.tmp$second[l,,] +
							 tmp9$second
		    }
		  }
		}

		fourth.tmp <- array(0,dim=c(k.size,k.size,block))

		tmp11 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  tmp12 <- matrix(0,r.size,block)

		  for(j1 in 1:d.size){
		    for(j2 in 1:d.size){

			tmp10 <- b1_b2(get_Y_e_V[j1,,],get_Y_e_V[j2,,],env)

			for(i1 in 1:d.size){
			  for(i2 in 1:d.size){

			    for(r in 1:r.size){
				tmp11[r,] <- get_x1_x2_e_f[l,i1,i2,r,] *
						 tmpY[i1,j1,] * tmpY[i2,j2,] *
						 tmp10
			    }

			    tmp12 <- tmp12 + tmp11
			  }
			}
		    }
		  }

		  tmp13 <- I_1(tmp12,env)/2

		  fourth.tmp[l,,] <- fourth.tmp[l,,] + tmp13
		}

		first <- third.tmp$first
		second <- second.tmp$first
		third <- first.tmp + third.tmp$second + fourth.tmp
		fourth <- second.tmp$second

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


##p.36(II)-7

	#first:l*k*k*t, second:l*k*t, third:l*t

	F_tilde2_7 <- function(tmpY, get_Y_e_V, get_Y_D, get_x_e_e_f, env){

		d.size <- env$d.size
		r.size <- env$r.size
		k.size <- env$k.size
		block <- env$block

		first.tmp <- array(0,dim=c(k.size,k.size,block))

		n1 <- length(get_x_e_e_f)
		n2 <- sum(get_x_e_e_f == 0)
		n3 <- length(get_Y_D)
		n4 <- sum(get_Y_D == 0)

		n <- (n1 - n2 != 0) * (n3 - n4 != 0)

		tmp1 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  tmp2 <- matrix(0,r.size,block)

		  for(i1 in 1:d.size){

		    for(r in 1:r.size){
			tmp1[r,] <- get_x_e_e_f[l,i1,r,] * get_Y_D[i1,]
		    }

		    tmp2 <- tmp2 + tmp1
		  }

		  tmp3 <- I_1(tmp2,env)/2

		  first.tmp[l,,] <- first.tmp[l,,] + tmp3
		}

		second.tmp <- list()
		second.tmp$first <- array(0,dim=c(k.size,k.size,k.size,block))
		second.tmp$second <- matrix(0,k.size,block)

		n <- n1 - n2

		tmp4 <- matrix(0,r.size,block)

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  for(j1 in 1:d.size){

		    tmp5 <- matrix(0,r.size,block)

		    for(i1 in 1:d.size){

			for(r in 1:r.size){
			  tmp4[r,] <- get_x_e_e_f[l,i1,r,] * tmpY[i1,j1,]
			}

			tmp5 <- tmp5 + tmp4
		    }

		    tmp6 <- I_12(tmp5/2,get_Y_e_V[j1,,],env)

		    second.tmp$first[l,,,] <- second.tmp$first[l,,,] +
							tmp6$first
		    second.tmp$second[l,] <- second.tmp$second[l,] +
						     tmp6$second
		  }
		}

		first <- second.tmp$first
		second <- first.tmp
		third <- second.tmp$second

		return(list(first=first,second=second,third=third))
	}


##p.37(II)-8

	#l*k*t

	F_tilde2_8 <- function(get_e_e_e_f,env){

		k.size <- env$k.size
		block <- env$block

		result <- array(0,dim=c(k.size,k.size,block))

		n1 <- length(get_e_e_e_f)
		n2 <- sum(get_e_e_e_f == 0)

		n <- n1 - n2

		for(l in 1:k.size){

		  if(n == 0){
		    break
		  }

		  result[l,,] <- I_1(get_e_e_e_f[l,,],env)/6
		}

		return(result)
	}


##p.6(1.19)

	#first:l*k*k*k*t, second:l*k*k*t, third:l*k*t, fourth:l*t

	F_tilde2 <- function(get_F_tilde2_1_1, get_F_tilde2_1_2, get_F_tilde2_2_1,
				   get_F_tilde2_2_2, get_F_tilde2_3_1, get_F_tilde2_3_2,
				   get_F_tilde2_4_1, get_F_tilde2_4_2, get_F_tilde2_5,
				   get_F_tilde2_6, get_F_tilde2_7, get_F_tilde2_8){

		result1_1 <- get_F_tilde2_1_1
		result1_2 <- get_F_tilde2_1_2
		result2_1 <- get_F_tilde2_2_1
		result2_2 <- get_F_tilde2_2_2
		result3_1 <- get_F_tilde2_3_1
		result3_2 <- get_F_tilde2_3_2
		result4_1 <- get_F_tilde2_4_1
		result4_2 <- get_F_tilde2_4_2
		result5 <- get_F_tilde2_5
		result6 <- get_F_tilde2_6
		result7 <- get_F_tilde2_7
		result8 <- get_F_tilde2_8

		first <- result1_1$first + result1_2$first + result2_1$first +
			   result2_2$first + result4_1$first + result4_2$first +
			   result5$first + result6$first

		second <- result1_1$second + result1_2$second + result2_1$second +
			    result2_2$second + result3_1$first + result3_2$first +
			    result4_1$second + result4_2$second + result5$second +
			    result6$second + result7$first

		third <- result1_1$third + result1_2$third + result2_1$third +
			   result2_2$third + result3_1$second + result3_2$second +
			   result4_1$third + result4_2$third + result5$third +
			   result6$third + result7$second + result8

		fourth <- result1_1$fourth + result1_2$fourth + result2_1$fourth +
			    result2_2$fourth + result3_1$third + result3_2$third +
			    result4_1$fourth + result4_2$fourth + result5$fourth +
			    result6$fourth + result7$third

		return(list(first=first,second=second,third=third,
				fourth=fourth))
	}


	F_tilde2_x <- function(l,x,get_F_tilde2,env){

		k.size <- env$k.size
		block <- env$block

		first <- array(get_F_tilde2$first[l,,,,],
				   dim=c(k.size,k.size,k.size,block))

		result1 <- 0

		for(k1 in 1:k.size){
		for(k2 in 1:k.size){
		    for(k3 in 1:k.size){
			result1 <- result1 + first[k1,k2,k3,block] *
				     x[k1] * x[k2] * x[k3]
		}
		}
		}


		second <- array(get_F_tilde2$second[l,,,],
				    dim=c(k.size,k.size,block))

		result2 <- 0

		for(k1 in 1:k.size){
		for(k2 in 1:k.size){
		    result2 <- result2 + second[k1,k2,block] *
				   x[k1] * x[k2]
		  }
		}


		third <- matrix(get_F_tilde2$third[l,,],
				    k.size,block)

		result3 <- 0

		for(k1 in 1:k.size){
		  result3 <- result3 + third[k1,block] * x[k1]
		}


		fourth <- get_F_tilde2$fourth[l,block]

		result <- result1 + result2 + result3 + fourth

		return(result)
	}


##p.38(III)

	#d*t

	x_rho <- function(X.t0,tmp.rho,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,block))

		assign(pars[1],0)

		de.rho <- list()

		for(i in 1:d.size){
		  de.rho[[i]] <- parse(text=deparse(D(tmp.rho,state[i])))
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i in 1:d.size){
		    result[i,t] <- eval(de.rho[[i]])
		  }
		}
		return(result)
	}


	#t

	e_rho <- function(X.t0,tmp.rho,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(block))

		assign(pars[1],0)

		tmp <- parse(text=deparse(D(tmp.rho,pars[1])))

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }
		  result[t] <- eval(tmp)
		}
           
		return(result)
	}


	#d*d*t

	x1_x2_rho <- function(X.t0,tmp.rho,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,d.size,block))

		assign(pars[1],0)

		dxdx.rho <- list()

		for(i1 in 1:d.size){
		  dxdx.rho[[i1]] <- list()

		  for(i2 in 1:d.size){
		    dxdx.rho[[i1]][[i2]] <- parse(text=deparse(D(D(tmp.rho,state[i2]),state[i1])))
		  }
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i1 in 1:d.size){
		    for(i2 in 1:d.size){
			result[i1,i2,t] <- eval(dxdx.rho[[i1]][[i2]])
		    }
		  }
		}

		return(result)
	}


	#d*t

	x_e_rho <- function(X.t0,tmp.rho,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- array(0,dim=c(d.size,block))

		assign(pars[1],0)

		dxde.rho <- list()

		for(i in 1:d.size){
		  dxde.rho[[i]] <- parse(text=deparse(D(D(tmp.rho,pars[1]),state[i])))
		}

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  for(i in 1:d.size){
		    result[i,t] <- eval(dxde.rho[[i]])
		  }
		}

		return(result)
	}


	#t

	e_e_rho <- function(X.t0,tmp.rho,env){

		d.size <- env$d.size
		state <- env$state
		pars <- env$pars
		block <- env$block
		my.range <- env$my.range

		result <- double(block)

		assign(pars[1],0)

		tmp <- parse(text=deparse(D(D(tmp.rho,pars[1]),pars[1])))

		for(t in 1:block){
		  for(d in 1:d.size){
		    assign(state[d],X.t0[my.range[t],d])
		  }

		  result[t] <- eval(tmp)
		}

		return(result)
	}


	#numeric

	scH_0_1 <- function(get_Y_D,get_e_rho,get_x_rho,env){

		d.size <- env$d.size
		block <- env$block

		result <- - I0(get_e_rho,env)[block]

		for(i in 1:d.size){
		  result <- result - I0(get_x_rho[i,] * get_Y_D[i,],env)[block]
		}

		return(result)
	}


	#d*t/d*r*t

	scH_1_1 <- function(tmpY,get_Y_e_V,get_x_rho,env){

		d.size <- env$d.size
		block <- env$block

		coef <- array(0,dim=c(d.size,block))

		for(j in 1:d.size){
		  tmp <- double(block)

		  for(t in 1:block){
		    tmp[t] <- get_x_rho[,t] %*% tmpY[,j,t]
		  }

#
#		  for(i in 1:d.size){
#		    tmp <- tmp + get_x_rho[i,] * tmpY[i,j,]
#		  }

		  coef[j,] <- I0(tmp,env)
		}

		integrand <- get_Y_e_V

		return(list(coef=coef,integrand=integrand))
	}


	#r*t

	scH_1_2 <- function(get_scH_1_1,env){

		r.size <- env$r.size
		block <- env$block

		result <- array(0,dim=c(r.size,block))

		coef <- get_scH_1_1$coef
		integrand <- get_scH_1_1$integrand

		for(t in 1:block){
		  for(r in 1:r.size){
		    result[r,t] <- integrand[,r,t] %*% coef[,t]
		  }
		}

		return(result)
	}


	F_tilde1H1 <- function(get_F_tilde1,get_scH_0_1,get_scH_1_1,get_scH_1_2,env){

		k.size <- env$k.size
		d.size <- env$d.size
		r.size <- env$r.size
		block <- env$block

		coef1 <- get_F_tilde1$result1[[1]][,block]
		coef2 <- get_F_tilde1$result2[[1]][,block]
		coef3 <- get_F_tilde1$result3[[1]][,block]
		coef4 <- get_F_tilde1$result4[[1]][,block]
		coef5 <- get_F_tilde1$result5[[1]][,block]
		coef6 <- get_F_tilde1$result6[[1]][,block]
		coef7 <- get_F_tilde1$result7[[1]][,block]
		coef8 <- array(get_F_tilde1$result8[[1]][,,block],dim=c(k.size,d.size))
		coef9 <- array(get_F_tilde1$result9[[1]][,,block],dim=c(k.size,d.size))
		coef10 <- get_F_tilde1$result10[[1]]
		coef11 <- array(get_F_tilde1$result11[[1]][,,block],dim=c(k.size,d.size))
		coef12 <- get_F_tilde1$result12[[1]]
		coef13 <- array(get_F_tilde1$result13[[1]][,,block],dim=c(k.size,d.size))
		coef14 <- get_F_tilde1$result14[[1]]
		coef15 <- array(get_F_tilde1$result15[[1]][,,block],dim=c(k.size,d.size))
		coef16 <- get_F_tilde1$result16[[1]]
		coef17 <- array(get_F_tilde1$result17[[1]][,,block],dim=c(k.size,d.size))
		coef18 <- array(get_F_tilde1$result18[[1]][,,block],dim=c(k.size,d.size))
		coef19 <- array(get_F_tilde1$result19[[1]][,,block],dim=c(k.size,d.size))
		coef20 <- array(get_F_tilde1$result20[[1]][,,block],dim=c(k.size,d.size))
		coef21 <- array(get_F_tilde1$result21[[1]][,,block],dim=c(k.size,d.size))
		coef22 <- array(get_F_tilde1$result22[[1]][,,,block],dim=c(k.size,d.size,d.size))
		coef23 <- get_F_tilde1$result23[[1]]
		coef24 <- array(get_F_tilde1$result24[[1]][,,,block],dim=c(k.size,d.size,d.size))
		coef25 <- get_F_tilde1$result25[[1]]
		coef26 <- array(get_F_tilde1$result26[[1]][,,block],dim=c(k.size,d.size))
		coef27 <- get_F_tilde1$result27[[1]]
		coef28 <- array(get_F_tilde1$result28[[1]][,,,block],dim=c(k.size,d.size,d.size))
		coef29 <- array(get_F_tilde1$result29[[1]][,,block],dim=c(k.size,d.size))
		coef30 <- array(get_F_tilde1$result30[[1]][,,,block],dim=c(k.size,d.size,d.size))

		integrand1 <- get_F_tilde1$result1[[2]]
		integrand2 <- get_F_tilde1$result2[[2]] 
		integrand3 <- get_F_tilde1$result3[[2]]
		integrand4 <- get_F_tilde1$result4[[2]]
		integrand5 <- get_F_tilde1$result5[[2]]
		integrand6 <- get_F_tilde1$result6[[2]]
		integrand7 <- get_F_tilde1$result7[[2]]
		integrand8 <- get_F_tilde1$result8[[2]]
		integrand9 <- get_F_tilde1$result9[[2]]
		integrand10 <- get_F_tilde1$result10[[2]]
		integrand11 <- get_F_tilde1$result11[[2]]
		integrand12 <- get_F_tilde1$result12[[2]]
		integrand13 <- get_F_tilde1$result13[[2]]
		integrand14 <- get_F_tilde1$result14[[2]]
		integrand15 <- get_F_tilde1$result15[[2]]
		integrand16 <- get_F_tilde1$result16[[2]]
		integrand17 <- get_F_tilde1$result17[[2]]
		integrand18 <- get_F_tilde1$result18[[2]]
		integrand19 <- get_F_tilde1$result19[[2]]
		integrand20 <- get_F_tilde1$result20[[2]]
		integrand21 <- get_F_tilde1$result21[[2]]
		integrand22 <- get_F_tilde1$result22[[2]]
		integrand23 <- get_F_tilde1$result23[[2]]
		integrand24 <- get_F_tilde1$result24[[2]]
		integrand25 <- get_F_tilde1$result25[[2]]
		integrand26 <- get_F_tilde1$result26[[2]]
		integrand27 <- get_F_tilde1$result27[[2]]
		integrand28 <- get_F_tilde1$result28[[2]]
		integrand29 <- get_F_tilde1$result29[[2]]
		integrand30 <- get_F_tilde1$result30[[2]]

		coef <- get_scH_1_1$coef[,block]
		integrand <- get_scH_1_1$integrand

##modified
		obj1 <- coef1 + coef2 + coef3 + coef4 + coef5 + coef6 + coef7
		obj2 <- coef10 * integrand10 + coef12 * integrand12 +
			  coef14 * integrand14 + coef16 * integrand16  
		obj3 <- array(0,dim=c(k.size,r.size,block))     

		result0 <- obj1 * get_scH_0_1         

		result1 <- array(0,dim=c(k.size,k.size,block))
		result3 <- list()

		for(l in 1:k.size){
##modified
		  tmp1 <-  coef10 * integrand10[l,,] + coef12 * integrand12[l,,] +
			     coef14 * integrand14[l,,] + coef16 * integrand16[l,,]

		  tmp1 <- obj2[l,,]

		  tmp2 <- array(0,dim=c(r.size,block))

		  for(j in 1:d.size){
		    tmp2 <- tmp2 + coef8[l,j] * integrand8[j,,] +
				coef9[l,j] * integrand9[j,,] +
				coef11[l,j] * integrand11[j,,] +
				coef13[l,j] * integrand13[j,,] +
				coef15[l,j] * integrand15[j,,] +
				coef17[l,j] * integrand17[j,,] +
				coef18[l,j] * integrand18[j,,] +
				coef19[l,j] * integrand19[j,,] +
				coef20[l,j] * integrand20[j,,] +
				coef21[l,j] * integrand21[j,,]

		    obj3[l,,] <- as.matrix(obj3[l,,] + coef[j] * integrand[j,,])
		  } 
                  
		  tmp3 <- get_scH_0_1 * (tmp1 + tmp2) +
			    obj1[l] * (obj3[l,,] + get_scH_1_2)

		  result1[l,,] <- I_1(array(tmp3,dim=c(r.size,block)),env)
		  result3[[l]] <- I_1_2(tmp1 + tmp2,obj3[l,,] + get_scH_1_2,env)

		}

		result2 <- list()
		result4 <- list()

		for(l in 1:k.size){

		  result2[[l]] <- list(first=array(0,dim=c(k.size,k.size,block)),second=double(block))
		  result4[[l]] <- list(first=array(0,dim=c(k.size,k.size,k.size,block)),second=array(0,dim=c(k.size,block)))

		  for(j1 in 1:d.size){
		    tmp23.1 <- I_12(get_scH_0_1 * integrand23[[1]][l,j1,,],integrand23[[2]][j1,,],env)
		    tmp25.1 <- I_12(get_scH_0_1 * integrand25[[1]][l,j1,,],integrand25[[2]][j1,,],env)
		    tmp27.1 <- I_12(get_scH_0_1 * integrand27[[1]][l,j1,,],integrand27[[2]][j1,,],env)

		    result2[[l]]$first <- result2[[l]]$first + coef23 * tmp23.1$first +
									     coef25 * tmp25.1$first +
									     coef27 * tmp27.1$first

		    result2[[l]]$second <- result2[[l]]$second + coef23 * tmp23.1$second +
									       coef25 * tmp25.1$second +
									       coef27 * tmp27.1$second

		    tmp23.2 <- I_1_23(get_scH_1_2 + obj3[l,,],integrand23[[1]][l,j1,,],integrand23[[2]][j1,,],env)
		    tmp25.2 <- I_1_23(get_scH_1_2 + obj3[l,,],integrand25[[1]][l,j1,,],integrand25[[2]][j1,,],env)
		    tmp27.2 <- I_1_23(get_scH_1_2 + obj3[l,,],integrand27[[1]][l,j1,,],integrand27[[2]][j1,,],env)

		    result4[[l]]$first <- result4[[l]]$first + coef23 * tmp23.2$first +
									     coef25 * tmp25.2$first +
									     coef27 * tmp27.2$first

		    result4[[l]]$second <- result4[[l]]$second + coef23 * tmp23.2$second +
									       coef25 * tmp25.2$second +
									       coef27 * tmp27.2$second

		    obj22 <- array(0,dim=c(r.size,block))
		    obj24 <- array(0,dim=c(r.size,block))
		    obj26 <- array(0,dim=c(r.size,block))
		    obj28 <- array(0,dim=c(r.size,block))
		    obj29 <- array(0,dim=c(r.size,block))
		    obj30 <- array(0,dim=c(r.size,block))
 
		    for(j2 in 1:d.size){
		      obj22 <- obj22 + coef22[l,j1,j2] * integrand22[[2]][j2,,]
		      obj24 <- obj24 + coef24[l,j1,j2] * integrand24[[2]][j2,,]
		      obj28 <- obj28 + coef28[l,j1,j2] * integrand28[[2]][j2,,]
		      obj30 <- obj30 + coef30[l,j1,j2] * integrand30[[2]][j2,,]
                                      
		      obj26 <- obj26 + coef26[l,j2] * integrand26[[1]][j1,j2,,]
		      obj29 <- obj29 + coef29[l,j2] * integrand29[[1]][j1,j2,,]
		    }

		    tmp22.1 <- I_12(get_scH_0_1 * integrand22[[1]][j1,,],array(obj22,dim=c(r.size,block)),env)
		    tmp24.1 <- I_12(get_scH_0_1 * integrand24[[1]][j1,,],array(obj24,dim=c(r.size,block)),env)
		    tmp28.1 <- I_12(get_scH_0_1 * integrand28[[1]][j1,,],array(obj28,dim=c(r.size,block)),env)
		    tmp30.1 <- I_12(get_scH_0_1 * integrand30[[1]][j1,,],array(obj30,dim=c(r.size,block)),env)
		    tmp26.1 <- I_12(get_scH_0_1 * array(obj26,dim=c(r.size,block)),integrand26[[2]][j1,,],env)
		    tmp29.1 <- I_12(get_scH_0_1 * array(obj29,dim=c(r.size,block)),integrand29[[2]][j1,,],env)

		    result2[[l]]$first <- result2[[l]]$first + tmp22.1$first +
									     tmp29.1$first +
									     tmp24.1$first +
									     tmp26.1$first +
									     tmp28.1$first +
									     tmp30.1$first

		    result2[[l]]$second <- result2[[l]]$second + tmp22.1$second +
									       tmp29.1$second +
									       tmp24.1$second +
									       tmp26.1$second +
									       tmp28.1$second +
									       tmp30.1$second
                       
		    tmp22.2 <- I_1_23(get_scH_1_2 + obj3[l,,],integrand22[[1]][j1,,],array(obj22,dim=c(r.size,block)),env)
		    tmp24.2 <- I_1_23(get_scH_1_2 + obj3[l,,],integrand24[[1]][j1,,],array(obj24,dim=c(r.size,block)),env)
		    tmp28.2 <- I_1_23(get_scH_1_2 + obj3[l,,],integrand28[[1]][j1,,],array(obj28,dim=c(r.size,block)),env)
		    tmp30.2 <- I_1_23(get_scH_1_2 + obj3[l,,],integrand30[[1]][j1,,],array(obj30,dim=c(r.size,block)),env)
		    tmp26.2 <- I_1_23(get_scH_1_2 + obj3[l,,],array(obj26,dim=c(r.size,block)),integrand26[[2]][j1,,],env)
		    tmp29.2 <- I_1_23(get_scH_1_2 + obj3[l,,],array(obj29,dim=c(r.size,block)),integrand29[[2]][j1,,],env)

		    result4[[l]]$first <- result4[[l]]$first + tmp22.2$first +
									     tmp29.2$first +
									     tmp24.2$first +
									     tmp26.2$first +
									     tmp28.2$first +
									     tmp30.2$first

		    result4[[l]]$second <- result4[[l]]$second + tmp22.2$second +
									       tmp29.2$second +
									       tmp24.2$second +
									       tmp26.2$second +
									       tmp28.2$second +
									       tmp30.2$second
		  }
		}

		return(list(result0=result0,result1=result1,result2=result2,
				result3=result3,result4=result4))
	}


	F_tilde1H1_x <- function(l,x,get_F_tilde1H1,env){

		block <- env$block

		result0 <- get_F_tilde1H1$result0
		result1 <- get_F_tilde1H1$result1
		result2 <- get_F_tilde1H1$result2
		result3 <- get_F_tilde1H1$result3
		result4 <- get_F_tilde1H1$result4

		result <- result0[l] + I_1_x(x,result1[l,,],env)[block] +
					     I_12_x(x,result2[[l]],env)[block] +
					     I_1_2_x(x,result3[[l]],env)[block] +
					     I_1_23_x(x,result4[[l]],env)[block]

		return(result)

	}  


##p.39(IV)-1

	#result1:numeric, result2:i*j/i*j;list, result3:numeric/list,
	#result4:i/i;list, result5:i/i*k*t, result6:numeric/k*t

	H2_1 <- function(get_scH_0_1,get_scH_1_1,get_scH_1_2,env){

		d.size <- env$d.size
		k.size <- env$k.size
		r.size <- env$r.size
		block <- env$block

		coef <- get_scH_1_1$coef
		integrand <- get_scH_1_1$integrand

		result1 <- get_scH_0_1^2  #mathcal{H}_{0;1}^2
		result2 <- list()         #mathcal{H}_{1;1}^2
		result3 <- list()         #mathcal{H}_{1;2}^2
		result4 <- list()         #2mathcal{H}_{0;1}mathcal{H}_{1;1}
		result5 <- list()         #2mathcal{H}_{1;1}mathcal{H}_{1;2}
		result6 <- list()         #2mathcal{H}_{1;2}mathcal{H}_{0;1}

		result2[[1]] <- matrix(coef[,block],d.size,1) %*% matrix(coef[,block],1,d.size)
		result3[[1]] <- 1
		result4[[1]] <- 2 * get_scH_0_1 * coef[,block]
		result5[[1]] <- 2 * coef[,block]
		result6[[1]] <- 2 * get_scH_0_1

		result2[[2]] <- array(list(),dim=c(d.size,d.size))
		result3[[2]] <- I_1_2(get_scH_1_2,get_scH_1_2,env)
		result4[[2]] <- array(0,dim=c(d.size,k.size,block))
		result5[[2]] <- list()
		result6[[2]] <- I_1(get_scH_1_2,env)

		for(i in 1:d.size){
		  for(j in 1:d.size){
		    result2[[2]][i,j][[1]] <- I_1_2(integrand[i,,],integrand[j,,],env)
		  }
		  result4[[2]][i,,] <- I_1(integrand[i,,],env)
		  result5[[2]][[i]] <- I_1_2(integrand[i,,],get_scH_1_2,env)
		} 

		return(list(result1=result1,result2=result2,result3=result3,
				result4=result4,result5=result5,result6=result6))

	}


##p.39(IV)-2
##remark:I0 is removed

	#result1:i1*i2*t/i1*i2*k*k*t & i1*i2*k*t & i1*i2*t,
	#result2:i1*t/i1*k*t & i1*t,
	#result3:t/numeric

	H2_2 <- function(get_D_b,get_D_c,get_x_x_rho,get_x_e_rho,get_e_e_rho){

		result1 <- list()
		result2 <- list()
		result3 <- list()

		result1[[1]] <- get_x_x_rho
		result2[[1]] <- 2 * get_x_e_rho
		result3[[1]] <- get_e_e_rho

		result1[[2]] <- get_D_b
		result2[[2]] <- get_D_c
		result3[[2]] <- 0

		return(list(result1=result1,result2=result2,result3=result3))
	}
         

##p.39(IV)-3
##remark:I0&H0_T is removed
     
	#coef:i*t, integrand:i*k*k*t & i*k*t & i*t                 

	H2_3 <- function(get_E0_bar,get_x_rho){
		return(list(coef=get_x_rho,integrand=get_E0_bar))
	}


##p.38(IV)
##remark:H0_T is removed

	H2_x <- function(x,get_H2_1,get_H2_2,get_H2_3,env){

		k.size <- env$k.size
		d.size <- env$d.size
		block <- env$block

		result1_1 <- get_H2_1$result1
		result2_1 <- get_H2_1$result2
		result3_1 <- get_H2_1$result3
		result4_1 <- get_H2_1$result4
		result5_1 <- get_H2_1$result5
		result6_1 <- get_H2_1$result6

##modified
		result1_2 <- get_H2_2$result1
		result2_2 <- get_H2_2$result2
#		result2_3 <- get_H2_2$result3
		result3_2 <- get_H2_2$result3

		H2_1_x <- result1_1 + result3_1[[1]] * I_1_2_x(x,result3_1[[2]],env) +
			    result6_1[[1]] * I_1_x(x,as.matrix(result6_1[[2]]),env)
           
#		H2_2_x <- I0(result2_3[[1]],env)[block]
		H2_2_x <- I0(result3_2[[1]],env)[block]
		H2_3_x <- 0

		first1_2 <- array(result1_2[[2]]$first,dim=c(d.size,d.size,k.size,k.size,block))
		second1_2 <- array(result1_2[[2]]$second,dim=c(d.size,d.size,k.size,block))
		third1_2 <- array(result1_2[[2]]$third,dim=c(d.size,d.size,block))

		first2_2 <- result2_2[[2]]$first
		second2_2 <- result2_2[[2]]$second

		H2_3coef <- get_H2_3$coef
		get_E0_bar_x <- as.matrix(E0_bar_x(x,get_H2_3$integrand,env))

		for(i in 1:d.size){
		  for(j in 1:d.size){
		    H2_1_x <- H2_1_x + result2_1[[1]][i,j] * I_1_2_x(x,result2_1[[2]][i,j][[1]],env)[block]

		    tmp <- list(first=first1_2[i,j,,,],second=third1_2[i,j,])
               
		    integrand1 <-  I_1_x(x,t(as.matrix(second1_2[i,j,,])),env) + I_12_x(x,tmp,env)           
               
		    H2_2_x <- H2_2_x + I0(result1_2[[1]][i,j,] * integrand1,env)[block]
		  }
 		  H2_1_x <- H2_1_x + result4_1[[1]][i] * I_1_2_x(x,result5_1[[2]][[i]],env)[block] +
					   result5_1[[1]][i] * I_1_x(x,t(as.matrix(result4_1[[2]][i,,])),env)[block]

		  integrand2 <- result2_2[[1]][i,] * (second2_2[i,]+I_1_x(x,t(as.matrix(first2_2[i,,])),env))

		  H2_2_x <- H2_2_x + I0(integrand2,env)[block]
            
		  H2_3_x <- H2_3_x + I0(H2_3coef[i,] * get_E0_bar_x[i,],env)[block]
		}

		result <- H2_1_x/2-H2_2_x/2-H2_3_x/2

		return(result)
	}


	ft_norm <- function(x,env){

		k.size <- env$k.size
		mu <- env$mu
		Sigma <- env$Sigma
		invSigma <- env$invSigma

		denominator <- (sqrt(2*pi))^k.size * sqrt(det(Sigma))
		numerator <- exp(-1/2 * t((x-mu)) %*% invSigma %*% (x-mu))

		result <- numerator/denominator

		return(result)
	}


#	p2 <- function(z,get_F_tilde1__2,get_F_tilde2,
#                     get_F_tilde1H1,get_H2_1,get_H2_2,get_H2_3,env){
#
#		k.size <- env$k.size
#		delta <- env$delta
#		mu <- env$mu
#		Sigma <- env$Sigma
#
#		z.tilde <- z - mu
#
#		first <- 0
#		second <- 0
#		third <- 0    #added
#
#		if(k.size == 1){
#
#		  first <- (F_tilde1__2_x(1,1,z.tilde+2*delta,get_F_tilde1__2,env) *
#				dnorm(z+2*delta,mean=mu,sd=sqrt(Sigma)) -
#				2 * F_tilde1__2_x(1,1,z.tilde+delta,get_F_tilde1__2,env) *
#				dnorm(z+delta,mean=mu,sd=sqrt(Sigma)) +
#				F_tilde1__2_x(1,1,z.tilde,get_F_tilde1__2,env) *
#				dnorm(z,mean=mu,sd=sqrt(Sigma)))/(delta)^2
#
#		  first <- first/2
#
#		  second <- (F_tilde2_x(1,z.tilde+delta,get_F_tilde2,env) *
#				 dnorm(z+delta,mean=mu,sd=sqrt(Sigma)) -
#				 F_tilde2_x(1,z.tilde,get_F_tilde2,env) *
#				 dnorm(z,mean=mu,sd=sqrt(Sigma)))/delta
#
##added:start
#		  obj1 <- F_tilde1H1_x(1,z.tilde+delta,get_F_tilde1H1,env)
#
#		  obj2 <- F_tilde1H1_x(1,z.tilde,get_F_tilde1H1,env)
#
#		  third <- (obj1 * dnorm(z+delta,mean=mu,sd=sqrt(Sigma)) -
#				obj2 * dnorm(z,mean=mu,sd=sqrt(Sigma)))/delta
#
#		  forth <- H2_x(z.tilde,get_H2_1,get_H2_2,get_H2_3,env) * #dnorm(z,mean=mu,sd=sqrt(Sigma))
#adde#d:end
#
#		}else{
#
#		  tmp1 <- matrix(0,k.size,k.size)
#
#		  for(l1 in 1:k.size){
#		    for(l2 in 1:k.size){
#
#			dif1 <- numeric(k.size)
#			dif2 <- numeric(k.size)
#			dif1[l1] <- dif1[l1] + delta
#			dif2[l2] <- dif2[l2] + delta
#
#			tmp1[l1,l2] <- (F_tilde1__2_x(l1,l2,z.tilde+dif1+dif2,get_F_tilde1__2,env) *
#					    ft_norm(z+dif1+dif2,env) -
#					    F_tilde1__2_x(l1,l2,z.tilde+dif1,get_F_tilde1__2,env) *
#					    ft_norm(z+dif1,env) -
#					    F_tilde1__2_x(l1,l2,z.tilde+dif2,get_F_tilde1__2,env) *
#					    ft_norm(z+dif2,env) +
#					    F_tilde1__2_x(l1,l2,z.tilde,get_F_tilde1__2,env) *
#					    ft_norm(z,env))/(delta)^2
#
#		    }
#		  }
#
#		  first <- sum(tmp1)/2
#
#		  tmp2 <- double(k.size)
#		  tmp3 <- double(k.size) #added
#
#		  for(l in 1:k.size){
#
#		    dif <- numeric(k.size)
#		    dif[l] <- dif[l] + delta
#
#		    tmp2[l] <- (F_tilde2_x(l,z.tilde+dif,get_F_tilde2,env) *
#				    ft_norm(z+dif,env) -
#				    F_tilde2_x(l,z.tilde,get_F_tilde2,env) *
#				    ft_norm(z,env))/delta
#
##added:start
#		    obj1 <- F_tilde1H1_x(l,z.tilde+dif,get_F_tilde1H1,env)
#
#		    obj2 <- F_tilde1H1_x(l,z.tilde,get_F_tilde1H1,env)
#
#		    tmp3[l] <- (obj1 * ft_norm(z+dif,env) -
#				    obj2 * ft_norm(z,env))/delta
##added:end
#
#		  }
#
#		  second <- sum(tmp2)
#		  third <- sum(tmp3)  #added
#		  forth <- H2_x(z.tilde,get_H2_1,get_H2_2,get_H2_3,env)*ft_norm(z,env) #added
#		}
#
#		result <- first - second -
#			    third + forth #added
#		return(result)
#	}
#
#
#	p2_z <- function(z){
#
#		if(k.size == 1){
#			zlen <- length(z)
#			result <- c()
#
#			for(i in 1:zlen){
#				result[i] <- p2(z[i],get_F_tilde1__2,get_F_tilde2,
#						    get_F_tilde1H1,get_H2_1,get_H2_2,get_H2_3,env)
#			}
#		}else{
#			result <- p2(z,get_F_tilde1__2,get_F_tilde2,
#					 get_F_tilde1H1,get_H2_1,get_H2_2,get_H2_3,env)
#		}
#
#		return(result)
#	}


#	d2.term <- function(){
#
#		gz_p2 <- function(z){
#
#			result <- H0 * G(z) * p2_z(z)
#			return(result)
#		}
#
#		if(k.size == 1){
#
#			ztmp <- seq(mu-7*sqrt(Sigma),mu+7*sqrt(Sigma),length=1000)
#			dt <- ztmp[2] - ztmp[1]
#
#			p2tmp <- gz_p2(ztmp)
#
#			result <- sum(p2tmp) * dt
#
#		}else if(2 <= k.size || k.size <= 20){
#
#			lambda <- eigen(Sigma)$values
#			matA <- eigen(Sigma)$vector
#
#			gz_p2 <- function(z){
#			  tmpz <- matA %*% z
#			  tmp <- H0 * G(tmpz) * p2_z(tmpz)	#det(matA) = 1
#			  return( tmp  )
#			}
#
#			my.x <- matrix(0,k.size,20^k.size)
#			dt <- 1
#
#			for(k in 1:k.size){
#				max <- 5 * sqrt(lambda[k])
#				min <- -5 * sqrt(lambda[k])
#				tmp.x <- seq(min,max,length=20)
#				dt <- dt * (tmp.x[2] - tmp.x[1])
#				my.x[k,] <- rep(tmp.x,each=20^(k.size-k),times=20^(k-1))
#			}
#
#			tmp <- 0
#
#			for(i in 1:20^k.size){
#				tmp <- tmp + gz_pi1(my.x[,i])
#			}
#
#			tmp <- tmp * dt
#
#		}else{
#			stop("length k is too long.")
#		}
#
#		return(result)
#	}


#