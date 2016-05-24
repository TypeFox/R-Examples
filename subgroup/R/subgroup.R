subgroup <-
function(data, overall.diff=NA, force.theoretical=FALSE, force.simulation=FALSE, design=FALSE, plots=TRUE)
{ 

#**************************************************************************
# Utility functions used in the calculations ******************************
#**************************************************************************

# Homoscedastic case --------------------------------------------------------
# Expectation of r-th order statistic ---------------------------------------
e.iid<-function(j, triangle, si.inv, overalleffect)
	{
		integrate(V.integrand.iid, -Inf, Inf, j=j, triangle=triangle, si.inv=si.inv, overalleffect=overalleffect)
	}

integrand.iid<-function(y, j, triangle, si.inv, overalleffect)
	{
		R<- length(si.inv)
		s.i<-si.inv[1]# IID case, so all si.inv are the same
		total<-0
		power<-c((j-1):(R-1))# Powering of the cumulative distribution function and the derivative
		derivative.component<-c(j:R)
		product<- s.i*dnorm(s.i*(y-overalleffect))*((pnorm(s.i*(y-overalleffect)))^power)
		product<-triangle*product*derivative.component
		total<-sum(product)
		y*total
	}

V.integrand.iid<-Vectorize(integrand.iid, vectorize.args="y")

# Function to work out the density of the range (Max - Min) ------------------
# Integrand of the marginal density of the range. ----------------------------
# Here x.min signifies the smallest order statistic, r the values of the range 
# and thus, x.min+r=the largest order statistic ------------------------------
integrand.marginal.density.iid<-function(x, r, overall, si) 
	{
		R*(R-1)*dnorm((x+r), overall, si[1])*dnorm(x, overall, si[1])*((pnorm((x+r),overall,si[1]) - pnorm(x,overall,si[1]))^(R-2))
	}

V.range.density.iid<-Vectorize(integrand.marginal.density.iid,vectorize.args=c("x"))

# Heteroscedastic case ------------------------------------------------------
# Theory based --------------------------------------------------------------
# ---------------------------------------------------------------------------
# Function that works out the expectation of all order statistics -----------
# Use the expectation of the last order statistic (routine above) -----------
# to work out the expectations of all the order statistics ------------------

# Expectation of largest  order statistic emax ------------------------------
# This is used repeatedly in the calculation of the expectation of ----------
# the r-th order statistic --------------------------------------------------
emaxR<-function(si.inv, overalleffect)
	{
		integrate(V.integrand.max,-Inf,Inf,si.inv, overalleffect)
	}
integrand.max<-function(y, si.inv, overalleffect)
	{
		R.temp<- length(si.inv)
		cum.mat<-matrix(pnorm(si.inv*(y-overalleffect)),R.temp,R.temp)
		diag(cum.mat)<-si.inv*dnorm(si.inv*(y-overalleffect))
		product<- apply(cum.mat,2,prod)
		y*sum(product)
	}
V.integrand.max<-Vectorize(integrand.max,vectorize.args="y")

expectation.heteroscedastic<-function(data, emaxR, si.inv, overalleffect) 
	{
		integral<-numeric() 
		for(j in c((ceiling(R/2)+1):R))
			{
			comb<-combn(R,j)
			integral.temp<-numeric()

			for(c in c(1:dim(comb)[2]))
				{
				sub.si.inv<-si.inv[comb[,c]]
				j.expectation<-emaxR(si.inv=sub.si.inv,overalleffect=overalleffect)[1] 
				integral.temp<-rbind(integral.temp, as.numeric(j.expectation))
				}
			integral.temp<-sum(integral.temp)
			integral<-rbind(integral, cbind(R, j, integral.temp))
			}

		# Work out the expectations using the integrals and the triangle matrix 
		expectation<-numeric()
		for(j in c(1:dim(integral)[1])) # Create triangular numbers
			{
			c.ir<-c((ceiling(R/2)+1):R)
			triangle<-as.matrix((-1)^(c(0:(dim(integral)[1]-j)))*choose((c.ir[j:length(c.ir)]-1), (c.ir[j]-1)), ncol=1)
			temp.integral<-matrix(integral[(j:dim(integral)[1]),],ncol=3)
			temp.integral<- temp.integral[,3]*triangle
			temp.exp<-sum(temp.integral)
			expectation<-rbind(expectation,as.numeric(temp.exp))
			}
		# Only the lower 1/2 of the expectations of the order statistics are processed.
		# So find the others using the fact that the expecations of the opposing order statistics add up to 2.

		lowerhalf<-(2*overalleffect-expectation[,1])
		# Adjustment for odd numbered subgroups 
		if((R%%2)==1){ lowerhalf=c(overalleffect, lowerhalf) } 
		lowerhalf<-matrix(lowerhalf, ncol=1)
		lowerhalf<-matrix(lowerhalf[order(lowerhalf),], ncol=1)
		fullexpectation<-rbind(lowerhalf, expectation)
		output<-cbind("Observed difference"=data[,1], "Observed SE of difference"=data[,2])
		output<-output[order(output[,1]), ] 
		output<-cbind(output,"Expected difference"=fullexpectation[,1], Order=seq(1:R))
		output
	}

# Function to work out the probability of number of subgroups that could have ---
# a treatment effect in favour of control by chance alone. ----------------------
greater.0 <-function(prob.G.0)
	{
		probG0.sum<-numeric()
		probG0.sum<-rbind(probG0.sum, cbind(0, prod(1-prob.G.0)))

		for(j in 1:R)
		{
			prob.comb<-numeric()
			probG0.comb<-combn(prob.G.0,j)
			for(i in 1:choose(R,j))
				{
				temp<-prob.G.0[is.na(pmatch(prob.G.0, probG0.comb[,i]))]  
				# These are the elements that are selected as not in favour of the control
				probL0<- (1-temp)  
				# Find the probability of not favouring control
				prob.comb<- rbind(prob.comb, (prod(probL0)*prod(probG0.comb[,i])))
				# Multiply the combinations in favour of control with those not in favour of control
				}
			probG0.sum<-rbind(probG0.sum, cbind("X, Number of subgroups favouring control=x"=j, "P(X=x)"=sum(prob.comb)))
			# Add the probability of all combinations of j out of R that favour control
		}
		probG0.sum
	}

# Function to work out the density of the range (Max - Min) -----------------------
# Integrand of the marginal density of the range. ---------------------------------
# Here x.min signifies the smallest order statistic, r the values of the range ----
# and thus, x.min+r=the largest order statistic -----------------------------------
integrand.marginal.density<-function(x, r, overall, si.rest, si.min.max.1, si.min.max.2) 
	{
		if(length(si.rest)>0)
			{
			prod.pnorm.rest<- apply(matrix((pnorm((x+r),overall,si.rest) - pnorm(x,overall,si.rest)), nrow=dim(combn(R,2))[2]),1,prod)
      		sum(
			(dnorm(x, overall, si.min.max.1)*dnorm((x+r), overall, si.min.max.2)  
			+ dnorm(x, overall, si.min.max.2)*dnorm((x+r), overall, si.min.max.1))*prod.pnorm.rest
			 )
			}
		else if(length(si.rest)==0)
			{
     			sum(
			(dnorm(x, overall, si.min.max.1)*dnorm((x+r), overall, si.min.max.2)  
			+ dnorm(x, overall, si.min.max.2)*dnorm((x+r), overall, si.min.max.1))
			 )
			}
	}
V.range.density<-Vectorize(integrand.marginal.density,vectorize.args=c("x"))

#**************************************************************************
# Function used in producing theoretical outputs in homoscedastic case ****
#**************************************************************************
output.homoscedastic<-function(data, e.iid, V.range.density.iid, overalleffect, si.inv, R, plots)
	{
	# Homoscedastic case: the standard errors are all the same ----------------

	# Expectation of the order statistics -------------------------------------
	expectation<-numeric() 
	for(j in c(1:R))
		{
		c.ir<-c(j:R)
		triangle<-as.matrix((-1)^(c.ir - c.ir[1])*choose((c.ir-1) ,(c.ir[1]-1))*choose(R, c.ir), ncol=1)
		j.expectation<-e.iid(j, triangle, si.inv, overalleffect)[1]
		expectation<-rbind(expectation, as.numeric(j.expectation))
		}
		expectation<-matrix(expectation[,1], ncol=1)
		order.stats.expectations<-cbind("Observed difference"=data[,1], "Observed SE of difference"=data[,2])
		order.stats.expectations<-order.stats.expectations[order(order.stats.expectations[,1]), ] 
		order.stats.expectations<-cbind(order.stats.expectations, "Expected difference"=expectation[,1], "Order"=seq(1:R))

		# Number of subgroups favouring control -----------------------------------
		standard.stat<-(0 - overall)/data[1,2]
		probG0<- (1-pnorm(standard.stat))
		probG0<- cbind("X, Number of subgroups favouring control=x"=c(0:R), "P(X=x)"=dbinom(0:R, R, probG0, log = FALSE))

	# Find the density function of the range. ---------------------------------
	# Integrate the marginal density over -Inf to Inf (the values that x.min can take)
	# for r in 0 to Inf as the range can only be positive.  The values of r can then be
	# plotted against the marginal density to dertermine the density function of the range.
	range.density<-numeric()
	epsilon<-1# Measures how much further to go in the density estimation, cumulative prob=1
	epsilon.count<-0# Extra parameters defined to help break the loop as required
	epsilon.old<-0
	cum.prob<-0# Measures how much of the cumulative probability has been accounted for 
	last.density<-0# Last density mass for calculating the AUC
	r=0
	exp.range<-max(order.stats.expectations[,3])-min(order.stats.expectations[,3]) # Expected range
	step=exp.range/100# Calculate probability density at these points
	observed.range<-max(data[,1])-min(data[,1]) # Observed range
	while(epsilon>0.000001 & epsilon.count<=10)# Cumulative probabilty of range distribution is nearly=1
		{ 
		current.density<-as.numeric(integrate(V.range.density.iid, lower=-Inf, upper=Inf, r=r, overall, si)[1])
		range.density<-rbind(range.density, cbind("Range=r"=r, "f(r)"=current.density))
		cum.prob<-cum.prob+(last.density+current.density)*step/2
		# Reset parameters for next iteration
		if(r==0) 
			{ 
			epsilon.old<-0 
			epsilon.count<-1  
			}
		else 
			{
			if(epsilon.old==epsilon) { epsilon.count=epsilon.count+1 }
			else { epsilon.count=0 }
			epsilon.old<-epsilon 
			} 
		r=r+step
		last.density<-current.density
		epsilon<-1-cum.prob
	}

	# Plots
	if(plots==TRUE)
		{
		par(mfrow=c(2,2), oma=c(0,0,3,0))  # Leave space in outer margin into which the main overall title will fit
		plots.subgroup(data=data, overall=overall, R=R, order.stats.expectations=order.stats.expectations, range.density=range.density, probG0=probG0, design=design)
		par(mfrow=c(1,1))
		}

	# Return values
	results<-list("expectations"=order.stats.expectations, "favourcontrol"=probG0, "rangedensity"=range.density, "overalldiff"=overall)
	return(results)
	} 

#**************************************************************************
# Function used in producing theoretical outputs in heteroscedastic case **
#**************************************************************************
output.heteroscedastic<-function(data, emaxR, V.range.density, overalleffect, si.inv, R, plots)
	{
	# Expectation of the order statistics -------------------------------------
	order.stats.expectations <- expectation.heteroscedastic(data=data, emaxR=emaxR, si.inv=si.inv, overalleffect=overall)

	# Number of subgroups favouring control -----------------------------------
	standard.stat<-(0 - overall)/data[,2]
	probG0<- (1-pnorm(standard.stat))
	probG0<- greater.0(prob.G.0=probG0)

	# Find the density function of the range. ------------------------------------
	# Integrate the marginal density over -Inf to Inf (the values that x.min can take)
	# for r in 0 to Inf as the range can only be positive.  The values of r can then be
	# plotted against the marginal density to dertermine the density function of the range.
	comb.min.max<-combn(R, 2)
	range.density<-numeric()
	epsilon<-1# Measures how much further to go in the density estimation, cumulative prob=1
	epsilon.count<-0# Extra parameters defined to help break the loop as required
	epsilon.old<-0
	cum.prob<-0# Measures how much of the cumulative probability has been accounted for 
	last.density<-0# Last density mass for calculating the AUC
	r=0
	exp.range<-max(order.stats.expectations[,3])-min(order.stats.expectations[,3]) # Expected range
	step=exp.range/100# Calculate probability density at these points
	observed.range<-max(data[,1])-min(data[,1]) # Observed range
	while(epsilon>0.000001 & epsilon.count<=10)# Cumulative probabilty of range distribution is nearly=1
		{ 
		si.mat<-numeric()
		for (i in c(1:dim(comb.min.max)[2]))
			{
			si.min.max<-si[comb.min.max[,i]]
			si.rest<-subset(c(1:R), c(1:R)!=as.vector(comb.min.max[1,i]) & c(1:R)!=as.vector(comb.min.max[2,i]), drop=TRUE)  
			si.rest<-si[si.rest]
			si.mat<-rbind(si.mat, matrix(c(i, c(si.min.max, si.rest)),nrow=1)) 
			}
		si.min.max.1<-si.mat[,2]
		si.min.max.2<-si.mat[,3]
		if(R>2){ si.rest<-si.mat[,4:(R+1)] }# More than 2 subgroups so si.mat is populated
			else {si.rest<-numeric() }# Only 2 subgroups so si.mat is undefined
		current.density<-as.numeric(integrate(V.range.density, lower=-Inf, upper=Inf, r=r, overall, si.rest, si.min.max.1, si.min.max.2)[1])
		range.density<-rbind(range.density, cbind("Range=r"=r, "f(r)"=current.density))
		cum.prob<-cum.prob+(last.density+current.density)*step/2
		# Reset parameters for next iteration
		if(r==0) 
			{ 
			epsilon.old<-0 
			epsilon.count<-1  
			}
		else 
			{
			if(epsilon.old==epsilon) { epsilon.count=epsilon.count+1 }
				else { epsilon.count=0 }
			epsilon.old<-epsilon 
			} 
		r=r+step
		last.density<-current.density
		epsilon<-1-cum.prob
		}

	# Plots
	if(plots==TRUE)
		{
		par(mfrow=c(2,2), oma=c(0,0,3,0))  # Leave space in outer margin to stick main title into
		plots.subgroup(data=data, overall=overall, R=R, order.stats.expectations=order.stats.expectations, range.density=range.density, probG0=probG0, design=design)
		par(mfrow=c(1,1))
		}

	# Return values
	results<-list("expectations"=order.stats.expectations, "favourcontrol"=probG0, "rangedensity"=range.density, "overalldiff"=overall)
	return(results)
	}

#**************************************************************************
# Function used in producing simulation outputs ***************************
#**************************************************************************
output.simulation<-function(data, overall, si.inv, R, plots)
	{
	# Simulate data (how many N per subgroup) ------------------------------
	# using each of the 2 methods of SE ------------------------------------
	nsim<-10000
	probG0<-numeric()
	range.density<-numeric()
	expected<-0
	samp.i.all<-NULL
	for (i in 1:nsim) 
		{
		samp.i<-matrix(sort(rnorm(R, overall, si))  , ncol=1)
		samp.i.all<-cbind(samp.i.all, samp.i)
		probG0<-rbind(probG0, cbind(length(samp.i[samp.i>0]), length(samp.i[samp.i>0])))
		range.density<-rbind(range.density, (max(samp.i) - min(samp.i)))
		}
	expected<-sort(rowMeans(samp.i.all))
	order.stats.expectations<-cbind("Observed difference"=data[,1], "Observed SE of difference"=data[,2])
	order.stats.expectations<-order.stats.expectations[order(order.stats.expectations[,1]), ] 
	order.stats.expectations<-cbind(order.stats.expectations,"Expected difference"=expected, Order=seq(1:R))

	# Number of subgroups favouring control -----------------------------------
	probG0<-as.data.frame(table(probG0[,2]))
	probG0<-cbind('X, Number of subgroups favouring control=x'=as.numeric(levels(probG0[,1])), 'P(X=x)'=probG0[,2])
	probG0[,2]<-probG0[,2]/nsim

	# Find the density function of the range. ---------------------------------
	range.density<-cbind(density(range.density[,1], from=0)$x, density(range.density[,1], from=0)$y)
	# Plots
	if(plots==TRUE)
		{
		par(mfrow=c(2,2), oma=c(0,0,3,0))  # Leave space in outer margin to stick main title into
		plots.subgroup(data=data, overall=overall, R=R, order.stats.expectations=order.stats.expectations, range.density=range.density, probG0=probG0, design=design)
		par(mfrow=c(1,1))
		}

	# Return values
	results<-list("expectations"=order.stats.expectations, "favourcontrol"=probG0, "rangedensity"=range.density, "overalldiff"=overall)
	return(results)
	}


#**************************************************************************
# Function used in plotting ***********************************************
#**************************************************************************
plots.subgroup<-function(data, overall, R, order.stats.expectations, range.density, probG0, design)
	{
	if(design==FALSE)# Observed summaries also presented in the plots at the analysis stage
		{
		boxplot(order.stats.expectations[,3], order.stats.expectations[,1], cex.axis=0.8, cex.lab=0.8, cex.main=0.8, range=0, names=c("Expected", "Observed"), ylab="Treatment difference", main="A")
		abline(h=overall, lty=2, col="red")

		lim <- ceiling(max(abs(order.stats.expectations[,3]), abs(data[,1])))
		plot(x=order.stats.expectations[,3], y=order.stats.expectations[,1], type="p", pch=20, cex.axis=0.8, cex.lab=0.8, cex.main=0.8, xlab="Expected treatment difference",ylab="Observed treatment difference", xlim=c(-lim, lim), ylim=c(-lim, lim), main="B")
		abline(a=0, b=1)
		abline(v=overall, lty=2, col="red")

		observed.range<-max(data[,1]) - min(data[,1])
		sub.range<-subset(range.density, range.density[, 1]>=observed.range, drop=TRUE)
		p.range<-0    # Contains the probability to be printed on the plot
		if(length(sub.range)==0) { p.range<-0 } # The observed value is outside the density of the range 
			else 
				{    # Calculate the probability that something more extreme than that which was observed could be observed
					for (i in (2:dim(sub.range)[1]))
					{
					x<-as.numeric(sub.range[i,2])+as.numeric(sub.range[(i-1),2])
					x<-x*(as.numeric(sub.range[i,1]) - as.numeric(sub.range[(i-1),1]))/2
					p.range<-p.range+x
					}
				}
		p.range<-format(round(p.range, 2), nsmall = 2)
		max.x.range<-max(c(range.density[dim(range.density)[1], 1], observed.range))
		plot(x=range.density[,1],y=range.density[,2], type="l", cex.axis=0.8, cex.lab=0.8, cex.main=0.8, xlab="Range of treatment differences", ylab="Probability density", xlim=c(0, max.x.range), main="C")
		abline(v=observed.range, lty=2, col="red")
		legend("topright", as.expression(c(bquote(paste(P[E], " = ",.(p.range))))))
		observed.G.0 <- length(subset(data[,1], data[,1]>0, drop=TRUE))

		p.G.0<-subset(probG0, probG0[,1]>=observed.G.0)
		p.G.0<-format(round(sum(p.G.0[,2]), 2), nsmall = 2) 

		plot(x=probG0[,1], y=probG0[,2], type="h", lwd=3, cex.axis=0.8, cex.lab=0.8, cex.main=0.8, xlab="Number of subgroups favouring control", ylab="Probability", xlim=c(0,(R+2)), main="D")
		abline(v=observed.G.0, lty=2, col="red")
		legend("topright", as.expression(c(bquote(paste(P[E], " = ",.(p.G.0))))))

		mtext(expression(bold("Observed and expected")), adj=0.5, side=3, outer=TRUE, cex=0.8) # The main title that is going to go above the 4 plots
		mtext(expression("Dotted line indicates observed overall treatment difference (A & B), range (C) and number of subgroups favouring control (D)"), padj=0.9, side=3, outer=TRUE, cex=0.6) # The main title that is going to go above the 4 plots
		mtext(expression(paste(P[E], " is the probability of obtaining an observation at least as extreme as that which was observed")), padj=1.9, side=3, outer=TRUE, cex=0.6) # The main title that is going to go above the 4 plots
		}
	else if(design==TRUE)
		{
		boxplot(order.stats.expectations[,3], cex.axis=0.8, cex.lab=0.8, cex.main=0.8, range=0, xlab="Expected", ylab="Treatment difference", main="A")
		plot(order.stats.expectations[,3], y=c(rep(0.5,dim(data)[1])), type="p", pch=124, cex.axis=0.8, cex.lab=0.8, cex.main=0.8, xlab="Expected treatment difference", ylab=" ", yaxt='n', main="B")
		abline(h=0.5, lty=1, col = 1)

		max.x.range<-max(range.density[, 1])
		plot(x=range.density[,1],y=range.density[,2], type="l", cex.axis=0.8, cex.lab=0.8, cex.main=0.8, xlab="Range of treatment differences", ylab="Probability density", xlim=c(0, max.x.range), main="C")
		plot(x=probG0[,1], y=probG0[,2], type="h", lwd=3, cex.axis=0.8, cex.lab=0.8, cex.main=0.8, xlab="Number of subgroups favouring control", ylab="Probability", xlim=c(0,(R+2)), main="D")

		mtext(expression(bold("Expected at design stage")), adj=0.5, side=3, outer=TRUE, cex=0.8) # The main title that is going to go above the 4 plots
		}
	} # End plots.subgroup


#**************************************************************************
# MAIN ROUTINE ************************************************************
#**************************************************************************
# Check that data that has been received is in correct format *************
#**************************************************************************
if(dim(data)[1]<2 || is.null(dim(data)[1]) || dim(data)[2]!=2 || is.null(dim(data)[2])) 
	{ 
	print("ERROR*****************************************")
	print("The data submitted must be a matrix of exactly 2 columns;  column 1=treatment differences, column 2=standard errors.")
	print("There should be at least 2 subgroups.  That is, the matrix should have at least 2 rows.")
	stop
	}
else if(min(data[,2]<=0))
	{
	print("ERROR*****************************************")
	print("Standard errors in the 2nd column of the data matrix must be > 0.")
	stop
	}
else
	{
	#**************************************************************************
	# Perform the calculations for the input data *****************************
	#**************************************************************************

	# Define values used in the calculations ----------------------------------
	R<-dim(data)[1]
	homoscedastic<-as.data.frame(table(data[,2]))

	si<-data[,2]# The standard errors
	si.inv<- si^(-1)# Inverse of the standard errors
	weight  <- 1/(si^2)# The weights

	# Set value of overall estimate if provided --------------------------------
	{if(exists("overall.diff")==TRUE & is.na(overall.diff)==FALSE) { overall <- overall.diff  } # The overall treatment difference
		else if(exists("overall.diff")==FALSE || is.na(overall.diff)==TRUE)
		{ overall <- sum(data[,1]*weight)/sum(weight) }}

	# Set default values and send messages to screen as required --------------
	if(dim(homoscedastic)[1]==1 & homoscedastic[1,2]==R) # Homoscedastic input
		{ 
		theoretical<-TRUE 
		simulation<-FALSE
 		if(force.simulation==TRUE & (exists("force.theoretical")==FALSE || force.theoretical==FALSE))
		# Simulation output is being requested even though theoretical output can be produced easily
			{
			theoretical<-FALSE
			simulation<-TRUE
			print("NOTE: As requested simulation based results are produced.")
			output.simulation(data=data, overall=overall, si.inv=si.inv, R=R, plots=plots)
			}
		else 
			{
			if(force.simulation==TRUE & force.theoretical==TRUE)
				{
				# Both theoretical and simulation outputs were requested.  Print error and break the loop. 
				print("ERROR: Both theory based and simulation based results requested. Please choose only one option.")
				stop()
				}
			output.homoscedastic(data=data, e.iid=e.iid, V.range.density.iid=V.range.density.iid, overalleffect=overall, si.inv=si.inv, R=R, plots=plots)
			}
		} # End homoscedastic 
	else if(R>20) 
		{ 
		theoretical<-FALSE
		simulation<-TRUE 
  		if(force.theoretical==TRUE & (exists("force.simulation")==FALSE || force.simulation==FALSE)) 
			{
			simulation<-FALSE
			theoretical<-TRUE 
			print("WARNING: Number of subgroups exceeds 20. Execution time for production of theory based outputs will be substantial.")
			output.heteroscedastic(data=data, emaxR=emaxR, V.range.density=V.range.density, overalleffect=overall, si.inv=si.inv, R=R, plots=plots)
			}
		else 
			{
			if(force.simulation==TRUE & force.theoretical==TRUE)
				{
				# Both theoretical and simulation outputs were requested.  Print error and break the loop. 
				print("ERROR: Both theory based and simulation based results requested. Please choose only one option.")
				stop()
				}
			else
				{ 
				print("NOTE: Only simulation based results are produced as the number of subgroups exceeds 20.") 
				output.simulation(data=data, overall=overall, si.inv=si.inv, R=R, plots=plots)
				}
			}
		} # End heteroscedastic R>20 
	else if(R<21) 
		{ 
		theoretical<-TRUE 
		simulation<-FALSE
		if(force.simulation==TRUE & (exists("force.theoretical")==FALSE || force.theoretical==FALSE))
		# Simulation output is being requested even though theoretical output can be produced easily
			{
			theoretical<-FALSE
			simulation<-TRUE
			print("NOTE: As requested simulation based results are produced.")
			output.simulation(data=data, overall=overall, si.inv=si.inv, R=R, plots=plots)
			}
		else 
			{
			if(force.simulation==TRUE & force.theoretical==TRUE)
				{
				# Both theoretical and simulation outputs were requested.  Print error and break the loop. 
				print("ERROR: Both theory based and simulation based results requested. Please choose only one option.")
				stop()
				}
			output.heteroscedastic(data=data, emaxR=emaxR, V.range.density=V.range.density, overalleffect=overall, si.inv=si.inv, R=R, plots=plots)
			}
		} # End heteroscedastic R<=20 
	} # End execution of valid data
}
