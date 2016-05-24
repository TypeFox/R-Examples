"icc" <-
function(ratings, model = c("oneway", "twoway"), type = c("consistency", "agreement"),
         unit = c("single", "average"), r0 = 0, conf.level = .95) {
	ratings <- as.matrix(na.omit(ratings))
  model <- match.arg(model)
  type  <- match.arg(type)
  unit  <- match.arg(unit)
  alpha <- 1-conf.level

  ns <- nrow(ratings)
	nr <- ncol(ratings)

  SStotal <- var(as.numeric(ratings))*(ns*nr-1)
  MSr <- var(apply(ratings,1,mean))*nr
	MSw <- sum(apply(ratings,1,var)/ns)
	MSc <- var(apply(ratings,2,mean))*ns
	MSe <- (SStotal-MSr*(ns-1)-MSc*(nr-1))/((ns-1)*(nr-1))

	#Single Score ICCs
	if (unit == "single") {
	    if (model=="oneway") {
			#Asendorpf & Wallbott, S. 245, ICu
			#Bartko (1966), [3]

		    icc.name <- "ICC(1)"
		    coeff  <- (MSr-MSw)/(MSr+(nr-1)*MSw)

		    Fvalue <- MSr/MSw*((1-r0)/(1+(nr-1)*r0))
		    df1    <- ns-1
		    df2    <- ns*(nr-1)
		    p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)

		    #confidence interval
		    FL <- (MSr/MSw)/qf(1-alpha/2, ns-1, ns*(nr-1))
		    FU <- (MSr/MSw)*qf(1-alpha/2, ns*(nr-1), ns-1)
		    lbound <- (FL-1)/(FL+(nr-1))
		    ubound <- (FU-1)/(FU+(nr-1))
	    }
	    else if (model=="twoway") {
			if (type == "consistency") {
		    #Asendorpf & Wallbott, S. 245, ICa
				#Bartko (1966), [21]
				#Shrout & Fleiss (1979), ICC(3,1)

		    icc.name <- "ICC(C,1)"
				coeff  <- (MSr-MSe)/(MSr+(nr-1)*MSe)

			  Fvalue <- MSr/MSe*((1-r0)/(1+(nr-1)*r0))
			  df1    <- ns-1
			  df2    <- (ns-1)*(nr-1)
			  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)

			  #confidence interval
			  FL <- (MSr/MSe)/qf(1-alpha/2, ns-1, (ns-1)*(nr-1))
			  FU <- (MSr/MSe)*qf(1-alpha/2, (ns-1)*(nr-1), ns-1)
			  lbound <- (FL-1)/(FL+(nr-1))
			  ubound <- (FU-1)/(FU+(nr-1))
			}
			else if (type == "agreement") {
				#Asendorpf & Wallbott, S. 246, ICa'
				#Bartko (1966), [15]
				#Shrout & Fleiss (1979), ICC(2,1)

				icc.name <- "ICC(A,1)"
				coeff  <- (MSr-MSe)/(MSr+(nr-1)*MSe+(nr/ns)*(MSc-MSe))

				a <- (nr*r0)/(ns*(1-r0))
				b <- 1+(nr*r0*(ns-1))/(ns*(1-r0))
				#v <- (a*MSc+b*MSe)^2/((a*MSc)^2/(nr-1)+(b*MSe)^2/((ns-1)*(nr-1)))

				Fvalue <- MSr/(a*MSc+b*MSe)

				a <- (nr*coeff)/(ns*(1-coeff))
				b <- 1+(nr*coeff*(ns-1))/(ns*(1-coeff))
				v <- (a*MSc+b*MSe)^2/((a*MSc)^2/(nr-1)+(b*MSe)^2/((ns-1)*(nr-1)))

			  df1     <- ns-1
			  df2     <- v
			  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)

			  #confidence interval (McGraw & Wong, 1996)
				FL <- qf(1-alpha/2, ns-1, v)
			  FU <- qf(1-alpha/2, v, ns-1)
			  lbound <- (ns*(MSr-FL*MSe))/(FL*(nr*MSc+(nr*ns-nr-ns)*MSe)+ns*MSr)
			  ubound <- (ns*(FU*MSr-MSe))/(nr*MSc+(nr*ns-nr-ns)*MSe+ns*FU*MSr)
			}
	  }
	}
	#Average Score ICCs
	else if (unit == "average") {
	    if (model=="oneway") {
			  #Asendorpf & Wallbott, S. 245, Ru
		    icc.name <- paste("ICC(",nr,")",sep="")
		    coeff  <- (MSr-MSw)/MSr

		    Fvalue <- MSr/MSw*(1-r0)
		    df1    <- ns-1
		    df2    <- ns*(nr-1)
		    p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)

		    #confidence interval
		    FL <- (MSr/MSw)/qf(1-alpha/2, ns-1, ns*(nr-1))
		    FU <- (MSr/MSw)*qf(1-alpha/2, ns*(nr-1), ns-1)
		    lbound <- 1-1/FL
		    ubound <- 1-1/FU
	    }
	    else if (model=="twoway") {
			if (type == "consistency") {
				#Asendorpf & Wallbott, S. 246, Ra
		    icc.name <- paste("ICC(C,",nr,")",sep="")
				coeff  <- (MSr-MSe)/MSr

			  Fvalue <- MSr/MSe*(1-r0)
			  df1    <- ns-1
			  df2    <- (ns-1)*(nr-1)
			  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)

			  #confidence interval
			  FL <- (MSr/MSe)/qf(1-alpha/2, ns-1, (ns-1)*(nr-1))
			  FU <- (MSr/MSe)*qf(1-alpha/2, (ns-1)*(nr-1), ns-1)
			  lbound <- 1-1/FL
			  ubound <- 1-1/FU
			}
			else if (type == "agreement") {
		    icc.name <- paste("ICC(A,",nr,")",sep="")
				coeff  <- (MSr-MSe)/(MSr+(MSc-MSe)/ns)

				a <- r0/(ns*(1-r0))
				b <- 1+(r0*(ns-1))/(ns*(1-r0))
				#v <- (a*MSc+b*MSe)^2/((a*MSc)^2/(nr-1)+(b*MSe)^2/((ns-1)*(nr-1)))

				Fvalue <- MSr/(a*MSc+b*MSe)

				a <- (nr*coeff)/(ns*(1-coeff))
				b <- 1+(nr*coeff*(ns-1))/(ns*(1-coeff))
				v <- (a*MSc+b*MSe)^2/((a*MSc)^2/(nr-1)+(b*MSe)^2/((ns-1)*(nr-1)))

			  df1    <- ns-1
			  df2    <- v
			  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)

			  #confidence interval (McGraw & Wong, 1996)
				FL <- qf(1-alpha/2, ns-1, v)
			  FU <- qf(1-alpha/2, v, ns-1)
			  lbound <- (ns*(MSr-FL*MSe))/(FL*(MSc-MSe)+ns*MSr)
			  ubound <- (ns*(FU*MSr-MSe))/(MSc-MSe+ns*FU*MSr)
			}
	  }
	}

  rval <- structure(list(subjects = ns, raters = nr,
                         model = model, type = type, unit = unit,
                         icc.name = icc.name, value = coeff,
                         r0 = r0, Fvalue = Fvalue, df1 = df1, df2 = df2, p.value = p.value,
                         conf.level = conf.level, lbound = lbound, ubound = ubound),
                    class="icclist")
  return(rval)
}

