"finn" <-
function(ratings, s.levels, model = c("oneway", "twoway")) {
  model <- match.arg(model)
	ratings <- as.matrix(na.omit(ratings))

	ns <- nrow(ratings)
	nr <- ncol(ratings)

  SStotal <- var(as.numeric(ratings))*(ns*nr-1)
  MSr <- var(apply(ratings,1,mean))*nr
	MSw <- sum(apply(ratings,1,var)/ns)
	MSc <- var(apply(ratings,2,mean))*ns
	MSe <- (SStotal-MSr*(ns-1)-MSc*(nr-1))/((ns-1)*(nr-1))

	MSexp <- 1/12*(s.levels^2-1)

	method <- paste("Finn-Coefficient (Model=",model,")",sep="")

	if (model=="oneway") {
		#Asendorpf & Wallbott, S. 245, Fu
		#Finn (1970)
		coeff  <- 1-(MSw/MSexp)
		Fvalue <- MSexp/MSw
	  df1    <- Inf
	  df2    <- ns*(nr-1)
	  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)
	}
	else {
		#Asendorpf & Wallbott, S. 246, Fa
		coeff  <- 1-(MSe/MSexp)
		Fvalue <- MSexp/MSe
	  df1    <- Inf
	  df2    <- ns*(nr-1)
	  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)
	}

  rval <- structure(list(method = method,
                         subjects = ns, raters = nr,
                         irr.name = "Finn", value = coeff,
                         stat.name = paste("F(Inf,",df2,")",sep=""), statistic = Fvalue, p.value = p.value),
                    class="irrlist")

  return(rval)
}

