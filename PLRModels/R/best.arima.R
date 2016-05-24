
#######################################################################################################################
# It gives us the orders p, q, P and Q of the best ARIMA (p, d, q)x(P, D, Q)_s according to one of the following criteria: AIC, AICC, BIC.
# It works with ARIMAs with constant term and d+D!=0, besides the ordinary cases.
#
# Arguments:
# x ==> vector or object of the class "time series", to which we want to fit an ARIMA
# order.max ==> vector of length 3: 
#                                  		order.max[1] ==> max.p
#                                  		order.max[2] ==> d
#                                  		order.max[3] ==> max.q
# seasonal$order.max ==> vector of length 3: 
#                                  			[1] ==> max.P
#                                  			[2] ==> D
#                                 			[3] ==> max.Q
#
# seasonal$period ==> period of the seasonal component
#
# include.mean ==> in reference to the inclusion or not of the mean/constant in the ARIMA. By default, TRUE. If d+D != 0, it does not allow mean/constant
# criterio ==> criterion to be chosen in the selection of the model: "AIC", "AICC" or "BIC". By default, BIC.
# dist.max.crit ==> the function will give the orders of all the ARMA models whose criterium function has a value which differs from the minimum at most dist.max.crit units 
#                   By default, 2.
# method ==> estimation method: it should be CSS-ML or ML, although it also allows CSS. By default, CSS-ML.
# Salida: orders and values of the criterion function for the selected ARMAs.
# 
########################################################################################################################

best.arima <- function(x=x, order.max=c(0,0,0), seasonal=list(order.max=c(0,0,0), period=1), include.mean=NULL, criterio=NULL, dist.max.crit=NULL, method=NULL)
{

if (is.null(include.mean)) include.mean <- TRUE

if (is.null(criterio)) criterio <- "BIC"

if (is.null(dist.max.crit)) dist.max.crit <- 2

p.max <- order.max[1]
d <- order.max[2]
q.max <- order.max[3]

P.max <- seasonal$order.max[1]
D <- seasonal$order.max[2]
Q.max <- seasonal$order.max[3]
if (is.ts(x)) period <- frequency(x) else period <- seasonal$period


num.x.perdidos <- d + period*D
T <- length(x) - num.x.perdidos


# We create a matrix which will contain all the combinations of the considered orders and the values of the criterium function for the corresponding ARIMAS
VALORES.CRITERIO <- matrix(0,(p.max+1)*(q.max+1)*(P.max+1)*(Q.max+1), 5)


# The penalty which, relating to the quantity of parameters, impose the criterium functions AIC or BIC depends on the factor defined below
if (criterio=="AIC")  factor <- 2
   else if (criterio=="BIC")  factor <- log(T)

fila <- 0

for (p in 0:p.max) 
for (q in 0:q.max)
for (P in 0:P.max)
for (Q in 0:Q.max)

   {
#optim.control=list(maxit=500)

	fila <- fila +1

	ajuste <- try(arima(x=x, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=period), include.mean=include.mean, method=method), silent=TRUE)

	if (class(ajuste)=="try-error") {
					VALORES.CRITERIO[fila, ] <- c(p, q, P, Q, NaN)	
 					next
						}


# The penalty which, relating to the quantity of parameters, impose the criterium function AICC depends on the factor defined below	
	if (criterio=="AICC")  factor <- 2*T/(T-length(ajuste$coef)-2)

	criterio.ajuste <- -2*ajuste$loglik + factor*(length(ajuste$coef)+1)

	VALORES.CRITERIO[fila, ] <- c(p, q, P, Q, criterio.ajuste)

   }


# We obtain the distance between each value of the criterium function and its minimum value
DISTANCIAS.AL.MINIMO <- VALORES.CRITERIO[,5] - min(VALORES.CRITERIO[,5], na.rm=TRUE)


# We create a vector with the values TRUE or FALSE depending on whether we want that the corresponding row to each corresponding vector is shown in the output
FILAS.OK <- rep(FALSE,length=(p.max+1)*(q.max+1)*(P.max+1)*(Q.max+1))

FILAS.OK[ DISTANCIAS.AL.MINIMO<= dist.max.crit] <- TRUE


# We just load the results (orders and value of the criterium function) that we want to show
VALORES.CRITERIO.OK <- VALORES.CRITERIO[FILAS.OK,]

if (!is.matrix(VALORES.CRITERIO.OK))  VALORES.CRITERIO.OK <- t(as.matrix(VALORES.CRITERIO.OK))


# We sort the VALORES.CRITERIO.OK from higher to lower according to the criterium function
DISTANCIA.MAXIMA.OK <- VALORES.CRITERIO.OK[order(VALORES.CRITERIO.OK[,5]),]

if (!is.matrix(DISTANCIA.MAXIMA.OK))  DISTANCIA.MAXIMA.OK <- t(as.matrix(DISTANCIA.MAXIMA.OK))


DISTANCIA.MAXIMA.OK <- data.frame(DISTANCIA.MAXIMA.OK)


if (criterio=="AIC") names(DISTANCIA.MAXIMA.OK) <- c("p", "q", "P", "Q", "AIC")

   else if (criterio=="AICC") names(DISTANCIA.MAXIMA.OK) <- c("p", "q", "P", "Q", "AICC")

      else names(DISTANCIA.MAXIMA.OK) <- c("p", "q", "P", "Q", "BIC")



return(DISTANCIA.MAXIMA.OK[,c(1*(p.max!=0), 2*(q.max!=0), 3*(P.max!=0), 4*(Q.max!=0), 5)])


}