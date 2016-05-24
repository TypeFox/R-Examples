cstrata.psa <- function(treatment, propensity, strata = NULL, 
	int = NULL, tree = FALSE, minsize = 2, graphic = TRUE, colors = c("dark blue", "dark green"), 	xlab = "Estimated Propensity Scores with Random Heights",  
	pch = c(16, 16))
	
{
#treatment == binary vector of 0s and 1s (necessarily? what if character, or 1, 2?)
#propensity == PS scores from some method or other.
#strata == either a vector of strata number for each row of covariate, or one number n in which case it is attempted to group rows by ps scores into n strata of size approximately 1/n.  
#This does not seem to work well in the case of few specific propensity values, as from a tree.
#int == a vector of cut points between 0 and 1 defining the subintervals (perhaps as suggested by loess.psa).  In either case these subintervals define strata, so strata can be of any size.  
#tree == logical, if unique ps scores are few, as from a recursively partitioned tree, then TRUE will force each ps value to define a stratum.  
#minsize == smallest allowable stratum-treatment size.  If violated, strata is removed.  Is this the best idea?
#colors = 2-ary color vector for the points
#xlab == x axis label
#ylab == y axis label
#pch == size of points  

######################################################## BEGIN A

if(tree & (!is.null(strata) | !is.null(int))){stop("If 'tree == TRUE' then 'strata' and 'int' must be NULL")}
if(tree){int <- sort(unique(propensity))}
if(is.null(strata) & is.null(int)){stop("One of 'strata' or 'int' must be defined.")}
if(!is.null(strata) & !is.null(int)){stop("Only one of 'strata' or 'int' may be defined.")}
if(!is.null(int)) int <- sort(unique(int))

# If strata is one number n, then rows are grouped into n similarly sized (length(treatment)/n) strata 
if(length(strata) == 1) 
    {strata <- findInterval(propensity, vec = quantile(propensity, seq(0,1,1/strata)), rightmost.closed = TRUE)
    }
    
#J: Defines sho (counts table) if strata are given explicitly 
if(!is.null(strata))
  {sho <- table(treatment,strata)
   psct <- strata
   if(length( unique(strata) )> min(26, nrow(treatment)/5) ) 
         {stop("The number of strata (before minsize adjustment) cannot exceed 26 or N/5")
         }
  }

#J: If strata is null, rely upon int
if(is.null(strata)){
	 if(is.null(int))
         {int = sort(unique(propensity))
	     }
   ed <- length(int)
   if(ed > min(26, nrow(treatment)/5))
         {stop("The number of strata (before minsize adjustment) cannot exceed 26 or N/5")
         }
#   if(ed == 1) 
#     {strat <- findInterval(propensity, vec = seq(0,1,1/int), rightmost.closed = TRUE)
#      sho = table(treatment,strat)  
#     }
   if(ed == 1){stop("'int' should be a vector of propensity scores defining the strata.")}
   if(ed > 1)
     {if(int[1] > 0 & int[ed] == 1){int <- c(0,int)}
      if(int[1] == 0 & int[ed] < 1){int <- c(int,1)}
      if(int[1] > 0 & int[ed] < 1){int <- c(0,int,1)}
      strat <- findInterval(propensity, vec = int, rightmost.closed = FALSE)
      sho = table(treatment, strat)  
     }                                                                               
   psct <- strat
  }

#This is to find the cutting points of propensity score between different strata
#It is named as cup finally
sn = sort(unique(psct))
alv = NULL
for(i in 1:length(sn))
   {alv[i*2-1] = min(propensity[psct == sn[i]])
    alv[i*2] = max(propensity[psct == sn[i]])
   }
alv <- sort(alv)
cup = NULL
for(j in 1:(length(sn) - 1))
   {cup[j] = (alv[j*2]+alv[j*2+1])/2
   }
   
######################################################## END A

######################################################## BEGIN B

shom = matrix(sho, nrow=2) 
dimnames(shom) = list( dimnames(sho)$treatment, dimnames(sho)$strat)

#J: Remove strata with number of observations < 'minsize'.
cutstr <- NULL 
for (i in 1:dim(shom)[2]) {if (min(shom[,i]) < minsize){cutstr <- c(cutstr, i)}  
                          } 
if (length(cutstr) > 0){shomn <- shom[,-cutstr]} else{shomn <- shom} 

#Define 'som' as matrix of strata counts for treatment/control and strata sums after minsize restriction is enforced. 
som = rbind(shomn, colSums(shomn)) 
rownames(som)[3] <- "Strata.Size" 
colnames(som) <- colnames(shomn)

######################################################## END B

######################################################## BEGIN C

if(graphic)
{
#the following codes are adapted from loess.psa, 
#minor changes are made to show sample size difference between treatment and control
#the y axis value in this graph is randomly generated.
        height = runif(length(treatment),1,100)
        htp <- data.frame(height, treatment, propensity)
        dimnames(htp)[2] <- list(c("height", "treatment", "propensity"))
        sut <- sort(unique(treatment))
if(!length(strata)==0){xlab<-"Strata as Input (random heights)"}
if(length(strata)==0){if(ed == 1){xlab <- paste(int, "Equally Spaced Strata (random heights)", sep = " ")}else{xlab <- "Strata Cutpoints as Input (random heights)"}}
par(mfrow=c(2,1))
    plot(htp$p, htp$h, type = "n", xlim =range(htp$p), xlab = "", 
         ylab = "", main=paste("Treatment =",sut[1]), yaxt="n"      )
    points(htp$p[treatment == sut[1]], htp$h[treatment == sut[1]], 
         col = colors[1], pch = pch[1])
    abline(v=cup)

    plot(htp$p, htp$h, type = "n", xlim =range(htp$p), xlab = xlab, 
         ylab = "", main=paste("Treatment =",sut[2]),  yaxt = "n"   )

    points(htp$p[treatment == sut[2]], htp$h[treatment == sut[2]], 
        col = colors[2], pch = pch[2])
    abline(v=cup)
par(mfrow=c(1,1))
}
######################################################## END C
out <- list(shom, som, psct)
names(out) <- c("Original.Strata", "Used.Strata", "strata")
return(out)
}
