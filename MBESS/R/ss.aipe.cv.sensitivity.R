ss.aipe.cv.sensitivity <- function(True.C.of.V=NULL, Estimated.C.of.V=NULL, width=NULL, degree.of.certainty=NULL, assurance=NULL, certainty=NULL, mean=100, Specified.N=NULL, conf.level=.95, G=1000, print.iter=TRUE)
{
if(!is.null(certainty)& is.null(degree.of.certainty)&is.null(assurance)) degree.of.certainty<-certainty
if (is.null(assurance) && !is.null (degree.of.certainty)& is.null(certainty)) assurance <-degree.of.certainty
if (!is.null(assurance) && is.null (degree.of.certainty)& is.null(certainty)) assurance -> degree.of.certainty

if(!is.null(assurance) && !is.null (degree.of.certainty) && assurance!=degree.of.certainty) 
stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")

if(!is.null(assurance) && !is.null (certainty) && assurance!=certainty) 
stop("The arguments 'assurance' and 'certainty' must have the same value.")

if(!is.null(degree.of.certainty) && !is.null (certainty) && degree.of.certainty!=certainty) 
stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")

if(is.null(Estimated.C.of.V))
{
if(is.null(Specified.N)) stop("Since you did not specify an \'Estimated.C.of.V\', \'Specified.N\' must be specified.")
N <- Specified.N
}

if(!is.null(Estimated.C.of.V))
{
if(!is.null(Specified.N))  stop("Since you specified an \'Estimated.C.of.V\', \'Specified.N\' should not be specified.")
N <- ss.aipe.cv(C.of.V=Estimated.C.of.V, mu=NULL, sigma=NULL, width=width, conf.level=conf.level, alpha.lower=NULL, alpha.upper=NULL, degree.of.certainty=degree.of.certainty, Suppress.Statement=TRUE)
}

CN <- c("Lower.Limit", "Upper.Limit", "CV", "Int.OK", "Width")
Results <- matrix(NA, G, length(CN))
colnames(Results) <- CN

for(i in 1:G)
{
if(print.iter==TRUE) cat(c(i),"\n")
X <- rnorm(N, mean=mean, sd=True.C.of.V*mean)
CI.for.CV <- ci.cv(data=X, conf.level=conf.level)

Results[i,1] <- CI.for.CV$Lower
Results[i,2] <- CI.for.CV$Upper
Results[i,3] <- CI.for.CV$C.of.V

Results[i,4] <- sum((Results[i,1] <= True.C.of.V) & (True.C.of.V <= Results[i,2]))
Results[i,5] <- Results[i,2] - Results[i,1]
}

# Observed coefficients of variation.
Obs.CV <- Results[,3]

Results <- as.data.frame(Results)

Summary.of.Results <- list(Mean.CV=mean(Obs.CV), Median.CV=median(Obs.CV), SD.CV=(var(Obs.CV))^.5, 
Mean.CI.width=mean(Results[,2]-Results[,1]), Median.CI.width=median(Results[,2]-Results[,1]), SD.CI.width=(var(Results[,2]-Results[,1]))^.5, 
Pct.CI.Less.w=mean((Results[,2]-Results[,1])<=width)*100, Pct.CI.Miss.Low=mean(c(True.C.of.V <= Results[,1]))*100, Pct.CI.Miss.High=mean(c(True.C.of.V >= Results[,2]))*100, Total.Type.I.Error=(mean((True.C.of.V <= Results[,1]) | (True.C.of.V >= Results[,2]))))

###################################################################################################
# Vector of specification values.
if(is.null(degree.of.certainty)) degree.of.certainty <- 0

if(is.null(Estimated.C.of.V)) MBESS.tmp <- NULL
if(!is.null(Estimated.C.of.V)) MBESS.tmp <- round(Estimated.C.of.V, 4)

Specifications <- list(Sample.Size=round(N), True.C.of.V=round(True.C.of.V, 4), Estimated.C.of.V=MBESS.tmp, 
conf.level=conf.level, desired.width=width, degree.of.certainty=degree.of.certainty, G=round(G))

return(list(Data.from.Simulation=Results, Specifications=Specifications, Summary.of.Results=Summary.of.Results))
}
