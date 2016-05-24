#roxygen embedded tags:

#' To  carry out a search partition analysis (SPAN)
#'
#' @param  formula  A formula of the standard form \code{ y ~ x + u + v + w....} 
#' giving the outcome \eqn{y} and predictor covariates \eqn{x, u, v, w....}. Operators other 
#' than \code{+} should not be used. A survival object is allowed for \eqn{y}. For example,
#' \code{Surv(time,death) ~ x + u + v + w....} in which case optimation is with respect to
#' log-rank chi-square survival differences
#' @param  data  A data frame  with the variables in the formula.
#' @param weight  A frequency weight attached to each row of data. Default, NA, indicates unit weight to each data row.
#' @param cc Indicates complete case analysis (default FALSE). If TRUE, a row of data is deleted if any one  
#' attribute is missing. Otherwise a case is only deleted if any attribute is missing in a Boolean combination, as evaluated during a search.
#' Default FALSE
#' @param  makepos  If TRUE, and an attribute is found to be negative, the direction of \eqn{x} is reversed.  
#' The rule for reversal is if  \eqn{mean of y|x=1 < mean of y|x=0}. When \code{y} is a survival object the rule for creversal is
#' if \eqn{ rate |x=1 < rate |x=0}  where \eqn{rate= case/person-time}.  Default is TRUE. 
#' @param beta Parameter controlling degree of complexity penalising. Zero for no complexity penalising. NA (default) 
#' or negative determines a value for beta automatically as 0.03 times the initial gradient of the compleity hull. 
#' @param size Defines the  upper allowable size parameters of a disjunctive normal form used in the initial iteration of a search. 
#' It is a list of length \eqn{q} defining \eqn{p_1,p_2,..p_q}. Default \code{c(2,2,1)} defines 
#' \eqn{p_1=2}, \eqn{p_2=2}, and \eqn{p_3=1}. 
#' @param gamma Parameter controlling balance of observations in  \eqn{A} and its complement  \eqn{!A}. 
#' Default is NA, corresponds to no balancing. Balancing multiplies either MSE reduction or log-rank by
#' \eqn{(P_A(1-P_A))^\gamma}  where \eqn{P_A} is proportion of data in \eqn{A} to make a new optimization criterion.
#' 
#' 
#' @return Object \code{spanr} with attributes: 
#' 
#' \code{A} Data frame of same length as input data that is a binary indicator of belonging to \eqn{A}.
#' 
#' \code{g} Data frame of same length as input data, columns indicating  belonging to  the
#' subgroups of \eqn{A} 
#' 
#' \code{h} Data frame of same length as input data, columns indicating  belonging to  the
#' subgroups of \eqn{!A} 
#' 
#'  
#' @export
#' @author Roger Marshall <rj.marshall@@auckland.ac.nz>, The University of Auckland, New Zealand
#' 
#' @details A function to search for an optimal Boolean combination partition. Optimization is with respect to
#' reduction in mean square error of \code{y} by split into partition \eqn{(A,!A)}, or if \code{y} 
#' is a survival object, with respect to log-rank chi-square for survival differences of \eqn{(A,!A)}. 
#' The Boolean expression for \eqn{A} is output in normal disjunctive form \eqn{A= g_1 | g_2 | g_3 | ...} and 
#' the Boolean expression for the complement \eqn{!A} is also output in normal disjunctive form 
#' \eqn{!A = h_1 | h_2 | h_3 | ...}.  Each element of the disjunctive forms, \eqn{g_i} of \eqn{A},  or \eqn{h_i} of \eqn{!A},  of the 
#' represents a subgroup.  Subgroups are returned data frames.  
#' 
#' If variables \code{x, u, v, w....} of the formula  are not coded binary, a pre-analysis is done to establish 
#' an optimal cut of the variable. This is done, again with respect to reduction in MSE, or log-rank for a survival formula, 
#' over values of the variable. If numeric, 
#' a dictotomy is made by above/below a cut, the possible cuts being unique values of the variable if there are 20 or fewer, 
#' otherwise at 20 equally spaced intervals. If factor variable, according each value of the factor. 
#' 
#' @import plyr stringr survival
#' 
#' @examples
#' ## 1. Simulate Bernoulli binary predictors x1, x2...x10, and outcome y
#' ## For (x1 x2 x3) | (x1 x4) | (x1 x9),  make y~N(11,0.5) and N(10,0.5) otherwise.
#' x <- matrix(data=rbinom(10000,1,0.5),nrow=1000,ncol=10)
#' colnames(x) <- paste("x", seq(1:10), sep = "")
#' P <- ifelse((x[,1]& x[,2] & x[,3])|(x[,1] & x[,4])|x[,9] & x[,1], 1,0)
#' y <- ifelse(P,rnorm(1000,11,0.5),rnorm(1000,10,0.5) )
#' d <- data.frame(cbind(y,x))
#' sp <- spanr(formula= y ~ x1 +x2+x3+x4+x5+x6+x7+x8+x9+x10,data=d,size=c(1,2,2),beta=NA)
#' ## 2. Survival analysis of pbc data
#' library(survival)
#' data(pbc) 
#' sp <-with(pbc, spanr(formula = Surv(time, status==2) ~ trt + age + sex + ascites 
#'                + hepato + spiders + edema + bili + chol + albumin 
#'                + copper  + ast + trig + platelet + protime + stage,
#'                  beta=NA,cc=TRUE,gamma=1)   )
#' test <- cbind(pbc,sp$A)
#' ##Kaplan-Meier curves of A  versus !A
#' x <- survfit(Surv(test$time,test$status==2) ~ test$A)
#' plot(x, col=c(1,2))


spanr<- function(formula,  weight=NA, data=NULL,
      cc=FALSE,makepos=TRUE,beta=NA,size=c(2,2,1),gamma=NA)
{  
  
 
#library(stringr)
#library(plyr) 
#library(survival)

  if(is.null(data)){
    
    YX <- model.frame(formula,na.action=NULL)}
    else
    {
    
  YX <- model.frame(formula,na.action=NULL,data=data)}

  nobs <- nrow(YX)
  if(is.na(weight)==TRUE) {weight <- matrix(data=1, ncol=1,nrow=nrow(YX))}


  
  ## use complete cases only  
ok <- complete.cases(YX)


if(sum(!ok) >0){message(paste(sum(!ok),"/",nrow(YX),"observations have missing items") )}
else
{message(paste("There are no missing observation out of ", nrow(YX)))}

  if(cc ==TRUE){ 
  
    
    YX <- YX[ok,]
    weight <- weight[ok]
    nobscc <- nrow(YX)
    if(nobscc<nobs)
      {
    print(paste("Complete cases, n reduced from",nobs,"to",nobscc))
 
    id <- c(seq(from=1,to=length(ok),1) )
    idcases <- id[ok==TRUE] }
    
    else
      
    {print(paste("Complete cases, no deletions. n is",nobs))}
  }

 ##browser()
  
 

Gtype <- "Reduction in MSE"
survival <- is.Surv(YX[,1])

 #  ensure size argument is decreasing sequence 
  size <- sort(size,decreasing=TRUE)
 # remove ..$ prefix from  colnames(XY) if present 
  dollarpos <- str_locate(colnames(YX),pattern="\\$")

  if(is.na(dollarpos[1])==FALSE)
  {
    ## ddply doesn't like $ characters in names, strip out with sub()
    ## replaces all occurrences of $ with nothing
    
    colnames(YX) <- sub(pattern="\\$",replacement="",x=colnames(YX))
    message("Stripping $ character from attribute strings")
    ## also strip out data name i.e xxxx  for xxxx$
    ## assumes that all xxxx sames??
    
##    dataname <- substr(colnames(YX),1,dollarpos[1]-1)
##    colnames(YX) <- sub(pattern=dataname,replacement="",x=colnames(YX))
    ##browser()
  }
  
  #  for reasons unknown column 1 name gets lost on Y
  nameY <- colnames(YX)[1]
 #  check whether is a factor
 
 if(is.factor(YX[,1])==TRUE) {
  message(paste("Error", nameY,"outcome var is not numeric. Cannot continue"))
   return()
 }

# check whether is a survival object
##Gtype <- "Reduction in MSE"
##survival <- is.Surv(YX[,1])

if(survival){print(paste(nameY,"is a survival outcome object"))
             TimeEvent <- as.matrix(YX[,1])
             Event <- TimeEvent[,2]
             Time  <- TimeEvent[,1]
             Gtype <- "Log-rank chi-sq"
             nameY <- "Time"
}

if(is.na(gamma)==FALSE & gamma !=0) {
  Gtype <- paste("Balanced", Gtype,"(gamma=",gamma,")")}
else
{gamma <- 0}
print(paste("Optimisation w.r.t G=",Gtype, "of", nameY))
  
                
   ##browser() 
   # check whether Y varies at all

if(survival) {uni <- unique(Time)}
else {uni <- unique(YX[,1])}
 
  if( length(which(is.na(uni)==FALSE)) ==1   ){

   stop(paste("Error",nameY,"outcome has no variability"))
    return()
   
 }
 #  separate the Y and X variables
  if(survival){
    
    Y<- as.matrix(Time)
    colnames(Y) <- "Time"
  }
 else
  {  
    Y <- as.data.frame(YX[,1])
  }
 

     
  names(Y) <- nameY
  X <- YX[,2:ncol(YX)]
  ##browser()
##  print(paste("Finding optimal cut, w.r.t ",nameY,", of the variables:"))
##  print(colnames(X))
  newdata <- matrix(nrow=nrow(X),ncol=0)
  MSE <- vector(length=ncol(X))


##------ loop over elements of X-----------------------------------------
print("Optimising cuts....")
  for(i in 1:(ncol(X))){
    
    
    # check whether Y varies at all
    
    uni <- unique(X[,i])
    
    if( length(which(is.na(uni)==FALSE)) ==1   ){
      
     warning(paste("Attribute -", colnames(X)[i],"- has no variability."))
   ##   return()
      
    }
    
   
    charvar <- is.factor(X[,i])
    if(charvar){
## replace any blanks in X data with underscore      
     X[,i] <- gsub(pattern=" ",replacement="_",x=X[,i])
    }  
    
  cuts <- sort(unique(X[,i]))
  ncuts <- length(cuts)
 #  is i-th strictly binary 0, 1 values ?
 binary <- FALSE
 if(charvar==FALSE)
{
   binary <- identical(as.integer(cuts),as.integer(c(0,1)))
      
 ## for non-character reduce number of cuts by one   
    ncuts <- length(cuts)-1
 }

     
    {if(ncuts>20){
## take a subsample of cuts, every ncut/20 th, start at second?
                 every <- floor(ncuts/21)+1
                 cuts <- cuts[seq(2,ncuts,every)]
                 ncuts <- length(cuts) - 1 
     
## or percentiles?? Test code:
#     x <- X[,i]
#     cuts <- quantile(x,probs=seq(0.02,0.98,0.02),type=4,rm.na=TRUE)
#    cuts <- as.vector(cuts)
#   ncuts=length(cuts)-1 

    }
    }
  
             
    
    testd <- matrix(ncol=ncuts,nrow=nrow(X))
    #  browser()
    ## if already binary
    if(ncuts>0){
      ## if strictly binary 0,1 retain same name
      if(binary){
        testd <- X[,i]
        testd <- as.matrix(testd)
        
        colnames(testd) <- paste(colnames(X)[i],"_eq_1",sep="")
        # can use span() to return mse
 ##browser()
        spcut <- span(d=testd,f=weight,Y=Y,makepos=FALSE,po=c(1),cc=cc,
                quietly=TRUE,gamma=gamma, survival=survival, Event=Event)
##        print(paste("Binary (0/1) variable ",colnames(testd), 
##                   "G=", round(spcut$mse,digits=7)))
        
      }
      else  ##of if binary
      {
        if(charvar)
        {
 #         library(stringr)
          ## replace " " signs with an "_"         
##          colnames(testd) <- gsub(pattern=" ",replacement="_",x=colnames(testd))
##cuts <- gsub(pattern=" ",replacement="_",x=cuts)

          colnames(testd) <- as.character(paste(colnames(X)[i],cuts[1:ncuts],sep="_eq_") )

#library(stringr)
## replace "-" signs with an "_"         
colnames(testd) <- gsub(pattern="-",replacement="_",x=colnames(testd))

##         browser()
        }
        else
        {
          #  use cut level >= as postfix to variable name, skip bottom cut
          ## deal with possibility of negative cut, as  - sign not allowed? 
          colnames(testd) <- as.character( paste( colnames(X)[i],round( cuts[1:(ncuts)],digits=7),sep="_gt_")  )
 
 #         library(stringr)
 ## replace "-" signs with an "m"         
          colnames(testd) <- gsub(pattern="-",replacement="m",x=colnames(testd))
          ## attempt to deal with  > and = in colnames using make.names. 
 ##  But simply replaces with dot .
 #colnames(testd) <- make.names(colnames(testd))
        }
      }
  #-----------------------------------------------------------------    
      for (j in 1:ncuts)
      {
        
        if(charvar) {
          # character variable, take 0 to be  first alpha-numeric order
          testd[,j] <- ifelse(X[,i] == cuts[j],1,0)  
          
        }
        else
        {testd[,j] <- ifelse(X[,i] > cuts[j],1,0)} 
      }  #for (j...
 #--------------------------------------------------------------------
 #browser()
      
      if(ncuts==1 ){
        newdata <- cbind(newdata,testd)
## i think binary must be FALSE to arrive here!
        if(binary==FALSE)
          
          {##print(paste("Binary variable created:",colnames(newdata)[i] ))
 
## this needs testing.....
spcut <- span(d=testd,f=weight,Y=Y,makepos=FALSE,po=c(1),cc=cc,
              quietly=TRUE,gamma=gamma,survival=survival, Event=Event)
##print(paste("Already 2-class variable ",colnames(testd), 
##            "with reduction in MSE", round(spcut$mse,digits=7)))

        } 
      }
      
      else
        ## invoke SPANR search for best cut  
      {

       if(charvar){print(paste("...over",ncuts, "cuts of factor", colnames(X)[i]))}
                  else
                  {print(paste("...over",ncuts, "cuts of numeric", colnames(X)[i]))}          

      
      for(j in 1:ncuts ) {
      if(charvar){  message(" ", cuts[j],appendLF=FALSE)}
      else
      {message(" ", signif(cuts[j],3) ,appendLF=FALSE)}
        #    message(" ",appendLF=FALSE)
      }
      message(" ")
      

        ##  carry out quiet search for optimum cutpoint
      #browser()
        spcut <- span(d=testd,f=weight,Y=Y,makepos=FALSE,po=c(1),cc=cc,
                      quietly=TRUE,gamma=gamma,survival=survival, Event=Event)
        
        newdata <- cbind(newdata,spcut$subgA[1])
      
      #browser()
 #     print(paste("Optimal", colnames(X)[i], "cut:" , colnames(newdata)[i],
 #     "G=", round(spcut$mse,digits=8) ))
       
      } 

 
     }  #ncuts>0
  MSE[i] <- spcut$mse

  }  #for(i in.....
##-----------------------------------------------------------------------------
## rank the attributes by descending MSE
rank <- order(MSE, decreasing=TRUE)
MSE <- MSE[rank]

newdata <- newdata[,rank]

print(paste("Optimal attributes ranked by G:"))

for(i in 1:(ncol(X))){
  message(paste("        ", colnames(newdata)[i],round(MSE[i],digits=8)))
}
print("Doing SPANR main search.....")
#  browser()
##test out truncating number of attributes
if(ncol(newdata) >10) {
  

  drop <- colnames(newdata)[11:ncol(X)]
  print(  "Limiting search over top 10 attributes: dropping:" )
  for(i in 1:(ncol(X)-10) ) {
  message(drop[i],appendLF=FALSE)
  message(" ",appendLF=FALSE)
}
message(" ")
#  message( list(paste(drop[1:(ncol(X)-10)],sep=" ") ), appendLF=FALSE) 
  newdata <- newdata[,1:10]
}

## automatically set beta here.
if(is.na(beta)){ beta <- 0.03*MSE[1]}
#-----------------------------------------------------------------------------
## carry out the SPAN search
test <- span(d=newdata,f=weight,Y=Y,makepos=makepos,po=size,beta=beta,
        quietly=FALSE,cc=cc,gamma=gamma,survival=survival, Event=Event)

#-------------------------------------------------------------------------------
  
if(cc ){
  if(nobs!=nobscc){
# need to return with data.frame of same length filled with NAs
 id <- as.matrix(id)
 colnames(id) <- "id"
 XX <- cbind(test$subgA,idcases)
 XX <- merge(x=id,y=XX, by.x="id",by.y="idcases",all.x=TRUE)
 test$subgA <- XX[,2:ncol(XX)]

 XX <- cbind(test$subgAc,idcases)
XX <- merge(x=id,y=XX, by.x="id",by.y="idcases",all.x=TRUE)
 test$subgAc <- XX[,2:ncol(XX)]
 
  }
}

## create output data frames g, h, A
qA <- ncol(test$subgA)-1
nmc <- colnames(test$subgA)
g <- as.data.frame(test$subgA[,1:qA])
colnames(g) <- nmc[1:qA]

#  check on balance of partition
a <- test$subgA[,qA+1]

a <- a[which(is.na(a)==FALSE)]
w <- weight[which(is.na(a)==FALSE)]
pA <- sum(a*w)/sum(w)

#if(pA < 0.05 | pA >0.95){
  message(paste("Balance  P(A):P(!A)= 1:",round( (1-pA)/pA,digits=3 ) ))
#}



A <- as.data.frame(test$subgA[,qA+1])
colnames(A) <- "A"

qAc <- ncol(test$subgAc)-1
nmc <- colnames(test$subgAc)
h <- as.data.frame(test$subgAc[,1:qAc])
colnames(h) <- nmc[1:qAc]

##if(is.null(test$warn) != TRUE){
##  print("There are notes/warning messages")
##}

output <- list(g=g, h=h, A=A)



  return(output)
  
}  ##end spanr
#---------------------------------------------------------------------------------------------
span <- function(d,  Y,  f=NA, cc=FALSE,
      makepos=TRUE,beta=NA,po=c(2,2,1),quietly=FALSE,gamma=0,
      survival=FALSE, Event=NA){
  
  
  d <- as.matrix(d)
  Y <- as.matrix(Y)
  yvar <- colnames(Y)
  f <- as.matrix(f)
  ## augment d (attribute) data with Event indicator if survival model
  ## then m-1 is actually number of attibutes is. m until .Fortran called 
  ## is the number of attributes + 1. Is reset m <- m-1  for .Fortran call
  ##  This is to ensure that collapsing and all other operation prior to 
  ##  .Fortran take account of event indicator too. 
  if(survival){d <- cbind(d,Event)}
  
  qo <- length(po)

  m <- ncol(d)
  ## m is size of attribute set, as dictated by size of binary input matrix d
  ##  but if survival model is is one plus the numnber of attributes
  ## check that Yvar is not one of the variables in d.
  ## stop if it is. 
  
  
  ##browser()
  for (i in 1:m){
  if(yvar==colnames(d)[i]){
   note <-paste("Y variable is same as an attribute: ",yvar)
   stop(note)
   return()
  }
  }
  
  
  ##-------------------------------------------------------------
  ## collapse data for turbo....
  nentry <- nrow(d)
  n <- nentry
  
  lenY <- length(unique(Y))
  
  ## do a collapse? When is it worthwhile? Try a criterion based on # possible
  ## combination of attribute (2**m) and number of unique values of Y
##  if(0.1*lenY*(2**m) <= n | TRUE) {
 ##  if(n >= 100
## Maybe always??
     if(TRUE)
      {
  
  dYf <- as.data.frame( cbind(d,Y,f) )

   
#  library(plyr)
  colnames(dYf)[ncol(dYf)] <- c("freq")
  if(quietly==FALSE) print("Collapsing data...")

  dYf <- ddply(dYf, colnames(dYf)[1:(ncol(dYf)-1)], function(xx)sum(xx$freq) )

  ## is possible repeated names in  attributes: ddply collapses over these 
  ## to so need to reset m. e.g. may happen if optimed cutpoints repeated
  
  m <- ncol(dYf)-2
  d_collapsed <- dYf[,1:m]
  Y <- dYf[,(ncol(dYf)-1)]
  
  f <- dYf[,ncol(dYf)]
  d_collapsed <- as.matrix(d_collapsed)
  Y <- as.matrix(Y)
  f <- as.matrix(f)

#   
  n <- nrow(d_collapsed)

note <- c(paste("Data frequency collapsed from ",nentry, " records to ",n))  
if(quietly==FALSE & n<nentry) print(note)

  }  
else   # of if collapsing data
  
 { 
  if(quietly==FALSE) print("No attempt to collapse data n<1000")

d_collapsed <- d
 }
#--------------------------------------------------------------------------
  warn <- NULL
  
## use complete cases only  
 if(cc ==TRUE){ 
  ok <- complete.cases(d_collapsed,Y)
 }

else{
  
  ## eliminate missing Y only

ok <- complete.cases(Y)

}

#browser()
  dcc <-  d_collapsed[ok,]
## possible when m=1 that dcc is not matrix, coerce into matrix
  dcc <- as.matrix(dcc)
  Y <- Y[ok] 
  f <- f[ok]
##browser()
ncc <- nrow(dcc)
##browser()
if(ncc < n){warn <- append(warn, paste("Complete cases reduced # records to =", ncc, "# observation =", sum(f,na.rm=TRUE)))
if(quietly==FALSE) print(warn)

}
 

   ##test for strictly binary  variables

  for (i in 1:m){
    test <-unique(dcc[,i])
    
    v <-colnames(d)[i] 
    if(length(test)>3){
     
     stop(paste("Attribute ",v, " has >3 states. Not binary."))
      return() }
    
    else
 
{  

       for(k in 1:length(test)){
         
         if(is.na(test[k])==FALSE) {
           
          if( test[k]!=0 & test[k] != 1){ 
          stop(paste("Attribute ",v, "  not coded 0 1. Quitting."))
           return()
         }
       }
       
  }
      }
}

  ## test whether positive attributes
  nonpos <- 0
v_nonpos <- NULL
  matt <- m
  if(survival) matt <- m-1
  for (i in 1:matt){
    
    
    if(survival){
##  compute rates  cases / person time to decide direction . Event counts in m th column of dcc
      meanYA  <- sum(f*dcc[,m]*dcc[,i],na.rm = TRUE)/sum(Y*f*dcc[,i],na.rm = TRUE)
      meanYAc <- sum(f*dcc[,m]*(1-dcc[,i]),na.rm = TRUE)/sum(Y*f*(1-dcc[,i]),na.rm = TRUE)
 ##     if(quietly==FALSE)  browser()
    }
    else
    {
    meanYA  <- sum(Y*f*dcc[,i],na.rm = TRUE)/sum(f*dcc[,i],na.rm = TRUE)
    meanYAc <- sum(Y*f*(1-dcc[,i]),na.rm = TRUE)/sum(f*(1-dcc[,i]),na.rm = TRUE)
    }
    
    if(is.na(meanYA) ==TRUE | is.na(meanYAc)==TRUE){

      v <-colnames(d)[i]
     stop(paste("Attribute  ", v," is redundant - all 0 or all 1. Quitting"))
     return()}

    if(meanYA < meanYAc){
      v <-colnames(d)[i]
##    if(!quietly) browser()
    warn <- append(warn, paste("Mean ",yvar, " greater for !", v, sep="") )

    nonpos <- nonpos +1 
 
 v_nonpos  <- append(v_nonpos,paste(v))
 if(makepos == TRUE){
   ##reverse the  direction, add ! to indicate  positive attibute is reversed
    colnames(d)[i] <- paste("!", colnames(d)[i],sep="")
    dcc[,i] <- 1-dcc[,i]
    d[,i]  <- 1-d[,i]
    
 }}
  
  }  # for i 1:matt
 
if(nonpos >0){
  
  v_nonpos <- paste(v_nonpos, collapse = ' ')
  
  if(quietly==FALSE){

 if(nonpos==1){
     
     if(makepos){print(paste("There is one non-positive attribute",
   v_nonpos,"Direction is reversed"))}
    
   else
   {print(paste("There is one non-positive attribute",v_nonpos))}
 
 }
 
  else # of if(nonpos=1) 
    
  {   if(makepos){ print(paste("There are ",nonpos,"non-positive attributes",v_nonpos
                                      ,"Directions are reversed"))}
      else 
{ print(paste("There are ",nonpos,"non-positive attributes",v_nonpos))}
}
}
}  # nonpos>0
#---------------------------------------------------------
## load DLL to  do SPAN analysis
##  dyn.load should be commented out in  CRAN package release

##  dyn.load("C:/spanr/q6000036408/Dll_spanr/Dll_spanr/x64/Debug/Dll_spanr.dll")
#---------------------------------------------------------
  ## invoke associated srd subroutine  to fit rectangle coordinates and % error
  
  
  combA <- rep(-1,100)
  combAc <- rep(-1,100)
  hc <- rep(-1,100)
  hg <- rep(-1,100)
  c_opt <- 0
  binw <- 0
   gfreq <- vector(length=2000)
   gfreq <- rep(0,2000)

##  make all NA s of  data -1, to handle by FORTRAN
##  if not complete cases analysis

if(cc==FALSE ) {dcc[is.na(dcc)] <- -1}

## invoke associated srd subroutine  to fit rectangle coordinates and % error
##  FORTRAN subroutine spanr(x,y,freq,n,m, combA, combAc,hc,hg,c_opt,beta, status)
## long is dimensioning device for array f() in spanr. Must be at least n*max number
##  of attributes (included added at each iteration). So I think is > n*(m+lc) and lc=20
##  setinside stanr.for. Be conservative, but keep an eye on.  
   long <- (m+10)*ncc
##if(!quietly) print("Doing SPANR search...")
##  negative beta for automatic setting, or if  NA
if(is.na(beta)){beta <- -1} 
  status <- 0
  npt <- 0
  gopt <- 0
## make another survival surv variable for safety. If survival  ensure 
## number of attributes is m-1, as last column of dcc is Event indicator
  surv <- survival 
  if(surv) m <- m-1

if(quietly==FALSE)print("Searching....")
  spanr <- .Fortran("spanr",dcc=as.integer(dcc),Y=as.single(Y),f=as.integer(f),ncc=as.integer(ncc),
                 m=as.integer(m),long=as.integer(long),combA=as.integer(combA),combAc=as.integer(combAc),
                 hc=as.single(hc),hg=as.single(hg), c_opt=as.integer(c_opt), beta=as.single(beta),
                 status=as.integer(status), npt=as.integer(npt), 
                 qo=as.integer(qo),po=as.integer(po),
                 gopt=as.single(gopt), gfreq=as.integer(gfreq), binw=as.single(binw),gamma=as.single(gamma),
                 surv=as.logical(surv), NAOK=TRUE )

# returned gopt is complexity penalised, de-penalise
#  returned may be MSE or log-rank, continue use "mse" generically
mse <- spanr$gopt+spanr$c_opt*spanr$beta   
status <-  "OK"
if(spanr$status !=0) status <- "Warnings"
if(quietly==FALSE) {print(paste("SPANR search complete. ", spanr$npt," partitions evaluated. Exit status ",status))
Gtype <- "reduction in MSE"
if(survival) Gtype <- "log-rank chi-sq"
if(gamma >0) Gtype <- paste("Balanced",Gtype)
                                    
print(paste("Exit", Gtype, round(mse,digits=8),"Exit beta ", round(spanr$beta,digits=8),
            "Optimal complexity",spanr$c_opt))
}
if(spanr$status !=0){
  warn <- append(warn, paste("Search status fail:",spanr$status ) )
if(spanr$status ==2) print("Status 2 fail: too large parameter beta")
## switch of what follows using quielty switch
quietly <- TRUE
}


##browser()
##if(!quietly) browser()
  
## need to trim off padded -1's in hc and hg
 spanr$hc <- spanr$hc[1:trim(spanr$hc)]
 spanr$hg <- spanr$hg[1:trim(spanr$hg)]


if(cc==FALSE ) {dcc[dcc==-1] <- NA}

##extract  q  and Boolean combination in readable form from  $combA and combAc

sizeA  <-  size(comb=spanr$combA, id=1, d)
sizeAc <-  size(comb=spanr$combAc,id=0, d)
##Boolean expressions in 2nd element
expA  <- sizeA[2]
expAc <- sizeAc[2]
partition <- paste("A=(", sizeA[2],")  !A=(",sizeAc[2],")")
## actual q, qAc in first element
## plot out the complexity hull 
##browser()
##if(quietly==FALSE) {print( paste("Optimal partition as a dnf is:"))
##message(paste("     A=(", sizeA[2],")") )
##print( paste("with dnf complement:"))
##message(paste("    !A=(", sizeAc[2],")") )}
if(quietly==FALSE){spanr.hull(spanr$hc,spanr$hg,spanr$c_opt,spanr$beta, 
                              partition,spanr$gfreq,spanr$binw, spanr$gamma,Gtype)
##print(paste("Plot of convex hull created. Optimal complexity =",spanr$c_opt))
}
qA <- as.numeric(sizeA[1])
qAc <- as.numeric(sizeAc[1])
# check on c_opt
if(qAc +qA -1 != spanr$c_opt) {
  warn <- append(warn, paste("Mean complexity check failed c_opt=",c_opt ) )}

## set up data matrix of  AND-ed combinations of each  subgroup in d.n.f
##  use the qA+1 th element for indicator of A itself



subgA <- subgroups(comb=spanr$combA, id=1, d,nentry,qA)
subgA <- as.data.frame(subgA)


#testing writing out subgroups more neatly
if(quietly==FALSE){
print( paste("Optimal partition as a dnf is:"))
message(" ")
message(paste("     A=(", colnames(subgA)[1],")") )
if(qA>1){for (i in 2:qA) {     message(paste("       (", colnames(subgA)[i]),")" ) }}
}

subgAc <- subgroups(comb=spanr$combAc, id=0, d,nentry,qAc)
subgAc <- as.data.frame(subgAc)
if(quietly==FALSE){
  message(" ")
  message(paste("    !A=(", colnames(subgAc)[1],")") )
  if(qAc>1){for (i in 2:qAc) {     message(paste("       (", colnames(subgAc)[i]),")" ) }}
  message(" ")}

output <- list(subgA=subgA, subgAc=subgAc, warn=warn, mse=mse)


return(output)

}
############################################################################
size <- function(comb, id, d){
  
  ## extracts Boolean expression associated with combination of the
  ## form comb=  0....0....0...0-1-1-1 that is created by FORTRAN spanr
  attributes <- colnames(d)

exp <- NULL
j <- 2
q <- 1
p <- 1


while (comb[j] != -1) {
  if(comb[j]==0){

    
    q <- q+ 1
    p <- 0

    
  }
  
  else
  { 
    
    if(p==0){
    exp  <- append(exp,")(")
    }
    p <- p+1
    
##    exp  <- append(exp,paste(attributes[comb[j]],"=",id, sep=""))

if(id ==1){
exp  <- append(exp,paste(attributes[comb[j]], sep=""))
}
else
{
  if(substr(attributes[comb[j]],1,1)=="!"){
## attribute has been switched, to become  !... need to remove the !     
    
    s <- substr(attributes[comb[j]],2,nchar(attributes[comb[j]]) )
    

    exp  <- append(exp,paste(s, sep=""))}
  else{
  exp  <- append(exp,paste("!",attributes[comb[j]], sep=""))}
}  

}
j <- j+1
}
##exit loop , q increment too much on last  0
q <- q -1
exp <- paste(exp, collapse = ' ')

output=c(q,exp)

return(output)
}
########################################################################

subgroups <- function(comb, id, d,n,q){
  
  ## creates data matrix  expression of the subgroups associated with
  ## a spanr analysis.  Data matrix has q+1 columns,  1,..q, are for each of
  ## q  subgroups. The q+1 th is of the overall (union) of subgroups. 
  

  subg <- matrix(nrow=n,ncol=q+1)
  subg_names <- matrix(nrow=1,ncol=q+1)
  subg[,q+1] <- 0

  j <- 2
  i <- 1
  p <- 1
  name=NULL
  prod <- rep(1,times=n)
  
  while (comb[j] != -1) {
    if(comb[j]==0){     #0  indicates an "or" in boolean expression <-> |
      
      subg[,i] <- prod
      # insert _  between elements of subgroup. This  is only necessary
      #  to avoid problem ddply in srd function
      subg_names[i] <-  paste(name, collapse = ' ')

      name <- NULL
      i <- i+ 1
      p <- 0
      prod <- rep(1,times=n)
      
    }
    
    else    #  of comb[j]==0, i.e. inside an & expression
    { 
  
      
      
      if(id==1)
        
        { prod <- prod*d[,comb[j]]
                 name <- append(name, colnames(d)[comb[j]]) 
      
      }
      else     
        
        { prod <- prod*(1-d[,comb[j]])
                 
                 
                 if(substr(colnames(d)[comb[j]],1,1)=="!"){
                   ## attribute has been switched, to become  !... need to remove the !     
                   
                   s <- substr(colnames(d)[comb[j]],2,nchar(colnames(d)[comb[j]]) )
                   
                   
                  name <- append(name,s)}
                 else{
                   name  <- append(name,paste("!",colnames(d)[comb[j]],sep=""))}
      }  
                 
      
    }
    j <- j+1
  }
  ##exit loop , q increment too much on last  0
  
  for (k in 1:n){
    subg[k,q+1] <-max(subg[k,])
  }

  if(id==1)
  
  {subg_names[q+1] <- "A"}
  else 
  
  {subg_names[q+1] <- "Ac"}

  colnames(subg) <- subg_names

  output=subg

  
  return(output)
}

#################################################################
trim <- function(x){
  j=1
  while(x[j]>=0)
    {j <- j+1 }
  return(j-1)
}

##################################################################

spanr.hull <- function(hc,hg,c_opt,beta, partition, gfreq,binw, gam,Gtype){
  
  
  
  #library(graphics)
## set up plotting position of frequency boxes 
## corresponding to gfreq(1:2000)
## centre frequency box, lower than midpoint of bin, i.e. at 0.8
## to avoid box beyond upper hull, impression of points greater than hull
  freqbin <- rep(seq(1:100)*binw - 0.8*binw,times=20)
  cfreq <- sort(rep(seq(1:20),times=100) )
  
  freqbin <- freqbin[which( gfreq!=0)]
  cfreq <-  cfreq[which( gfreq!=0)]
  gfreq <- gfreq[which( gfreq!=0)]
  #browser()
  hc[length(hc)]
  maxhg <- hg[length(hc)]
  ##browser()
#  plot.window(xlim=c(0,hc[length(hc)]), ylim=c(0,hg[length(hc)]))
par(bg = "#FFFFFF")
plot.new()     
plot.window(xlim=c(0,max(cfreq)), ylim=c(0,hg[length(hc)]))

  axis(1)
  axis(2)
  
  title(main="Convex hull and pattern of search") 
   if(gam != 0){
     part2 <- as.character(round(gam,digits=1))
   
     part1 <- paste(Gtype," ")
## josh's code!
     title(ylab = bquote(.(part1) ~ gamma ~ "=" ~ .(part2)))
    }
    else{title(ylab=Gtype)}
  
  title(xlab="Complexity c")
  
#  add frequency of search points
gfmax <- max(gfreq)
#  draw frequency boxes w area proportional to frequency 
points(cfreq,freqbin, col="grey", pch=15, cex=3.5*sqrt(gfreq/gfmax)+0.4 )
#  draw upper hull points
  points(hc,hg, cex=1.4, pch=19, bg="black")
#  join up points  
  lines(hc,hg,lwd=2)
  w <- which( hc==c_opt)
#  optimal point shown in red
  points(hc[w],hg[w],col="red",cex=1.2, pch=21, bg="red")
  
  int <- hg[w]-beta*c_opt
  abline(a=int,b=beta,col="blue",lty=2)
  
  text(0.8,int+beta*0.8,expression( paste("slope ",beta)  ) )
  text(0.8,0.95*(int+beta*0.8),paste(round(beta,digits=6) ), cex=0.7)
  #arrows(hc[w],hg[w],hc[w],0,col="gray",length = 0.1, angle = 25)
  text(hc[w],0,paste("Optimal at c= ",c_opt), cex=0.8)
  segments(hc[w],-0.05*hg[w],hc[w],0.05*hg[w],col="black")
  ##text(0,0.05*maxhg,paste(partition),cex=0.8,adj=c(0,0),col="blue")
  box()
  ##browser()
}
