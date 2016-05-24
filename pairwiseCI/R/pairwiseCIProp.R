
"pairwiseCIProp" <-
function(formula, data, alternative="two.sided",
 conf.level=0.95, control=NULL, method=NULL, ...)

{

arglist<-list(...)

arglist$alternative <- alternative
arglist$conf.level <- conf.level


if(method %in% c("Prop.ratio", "Prop.or", "Quasibin.ratio", "ODbin.ratio", "Betabin.ratio")) 
 {sepcompname <- "/"}
else
 {sepcompname <- "-"}

# check the arguments

mf <- model.frame(data, formula=formula)

if(!is.numeric(mf[,1]))
 {stop("response variable must be numeric")}

if(!is.factor(mf[,2]))
 {mf[,2]<-as.factor(mf[,2])}

if(any(as.integer(floor(mf[,1])) - as.numeric(mf[,1])<0))
 {warning("At least one value in the response is not an integer!")}

# # # Split up into pairwise comparisons

grouplist <- split(x=mf, f=mf[,2], drop=TRUE  )
datalist <-split(x=mf[,1, drop=FALSE], f=mf[,2], drop=TRUE  )

levelnames <- names(grouplist)

k=length(grouplist)


# # # ALL pairs comparisons

if(is.null(control))
{
 {
 CM <- c()
 rnames <- c()
 kindx <- 1:k
        for (i in 1:(k - 1))
         {
            for (j in (i + 1):k) 
             {
                CM <- rbind(CM, as.numeric(kindx == j) - as.numeric(kindx == 
                  i))
                #for(i in 1:length(CM)) { if(CM[i]==-1){CM[i]<-1} }
                rnames <- c(rnames, paste(levelnames[j], sepcompname, levelnames[i], 
                  sep = ""))
             }
         }
 } 
 


conf.int <- list()

if (length(grouplist) != ncol(CM))
 {stop("length of group list is not equal to length of contrasts")}


compnames <- numeric(length=nrow(CM))
lower <- numeric(length=nrow(CM))
upper <- numeric(length=nrow(CM)) 
estimate <- numeric(length=nrow(CM)) 


for (a in 1:nrow(CM))
  {Cvec <- as.logical(CM[a,])
  groupnames <- levelnames[Cvec]
  data <- datalist[Cvec]

  arglist$x <- data[[groupnames[2] ]]
  arglist$y <- data[[groupnames[1] ]]

  temp <- do.call(what=method, args=arglist)

  compnames[a] <- paste(groupnames[c(2,1)], collapse=sepcompname)
  lower[a] <- temp$conf.int[[1]]
  upper[a] <- temp$conf.int[[2]]
  estimate[a] <- temp$estimate

  }
}

#end all-pairs


# # # comparisons to a CONTROL:

if(is.null(control) == FALSE)
{

if(all(levelnames != control))
 {stop(paste("there is no control named",control, "among the levels of the grouping variable"))}

# extract the control group 
cgroup <- datalist[[control]]

# extract the list of treatments
treatlist <- datalist[-which(levelnames==control)]
treatnames <- levelnames[-which(levelnames==control)]


idmat <- diag(x=1, nrow=(k-1))

compnames <- numeric(length=length(treatlist))
lower <- numeric(length=length(treatlist))
upper <- numeric(length=length(treatlist)) 
estimate <- numeric(length=length(treatlist)) 



for (a in 1:(k-1))
  {Cvec <- as.logical(idmat[a,])
  groupnames <- treatnames[Cvec]

  arglist$x <- treatlist[Cvec][[1]]
  arglist$y <- cgroup

  temp <- do.call(what=method, args=arglist)

  compnames[a] <- paste(groupnames,sepcompname, control)
  lower[a] <- temp$conf.int[[1]]
  upper[a] <- temp$conf.int[[2]]
  estimate[a] <- temp$estimate
  
  }


}
# end comparisons to a CONTROL

names(estimate)<-compnames

names(lower)<-compnames

names(upper)<-compnames


return(
list(
estimate=estimate,
lower=lower,
upper=upper,
compnames=compnames,
method=attr(temp$conf.int, "methodname")
))
 
}

