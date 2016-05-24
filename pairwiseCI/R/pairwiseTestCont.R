
"pairwiseTestCont" <-
function(formula, data,
 control=NULL, method, ...)

{
arglist<-list(...)


# # # check the arguments

# # # Split up into pairwise comparisons

mf <- model.frame(data, formula=formula)

if(!is.numeric(mf[,1]))
 {stop("response variable must be numeric")}

if( method=="Prop.test" && !is.matrix(mf[,1]) )
 {stop("for proportions, the response must be specified like 'cbind(success,failure)'") }

if(!is.factor(mf[,2]))
 {mf[,2]<-as.factor(mf[,2])}

grouplist <- split(x=mf, f=mf[,2], drop=TRUE  )
datalist <-split(x=mf[,1], f=mf[,2], drop=TRUE  )

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
                rnames <- c(rnames, paste(levelnames[j], "-", levelnames[i], 
                  sep = ""))
             }
         }
 } 


conf.int <- list()

if (length(grouplist) != ncol(CM))
 {stop("length of group list is not equal to length of contrasts")}

p.value <- numeric(length=nrow(CM))
compnames <- character(length=nrow(CM))
groupx <- character(length=nrow(CM))
groupy <- character(length=nrow(CM))

for (a in 1:nrow(CM))
  {Cvec <- as.logical(CM[a,])
  groupnames <- levelnames[Cvec]
  data <- datalist[Cvec]

  compnames[a] <- paste(groupnames[c(2,1)], collapse="-")
  groupx[a] <- groupnames[1]
  groupy[a] <- groupnames[2]

  arglist$x <- data[[groupnames[2] ]]
  arglist$y <- data[[groupnames[1] ]]

  temp <- do.call(what=method, args=arglist)
  p.value[a] <- temp$p.value

  }
}

#end all-pairs


# # # comparisons to a CONTROL:

if(is.null(control) == FALSE)
{

if(all(levelnames != control))
 {stop(paste("there is no control",control, "among the levels of the grouping variable"))}

# extract the control group 
cgroup <- datalist[[control]]

# extract the list of treatments
treatlist <- datalist[-which(levelnames==control)]
treatnames <- levelnames[-which(levelnames==control)]


idmat <- diag(x=1, nrow=(k-1))

groupx <- character(length=length(treatlist))
groupy <- character(length=length(treatlist))
compnames <- numeric(length=length(treatlist))
p.value <- numeric(length=length(treatlist))
 
for (a in 1:(k-1))
  {Cvec <- as.logical(idmat[a,])
  groupnames <- treatnames[Cvec]
  groupx[a] <- groupnames
  groupy[a] <- control
  compnames[a] <- paste(groupnames,"-", control)

  arglist$x <- treatlist[Cvec][[1]]
  arglist$y <- cgroup
  temp <- do.call(what=method, args=arglist)
  p.value[a] <- temp$p.value
  }


}
# end comparisons to a CONTROL

names(p.value)<-compnames


out<-data.frame(
p.value=p.value,
compnames=compnames,
groupx=groupx,
groupy=groupy
)
attr(out, "method")<-temp$method

return(out)
 
}

