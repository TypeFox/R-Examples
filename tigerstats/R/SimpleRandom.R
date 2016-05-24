#' @title Simple Random Sample 

#' @description An app to investigate the visual and numerical differences between a
#' sample and a population.  A sample is drawn from the input population and then a 
#' variable of choice is selected by the user.  If a categorical variable is chosen, the user sees a bar
#' chart with red bars designating the population and blue bars designating the sample.  
#' Simultaneously, a summary table (of percents) is output to the console for both the
#' population and the sample.  If a numerical variable is chose, the kernel density 
#' plot for the population is plotted in red and the histogram for each new sample 
#' is plotted in blue.  Simultaneously, the summary information for minimum, maximum, quartiles, 
#' median, mean, and standard deviation are output to the console for both the population
#' and the sample.  The size of the sample can be changed to explore how this affects 
#' statistics and the plots.
#' 
#' @rdname SimpleRandom
#' @usage SimpleRandom()
#' @return Graphical and numerical output
#' @export
#' @author Rebekah Robinson \email{rebekah_robinson@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' if (require(manipulate)) SimpleRandom()
#' }
SimpleRandom<-function(){
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  manipulate(
    n=slider(5,10000,initial=100,label="Sample Size n"),
    type=picker("sex", "math", "cappun", "income", "height", "idealheight", "diff", "kkardashtemp",label="Variable"),
{
  data=popsamp(n,imagpop)
  if(type=="sex"){
    poptab=rowPerc(xtabs(~sex,imagpop))
    samptab=rowPerc(xtabs(~sex,data))
    ymax=max(poptab[1]/100,poptab[2]/100)
    barplot(cbind(male=c(poptab[2]/100,samptab[2]/100),female=c(poptab[1]/100,samptab[1]/100)),
            beside=TRUE,col=c("red","lightblue"),main="Sex",ylim=c(0,ymax+0.2),legend=(c("Population","Sample")))
    cat("\n Sample Summary for sample size", n,"\n")
    print(round(samptab,2))
    cat("\n Population Summary \n")
    print(round(poptab,2))
  }
  if(type=="math"){
    poptab=rowPerc(xtabs(~math,imagpop))
    samptab=rowPerc(xtabs(~math,data))
    ymax=max(poptab[1]/100,poptab[2]/100)
    barplot(cbind(No=c(poptab[1]/100,samptab[1]/100),Yes=c(poptab[2]/100,samptab[2]/100)),
            beside=TRUE,col=c("red","lightblue"),main="Math Major",ylim=c(0,ymax+0.2),legend=(c("Population","Sample")))
    cat("\n Sample Summary for sample size", n,"\n")
    print(round(samptab,2))
    cat("\n Population Summary \n")
    print(round(poptab,2))
  }
  if(type=="income"){
    pop.dens = density(imagpop$income)
    max.dens = pop.dens$y[which.max(pop.dens$y)]
    ymax = 1.5 * max.dens
    max=imagpop$income[which.max(imagpop$income)]
    min=imagpop$income[which.min(imagpop$income)]
    hist(data$income,freq=FALSE,breaks=sqrt(n),col="lightblue",xlim = c(min(imagpop$income), max(imagpop$income)), 
         ylim = c(0, ymax),xlab="Income (to the nearest $100)", main="Income")
    lines(density(imagpop$income,from=min,to=max,bw=(max-min)/sqrt(10000)),col="red")
    pop.var=as.vector(imagpop$income)
    pop.quant=round(quantile(pop.var),2)
    pop.mean=round(mean(pop.var),2)
    pop.sd=round(sd(pop.var),2)
    pop.stats=data.frame(pop.quant[1],pop.quant[2],pop.quant[3],pop.quant[4],pop.quant[5],pop.mean,pop.sd)
    rownames(pop.stats)=""
    colnames(pop.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    samp.var=as.vector(data$income)
    samp.quant=round(quantile(samp.var),2)
    samp.mean=round(mean(samp.var),2)
    samp.sd=round(sd(samp.var),2)
    samp.stats=data.frame(samp.quant[1],samp.quant[2],samp.quant[3],samp.quant[4],samp.quant[5],samp.mean,samp.sd)
    rownames(samp.stats)=""
    colnames(samp.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    cat("\n Sample Summary for sample size", n,"\n")
    print(samp.stats)
    cat("\n Population Summary \n")
    print(pop.stats)
  }
  
  if(type=="cappun"){
    poptab=rowPerc(xtabs(~cappun,imagpop))
    samptab=rowPerc(xtabs(~cappun,data))
    ymax=max(poptab[1]/100,poptab[2]/100)
    barplot(cbind(Favor=c(poptab[1]/100,samptab[1]/100),Oppose=c(poptab[2]/100,samptab[2]/100)),
            beside=TRUE,col=c("red","lightblue"),main="Opinion on Capital Punishment",ylim=c(0,ymax+0.2),legend=(c("Population","Sample")))
    cat("\n Sample Summary for sample size", n,"\n")
    print(round(samptab,2))
    cat("\n Population Summary \n")
    print(round(poptab,2))
  }
  if(type=="height"){
    pop.dens = density(imagpop$height)
    max.dens = pop.dens$y[which.max(pop.dens$y)]
    min.dens = pop.dens$y[which.max(pop.dens$y)]
    ymax = 1.5 * max.dens
    max=imagpop$height[which.max(imagpop$height)]
    min=imagpop$height[which.min(imagpop$height)]
    hist(data$height,freq=FALSE,breaks=sqrt(n),col="lightblue",xlim = c(min(imagpop$height), max(imagpop$height)), 
         ylim = c(0, ymax),xlab="Height (in inches)", main="Height")
    lines(density(imagpop$height,from=min,to=max,bw=(max-min)/sqrt(10000)),col="red")
    pop.var=as.vector(imagpop$height)
    pop.quant=round(quantile(pop.var),2)
    pop.mean=round(mean(pop.var),2)
    pop.sd=round(sd(pop.var),2)
    pop.stats=data.frame(pop.quant[1],pop.quant[2],pop.quant[3],pop.quant[4],pop.quant[5],pop.mean,pop.sd)
    rownames(pop.stats)=""
    colnames(pop.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    samp.var=as.vector(data$height)
    samp.quant=round(quantile(samp.var),2)
    samp.mean=round(mean(samp.var),2)
    samp.sd=round(sd(samp.var),2)
    samp.stats=data.frame(samp.quant[1],samp.quant[2],samp.quant[3],samp.quant[4],samp.quant[5],samp.mean,samp.sd)
    rownames(samp.stats)=""
    colnames(samp.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    cat("\n Sample Summary for sample size", n,"\n")
    print(samp.stats)
    cat("\n Population Summary \n")
    print(pop.stats)
  }
  if(type=="idealheight"){
    pop.dens = density(imagpop$idealheight)
    max.dens = pop.dens$y[which.max(pop.dens$y)]
    ymax = 1.5 * max.dens
    max=imagpop$idealheight[which.max(imagpop$idealheight)]
    min=imagpop$idealheight[which.min(imagpop$idealheight)]
    hist(data$idealheight,freq=FALSE,col="lightblue",breaks=seq(min,max,by=1),xlim = c(min(imagpop$idealheight), max(imagpop$idealheight)), 
         ylim = c(0, ymax),xlab="Ideal Height (in inches)", main="Ideal Height")
    lines(density(imagpop$idealheight,from=min,to=max,bw=.5),col="red")
    pop.var=as.vector(imagpop$idealheight)
    pop.quant=round(quantile(pop.var),2)
    pop.mean=round(mean(pop.var),2)
    pop.sd=round(sd(pop.var),2)
    pop.stats=data.frame(pop.quant[1],pop.quant[2],pop.quant[3],pop.quant[4],pop.quant[5],pop.mean,pop.sd)
    rownames(pop.stats)=""
    colnames(pop.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    samp.var=as.vector(data$idealheight)
    samp.quant=round(quantile(samp.var),2)
    samp.mean=round(mean(samp.var),2)
    samp.sd=round(sd(samp.var),2)
    samp.stats=data.frame(samp.quant[1],samp.quant[2],samp.quant[3],samp.quant[4],samp.quant[5],samp.mean,samp.sd)
    rownames(samp.stats)=""
    colnames(samp.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    cat("\n Sample Summary for sample size", n,"\n")
    print(samp.stats)
    cat("\n Population Summary \n")
    print(pop.stats)
  }
  if(type=="diff"){
    pop.dens = density(imagpop$diff)
    max.dens = pop.dens$y[which.max(pop.dens$y)]
    ymax = 1.5 * max.dens
    max=imagpop$diff[which.max(imagpop$diff)]
    min=imagpop$diff[which.min(imagpop$diff)]
    hist(data$diff,freq=FALSE,col="lightblue",xlim = c(min(imagpop$diff), max(imagpop$diff)), 
         ylim = c(0, ymax),xlab="Difference (in inches)", breaks=seq(min,max,by=0.1), main="Difference between Ideal Height and Actual Height")
    lines(density(imagpop$diff,from=min,to=max,bw=.15),col="red")
    pop.var=as.vector(imagpop$diff)
    pop.quant=round(quantile(pop.var),2)
    pop.mean=round(mean(pop.var),2)
    pop.sd=round(sd(pop.var),2)
    pop.stats=data.frame(pop.quant[1],pop.quant[2],pop.quant[3],pop.quant[4],pop.quant[5],pop.mean,pop.sd)
    rownames(pop.stats)=""
    colnames(pop.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    samp.var=as.vector(data$diff)
    samp.quant=round(quantile(samp.var),2)
    samp.mean=round(mean(samp.var),2)
    samp.sd=round(sd(samp.var),2)
    samp.stats=data.frame(samp.quant[1],samp.quant[2],samp.quant[3],samp.quant[4],samp.quant[5],samp.mean,samp.sd)
    rownames(samp.stats)=""
    colnames(samp.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    cat("\n Sample Summary for sample size", n,"\n")
    print(samp.stats)
    cat("\n Population Summary \n")
    print(pop.stats)
  }
  if(type=="kkardashtemp"){
    pop.dens = density(imagpop$kkardashtemp)
    max.dens = pop.dens$y[which.max(pop.dens$y)]
    ymax = 3 * max.dens
    max=imagpop$kkardashtemp[which.max(imagpop$kkardashtemp)]
    min=imagpop$kkardashtemp[which.min(imagpop$kkardashtemp)]
    hist(data$kkardashtemp,freq=FALSE,breaks=sqrt(n),col="lightblue",xlim = c(min(imagpop$kkardashtemp)-10, max(imagpop$kkardashtemp)+10), 
         ylim = c(0, ymax),xlab="Feelings (0-100 scale)", main="Feelings about Kim Kardashian")
    lines(density(imagpop$kkardashtemp,from=0,to=100,bw=0.5),col="red")
    pop.var=as.vector(imagpop$kkardashtemp)
    pop.quant=round(quantile(pop.var),2)
    pop.mean=round(mean(pop.var),2)
    pop.sd=round(sd(pop.var),2)
    pop.stats=data.frame(pop.quant[1],pop.quant[2],pop.quant[3],pop.quant[4],pop.quant[5],pop.mean,pop.sd)
    rownames(pop.stats)=""
    colnames(pop.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    samp.var=as.vector(data$kkardashtemp)
    samp.quant=round(quantile(samp.var),2)
    samp.mean=round(mean(samp.var),2)
    samp.sd=round(sd(samp.var),2)
    samp.stats=data.frame(samp.quant[1],samp.quant[2],samp.quant[3],samp.quant[4],samp.quant[5],samp.mean,samp.sd)
    rownames(samp.stats)=""
    colnames(samp.stats)=c("min", "Q1", "median", "Q3", "max", "mean", "stdev")
    cat("\n Sample Summary for sample size", n,"\n")
    print(samp.stats)
    cat("\n Population Summary \n")
    print(pop.stats)
  }
}
  )}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("picker","type","imagpop"))
