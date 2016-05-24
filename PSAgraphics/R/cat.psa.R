cat.psa<-function(categorical, treatment = NULL, strata = NULL,
    catnames = NULL, catcol = "terrain.colors", width = .25, barlab = c("A","B"), barnames = NULL,
    rtmar = 1.5, balance = FALSE, B = 1000, tbl = TRUE, cex.leg = 1, ...)
{

#categorical should be a categorical variable with levels 1:n
#treatment should be the binary treatment variable from a PSA
#strata should be the numerical strata that subjects belong to (from a logistic PSA, usually 1:5)
#width will be width of bars

#If "balance" is TRUE, then two different heuristic bootstrap distributions 
#and balance statistics are found. First a histogram of balance measures is 
#given from randomly generated strata, and compared with the measure for the 
#actual strata. Second, within strata, the covariate balance between treatments 
#is bootstrapped and measured (in a standard way[I think - need to talk to Bob]), 
#and then a "p-value" is calculated and given in the bar chart at the bottom of 
#each stratum in red.

#Ok, more on the balance statistics calculated:  For the first, randomly 
#generated strata are generated. In a given stratum, proportions in each 
#treatment are calculated, and then the absolute difference is between 
#these proportions in each category is found, and summed across categories.
#This value is summed across strata for the final balance statistic.  The
#same statistic from the original is calculated and compared to the 
#randomly generated values via a histogram and ranking.  In the second 
#measure, the strata are considered fixed, and the covariate values in a 
#given statum are randomly sampled without replacement and assigned to 
#the two treatments.  Absolute differences in size of categories are 
#found, and summed across categories, to be compared with the true 
#covariate distribution. 

#If "categorical" has three columns, treat as c, t, s.
if(dim(as.data.frame(categorical))[2]==3){ treatment   <- categorical[,2]
                                           strata      <- categorical[,3]
                                           categorical <- categorical[,1]
                                         }
cat.f<-as.factor(categorical)
cat.levels<-levels(cat.f)

 table.cts<-table(categorical,treatment,strata)
 cat.dim<-dim(table.cts)[1]
 strata.dim<-dim(table.cts)[3]
 if(is.null(catnames)&!is.factor(categorical)){catnames<-sort(unique(categorical))}
 if(is.null(catnames)){if(is.factor(categorical)){catnames<-levels(categorical)}else{catnames<-c(1:cat.dim)}}
 if(catcol[1]=="terrain.colors"){catcol<-terrain.colors(cat.dim)}

#Creating the x coordinates of the rectangles
 minlim<-1:strata.dim-width
 midlim<-1:strata.dim
 maxlim<-1:strata.dim+width

x.left<-as.vector(matrix(rep(minlim,cat.dim),nrow=cat.dim,byrow=TRUE))
x.mid<-as.vector(matrix(rep(midlim,cat.dim),nrow=cat.dim,byrow=TRUE))
x.right<-as.vector(matrix(rep(maxlim,cat.dim),nrow=cat.dim,byrow=TRUE))

#x.l/x.r have left x coordinates of all 0/1 rectangles for each stratum
x.l<-c(x.left,x.mid)
x.r<-c(x.mid,x.right)

#Want to create the y coordinates of the rectangles.
#Finding the numbers of 0/1 by stratum in sum.strata.01

sum.strata.01<-NULL
for(j in 1:strata.dim) sum.strata.01<-rbind(sum.strata.01,apply(table.cts[,,j],2,sum))

#Calculating the vectors of bottom and top y values for the rectangles in each treatment
 y.b.0<-NULL
 y.t.0<-NULL
 y.b.1<-NULL
 y.t.1<-NULL
 prop.zo<-NULL
  for(j in 1:strata.dim) {
    # prop.i is a vector of the proportion of 0/1 by categorical level
    # level for stratum j
        prop.0<-table.cts[,1,j]/sum.strata.01[j,1]
        prop.1<-table.cts[,2,j]/sum.strata.01[j,2]
        C<-prop.0
        T<-prop.1
        prop.zo<-cbind(prop.zo,C)
        prop.zo<-cbind(prop.zo,T)       
    # prop.i.sum is a vector of the cumulative proportions
    prop.0.sum<-NULL
        prop.1.sum<-NULL
    for(i in 1:cat.dim) prop.0.sum<-cbind(prop.0.sum,sum(prop.0[1:i]))
        for(i in 1:cat.dim) prop.1.sum<-cbind(prop.1.sum,sum(prop.1[1:i]))
    # The bottom y's should start at 0 and end at the penultimate entry of prop.i.sum.
    # The top y's should start with first entry of prop.i.sum and end at 1, which is last entry of prop.i.sum
    y.b.0<-c(y.b.0,0,prop.0.sum[-cat.dim])
    y.t.0<-c(y.t.0,prop.0.sum)
        y.b.1<-c(y.b.1,0,prop.1.sum[-cat.dim])
    y.t.1<-c(y.t.1,prop.1.sum)
                         }
 y.b<-c(y.b.0,y.b.1)
 y.t<-c(y.t.0,y.t.1)

if(balance){bal.cs<-bal.cs.psa(categorical,treatment,strata,B=B)
            cat("Histogram of Random Strata Balance. Press <enter> for next chart...")
            readline()
            bal.fe<-bal.fe.psa(categorical,treatment,strata)
           } 

#Creating Graph
 xlimits<-c(.5,strata.dim+rtmar)
 
 plot(c(.5,2), type="n", log = "", axes=FALSE, xlim = xlimits, ylim = c(-.05,1.05),...)
 rect(x.l,y.b,x.r,y.t,col=catcol,lty="solid")
 axis(1,at=1:strata.dim, labels=sort(unique(strata)))
 axis(2,at=c(0,.5,1)) 
 
 if(is.null(barnames)){barnames<-unlist(dimnames(table.cts)[2])}
 
 for(i in 1:strata.dim) {text  (i-.125,0, barlab[1], cex=.7, pos=1)}
 for(i in 1:strata.dim) {text  (i+.125,0, barlab[2], cex=.7, pos=1)}
 table.ts<-table(treatment,strata)
 for(i in 1:strata.dim) {text (i-width/2,1,table.ts[1,i],cex=.7,pos=3,col=4)}
 for(i in 1:strata.dim) {text (i+width/2,1,table.ts[2,i],cex=.7,pos=3,col=4)}
 if(balance){lnd<-c("No. Obs.","F.E. p-val",barnames)
             ppch<-c("#","#",barlab)
             ccol<-c(4,2,1,1)
             for (i in 1:strata.dim) {text  (i,-.067,round(bal.fe[i],2),col=2,cex=.7)}}else{
             lnd<-c("No. Obs.",barnames)
             ppch<-c("#",barlab)
             ccol<-c(4,1,1)}
 legend(strata.dim+.4,.3,legend=lnd,pch=ppch,col=ccol,bty="n",cex = cex.leg)
 legend(strata.dim+.4,.97,c("Levels",catnames),pch=c(0,rep(15,cat.dim)),col=c(0,catcol),bty="n", cex = cex.leg)
            
##Output a table of percents for catagorical levels within each treatment, stratum.
 if(tbl){
    Levels<-NULL
    numprop<-NULL
    tns<-NULL
    t.names<-unlist(dimnames(table.cts)[2])
    s.names<-unlist(dimnames(table.cts)[3])
    for(i in 1:strata.dim){
    for(j in 1:2){
    tns<-c(tns,paste(t.names[j],":",s.names[i],sep=""))
    }}
    strata.names<-unlist(dimnames(table.cts)[3])
    treatment.stratum.proportions<-round(prop.zo,3)
    colnames(treatment.stratum.proportions)<-tns
    out<-list(treatment.stratum.proportions)
    names(out)<-c(paste("treatment",":","stratum",".proportions",sep=""))
    return(out)
         }
}



