#' @name zooaRch-package
#' @docType package
#' @aliases zooaRch
#' @title Analytical Tools for Zooarchaeological Data
#' @author Erik Otarola-Castillo, Jesse Wolfhagen, and Max D. Price
#'
#' @description Functions in this package allow users to import, manipulate, analyze and visualize
#' zooarchaeological data - the faunal remains recovered from archaeological sites. 
NULL

#' @import stats
NULL

#' @import graphics
NULL

#' @import ggplot2
NULL

#' Survival data from Marj Rabba, using Payne's (1973) age classes
#'
#' @name marjRab
#' @docType data
#' @author M.D. Price, M. Buckley, Y.M. Rowan, and M. Kersel
#' @references Payne, S. 1973  Kill - off Patterns in Sheep and Goats: The Mandibles from Asvan Kale. Anatolian Studies 23:281 - 303.
#' @references Price, M.D., Buckley, M., Rowan, Y.M., Kersel, M., 2013. Animal Management Strategies during the Chalcolithic in the Lower Galilee: New Data from Marj Rabba, Paleorient 39, 183 - 200.
#' @keywords datasets
NULL

#' Fusion Survival data from Marj Rabba
#'
#' @docType data
#' @author M.D. Price, M. Buckley, Y.M. Rowan, and M. Kersel
#' @references Price, M.D., Buckley, M., Rowan, Y.M., Kersel, M., 2013. Animal Management Strategies during the Chalcolithic in the Lower Galilee: New Data from Marj Rabba, Paleorient 39, 183 - 200.
#' @keywords datasets
#' @name marjRab.fuse
NULL

#' Fusion Survival data for cattle remains from the 
#' Winslow site, a colonial period farm near Boston, MA.
#'
#' @name winslow.fuse
#' @docType data
#' @author D. Landon
#' @references Landon, David B. 1993	Feeding Colonial Boston: A Zooarchaeological Study. Historical Archaeology 30:i-vii, 1-153
#' @keywords datasets
NULL

#' Bison survival data from Speth 1983
#'
#' @name speth83
#' @docType data
#' @author John D. Speth
#' @references Speth, J. D. 1983  Bison Kills and Bone Counts: 
#' Decision Making by Ancient Hunters. University of Chicago Press, London.
#' @keywords datasets
NULL

#' General survival analysis
#'
#' A general function to perform survival analysis of zooarchaeological data
#'
#' The function constructs Kaplan-Meier Estimator (KME) Confidence Intervals
#' Using Dental Eruption Wear Data
#'
#' @param SurviveData This function inputs datasets composed of three columns. 
#' The first column denotes the genus; the second is the age class (this MUST be numeric) 
#' if data contains nominal age classes (e.g., 'A', 'B', 'C', etc.) these data must be converted
#' to numbers (e.g., A = 1, B = 2, etc.).
#' @param labels Character value indicating wether age class labels wishing to be displayed.
#' @param models A numerical value (1-5) indiacting the models to compare the data to. Currently surv.func
#' makes use of 5 survival models: 1) Security (ref); 2) Milk (ref); 3) Wool (ref); 
#' 4) Catastrophic (Stiner 1990); and 5) Attritional (Stiner 1990). More models will be added soon. 
#' An option to include user's own model will also be available.
#' @param ci Numerical value indicating desired CI level (e.g., 90, 95, 99). Defaults to 95.
#' @param usermod numeric list (see help(list)) user-specified survivorship model. See example 3 below. 
#' Data must be entered as a list, else user will receive error.
#' @param plot A logical value indicating wether user wishes an output plot. Default = TRUE.
#' @param iter A numeric value indicating the number of bootstrap iterations. Defaults to 1000.
#' @keywords analysis
#' @export
#' @author Jesse WolfHagen and Erik Otarola-Castillo.
#' @return Function returns a matrix with the following components 
#'   \item{Lower and Upper CI}{typically the 97.5 and 2.5 percentile markers}
#'   \item{Point Value}{the y value on the survivorship curve}
#' @references Klein, R.G., Cruz-Uribe, K., 1983. The Analysis of Animal Bones from Archaeological Sites, University of Chicago Press, Chicago.
#' @references Stiner, M. C. 1990	The Use of Mortality Patterns in Archaeological Studies 
#' of Homonid Predatory Adaptations. Journal of Anthropological Archaeology 9:305 - 351.
#' @references Lyman, R.L., 1994. Vertebrate Taphonomy, Cambridge University Press, Cambridge.
#' @references Zeder, M.A., 2006. Reconciling Rates of Long Bone Fusion and Tooth Eruption in Sheep (Ovis) and Goat (Capra), in: Ruscillo, D. (Ed.), Recent Advances in Ageing and Sexing Animal Bones, Oxbow Books, Oxford. 
#' @references Twiss, K.C., 2008. An Assessment of the Archaeological Applicability of Faunal Ageing Methods Based on Dental Wear, International Journal of Osteoarchaeology 18, 329-351.
#' @examples
#' # Example 1: Payne 1973
#'  data(marjRab)
#' 
#' # Inspect data structure  
#'  head(marjRab)
#'  
#' # create age-class labels: Payne 1973 uses ageclasses A-I
#'  Labels <-c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')
#'  surv.func(SurviveData=marjRab,labels=Labels, models=1:3, ci=95, plot=TRUE, iter=1000)
#'  
#' # Example 2: Garnsey Site Bison Data (from Speth 1983)
#'  data(speth83)
#' 
#' # Inspect data structure  
#'  head(speth83)
#' 
#' # create age-class labels using the 13 age classes of Speth's (1983) scheme.
#'  Labels <-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
#' 
#' # Use the catastrophic and attritional mortality curves (after Stiner 1990).
#'  surv.func(SurviveData=speth83,labels=Labels, models=4:5, ci=95, plot=TRUE, iter=1000)
#'  
#' # Example 3: marjRab, input user-defined survivorship models.
#' data(marjRab)
#' # extract age classes from marjRab
#' age<-unique(marjRab$Ageclass)
#' age
#' # model survivorship using an exponential decay function 
#' # with parameter b: survivorship = age^(-1/b)
#' # surv 1: b= .95 
#' surv1<-c(1.00, 0.48, 0.31, 0.23, 0.18, 0.15, 0.12, 0.11, 0.09)
#' # surv 2: b= .35
#' surv2<-c(1.00, 0.13, 0.04, 0.01, 0.01, 0.005, 0.003, 0.002, 0.001)
#'
#' plot(age, surv1,type='l',xlim=range(age),ylim=c(0,1))
#' lines(age, surv2,col='red',)
#' 
#' # usermods in surv.func must be a list (if not a list, then user will receive error message)
#' mods<-list(surv1=surv1,surv2=surv2)
#'
#' surv.func(marjRab,models=NULL,usermod=mods)
surv.func<-function(SurviveData,labels=NULL, models=NULL, ci=95, plot=TRUE, iter=1000,usermod=NULL){
  data<-SurviveData$Ageclass 
  N.ages<-max(data)
  if(!is.null(labels)){
    Labels.ageclass <-labels
  }  
  if(is.null(labels)){
    Labels.ageclass <-1:N.ages
  }  
  if(!is.null(usermod) & !is.list(usermod)){
    stop("usmod MUST be entered as a list")
  }  
  survivorcurve.Eq4 <- function(data){
    vector <- rep(1, N.ages+1)   
    for(i in 1:N.ages+1) {
      vector[i] <- vector[i-1]*(sum(data >= (i-1)) - sum(data == (i-1)))/sum(data >= (i-1)) 
    }
    vector[is.na(vector)] <- 0
    round(vector,4) 
  }
  survivorcurve.Eq4(data)
  survive.matrix <- matrix(NA, ncol = N.ages+1, nrow = iter)
  survive.matrix[1,] <- survivorcurve.Eq4(data)
  for(i in 2:iter){
    bootstrap <- sample(1:length(data), length(data), replace = TRUE)
    survive.matrix[i,] <- survivorcurve.Eq4(data[bootstrap])
  }  
  ci<-ci/100
  upCI<-apply(survive.matrix[,-1], MARGIN = 2, FUN = quantile, 
              probs = ci+((1-ci)/2))
  loCI<-apply(survive.matrix[,-1], MARGIN = 2, FUN = quantile, 
              probs = ((1-ci)/2))
  if (plot==TRUE){
    plot(x = (1:N.ages), y = survive.matrix[1,-1], type = "n", xlab = "Age Class", 
         ylab = "Proportion Survived", axes = FALSE, ylim = c(0,1))
    polygon(c(1:N.ages,rev(1:N.ages)),c(loCI,rev(upCI)),col="light gray",border=NA)
    ###         Format the axes
    axis(side = 1, at = 1:N.ages, labels = Labels.ageclass)
    axis(side = 2)   
    lines(x = (1:N.ages), y = survive.matrix[1,-1], lwd = 2)
    ###         Plot the lower confidence interval (2.5%)
    lines(x = 1:N.ages, y = loCI, lty = "dashed", ylim = c(0,1))    
    ###         Plot the upper confidence interval (97.5%)
    lines(x = 1:N.ages, y = upCI, lty = "dashed", ylim = c(0,1))
    
    if(!is.null(models) & is.null(usermod)){
      surv.model<-list(Security = c(.904, .904, .645, .38, .25, .239, .171, .118, 0),
                       Milk = c(.47, .42, .39, .35, .28, .23, .19, .1, 0),
                       Wool = c(.85, .75, .65, .63, .57, .50, .43, .20, 0),
                       Catastrophic = c(.86, .73, .60, .49, .39, .31, .23, .16, .11, .07, .03, .01, 0),
                       Attritional = c(.76, .53, .48, .46, .45, .42, .39, .34, .29, .20, .11, .03, 0))
      for(j in 1:length(models)){ 
        i<-models[j]
        lines(x = 1:N.ages, y = surv.model[[i]], col = j+1, lwd = 2, lty = j)
      }    
      legend(x = "topright", cex = .75, lwd = 2,lty = c(1,2,1:length(models)), 
             col = c(1,1,2:(length(models)+1)),
             legend = c("Survivorship", paste(ci,"% Confidence Interval",sep=""), names(surv.model)[models]))
    }    
    if(!is.null(usermod)){
      surv.model<-usermod
      models<-1:length(surv.model)
      for(j in 1:length(models)){ 
        i<-models[j]
        lines(x = 1:N.ages, y = surv.model[[i]], col = j+1, lwd = 2, lty = j)
      }    
      legend(x = "topright", cex = .75, lwd = 2,lty = c(1,2,1:length(models)), 
             col = c(1,1,2:(length(models)+1)),
             legend = c("Survivorship", paste(ci,"% Confidence Interval",sep=""), names(surv.model)[models]))
    } 
  }  
  data.UpperCI <- upCI
  data.PointValue <- survive.matrix[1,-1]
  data.LowerCI <- loCI
  Taxon <- rep(unique(SurviveData$Genus), times = N.ages)
  Output.Matrix <- data.frame(Taxon = Taxon, AgeClassLabs = Labels.ageclass, LowerCI = data.LowerCI, 
                              PointValue = data.PointValue, UpperCI = data.UpperCI)
  return(Output.Matrix)
}

#' Epiphyseal fusion survival analysis
#'
#' A general function to perform survival analysis of zooarchaeological epiphyseal fusion data.
#'
#' The function constructs Confidence Intervals based off bootstraps of percent Fused values
#'
#' @param data This function inputs a dataframe composed of three columns, names must be 
#' 'Identification', 'Element', 'Fusion'. The first column denotes the arbitrary ID # and can be left blank if desired;
#' the second is the element name (differentiate proximal and distal as needed);
#' the third is the state of fusion. It must be 'Fused', 'Fusing', or 'Unfused' 
#' (NOTE: elements coded as 'Fusing' will be counted as 'Fused').
#' @param ci Numerical value indicating desired CI level (e.g., 90, 95, 99). Defaults to 95.
#' @param plotci A logical value indicating wether user wishes an output plot. Default = TRUE.
#' @param iter A numeric value indicating the number of bootstrap iterations. Defaults to 100.
#' @param plot.title A character value providing a title for the plot. Default is NULL.
#' @keywords analysis
#' @export
#' @author Jesse Wolfhagen, Max Price, and Erik Otarola-Castillo.
#' @return Function returns a matrix with the following components 
#'   \item{Lower and Upper CI}{typically the 97.5 and 2.5 percentile markers}
#'   \item{Point Value}{the y value on the percent Fused survivorship curve} 
#' @references Klein, R.G., Cruz-Uribe, K., 1983. The Analysis of Animal Bones from Archaeological Sites, University of Chicago Press, Chicago.
#' @references Lyman, R.L., 1994. Vertebrate Taphonomy, Cambridge University Press, Cambridge.
#' @references Zeder, M.A., 2006. Reconciling Rates of Long Bone Fusion and Tooth Eruption in Sheep (Ovis) and Goat (Capra), in: Ruscillo, D. (Ed.), Recent Advances in Ageing and Sexing Animal Bones, Oxbow Books, Oxford. 
#' @references Twiss, K.C., 2008. An Assessment of the Archaeological Applicability of Faunal Ageing Methods Based on Dental Wear, International Journal of Osteoarchaeology 18, 329-351.
#' @references Price, M.D., Buckley, M., Rowan, Y.M., Kersel, M., 2013. Animal Management Strategies during the Chalcolithic in the Lower Galilee: New Data from Marj Rabba, Paleorient 39, 183-200. 
#' @examples
#' # Example 1
#' # fusedat<-data(marjRab.fuse)
#' # test<-fuse.func(fusedat, iter=100, plotci=TRUE, plot.title="Fusion Example") 
#' # send the following into the console as you are prompted
#' # interactively
#' # 5
#' # 2
#' # 1
#' # 1
#' # 1
#' # 1
#' # Px.Humerus
#' # Ds.Humerus
#' # Calcaneus
#' # Ds.Tibia
#' # Px.Femur
#' # Phalanx1
fuse.func<-function(data,iter=1000,ci=95,plotci=TRUE,plot.title=NULL){
  cat(paste("Enter number of fusion groups"), "\n")
  ans<-readLines(n = 1)
  ans <- as.numeric(ans)
  fu.grps<-LETTERS[1:ans]
  ske.n<-numeric(length(fu.grps))
  for(i in 1:length(ske.n)){
    cat(paste("Enter number of skeletal elements for fusion group"), fu.grps[i],"\n")
    ans<-readLines(n = 1)
    ske.n[i] <- as.numeric(ans)
  }
  ele.list<-as.list(rep(NA,length(ske.n)))
  names(ele.list)<-fu.grps
  for(i in 1: length(ske.n)){
    cat(paste("Enter the", ske.n[i],"names of skeletal elements for fusion group"), fu.grps[i],"then press enter","\n")
    ele.list[[i]]<-readLines(n = ske.n[i])
  }  
  pctfuse<-function(dat){
    pct.ufu<-n<-numeric(length(ele.list))
    names(pct.ufu)<-fu.grps
    wh<-function(it){which(dat$Element==ele.list[[i]][it])} 
    for (i in 1:length(pct.ufu)){
      tab<-table(dat$Fusion[unlist(lapply(1:ske.n[i],wh))])
      fu<-tab[which(names(table(dat$Fusion[unlist(lapply(1:ske.n[i],wh))] ))=="Fused")]
      ufu<-tab[which(names(table(dat$Fusion[unlist(lapply(1:ske.n[i],wh))] ))=="Unfused")]
      if (is.nan(fu/(fu+ufu))){
        pct.ufu[i]<-0} else {pct.ufu[i]<-fu/(fu+ufu)}
      n[i]<-(fu+ufu)    
    }
    return(list(pct.ufu,n))
  }  
  boot <- matrix(NA, ncol = length(ele.list), nrow = iter) 
  boot[1,] <- pctfuse(data)[[1]]
  for(i in 2:iter){
    data.boot<-data[sample(1:dim(data)[1],dim(data)[1],replace=T),]
    boot[i,] <- pctfuse(data.boot)[[1]]
  }
  ci<-ci/100
  upCI<-ci+((1-ci)/2)
  loCI<-((1-ci)/2)
  quantilematrix <- matrix(NA, ncol = 2, nrow = length(fu.grps))
  for(i in 1:ncol(boot)){
    quantilematrix[i,] <- quantile(boot[,i], probs = c(loCI, upCI), na.rm = T)
  }
  Fusion.groups<-Data<-LowerCI<-UpperCI<-Count<-NULL
  outputtable <- data.frame(Fusion.groups = fu.grps, 
                            Data = round(boot[1,],2), 
                            LowerCI = round(quantilematrix[,1],2), UpperCI = round(quantilematrix[,2],2), 
                            Count = pctfuse(data)[[2]])
  
  ciplot<-ggplot(outputtable, aes(x = Fusion.groups, y = Data))+
    geom_point(size = 3)+
    geom_errorbar(aes(ymax = UpperCI, ymin = LowerCI), width = 0.2)+
    geom_text(aes(x = Fusion.groups, y = rep(1.05, length(Fusion.groups)), label = Count))+
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, color = "black"),
          axis.title = element_text(color = "black", size = 20),
          axis.text = element_text(color = "black", size = 15),
          axis.ticks = element_line(color = "black", size = 0.75),
          plot.title = element_text(color = "black", size = 24))+
    ylab("%Fused")+xlab("Fusion Group")+
    ggtitle(plot.title)
  if(plotci==TRUE){print(ciplot)}
  return(list(Output = outputtable, Bootstrap.Data = boot))
}

#' Analysis of Mortality Profiles
#'
#' This is a function used to conduct mortality analyses of zooarchaeological data
#'
#' This function plots mortality profiles, along with Confidence Intervals
#' using dental eruption and wear data. Optionally, plotted mortality profiles can
#' be compared to idealized models of mortality.
#'
#' @param mortData is the age-at-death dataset. This function inputs datasets composed of three columns. 
#' The first column denotes the genus; the second is the age class (this MUST be numeric) 
#' if data contains nominal age classes (e.g., 'A', 'B', 'C', etc.) these data must be converted
#' to numbers (e.g., A = 1, B = 2, etc.).
#' @param labels Character value indicating wether age class labels wishing to be displayed.
#' @param models A numerical value (1-5) indiacting the models to compare the data to. Currently mort.func
#' makes use of 5 mortality models: 1) Security (ref); 2) Milk (ref); 3) Wool (ref); 
#' 4) Catastrophic (Stiner 1990); and 5) Attritional (Stiner 1990). More models will be added soon. 
#' An option to include user's own model will also be available.
#' @param ci Numerical value indicating desired CI level (e.g., 90, 95, 99). Defaults to 95.
#' @param usermod numeric list (see help(list)) user-specified mortality model. See example 3 below. 
#' Data must be entered as a list, else user will receive error.
#' @param plot A logical value indicating wether user wishes an output plot. Default = TRUE.
#' @param iter A numeric value indicating the number of bootstrap iterations. Defaults to 1000.
#' @param lsize A numeric value indicating the vertical spacing value in a legend.
#' @keywords analysis
#' @export
#' @author Erik Otarola-Castillo.
#' @return Function returns a matrix with the following components 
#'   \item{Lower and Upper CI}{typically the 97.5 and 2.5 percentile markers}
#'   \item{Point Value}{the y value on the mortality profile}
#' @references Klein, R.G., Cruz-Uribe, K., 1983. The Analysis of Animal Bones from Archaeological Sites, University of Chicago Press, Chicago.
#' @references Stiner, M. C. 1990  The Use of Mortality Patterns in Archaeological Studies 
#' of Hominid Predatory Adaptations. Journal of Anthropological Archaeology 9:305 - 351.
#' @references Lyman, R.L., 1994. Vertebrate Taphonomy, Cambridge University Press, Cambridge.
#' @references Voorhies, M. R., 1969  Taphonomy and Population Dynamics of an Early Pliocene Vertebrate Fauna, Knox County, Nebraska. University of Wymonig Press. Contributions to Geology, Special Paper No. 1, Laramie (WY).
#' @references Reitz, E. and E. Wing 2008	Zooarchaeology. Second Edition. Cambridge University Press, Cambridge.
#' @examples
#' # Example 1: Payne 1973
#'  data(marjRab)
#' 
#' # Inspect data structure  
#'  head(marjRab)
#'  
#' # create age-class labels: Payne 1973 uses ageclasses A-I
#'  Labels <-c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')
#'  mort.func(mortData=marjRab,labels=Labels, models=1:3, ci=95, plot=TRUE, iter=1000)
#'  
#' # Example 2: Garnsey Site Bison Data (from Speth 1983)
#'  data(speth83)
#' 
#' # Inspect data structure  
#'  head(speth83)
#' 
#' # create age-class labels using the 13 age classes of Speth's (1983) scheme.
#'  Labels <-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
#' 
#' # Use the catastrophic and attritional mortality curves (after Stiner 1990).
#'  mort.func(mortData=speth83,labels=Labels, models=4:5, ci=95, plot=TRUE, iter=1000)
mort.func<-function(mortData,labels=NULL, models=NULL, ci=95, plot=TRUE, iter=1000,usermod=NULL,lsize=.1){
  data<-mortData$Ageclass 
  N.ages<-max(data)
  if(!is.null(labels)){
    Labels.ageclass <-labels
  }  
  if(is.null(labels)){
    Labels.ageclass <-1:N.ages
  }  
  if(!is.null(usermod) & !is.list(usermod)){
    stop("usmod MUST be entered as a list")
  }  
  mortprof <- function(data){
    vector <- rep(NA, N.ages)   
    for(i in 1:N.ages) {
      vector[i] <- sum(data==i)/length(data)
    }
    vector[is.na(vector)] <- 0
    round(vector,4) 
  }  
  mortprof(data)
  mortality.matrix <- matrix(NA, ncol = N.ages, nrow = iter)
  mortality.matrix[1,] <- mortprof(data)
  for(i in 2:iter){
    bootstrap <- sample(1:length(data), length(data), replace = TRUE)
    mortality.matrix[i,] <- unlist(mortprof(data[bootstrap]))
  }  
  ci<-ci/100  
  upCI<-apply(mortality.matrix[,], MARGIN = 2, FUN = quantile, 
              probs = ci+((1-ci)/2))
  loCI<-lo<-apply(mortality.matrix[,], MARGIN = 2, FUN = quantile, 
                  probs = ((1-ci)/2))
  
  if (plot==TRUE){
    ylab="Mortality"
    xlab = "Age Class"
    if(is.null(models) & is.null(usermod)){
      bar<-barplot(mortality.matrix[1,],ylim=c(0,(max(upCI)+.1)),
                  names=Labels.ageclass,ylab=ylab,
                  xlab=xlab,beside=T)
      g<-(max(bar)-min(bar))/110
      for (i in 1:length(bar))         {
        lines(c(bar[i],bar[i]),c(upCI[i],loCI[i]))
        lines(c(bar[i]-g,bar[i]+g),c(upCI[i],upCI[i]))
        lines(c(bar[i]-g,bar[i]+g),c(loCI[i],loCI[i]))
      }
    }
    
    if(!is.null(models) & is.null(usermod)){      
      mort.model<-list(Security = c(.096, 0, .259, .265, .13, .011, .068, .053, .118),
                       Milk = c(.53, .05, .03, .04, .07, .05, .04, .09, .1),
                       Wool = c(.15, .10, .10, .02, .06, .07, .07, .23, .20),
                       Catastrophic = c(.14, .13, .13, .11, .10, .08, .08, .07, .05, .04, .04, .02, .01),
                       Attritional = c(.24, .23, .05, .02, .01, .03, .03, .05, .05, .09, .09, .08, .03))
      dat<-rbind(mortality.matrix[1,], 
                 do.call(rbind,lapply(mort.model[models],matrix,
                                      ncol=length(mortality.matrix[1,]),
                                      byrow=TRUE)))      
      par(mfrow=c(length(models),1))
      ylab="Mortality"
      xlab = "Age Class"      
      for(j in 2:dim(dat)[1]){
        bar<-barplot(rbind(mortality.matrix[1,],dat[j,]),ylim=c(0,(max(c(upCI,dat))+.1)),names=Labels.ageclass,ylab=ylab,
                    xlab=xlab,beside=T,col=c("gray","black") )
        bar<-bar[1,]
        g<-(max(bar)-min(bar))/110        
        for (i in 1:length(bar))         {
          lines(c(bar[i],bar[i]),c(upCI[i],loCI[i]))
          lines(c(bar[i]-g,bar[i]+g),c(upCI[i],upCI[i]))
          lines(c(bar[i]-g,bar[i]+g),c(loCI[i],loCI[i]))
        } 
        #lines(x = bar, y = dat[j,],  lwd = 2,col="red")
        legend(x = "topright",
               fill=c(NA,"gray","black"),border=c(NA,1,1),lty=c(1,0,0),
               ncol=3,cex=.75,y.intersp=lsize,
               legend = c(paste(ci*100,"% CI",sep=""), "Observed",names(mort.model)[models][j-1]))
      }
      
      par(mfrow=c(1,1))
      
    }    
    if(!is.null(usermod)){
      mort.model<-usermod
      models<-1:length(mort.model)
      dat<-rbind(mortality.matrix[1,], 
                 do.call(rbind,lapply(mort.model[models],matrix,
                                      ncol=length(mortality.matrix[1,]),
                                      byrow=TRUE)))      
      par(mfrow=c(length(models),1))
      ylab="Mortality"
      xlab = "Age Class"      
      for(j in 2:dim(dat)[1]){
        bar<-barplot(rbind(mortality.matrix[1,],dat[j,]),ylim=c(0,(max(c(upCI,dat))+.1)),names=Labels.ageclass,ylab=ylab,
                    xlab=xlab,beside=T,col=c("gray","black") )
        bar<-bar[1,]
        g<-(max(bar)-min(bar))/110        
        for (i in 1:length(bar))         {
          lines(c(bar[i],bar[i]),c(upCI[i],loCI[i]))
          lines(c(bar[i]-g,bar[i]+g),c(upCI[i],upCI[i]))
          lines(c(bar[i]-g,bar[i]+g),c(loCI[i],loCI[i]))
        } 
        #lines(x = bar, y = dat[j,],  lwd = 2,col="red")
        legend(x = "topright",
               fill=c(NA,"gray","black"),border=c(NA,1,1),lty=c(1,0,0),
               ncol=3,cex=.75,y.intersp=lsize,
               legend = c(paste(ci*100,"% CI",sep=""), "Observed",names(mort.model)[models][j-1]))
      }
      
      par(mfrow=c(1,1))
    } 
  }  
  data.LowerCI <-loCI
  data.PointValue <- mortality.matrix[1,]
  data.UpperCI <- upCI
  Taxon <- rep(unique(mortData$Genus), times = N.ages)
  Output.Matrix <- data.frame(Taxon = Taxon, AgeClassLabs = Labels.ageclass, LowerCI = data.LowerCI, 
                              PointValue = data.PointValue, UpperCI = data.UpperCI)
  return(Output.Matrix)
}