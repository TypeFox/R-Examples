#####################################################################################
#' Plots the CI estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of the Confidence Interval using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} and \code{h}
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05;h=2
#' PlotciAAll(n,alp,h)
#' @export
#8. Plot all methods
PlotciAAll<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  Abberation=ID=Value=method=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciAAll(n,alp,h)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,4)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,5)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,4)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted methods") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit,
                                           color= method),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red", "black", "orange","brown")) +
      ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red", "orange")) +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted methods") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit, color= method),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plots the CI estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) grouped by x value
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adjustment
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} and \code{h} grouped by x for 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05; h=2
#' PlotciAAllg(n,alp,h)
#' @export
#9.All methods plots with grouping
PlotciAAllg<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || is.integer(h)) stop("'h' has to be an integer greater than or equal to 0")
  Abberation=ID=Value=method=val1=val2=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciAAll(n,alp,h)
  nss= ss1[order(ss1$x, (ss1$UpperLimit-ss1$LowerLimit)),]
  id=1:nrow(ss1)
  ss= data.frame(ID=id,nss)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,4)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,5)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,4)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if((max(as.numeric(unique(ss$method)))-nrow(ss))==0){
    if(nrow(ldf)>0){
      oo= ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit,
                                             color= method),
                                size = 0.5)+
        ggplot2::geom_point(data=ldf,
                            ggplot2::aes(x=Value, y=ID,
                                         group = Abberation,shape=Abberation),
                            size = 4, fill = "red") +
        ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red", "black", "orange","brown")) +
        ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red", "orange")) +
        ggplot2::scale_shape_manual(values=c(21,22,23))
    }
    else {
      oo=  ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit, color= method),
                                size = 0.5)
    }
    oo
  }
  else {

    ff= data.frame(val1=seq(0.5,max(ss$ID),by=(max(ss$ID)/(max(ss$x)+1))),val2=(0:max(ss$x)))

    if(nrow(ldf)>0){
      oo= ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit,
                                             color= method),
                                size = 0.5)+
        ggplot2::geom_point(data=ldf,
                            ggplot2::aes(x=Value, y=ID,
                                         group = Abberation,shape=Abberation),
                            size = 4, fill = "red") +
        ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red", "black", "orange","brown")) +
        ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red", "orange")) +
        ggplot2::scale_shape_manual(values=c(21,22,23))    +
        ggplot2::geom_hline(ggplot2::aes(yintercept=val1),data=ff) +
        ggplot2::geom_text(ggplot2::aes(0,val1,label = paste("x=", sep="", val2),hjust=1.1, vjust = -1), data=ff)
    }
    else {
      oo=  ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit, color= method),
                                size = 0.5) +
        ggplot2::geom_hline(ggplot2::aes(yintercept=val1),data=ff) +
        ggplot2::geom_text(ggplot2::aes(0,val1,label = paste("x=", sep="", val2),hjust=1.1, vjust = -1), data=ff)
    }
    oo
  }
}

#####################################################################################
#' Plot of CI estimation of adjusted Wald
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of the Confidence Interval using adjusted Wald method for \code{n} given \code{alp} and \code{h}
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05;h=2
#' PlotciAWD(n,alp,h)
#' @export
#8. Plot all methods
PlotciAWD<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  WaldCI.df    = ciAWD(n,alp,h)
  ss1 = data.frame(x=WaldCI.df$x, LowerLimit = WaldCI.df$LAWD, UpperLimit = WaldCI.df$UAWD, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,4)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,3)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Wald method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Wald method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plot of CI estimation of adjusted ArcSine
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of the Confidence Interval using adjusted ArcSine method for \code{n} given \code{alp} and \code{h}
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05;h=2
#' PlotciAAS(n,alp,h)
#' @export
#8. Plot all methods
PlotciAAS<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ArcSineCI.df = ciAAS(n,alp,h)

  ss1= data.frame( x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LAAS, UpperLimit = ArcSineCI.df$UAAS, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,4)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,3)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted ArcSine method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted ArcSine method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plot of CI estimation of adjusted Likelihood Ratio
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of the Confidence Interval using adjusted Likelihood Ratio method for \code{n} given \code{alp} and \code{h}
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05;h=2
#' PlotciALR(n,alp,h)
#' @export
#8. Plot all methods
PlotciALR<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  LRCI.df      = ciALR(n,alp,round(h,0))		#h must be +ve integer

  ss1 = data.frame( x=LRCI.df$x, LowerLimit = LRCI.df$LALR, UpperLimit = LRCI.df$UALR, LowerAbb = LRCI.df$LABB, UpperAbb = LRCI.df$UABB, ZWI = LRCI.df$ZWI)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,4)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,3)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Likelihood Ratio method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Likelihood Ratio method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plot of CI estimation of adjusted Score
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of the Confidence Interval using adjusted Score method for \code{n} given \code{alp} and \code{h}
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05;h=2
#' PlotciASC(n,alp,h)
#' @export
PlotciASC<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ScoreCI.df   = ciASC(n,alp,h)

  ss1 = data.frame( x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LASC, UpperLimit = ScoreCI.df$UASC, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,4)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,3)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Score method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Score method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plot of CI estimation of adjusted Logit Wald
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of the Confidence Interval using adjusted Logit Wald method for \code{n} given \code{alp} and \code{h}
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05;h=2
#' PlotciALT(n,alp,h)
#' @export
#8. Plot all methods
PlotciALT<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  WaldLCI.df   = ciALT(n,alp,h)
  ss1 = data.frame( x=WaldLCI.df$x, LowerLimit = WaldLCI.df$LALT, UpperLimit = WaldLCI.df$UALT, LowerAbb = WaldLCI.df$LABB, UpperAbb = WaldLCI.df$UABB, ZWI = WaldLCI.df$ZWI)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,4)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,3)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Logit Wald method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Logit Wald method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plot of CI estimation of adjusted Wald-T
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of the Confidence Interval using adjusted Wald-T method for \code{n} given \code{alp} and \code{h}
#' @family Adjusted methods of CI estimation
#' @examples
#' n=5; alp=0.05;h=2
#' PlotciATW(n,alp,h)
#' @export
PlotciATW<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  AdWaldCI.df  = ciATW(n,alp,h)
  ss1 = data.frame( x=AdWaldCI.df$x, LowerLimit = AdWaldCI.df$LATW, UpperLimit = AdWaldCI.df$UATW, LowerAbb = AdWaldCI.df$LABB, UpperAbb = AdWaldCI.df$UABB, ZWI = AdWaldCI.df$ZWI)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,4)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,3)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Wald-T method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted Wald-T method") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

