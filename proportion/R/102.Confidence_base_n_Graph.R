#####################################################################################
#' Plots the CI estimation of exact methods
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P},
#' The input can also be a range of values between 0 and 1.
#' @details  The plot of Confidence Interval of \code{n} given \code{alp}.
#' @family Basic methods of CI estimation
#' @examples
#' \dontrun{
#' n=5; alp=0.05; e=0.5 #Mid-p
#' PlotciEX(n,alp,e)
#' n=5; alp=0.05;e=1 #Clopper-Pearson
#' PlotciEX(n,alp,e)
#' n=5; alp=0.05;e=c(0.05,0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' PlotciEX(n,alp,e)
#' }
#' @export
PlotciEX<-function(n,alp,e)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10 ) stop("Plot of only 10 intervals of 'e' is possible")
  Abberation=ID=Value=UEX=LEX=LABB=UABB=LowerLimit=UpperLimit=ZWI=NULL

  ss1=ciEX(n,alp,e)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)
  ss$e = as.factor(ss$e)

  ll=subset(ss, LABB=="YES")
  ul=subset(ss, UABB=="YES")
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
      ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UEX,y = ID,
                                           xmin = LEX,
                                           xmax = UEX,
                                           color= e),
                              size=0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Exact method") +
      ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red",
                                          "black", "orange","brown","chartreuse4",
                                           "blueviolet" , "deeppink", "darksalmon", "tan1" )) +
      ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red",
                                            "orange","chartreuse4",
                                            "blueviolet" , "deeppink", "darksalmon", "tan1")) +
      ggplot2::scale_shape_manual(values=c(21,22,23))                  # Change shapes
  }
  else {
    ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Exact method") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UEX,y = ID,
                                           xmin = LEX,
                                           xmax = UEX,
                                           color= e),
                              size = 0.5) +
      ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red",
                                          "black", "orange","brown","chartreuse4",
                                          "blueviolet" , "deeppink", "darksalmon", "tan1" )) +
      ggplot2::scale_colour_manual(values=c("red", "black", "blue", "cyan4", "orange",
                                            "deeppink","chartreuse4",
                                            "blueviolet" , "brown", "darksalmon", "tan1"))

  }

}

#####################################################################################
#' Plots the CI estimation of Bayesian method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Shape parameter 1 for prior Beta distribution in Bayesian model. Can also be a vector of priors.
#' @param b - Shape parameter 2 for prior Beta distribution in Bayesian model. Can also be a vector of priors.
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} using Bayesian method
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05; a=0.5;b=0.5;
#' PlotciBA(n,alp,a,b)
#' n=5; alp=0.05; a=c(0.5,2,1,1,2,0.5);b=c(0.5,2,1,1,2,0.5)
#' PlotciBA(n,alp,a,b)
#' @export
PlotciBA<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  if(length(a)>1 || length(b)>1){
    if ((class(a) != "integer") & (class(a) != "numeric") || any(a < 0)) stop("'a' has to be a set of positive numeric vectors")
    if ((class(b) != "integer") & (class(b) != "numeric") || any(b < 0)) stop("'b' has to be a set of positive numeric vectors")
    if (length(a) <  n || length(b) <  n ) stop("'a' and 'b' vectors have to be equal to length n")
    PlotciBAD(n,alp,a,b)
    }
  else{
    if ((class(a) != "integer") & (class(a) != "numeric") || a<0  ) stop("'a' has to be greater than or equal to 0")
    if ((class(b) != "integer") & (class(b) != "numeric") || b<0  ) stop("'b' has to be greater than or equal to 0")
    UB=LB=Abberation=ID=method=Value=LowerLimit=UpperLimit=ZWI=NULL

  ss1=ciBA(n,alp,a,b)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)
  dfq=ss[,c(1,4,5)]
  dfh=ss[,c(1,6,7)]
  df1=data.frame(ID=dfq$ID, LB=dfq$LBAQ, UB=dfq$UBAQ, method="Quantile")
  df2=data.frame(ID=dfh$ID, LB=dfh$LBAH, UB=dfh$UBAH, method="HPD")

  vs=rbind(df1,df2)

    ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= vs,
                              ggplot2::aes(x = UB,y = ID,
                                           xmin = LB,
                                           xmax = UB,
                                           color=method),
                              size=0.5)+
    ggplot2::labs(title = "Bayesian methods - Quantile and HPD") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::scale_colour_manual(values=c("blue", "red"))
  }
}
#######################################################################################
PlotciBAD<-function(n,alp,a,b)
{
  ss1=ciBAD(n,alp,a,b)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)
  dfq=ss[,c(1,4,5)]
  dfh=ss[,c(1,6,7)]
  df1=data.frame(ID=dfq$ID, LB=dfq$LBAQ, UB=dfq$UBAQ, method="Quantile")
  df2=data.frame(ID=dfh$ID, LB=dfh$LBAH, UB=dfh$UBAH, method="HPD")
  UB=LB=Abberation=ID=method=Value=LowerLimit=UpperLimit=ZWI=NULL

  vs=rbind(df1,df2)

  ggplot2::ggplot()+
    ggplot2::geom_errorbarh(data= vs,
                            ggplot2::aes(x = UB,y = ID,
                                         xmin = LB,
                                         xmax = UB,
                                         color=method),
                            size=0.5)+
    ggplot2::labs(title = "Bayesian methods  - Quantile and HPD") +
    ggplot2::labs(x = "Lower and Upper limits") +
    ggplot2::scale_colour_manual(values=c("blue", "red"))

}

#####################################################################################
#' Plots the CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) methods
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}.
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05;
#' PlotciAll(n,alp)
#' @export
PlotciAll<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  Abberation=ID=method=Value=LABB=UABB=LowerLimit=UpperLimit=ZWI=NULL

  ss1=ciAll(n,alp)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LABB=="YES")
  ul=subset(ss, UABB=="YES")
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
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit,
                                           color= method),
                              size=0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(title = "Basic methods  of CI estimation") +
      ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red", "black", "orange","brown")) +
      ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red", "orange")) +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(title = "Basic methods  of CI estimation") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit, color= method),
                              size=0.5)
  }
  oo
}


#####################################################################################
#' Plots the CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) grouped by x value
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}.
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05;
#' PlotciAllg(n,alp)
#' @export
PlotciAllg<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  Abberation=ID=method=Value=val1=val2=LABB=UABB=LowerLimit=UpperLimit=ZWI=NULL

  ss1=ciAll(n,alp)
  nss= ss1[order(ss1$x, (ss1$UpperLimit-ss1$LowerLimit)),]
  id=1:nrow(ss1)
  ss= data.frame(ID=id,nss)

  ll=subset(ss, LABB=="YES")
  ul=subset(ss, UABB=="YES")
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
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::labs(title = "Confidence interval for adjusted methods sorted by x") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit,
                                             color= method),
                                size=0.5)+
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
        ggplot2::labs(title = "Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit, color= method),
                                size=0.5)
    }
    oo
  }
  else {

    ff= data.frame(val1=seq(0.5,max(ss$ID),by=(max(ss$ID)/(max(ss$x)+1))),val2=(0:max(ss$x)))

    if(nrow(ldf)>0){
      oo= ggplot2::ggplot()+
        ggplot2::labs(title = "Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit,
                                             color= method),
                                size=0.5)+
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
        ggplot2::labs(title = "Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit, color= method),
                                size=0.5) +
        ggplot2::geom_hline(ggplot2::aes(yintercept=val1),data=ff) +
        ggplot2::geom_text(ggplot2::aes(0,val1,label = paste("x=", sep="", val2),hjust=1.1, vjust = -1), data=ff)
    }
    oo
  }
}

#####################################################################################
#' Plots the CI estimation of  Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} using Wald method
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05
#' PlotciWD(n,alp)
#' @export
PlotciWD<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  WaldCI.df = ciWD(n,alp)
  ss1 = data.frame( x=WaldCI.df$x, LowerLimit = WaldCI.df$LWD, UpperLimit = WaldCI.df$UWD, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)

  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")
  }
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

  if(nrow(ldf)>0)
  {
    ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Wald method")  +
      ggplot2::scale_shape_manual(values=c(21,22,23))                  # Change shapes
  }
  else {
    ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Wald method") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)
  }
}

#####################################################################################
#' Plots the CI estimation of ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} using ArcSine method
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05
#' PlotciAS(n,alp)
#' @export
PlotciAS<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ArcSineCI.df = ciAS(n,alp)
  ss1 = data.frame(x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LAS, UpperLimit = ArcSineCI.df$UAS, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)


  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")
  }
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

  if(nrow(ldf)>0)
  {
    ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5) +
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - ArcSine method")  +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - ArcSine method") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)
  }
}

#####################################################################################
#' Plots the CI estimation of Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} using Likelihood Ratio method
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05
#' PlotciLR(n,alp)
#' @export
PlotciLR<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

    LRCI.df = ciLR(n,alp)
  ss1 = data.frame( x=LRCI.df$x, LowerLimit = LRCI.df$LLR, UpperLimit = LRCI.df$ULR, LowerAbb = LRCI.df$LABB, UpperAbb = LRCI.df$UABB, ZWI = LRCI.df$ZWI)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")
  }
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

  if(nrow(ldf)>0)
  {
    ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Likelihood Ratio method")  +
      ggplot2::scale_shape_manual(values=c(21,22,23))                  # Change shapes
  }
  else {
    ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Likelihood Ratio method") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)
  }
}

#####################################################################################
#' Plots the CI estimation of Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} using Score method
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05
#' PlotciSC(n,alp)
#' @export
PlotciSC<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ScoreCI.df = ciSC(n,alp)
  ss1 = data.frame( x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LSC, UpperLimit = ScoreCI.df$USC, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)

  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")
  }
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

  if(nrow(ldf)>0)
  {
    ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Score method")  +
      ggplot2::scale_shape_manual(values=c(21,22,23))                  # Change shapes
  }
  else {
    ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Score method") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)
  }
}

#####################################################################################
#' Plots the CI estimation of Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} using Wald-T method
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05
#' PlotciTW(n,alp)
#' @export
PlotciTW<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  WaldTCI.df = ciTW(n,alp)
  ss1 = data.frame( x=WaldTCI.df$x, LowerLimit = WaldTCI.df$LTW, UpperLimit = WaldTCI.df$UTW, LowerAbb = WaldTCI.df$LABB, UpperAbb = WaldTCI.df$UABB, ZWI = WaldTCI.df$ZWI)

  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")
  }
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

  if(nrow(ldf)>0)
  {
    ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Wald-T method")  +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Wald-T method") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)
  }
}

#####################################################################################
#' Plots the CI estimation of  Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The plot of Confidence Interval of \code{n} given \code{alp} using Logit Wald method
#' @family Basic methods of CI estimation
#' @examples
#' n=5; alp=0.05
#' PlotciLT(n,alp)
#' @export
PlotciLT<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  LogitWald.df = ciLT(n,alp)
  ss1 = data.frame(x=LogitWald.df$x, LowerLimit = LogitWald.df$LLT, UpperLimit = LogitWald.df$ULT, LowerAbb = LogitWald.df$LABB, UpperAbb = LogitWald.df$UABB, ZWI = LogitWald.df$ZWI)

  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,3)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")
  }
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

  if(nrow(ldf)>0)
  {
    ggplot2::ggplot()+
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Logit Wald method")  +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::labs(title = "Confidence Interval - Logit Wald method") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size=0.5)
  }
}

