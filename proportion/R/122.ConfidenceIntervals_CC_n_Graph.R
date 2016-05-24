#####################################################################################
#' Plots the CI estimation of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) given n, alp and c
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @family Continuity correction methods of CI estimation
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotciCAll(n,alp,c)
#' @export
#7. Plot all CC methods
PlotciCAll<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  Abberation=ID=Value=method=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciCAll(n,alp,c)
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
      ggplot2::ggtitle("Confidence interval for continuity corrected methods") +
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
      ggplot2::ggtitle("Confidence interval for continuity corrected methods") +
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
#' Plots the CI estimation of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) grouped by x value
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) grouped by x for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @family Continuity correction methods of CI estimation
#' @examples
#' n=5; alp=0.05; c=1/(2*n)
#' PlotciCAllg(n,alp,c)
#' @export
#8.All methods plots with grouping
PlotciCAllg<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  Abberation=ID=Value=method=val1=val2=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciCAll(n,alp,c)
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
        ggplot2::ggtitle("Confidence interval for continuity corrected methods sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for Continuity corrected methods sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for Continuity corrected methods sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for Continuity corrected methods sorted by x") +
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
#' Plots the CI estimation of  continuity corrected Wald method given n, alp and c
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval for continuity corrected Wald method
#'  for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @family Continuity correction methods of CI estimation
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotciCWD(n,alp,c)
#' @export
#7. Plot all CC methods
PlotciCWD<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (c<0 || length(c)>1) stop("'c' has to be positive")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  WaldCI.df    = ciCWD(n,alp,c)

  ss1 = data.frame( x=WaldCI.df$x, LowerLimit = WaldCI.df$LCW, UpperLimit = WaldCI.df$UCW, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Wald method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Wald method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plots the CI estimation of  continuity corrected ArcSine method given n, alp and c
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval for continuity corrected ArcSine method
#'  for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @family Continuity correction methods of CI estimation
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotciCAS(n,alp,c)
#' @export
PlotciCAS<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (c<0 || length(c)>1) stop("'c' has to be positive")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ArcSineCI.df = ciCAS(n,alp,c)
  ss1 = data.frame(x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LCA, UpperLimit = ArcSineCI.df$UCA, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected ArcSine method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected ArcSine method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plots the CI estimation of  continuity corrected Score method given n, alp and c
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval for continuity corrected Score method
#'  for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @family Continuity correction methods of CI estimation
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotciCSC(n,alp,c)
#' @export
#7. Plot all CC methods
PlotciCSC<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (c<0 || length(c)>1) stop("'c' has to be positive")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ScoreCI.df   = ciCSC(n,alp,c)
  ss1 = data.frame(x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LCS, UpperLimit = ScoreCI.df$UCS, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Score method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Score method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plots the CI estimation of  continuity corrected Logit Wald method given n, alp and c
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval for continuity corrected Logit Wald method
#'  for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @family Continuity correction methods of CI estimation
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotciCLT(n,alp,c)
#' @export
#7. Plot all CC methods
PlotciCLT<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  WaldLCI.df   = ciCLT(n,alp,c)
  ss1= data.frame(x=WaldLCI.df$x, LowerLimit = WaldLCI.df$LCLT, UpperLimit = WaldLCI.df$UCLT, LowerAbb = WaldLCI.df$LABB, UpperAbb = WaldLCI.df$UABB, ZWI = WaldLCI.df$ZWI)
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Logit Wald method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Logit Wald method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

#####################################################################################
#' Plots the CI estimation of  continuity corrected Wald-T method given n, alp and c
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval for continuity corrected Wald-T method
#'  for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @family Continuity correction methods of CI estimation
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotciCTW(n,alp,c)
#' @export
#7. Plot all CC methods
PlotciCTW<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  Abberation=ID=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  WaldTCI.df   = ciCTW(n,alp,c)

  ss1 = data.frame( x=WaldTCI.df$x, LowerLimit = WaldTCI.df$LCTW, UpperLimit = WaldTCI.df$UCTW, LowerAbb = WaldTCI.df$LABB, UpperAbb = WaldTCI.df$UABB, ZWI = WaldTCI.df$ZWI)
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Wald-T method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
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
      ggplot2::labs(title = "Confidence Interval - Continuity corrected Wald-T method") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "ID") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit),
                              size = 0.5)
  }
  oo
}

