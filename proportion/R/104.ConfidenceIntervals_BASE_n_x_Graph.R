#####################################################################################
#' Plots the CI estimation of the exact method
#' @param n - Number of trials
#' @param x - Number of sucess
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' @details  Plot of the Confidence interval for exact method
#' @family  Base methods of CI estimation given x & n
#' @examples
#' x=5; n=5; alp=0.05;e=0.5
#' PlotciEXx(x,n,alp,e) #Mid-p
#' x=5; n=5; alp=0.05;e=1 #Clopper Pearson
#' PlotciEXx(x,n,alp,e)
#' x=5; n=5; alp=0.05;e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' PlotciEXx(x,n,alp,e)
#' @export
#10. Plot all methods
PlotciEXx<-function(x,n,alp,e)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (e>1 || e<0) stop("'e' has to be between 0 and 1")
  if (length(e)>10 ) stop("Plot of only 10 interavals of 'e' is possible")
  Abberation=ID=method=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciEXx(x,n,alp,e)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,x=x,LowerLimit=ss1$LEXx,UpperLimit=ss1$UEXx,
                 LowerAbb=ss1$LABB,UpperAbb=ss1$UABB,ZWI=ss1$ZWI,e=ss1$e)
  ss$e = as.factor(ss$e)

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
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(title = "Exact method given x") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit,
                                           color= e),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red",
                                          "black", "orange","brown","chartreuse4",
                                          "blueviolet" , "deeppink", "darksalmon", "tan1" )) +
      ggplot2::scale_colour_manual(values=c("red", "black", "blue", "cyan4", "orange",
                                            "deeppink","chartreuse4",
                                            "blueviolet" , "brown", "darksalmon", "tan1")) +
      ggplot2::scale_shape_manual(values=c(21,22,23))                  # Change shapes
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::labs(y = "x values") +
      ggplot2::labs(title = "Exact method given x") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit, color= e),
                              size = 0.5)  +
      ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red",
                                          "black", "orange","brown","chartreuse4",
                                          "blueviolet" , "deeppink", "darksalmon", "tan1" )) +
      ggplot2::scale_colour_manual(values=c("red", "black", "blue", "cyan4", "orange",
                                            "deeppink","chartreuse4",
                                            "blueviolet" , "brown", "darksalmon", "tan1"))
  }
  oo
}

#####################################################################################
#' Plots the CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param x - Number of sucess
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Plots of the Confidence Intervals of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} and \code{x}
#' @family Base methods of CI estimation given x & n
#' @examples
#' x=5; n=5; alp=0.05;
#' PlotciAllx(x,n,alp)
#' @export
#10. Plot all methods
PlotciAllx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=method=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciAllx(x,n,alp)
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
      ggplot2::ggtitle("Confidence interval for base methods") +
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
      ggplot2::ggtitle("Confidence interval for base methods") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit, color= method),
                              size = 0.5)
  }
  oo
}
#############################################
#' Plots the CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) grouped by x value
#' @param x - Number of sucess
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Plots of the Confidence Interval of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} and \code{x} grouped by x
#' @family Base methods of CI estimation given x & n
#' @examples
#' x=5; n=5; alp=0.05;
#' PlotciAllxg(x,n,alp)
#' @export
#12.All methods plots with grouping
PlotciAllxg<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  Abberation=ID=method=Value=val1=val2=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciAllx(x,n,alp)
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
        ggplot2::ggtitle("Confidence interval for Base methods of CI estimation sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for Base methods of CI estimation sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for Base methods of CI estimation sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for Base methods of CI estimation sorted by x") +
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
