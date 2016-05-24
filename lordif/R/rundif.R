rundif <-
function(item,resp,theta,gr,criterion,alpha,beta.change,pseudo.R2,R2.change,wt) {
    ncat<-apply(resp,2,max,na.rm=T) 
    ni<-length(item)
    beta12<-numeric(ni)
    chi12<-numeric(ni)
    chi13<-numeric(ni)
    chi23<-numeric(ni)
    df12<-numeric(ni)
    df13<-numeric(ni)
    df23<-numeric(ni)
    pseudo12.CoxSnell<-numeric(ni)
    pseudo13.CoxSnell<-numeric(ni)
    pseudo23.CoxSnell<-numeric(ni)
    pseudo12.Nagelkerke<-numeric(ni)
    pseudo13.Nagelkerke<-numeric(ni)
    pseudo23.Nagelkerke<-numeric(ni)
    pseudo12.McFadden<-numeric(ni)
    pseudo13.McFadden<-numeric(ni)
    pseudo23.McFadden<-numeric(ni)
    flag.post<-logical(ni)
    for (i in 1:ni) {
      output<-try(runolr(resp[,i],theta,as.factor(gr),wt),silent=T)
      if (class(output)=="try-error") output<-runolr(resp[,i],log((theta-min(theta)+0.01)/(max(theta)-theta+.01)),as.factor(gr),wt)
      if (exists("output")) {
        beta12[i]<-output$beta12
        chi12[i]<-output$chi12
        chi13[i]<-output$chi13
        chi23[i]<-output$chi23
        df12[i]<-output$df12
        df13[i]<-output$df13
        df23[i]<-output$df23
        pseudo12.CoxSnell[i]<-output$pseudo12.CoxSnell
        pseudo13.CoxSnell[i]<-output$pseudo13.CoxSnell
        pseudo23.CoxSnell[i]<-output$pseudo23.CoxSnell
        pseudo12.Nagelkerke[i]<-output$pseudo12.Nagelkerke
        pseudo13.Nagelkerke[i]<-output$pseudo13.Nagelkerke
        pseudo23.Nagelkerke[i]<-output$pseudo23.Nagelkerke
        pseudo12.McFadden[i]<-output$pseudo12.McFadden
        pseudo13.McFadden[i]<-output$pseudo13.McFadden
        pseudo23.McFadden[i]<-output$pseudo23.McFadden
      }
    }
    if (toupper(criterion)=="CHISQR") flag.post<-(chi12<=alpha | chi13<=alpha | chi23<=alpha)
    else if (toupper(criterion)=="BETA") flag.post<-beta12>=beta.change
    else if (toupper(criterion)=="R2") {
      if (toupper(pseudo.R2)=="MCFADDEN") flag.post<-pseudo13.McFadden>=R2.change
      else if (toupper(pseudo.R2)=="NAGELKERKE") flag.post<-pseudo13.Nagelkerke>=R2.change
      else if (toupper(pseudo.R2)=="COXSNELL") flag.post<-pseudo13.CoxSnell>=R2.change
      else {
        warning("invalid pseudo R^2 is selected: McFadden will be used instead")
        flag.post<-(chi12<=alpha | chi13<=alpha | chi23<=alpha)
      }
    } else {
      warning("invalid flagging criterion is selected: \"Chisqr\" will be used instead")
      flag.post<-chi13<=alpha
    }
    stats<-data.frame(item,ncat,chi12,chi13,chi23,beta12,
                      pseudo12.McFadden,pseudo13.McFadden,pseudo23.McFadden,
                      pseudo12.Nagelkerke,pseudo13.Nagelkerke,pseudo23.Nagelkerke,
                      pseudo12.CoxSnell,pseudo13.CoxSnell,pseudo23.CoxSnell,df12,df13,df23)
    row.names(stats)<-NULL
    return(list(stats=stats,flag=flag.post))
  }
