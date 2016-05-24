tdt<-function(sample){
  dec=3
  idx.case=grep(TRUE,(sample[,1]==1))
  idx.ctrl=grep(TRUE,(sample[,1]==0))

  case=sample[idx.case,]
  ctrl=sample[idx.ctrl,]

  #----------------------------------------
  # Compute diagonal elements of TDT table
  # for case-trios and control-trios, find
  # TDT p-values, and adjusted TDT p-value
  #----------------------------------------
  if (length(idx.case)==0){
    cat('There are only control-trios in this sample.',fill=T)

    b2=tdt.cnt(ctrl)[5]
    c2=tdt.cnt(ctrl)[6]
    r2=b2/(b2+c2)
    t.tdt.ctrl=(b2-c2)^2/(b2+c2)
    p.tdt.ctrl=1-pchisq(t.tdt.ctrl,1)

    out=data.frame(rbind(c(NA,t.tdt.ctrl,NA),
                         c(NA,p.tdt.ctrl,NA)))
  }
  else if (length(idx.ctrl)==0){
    cat('There are only case-trios in this sample.',fill=T)

    b1=tdt.cnt(case)[5]
    c1=tdt.cnt(case)[6]
    t.tdt.case=(b1-c1)^2/(b1+c1)
    p.tdt.case=1-pchisq(t.tdt.case,1)

    out=data.frame(rbind(c(t.tdt.case,NA,NA),
                         c(p.tdt.case,NA,NA)))
  } else {
    b1=tdt.cnt(case)[5]
    c1=tdt.cnt(case)[6]
    t.tdt.case=(b1-c1)^2/(b1+c1)
    p.tdt.case=1-pchisq(t.tdt.case,1)

    b2=tdt.cnt(ctrl)[5]
    c2=tdt.cnt(ctrl)[6]
    r2=b2/(b2+c2)
    t.tdt.ctrl=(b2-c2)^2/(b2+c2)
    p.tdt.ctrl=1-pchisq(t.tdt.ctrl,1)

    t.tdt.adj=((1-r2)*b1-r2*c1)^2/(r2*(1-r2)*(b1+c1))
    p.tdt.adj=1-pchisq(t.tdt.adj,1)

    out=data.frame(rbind(c(t.tdt.case,t.tdt.ctrl,t.tdt.adj),
                         c(p.tdt.case,p.tdt.ctrl,p.tdt.adj)))
  }
  out=signif(out,dec)
  rownames(out)=c("Statistics","p-value")
  names(out)=c("Case","Control","Adjusted-Case")
  out
}

