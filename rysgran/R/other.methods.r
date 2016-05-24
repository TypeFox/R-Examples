other.methods<-function(clz, vbz, method, output = "phi"){
  if(output == "metric"){
    r <- length(clz)
    soma <- sum(vbz)
    vfz <- vbz*100/soma 
    pz <- vfz
    pz[1] <- vfz[1];
    for (i in 2:r) pz[i] <- pz[i-1] + vfz[i]
    
    p<- quantile(pz, probs = seq(0, 1, 0.0005), na.rm = FALSE,names = TRUE, type = 8)
    c<- quantile(clz, probs = seq(0, 1, 0.0005), na.rm = FALSE,names = TRUE, type = 8)
    c<-((1/2^c) * 1000)

    PT03<-as.numeric((c[p>=3])[1])
    PT05<-as.numeric((c[p>=5])[1])
    PT10<-as.numeric((c[p>=10])[1])
    PT15<-as.numeric((c[p>=15])[1])
    PT16<-as.numeric((c[p>=16])[1])
    PT20<-as.numeric((c[p>=20])[1])
    PT25<-as.numeric((c[p>=25])[1])
    PT30<-as.numeric((c[p>=30])[1])
    PT35<-as.numeric((c[p>=35])[1])
    PT45<-as.numeric((c[p>=45])[1])
    PT50<-as.numeric((c[p>=50])[1])
    PT55<-as.numeric((c[p>=55])[1])
    PT65<-as.numeric((c[p>=65])[1])
    PT70<-as.numeric((c[p>=70])[1])
    PT75<-as.numeric((c[p>=75])[1])
    PT80<-as.numeric((c[p>=80])[1])
    PT84<-as.numeric((c[p>=84])[1])
    PT85<-as.numeric((c[p>=85])[1])
    PT90<-as.numeric((c[p>=90])[1])
    PT95<-as.numeric((c[p>=95])[1])
    PT97<-as.numeric((c[p>=97])[1])
    
    if (method=="folk")
      {
      mean<-exp( (log(PT16)+log(PT50)+log(PT84))/3 )
      sort<-exp( ((log(PT16)-log(PT84))/4)+((log(PT05)-log(PT95))/6.6) )
      }
    if (method=="mcA")
      {
      mean<-exp((log(PT10)+log(PT30)+log(PT50)+log(PT70)+log(PT90))/5)
      sort<-exp((log(PT15)+log(PT05)-log(PT95)-log(PT85))/5.4)
      }
    if (method=="mcB")
      {
      mean<-exp((log(PT05)+log(PT15)+log(PT25)+log(PT35)+log(PT45)+log(PT55)+log(PT65)+log(PT75)+log(PT85)+log(PT95))/10)
      sort<-exp((log(PT30)+log(PT20)+log(PT10)+log(PT03)-log(PT70)-log(PT80)-log(PT90)-log(PT97))/9.1)
      }
    if (method=="trask")
      {
      mean<-exp(log(PT50))
      sort<-exp((log(PT25)-log(PT75))/1.35)
      }
    if (method=="otto")
      {
      mean<-exp((log(PT16)+log(PT84))/2)
      sort<-exp((log(PT16)-log(PT84))/2)
      }
    median<-exp(log(PT50))
    skew<-(((log(PT16)+log(PT84)-2*log(PT50))/(2*(log(PT16)-log(PT84))))+((log(PT05)+log(PT95)-2*log(PT50))/(2*(log(PT05)-log(PT95)))))
    kurt<-(log(PT05)-log(PT95))/(2.44*(log(PT25)-log(PT75)))
    }
  if (output=="phi"){
    r <- length(clz)
    soma <- sum(vbz)
    vfz <- vbz*100/soma 
    pz <- vfz
    pz[1] <- vfz[1];
    for (i in 2:r) pz[i] <- pz[i-1] + vfz[i]
    
    p<- quantile(pz, probs = seq(0, 1, 0.0005), na.rm = FALSE,names = FALSE, type = 8)
    c<- quantile(clz, probs = seq(0, 1, 0.0005), na.rm = FALSE,names = FALSE, type = 8)
    
    PT03<-as.numeric((c[p>=3])[1])
    PT05<-as.numeric((c[p>=5])[1])
    PT10<-as.numeric((c[p>=10])[1])
    PT15<-as.numeric((c[p>=15])[1])
    PT16<-as.numeric((c[p>=16])[1])
    PT20<-as.numeric((c[p>=20])[1])
    PT25<-as.numeric((c[p>=25])[1])
    PT30<-as.numeric((c[p>=30])[1])
    PT35<-as.numeric((c[p>=35])[1])
    PT45<-as.numeric((c[p>=45])[1])
    PT50<-as.numeric((c[p>=50])[1])
    PT55<-as.numeric((c[p>=55])[1])
    PT65<-as.numeric((c[p>=65])[1])
    PT70<-as.numeric((c[p>=70])[1])
    PT75<-as.numeric((c[p>=75])[1])
    PT80<-as.numeric((c[p>=80])[1])
    PT84<-as.numeric((c[p>=84])[1])
    PT85<-as.numeric((c[p>=85])[1])
    PT90<-as.numeric((c[p>=90])[1])
    PT95<-as.numeric((c[p>=95])[1])
    PT97<-as.numeric((c[p>=97])[1])
    
    if (method=="folk")
      {
      mean<-(PT16+PT50+PT84)/3
      sort<-((PT84-PT16)/4)+((PT95-PT05)/6.6)
      }
    if (method=="mcA")
      {
      mean<-((PT10+PT30+PT50+PT70+PT90)/5)
      sort<-((PT85+PT95-PT05-PT15)/5.4)
      }
    if (method=="mcB")
      {
      mean<-(PT05+PT15+PT25+PT35+PT45+PT55+PT65+PT75+PT85+PT95)/10
      sort<-(PT70+PT80+PT90+PT97-PT03-PT10-PT20-PT30)/9.1
      }
    if (method=="trask")
      {
      mean<-PT50
      sort<-(PT75-PT25)/1.35
      }
    if (method=="otto")
      {
      mean<-(PT16+PT84)/2
      sort<-(PT84-PT16)/2
      }
    median<-PT50
    skew<-((PT16+PT84-2*PT50)/(2*(PT84-PT16)))+((PT05+PT95-2*PT50)/(2*(PT95-PT05)))
    kurt<-(PT95-PT05)/(2.44*(PT75-PT25))
    }
  table<-list(mean=mean, median=median, sort=sort, skew=skew, kurt=kurt)
  return(table)
  }
