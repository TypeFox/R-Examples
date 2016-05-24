SEISNtime<-function(GH)
{
  rdate = recdatel(GH$info)
  
  yd = yeardate(rdate$yr,rdate$jd, rdate$hr, rdate$mi,rdate$sec      )
  
  w1 = which.min(yd)

  ymax = yeardate(rdate$yr,rdate$jd, rdate$hr, rdate$mi,rdate$sec +   GH$info$dt* GH$info$n  )
  
  w2 = which.max(ymax)


  rdate2 = recdate( rdate$jd[w2], rdate$hr[w2], rdate$mi[w2],rdate$sec[w2] +   GH$info$dt[w2]* GH$info$n[w2], rdate$yr[w2])
  
  rdate1 =  recdate( rdate$jd[w1], rdate$hr[w1], rdate$mi[w1],  rdate$sec[w1],  rdate$yr[w1]    )


  cout = rbind(data.frame(rdate1) , data.frame(rdate2)  )

  sdf = secdifL(rdate1 ,  rdate2   )


  
  ##  result =  list( yr=cout$yr , jd =cout$jd , hr =cout$hr  , mi = cout$mi , sec = cout$sec , w1=c(w1, w2) , sdif = sdf )

  result = list(start=rdate1, end=rdate2, sdif = sdf ) 
  
  return(result  )
}

