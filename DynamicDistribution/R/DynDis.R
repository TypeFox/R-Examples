#' @export DynDis
DynDis <-
function(name,par_matrix,total=c(100,100,100),choice='cdf',interval=0.05,const_par=c(NULL,NULL,NULL))
{                                                                                                   
  
  par_no<-ncol(par_matrix) 
  if(par_no>=4||par_no==0)
    stop("Error! The number of parameters must be 1,2 or 3!")
  case<-c('cdf','Mean','Variance','Mode','Skewness','Kurtosis')
  a<-which(case==choice)
  if(par_no==1)
  {
    if(identical(as.character(substitute(name)),"Dis_Uniform")!=T){
      total[1]=total[1]
      para1<-seq(par_matrix[1,1],par_matrix[2,1],length.out=total[1])
    } else {
      total[1]=1+round(par_matrix[2,1])-round(par_matrix[1,1])
      para1<-seq(round(par_matrix[1,1]),round(par_matrix[2,1]),length.out=total[1]) 
      
    }  
    
    
    para1<-t(as.matrix(para1))
    k_y_range<-name(para=para1)
    k_range<-k_y_range[1]:k_y_range[2]
    y_max<-k_y_range[3]
    
    paramatrix<-matrix(nrow=total,ncol=3)
    for(i in 1:total[1])
    {
      num_char<-name(k_range,as.matrix(para1[i]))
      dev.hold()
      par(mfrow=c(1,2),oma=c(1,1,1,1))
      
      plot(k_range,num_char$density,type='h',xlim=c(min(k_range),max(k_range)),ylim=c(0,y_max),
           xlab='k',ylab='p(k)',main='Probability Mass Function',font.lab=9,font.main=7,
           cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
      points(k_range,num_char$density)
      legend("topright",bty='n',paste(paste(num_char[9],'='),round(as.numeric(num_char[8]),digits=3)))
      
      if(a==1)
      {
        plot(0,type="n",xlim=c(min(k_range),max(k_range)+1),ylim=c(0,1),xlab='k',ylab='F(k)',
             main='Cumulative Distribution Function',font.lab=9,font.main=7,cex.axis=0.8,
             cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        segments(k_range,num_char[[a+1]],k_range+1,num_char[[a+1]])  
        legend("topleft",bty='n',paste(paste(num_char[9],'='),round(as.numeric(num_char[8]),digits=3)))
      }
      
      if(a!=1)
      {
        if(num_char[[a+1]]==Inf)
        {
          paramatrix[i,1]<-para1[i];
          paramatrix[i,2]<-0
          paramatrix[i,3]<-4
        }
        else
        {
          paramatrix[i,1]<-para1[i];
          paramatrix[i,2]<-num_char[[a+1]]
          paramatrix[i,3]<-21
        }
        plot(paramatrix[1:i,1],paramatrix[1:i,2],type='p',pch=paramatrix[1:i,3],bg="grey",xlim=c(min(para1),max(para1)),
             xlab=num_char[9],ylab=case[a],main=case[a],font.lab=7,font.main=7,cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,
             panel.first=grid())
        
        if(num_char[[a+1]]==Inf)
          leg=0
        else leg=paramatrix[i,2]
        
        legend("topleft",bty='n',paste(paste(case[a],'='),round(leg,digits=3)))
      }
      
      title(main=num_char[10],outer=T,font.main=9,cex.main=1.6)
      
      if(dev.interactive())
      {
        dev.flush()
        Sys.sleep(interval)
      }
      
    } 
    
  }
  
  if(par_no==2)        
  {
    
    total[1]=1+round(par_matrix[2,1])-round(par_matrix[1,1])
    total[2]=total[2]
    para1<-seq(round(par_matrix[1,1]),round(par_matrix[2,1]),length.out=total[1])
    para2<-seq(par_matrix[1,2],par_matrix[2,2],length.out=total[2])
    
    para1_no=length(para1)
    para2_no=length(para2)
    if(para1_no>para2_no){
      para2<-c(para2,seq(par_matrix[2,2],par_matrix[2,2],length.out=para1_no-para2_no))
    }else{
      para1<-c(para1,seq(par_matrix[2,1],par_matrix[2,1],length.out=para2_no-para1_no))}
    
    Para<-t(cbind(para1,para2))
    k_y_range<-name(para=Para,const_par=const_par)
    k_range1<-k_y_range[1]:k_y_range[2]   
    k_range2<-k_y_range[4]:k_y_range[5]    
    y_max1<-k_y_range[3]
    y_max2<-k_y_range[6]
    
    paramatrix1<-matrix(nrow=max(total[1],total[2]),ncol=3)
    paramatrix2<-matrix(nrow=max(total[1],total[2]),ncol=3)
    
    
    for(j in 1:max(total[1],total[2]))
    {
      dev.hold()
      num_char<-name(k_range1,matrix(c(para1[j],const_par[2])),const_par=const_par)
      par(mfrow=c(2,2),oma=c(1,1,1,1)) 
      
      plot(k_range1,num_char$density,type='h',xlim=c(min(k_range1),max(k_range1)),ylim=c(0,y_max1),
           xlab='k',ylab='p(k)',main='Probability Mass Function',font.lab=9,font.main=7,
           cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
      points(k_range1,num_char$density)
      legend("topright",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3))))
      
      if(a==1)
      {
        plot(0,type='n',xlim=c(min(k_range1),max(k_range1)),ylim=c(0,1),
             xlab='k',ylab='F(k)',main='Cumulative Distribution Function',font.lab=9,font.main=7,
             cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        segments(k_range1,num_char[[a+1]],k_range1+1,num_char[[a+1]])  
        legend("topleft",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3))))
      }
      if(a!=1)  
      {
        if(num_char[[a+1]]==Inf)
        {
          paramatrix1[j,1]<-para1[j];
          paramatrix1[j,2]<-0
          paramatrix1[j,3]<-4
        }
        else
        {
          paramatrix1[j,1]<-para1[j];
          paramatrix1[j,2]<-num_char[[a+1]]
          paramatrix1[j,3]<-21
        }
        plot(paramatrix1[1:j,1],paramatrix1[1:j,2],pch=paramatrix1[1:j,3],bg="grey",type='p',
             xlim=c(min(para1),max(para1)),xlab=num_char[9],ylab=case[a],main=case[a],font.lab=7,
             font.main=7,cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        
        if(num_char[[a+1]]==Inf)
          leg=0
        else leg=paramatrix1[j,2]
        
        legend("topleft",bty='n',paste(paste(case[a],'='),round(leg,digits=3)))
      }
      
      
      num_char<-name(k_range2,as.matrix(c(const_par[1],para2[j])),const_par=const_par)
      dev.hold()
      plot(k_range2,num_char$density,type='h',xlim=c(min(k_range2),max(k_range2)),ylim=c(0,y_max2),
           xlab='k',ylab='p(k)',main='Probability Mass Function',font.lab=9,font.main=7,
           cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
      points(k_range2,num_char$density)
      legend("topright",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3))))
      if(a==1)
      {
        plot(0,type='n',xlim=c(min(k_range2),max(k_range2)),ylim=c(0,1),
             xlab='k',ylab='F(k)',main='Cumulative Distribution Function',font.lab=9,font.main=7,
             cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid()) 
        segments(k_range2,num_char[[a+1]],k_range2+1,num_char[[a+1]])   
        legend("topleft",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3)))) 
      }
      if(a!=1)
      {
        if(num_char[[a+1]]==Inf) 
        {
          paramatrix2[j,1]<-para2[j];
          paramatrix2[j,2]<-0
          paramatrix2[j,3]<-4
        }
        else
        {
          paramatrix2[j,1]<-para2[j];
          paramatrix2[j,2]<-num_char[[a+1]]
          paramatrix2[j,3]<-21 
        }
        plot(paramatrix2[1:j,1],paramatrix2[1:j,2],pch=paramatrix2[1:j,3],bg="grey",type='p',
             xlim=c(min(para2),max(para2)),xlab=num_char[11],ylab=case[a],main=case[a],font.lab=7,
             font.main=7,cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        
        if(num_char[[a+1]]==Inf)
          leg=0
        else leg=paramatrix2[j,2]
        
        legend("topleft",bty='n',paste(paste(case[a],'='),round(leg,digits=3)))
      }
      
      title(main=num_char[12],outer=T,font.main=9,cex.main=1.6)
      
      if(dev.interactive())
      {
        dev.flush()
        Sys.sleep(interval)
      }
      
    }
  }
  if(par_no==3)        
  {
    total[1]=1+round(par_matrix[2,1])-round(par_matrix[1,1])
    total[2]=1+round(par_matrix[2,2])-round(par_matrix[1,2])
    total[3]=1+round(par_matrix[2,3])-round(par_matrix[1,3])
    para1<-seq(round(par_matrix[1,1]),round(par_matrix[2,1]),length.out=total[1])
    para2<-seq(round(par_matrix[1,2]),round(par_matrix[2,2]),length.out=total[2])
    para3<-seq(round(par_matrix[1,3]),round(par_matrix[2,3]),length.out=total[3])
    
    para1_no=length(para1)
    para2_no=length(para2)
    para3_no=length(para3)
    
    if(para1_no>para2_no){
      if(para2_no>para3_no){
        para2<-c(para2,seq(par_matrix[2,2],par_matrix[2,2],length.out=para1_no-para2_no))
        para3<-c(para3,seq(par_matrix[2,3],par_matrix[2,3],length.out=para1_no-para3_no))
      }
      else{
        if(para1_no>para3_no){
          para2<-c(para2,seq(par_matrix[2,2],par_matrix[2,2],length.out=para1_no-para2_no))
          para3<-c(para3,seq(par_matrix[2,3],par_matrix[2,3],length.out=para1_no-para3_no))
        }else{
          para2<-c(para2,seq(par_matrix[2,2],par_matrix[2,2],length.out=para3_no-para2_no))
          para1<-c(para1,seq(par_matrix[2,1],par_matrix[2,1],length.out=para3_no-para1_no))
        }
      }
    }else{
      if(para1_no>para3_no){
        para1<-c(para1,seq(par_matrix[2,1],par_matrix[2,1],length.out=para2_no-para1_no))
        para3<-c(para3,seq(par_matrix[2,3],par_matrix[2,3],length.out=para2_no-para3_no))
      }else{
        if(para2_no>para3_no){
          para1<-c(para1,seq(par_matrix[2,1],par_matrix[2,1],length.out=para2_no-para1_no))
          para3<-c(para3,seq(par_matrix[2,3],par_matrix[2,3],length.out=para2_no-para3_no))
        }else{
          para1<-c(para1,seq(par_matrix[2,1],par_matrix[2,1],length.out=para3_no-para1_no))
          para2<-c(para2,seq(par_matrix[2,2],par_matrix[2,2],length.out=para3_no-para2_no))
        }
      }
    }
    
    
    Para<-rbind(para1,para2,para3)
    k_y_range<-name(para=Para,const_par=const_par)
    k_range1<-k_y_range[1]:k_y_range[2]   
    k_range2<-k_y_range[4]:k_y_range[5]  
    k_range3<-k_y_range[7]:k_y_range[8]  
    y_max1<-k_y_range[3]
    y_max2<-k_y_range[6]
    y_max3<-k_y_range[9]
    
    paramatrix1<-matrix(nrow=max(total[1],total[2],total[3]),ncol=3)
    paramatrix2<-matrix(nrow=max(total[1],total[2],total[3]),ncol=3)
    paramatrix3<-matrix(nrow=max(total[1],total[2],total[3]),ncol=3)    
    
    for(j in 1:max(total[1],total[2],total[3]))
    {
      dev.hold()
      num_char<-name(k_range1,matrix(c(para1[j],const_par[2],const_par[3])),const_par=const_par)
      par(mfrow=c(3,2),oma=c(1,1,1,1)) 
      
      plot(k_range1,num_char$density,type='h',xlim=c(min(k_range1),max(k_range1)),ylim=c(0,y_max1),
           xlab='k',ylab='p(k)',main='Probability Mass Function',font.lab=9,font.main=7,
           cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
      points(k_range1,num_char$density)
      legend("topright",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'='),paste(num_char[13],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3),round(as.numeric(num_char[12]),digits=3))))
      
      if(a==1)
      {
        plot(0,type='n',xlim=c(min(k_range1),max(k_range1)),ylim=c(0,1),
             xlab='k',ylab='F(k)',main='Cumulative Distribution Function',font.lab=9,font.main=7,
             cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        segments(k_range1,num_char[[a+1]],k_range1+1,num_char[[a+1]])  
        legend("topleft",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'='),paste(num_char[13],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3),round(as.numeric(num_char[12]),digits=3))))
      }
      if(a!=1)  
      {
        if(num_char[[a+1]]==Inf)
        {
          paramatrix1[j,1]<-para1[j];
          paramatrix1[j,2]<-0
          paramatrix1[j,3]<-4
        }
        else
        {
          paramatrix1[j,1]<-para1[j];
          paramatrix1[j,2]<-num_char[[a+1]]
          paramatrix1[j,3]<-21
        }
        plot(paramatrix1[1:j,1],paramatrix1[1:j,2],pch=paramatrix1[1:j,3],bg="grey",type='p',
             xlim=c(min(para1),max(para1)),xlab=num_char[9],ylab=case[a],main=case[a],font.lab=7,
             font.main=7,cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        
        if(num_char[[a+1]]==Inf)
          leg=0
        else leg=paramatrix1[j,2]
        
        legend("topleft",bty='n',paste(paste(case[a],'='),round(leg,digits=3)))
      }
      
      
      num_char<-name(k_range2,as.matrix(c(const_par[1],para2[j],const_par[3])),const_par=const_par)
      dev.hold()
      plot(k_range2,num_char$density,type='h',xlim=c(min(k_range2),max(k_range2)),ylim=c(0,y_max2),
           xlab='k',ylab='p(k)',main='Probability Mass Function',font.lab=9,font.main=7,
           cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
      points(k_range2,num_char$density)
      legend("topright",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'='),paste(num_char[13],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3))))
      if(a==1)
      {
        plot(0,type='n',xlim=c(min(k_range2),max(k_range2)),ylim=c(0,1),
             xlab='k',ylab='F(k)',main='Cumulative Distribution Function',font.lab=9,font.main=7,
             cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid()) 
        segments(k_range2,num_char[[a+1]],k_range2+1,num_char[[a+1]])   
        legend("topleft",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'='),paste(num_char[13],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3)))) 
      }
      if(a!=1)
      {
        if(num_char[[a+1]]==Inf) 
        {
          paramatrix2[j,1]<-para2[j];
          paramatrix2[j,2]<-0
          paramatrix2[j,3]<-4
        }
        else
        {
          paramatrix2[j,1]<-para2[j];
          paramatrix2[j,2]<-num_char[[a+1]]
          paramatrix2[j,3]<-21 
        }
        plot(paramatrix2[1:j,1],paramatrix2[1:j,2],pch=paramatrix2[1:j,3],bg="grey",type='p',
             xlim=c(min(para2),max(para2)),xlab=num_char[11],ylab=case[a],main=case[a],font.lab=7,
             font.main=7,cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        
        if(num_char[[a+1]]==Inf)
          leg=0
        else leg=paramatrix2[j,2]
        
        legend("topleft",bty='n',paste(paste(case[a],'='),round(leg,digits=3)))
      }
      
      num_char<-name(k_range3,as.matrix(c(const_par[1],const_par[2],para3[j])),const_par=const_par)
      dev.hold()
      plot(k_range3,num_char$density,type='h',xlim=c(min(k_range3),max(k_range3)),ylim=c(0,y_max3),
           xlab='k',ylab='p(k)',main='Probability Mass Function',font.lab=9,font.main=7,
           cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
      points(k_range3,num_char$density)
      legend("topright",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'='),paste(num_char[13],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3),round(as.numeric(num_char[12]),digits=3))))
      if(a==1)
      {
        plot(0,type='n',xlim=c(min(k_range3),max(k_range3)),ylim=c(0,1),
             xlab='k',ylab='F(k)',main='Cumulative Distribution Function',font.lab=9,font.main=7,
             cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid()) 
        segments(k_range3,num_char[[a+1]],k_range3+1,num_char[[a+1]])   
        legend("topleft",bty='n',paste(c(paste(num_char[9],'='),paste(num_char[11],'='),paste(num_char[13],'=')),c(round(as.numeric(num_char[8]),digits=3),round(as.numeric(num_char[10]),digits=3),round(as.numeric(num_char[12]),digits=3)))) 
      }
      if(a!=1)
      {
        if(num_char[[a+1]]==Inf) 
        {
          paramatrix3[j,1]<-para3[j];
          paramatrix3[j,2]<-0
          paramatrix3[j,3]<-4
        }
        else
        {
          paramatrix3[j,1]<-para3[j];
          paramatrix3[j,2]<-num_char[[a+1]]
          paramatrix3[j,3]<-21 
        }
        plot(paramatrix3[1:j,1],paramatrix3[1:j,2],pch=paramatrix3[1:j,3],bg="grey",type='p',
             xlim=c(min(para3),max(para3)),xlab=num_char[13],ylab=case[a],main=case[a],font.lab=7,
             font.main=7,cex.axis=0.8,cex.lab=1.2,cex.main=1.3,tcl=0.2,panel.first=grid())
        
        if(num_char[[a+1]]==Inf)
          leg=0
        else leg=paramatrix3[j,2]
        
        legend("topleft",bty='n',paste(paste(case[a],'='),round(leg,digits=3)))
      }
      title(main=num_char[[length(num_char)]],outer=T,font.main=9,cex.main=1.6)
      
      if(dev.interactive())
      {
        dev.flush()
        Sys.sleep(interval)
      }
      
    }
  }
  
}
