table.nw <-
function(tdcor_out,outfile="TDCor_output.txt",n=100,...)

{ 

  Lst<-list(...)

  if (!is.null(Lst$stimulus))

  {stimulus=Lst$stimulus}else

  {stimulus="Stimulus"}



  nw=data.frame(tdcor_out$network)

  nw_bool=ceiling(abs(nw)/max(abs(nw)))*sign(nw)

  tbl=data.frame(matrix(0,6,0))



iin=as.vector(which(apply(abs(nw),2,sum)>0))



for (i in iin)

{

m1=matrix(0,6,length(names(nw)[nw[,i]!=0]))

m1[1,]= rep(rownames(nw)[i],dim(m1)[2])

 m1[2,]= as.numeric(nw_bool[nw[,i]!=0,i])

 m1[3,]= names(nw)[nw[,i]!=0]

 m1[4,]= round(as.numeric(abs(nw)[nw[,i]!=0,i])/n*100,1)

 m1[5,]= as.numeric(tdcor_out$ID[nw[,i]!=0,i])

 m1[6,]= round(as.numeric(tdcor_out$delay[nw[,i]!=0,i]),1)

 tbl<-data.frame(tbl,m1)

}



if (length(which(tdcor_out$EP>0))>0)

{for (i in which(tdcor_out$EP>0))

{

m1=matrix(0,6,1)

m1[1,]= rep(stimulus,dim(m1)[2])

m1[2,]= rep(1,dim(m1)[2])

 m1[3,]= names(nw)[i]

m1[4,]= round(tdcor_out$EP[i]/n*100,1)

m1[5,]= rep(0,dim(m1)[2])

m1[6,]= rep(0,dim(m1)[2])

 tbl<-data.frame(tbl,m1)

}

}



 tbl=data.frame(t(tbl)); names(tbl)=c("TF1","inter.","TF2","bootstrap","ID","delay"); rownames(tbl)=seq(1,dim(tbl)[1])

 write.table(tbl,outfile,col.names=F,row.names=F,sep="\t")



 tbl_param=data.frame(matrix(0,31,2))

 rownames(tbl_param)=c("Version","tol","delayspan","thr_cor","delaymin","delaymax","delay","thr_ind1","thr_ind2","thr_overlap","thrpTPI",

"thrpDPI","thr_isr","thr_bool_EP","MinTarNumber","MinProp","MaxEPNumber","regmax","n0 n1",

"TPI_N","TPI_ks","TPI_kd","TPI_delta","TPI_delay","TPI_noise",

"DPI_N","DPI_ks","DPI_kd","DPI_delta","DPI_delay","DPI_noise")



 tbl_param["Version",]=as.character(c(tdcor_out$input$version,""))

 tbl_param["tol",]=as.character(c(tdcor_out$input$tol,""))

 tbl_param["delayspan",]=as.character(c(tdcor_out$input$delayspan,""))

 tbl_param["thr_cor",]=as.character(tdcor_out$input$thr_cor)

 tbl_param["delaymin",]=as.character(tdcor_out$input$delaymin)

 tbl_param["delaymax",]=as.character(tdcor_out$input$delaymax)

 tbl_param["delay",]=as.character(c(tdcor_out$input$delay,""))

 tbl_param["thr_ind1",]=as.character(tdcor_out$input$thr_ind1)

 tbl_param["thr_ind2",]=as.character(tdcor_out$input$thr_ind2)

 tbl_param["thr_overlap",]=as.character(tdcor_out$input$thr_overlap)

 tbl_param["thrpTPI",]=as.character(tdcor_out$input$thrpTPI)

 tbl_param["thrpDPI",]=as.character(tdcor_out$input$thrpDPI)

 tbl_param["thr_isr",]=as.character(tdcor_out$input$thr_isr)

 tbl_param["thr_bool_EP",]=rep(tdcor_out$input$thr_bool_EP,2)

 tbl_param["MinTarNumber",]=as.character(c(tdcor_out$input$MinTarNumber,""))

 tbl_param["MinProp",]=as.character(c(tdcor_out$input$MinProp,""))

 tbl_param["MaxEPNumber",]=as.character(c(tdcor_out$input$MaxEPNumber,""))

 tbl_param["regmax",]=as.character(c(tdcor_out$input$regmax,""))

 tbl_param["n0 n1",]=as.character(c(tdcor_out$input$n0,tdcor_out$input$n1))

 tbl_param["TPI_N",]=as.character(c(tdcor_out$input$TPI_input$N,""))

 tbl_param["TPI_ks",]=as.character(tdcor_out$input$TPI_input$ks_int)

 tbl_param["TPI_kd",]=as.character(tdcor_out$input$TPI_input$ks_int)

 tbl_param["TPI_delta",]=as.character(tdcor_out$input$TPI_input$delta_int)

 tbl_param["TPI_delay",]=as.character(tdcor_out$input$TPI_input$delay)

 tbl_param["TPI_noise",]=as.character(tdcor_out$input$TPI_input$noise)

 tbl_param["DPI_N",]=as.character(c(tdcor_out$input$DPI_input$N,""))

 tbl_param["DPI_ks",]=as.character(tdcor_out$input$DPI_input$ks_int)

 tbl_param["DPI_kd",]=as.character(tdcor_out$input$DPI_input$ks_int)

 tbl_param["DPI_delta",]=as.character(tdcor_out$input$DPI_input$delta_int)

 tbl_param["DPI_delay",]=as.character(tdcor_out$input$DPI_input$delay)

 tbl_param["DPI_noise",]=as.character(tdcor_out$input$DPI_input$noise)



 write.table(tbl_param,paste("Param-",outfile,sep=""),col.names=F,sep="\t")



return(tbl)}
