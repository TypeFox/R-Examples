orddom_f <- function(x,y, ... ,outputfile="orddom_csv.txt",quotechar=TRUE,decimalpt=".",separator="\t",notavailable="NA",endofline="\n") {
 filename<-outputfile
 outputfile<-""
 da<-orddom(x,y, ...)
 argum<-list(...)
  if(is.null(argum$paired)){argum[["paired"]]<-FALSE}
  if(is.null(argum$studdist)){argum[["studdist"]]<-TRUE}
  if(is.null(argum$symmetric)){argum[["symmetric"]]<-FALSE}
 if(is.matrix(da)){
 #Write header if file exists
 if (file.exists(filename)==FALSE){
 sink(file=filename, append=FALSE, type="output", split=FALSE) 
 #write header here
 header<-c("title","var1_x","var2_y","type","distribz","symm_ci")
 header<-c(header,"n_x","n_y","1-alpha","h1_tails")
 header<-c(header,"n_ygtx_w","n_yltx_w","n_yeqx_w","ps_ygtx_w","ps_xgty_w","a_ygtx_w","a_xgty_w")
 header<-c(header,"delta_dw","sd_dw","z_dw","p_dw","ci_lo_dw","ci_hi_dw","cohd_dw","cd_dw_ci_lo","cd_dw_ci_hi","NNT_dw")
 header<-c(header,"n_ygtx_b","n_yltx_b","n_yeqx_b","ps_ygtx_b","ps_xgty_b","a_ygtx_b","a_xgty_b")
 header<-c(header,"delta_db","sd_db","z_db","p_db","ci_lo_db","ci_hi_db","coh_d_db","cd_db_ci_lo","cd_db_ci_hi","NNT_db")
 header<-c(header,"n_ygtx_c","n_yltx_c","n_yeqx_c","ps_ygtx_c","ps_xgty_c","a_ygtx_c","a_xgty_c")
 header<-c(header,"delta_dc","sd_dc","z_dc","p_dc","ci_lo_dc","ci_hi_dc","coh_d_dc","cd_dc_ci_lo","cd_dc_ci_hi")
 header<-c(header,"M_x","M_y","var_x","var_y")
 header<-c(header,"delta_M","sd_t","se_t","t","p_t","ci_lo_t","ci_hi_t","coh_d_t","cd_t_ci_lo","cd_t_ci_hi","NNT_coh_d_t")
 header<-c(header,"var_d_i","var_dj_","var_dij","cov_di_dj","cov_dih_dhi","cov_db_dw","df_d","df_t")
 write.table(t(header),row.names=FALSE,col.names=FALSE,quote=quotechar,dec=decimalpt,sep=separator,eol=endofline,na=notavailable)
 }
 else {sink(file=filename, append=TRUE, type="output", split=FALSE)}
 #compile and write relevant information
if(argum$paired==FALSE){#for independent data ----------------------------------------------------------
 if(da["type_title",1]==da["type_title",2]){#no description->"NA"
   datatable1<-c("NA") } else {datatable1<-c(da["type_title",2])}
 datatable1<-c(datatable1,da["var1_X",1],da["var2_Y",1])
 datatable2<-c(da["n in X",1],da["n in Y",1],da["1-alpha",1],da["H1 tails p/CI",1])
 datatable2<-c(datatable2,NA,NA,NA,NA,NA,NA,NA)
 datatable2<-c(datatable2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
 datatable2<-c(datatable2,da["N #Y>X",1],da["N #Y<X",1],da["N #Y=X",1],da["PS Y>X",1],da["PS X>Y",1],da["A Y>X",1],da["A X>Y",1])
 datatable2<-c(datatable2,da["delta",1],da["s delta",1],da["z/t score",1],da["p",1],da["CI low",1],da["CI high",1],da["Cohen's d",1],da["d CI low",1],da["d CI high",1],da["NNT",1])
 datatable2<-c(datatable2,NA,NA,NA,NA,NA,NA,NA)
 datatable2<-c(datatable2,NA,NA,NA,NA,NA,NA,NA,NA,NA)
 datatable2<-c(datatable2,mean(return1colmatrix(x)),mean(return1colmatrix(y)),var(return1colmatrix(x)),var(return1colmatrix(y)))
 datatable2<-c(datatable2,da["delta",2],da["s delta",2],da["se delta",2],da["z/t score",2],da["p",2],da["CI low",2],da["CI high",2],da["Cohen's d",2],da["d CI low",2],da["d CI high",2],da["NNT",2])
 datatable2<-c(datatable2,da["var d.i",1],da["var dj.",1],da["var dij",1],NA,NA,NA,da["df",1],da["df",2])
 }else{#for paired data ----------------------------------------------------------
 if(da["type_title",1]==da["type_title",4]){#no description->"NA"
   datatable1<-c("NA") } else {datatable1<-c(da["type_title",4])}
 datatable1<-c(datatable1,da["var1_X_pre",1],da["var2_Y_post",1])
 datatable2<-c(length(return1colmatrix(x)),length(return1colmatrix(y)),da["1-alpha",1],da["H1 tails p/CI",1])
 for (i in 1:3){
  datatable2<-c(datatable2,da["N #Y>X",i],da["N #Y<X",i],da["N #Y=X",i],da["PS Y>X",i],da["PS X>Y",i],da["A Y>X",i],da["A X>Y",i])
  datatable2<-c(datatable2,da["delta",i],da["s delta",i],da["z/t score",i],da["p",i],da["CI low",i],da["CI high",i],da["Cohen's d",i],da["d CI low",i],da["d CI high",i])
  if(i<=2){datatable2<-c(datatable2,da["NNT",i])}}
 datatable2<-c(datatable2,mean(return1colmatrix(x)),mean(return1colmatrix(y)),var(return1colmatrix(x)),var(return1colmatrix(y)))
 datatable2<-c(datatable2,da["delta",4],da["s delta",4],NA,da["z/t score",4],da["p",4],da["CI low",4],da["CI high",4],da["Cohen's d",4],da["d CI low",4],da["d CI high",4],da["NNT",4]) 
 datatable2<-c(datatable2,da["var d.i",3],da["var dj.",3],da["var dij",3],da["cov(di,dj)",3],da["cov(dih,dhi)",3],da["cov(db,dw)",3],da["df",1],da["df",2])}
 datatable1<-c(datatable1,da["type_title",1])
 if(argum$studdist==TRUE) {datatable1<-c(datatable1,"student")} else {datatable3<-c(datatable1,"normal")}
 if(argum$symmetric==FALSE) {datatable1<-c(datatable1,"FALSE")} else {datatable1<-c(datatable3,"TRUE")}
 write.table(t(datatable1),row.names=FALSE,col.names=FALSE,quote=quotechar,dec=decimalpt,sep=separator,eol=separator,na=notavailable)
 write.table(t(as.numeric(datatable2)),row.names=FALSE,col.names=FALSE,quote=FALSE,dec=decimalpt,sep=separator,eol=endofline,na=notavailable)
 sink()
 closeAllConnections()
 return(cat("Relevant information written to file '",filename,"'.\n",sep=""))
  } else {return(cat("Error. Confirm that the basic orddom function yields any results using the given data sets."))}
  }
  
 
 
 
 
  
 