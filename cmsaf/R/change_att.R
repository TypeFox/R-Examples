change_att <-
function(var,infile,v_name="No_change",s_name="No_change",l_name="No_change",u_name="No_change",F_val="No_change",m_val="No_change",val_prec="double"){

# define standard names of variables and dimensions

   info = "Created with the CM SAF R toolbox." 
  
# open infile

  cat("open infile", "\n")

  id <- nc_open(infile, write=TRUE)

  # get information about variables
	
  varnames <- names(id$var)

   if (var %in% varnames){
     if (s_name!="No_change"){
      ncatt_put(id,var,"standard_name",s_name,prec="text")
      cat("standard_name changed to ",s_name,sep="", "\n")}
     if (l_name!="No_change"){
      ncatt_put(id,var,"long_name",l_name,prec="text")
      cat("long_name changed to ",l_name,sep="", "\n")}
     if (u_name!="No_change"){
      ncatt_put(id,var,"units",u_name,prec="text")
      cat("units changed to ",u_name,sep="", "\n")}
     if (F_val!="No_change"){
      ncatt_put(id,var,"_FillValue",F_val,prec=val_prec)
      cat("FillValue changed to ",F_val,sep="", "\n")}
     if (m_val!="No_change"){
      ncvar_change_missval( id, var, m_val )
      cat("missing_value changed to ",m_val,sep="", "\n")}
    if (v_name!="No_change"){
      ncvar_rename( id, old_varname=var, new_varname=v_name)
      cat("variable name changed to ",v_name,sep="", "\n")}

      ncatt_put(id,0,"Info",info,prec="text")

      nc_close(id)

   }else{
      nc_close(id)
      stop(cat(paste("Variable ",var," not found! File contains: ",varnames,sep="")),"\n")}
}