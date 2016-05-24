#retreive desired parameter values specified in param_list from wasim-parameter file wasim_param_file

#Till 2.4.08

extract_wasim_params=function(wasim_param_file,param_list)
{
  params=data.frame(a=3)
  content=scan(wasim_param_file,what="character",sep="\n",blank.lines.skip=TRUE,comment.char="#",strip.white=TRUE,quiet=TRUE)

  if (any(param_list %in% c("t0r","t0","c0")))  #extract snow parameters
  {
    start_line=which(content=="[snow_model]")  #find [snow_model] section
    params$t0r=content[start_line+5]
    params$t0=content[start_line+6]
    params$c0=content[start_line+9]
    params$c1=content[start_line+10]
  }

  if (any(param_list %in% c("m","tkorr","kkorr","kd","hmax","kh","cmelt")))  #extract soilmodel parameters
  {
    start_line=which(content=="[soil_model]")  #find [soil model] section
    params$m    =content[start_line+54]
    params$tkorr=content[start_line+55]
    params$kkorr=content[start_line+56]
    params$kd   =content[start_line+57]
    params$hmax =content[start_line+58]
    params$kh   =content[start_line+59]
    params$cmelt=content[start_line+67]
  }


  requested_pars=names(params) %in% param_list  #find reqested parameters that have been extracted
  params=params[,requested_pars]                #only return requested pars
  return(params)
}
