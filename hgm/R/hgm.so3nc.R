# $OpenXM: OpenXM/src/R/r-packages/hgm/R/hgm.so3nc.R,v 1.5 2014/03/22 00:15:40 takayama Exp $
"hgm.tk.so3nc" <-
function(a=1,b=2,c=3,t0=0.0,q=1,deg=0) { 
#  library.dynam("mypkg",package="mypkg",lib.loc="./");
  if (!is.loaded("hgm")) library.dynam("hgm",package="hgm",lib.loc=NULL);
  .C("so3_main",as.double(a),as.double(b),as.double(c),
     as.double(t0),
     as.integer(q),as.integer(deg),
     retv=as.double(0),PACKAGE="hgm")$retv
}
