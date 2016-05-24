`GetRakeSense` <-
function(uaz, upl, vaz,  vpl , paz,  ppl,  taz, tpl)
{
    A1 = mc2cart(uaz, upl);
    A2 = mc2cart(vaz, vpl);
  x1=  A1$x;
  y1 = A1$y;
  z1 = A1$z;

  x2 = A2$x;
  y2 = A2$y;
  z2 = A2$z;

    x3 = x1+x2;
    y3 = y1+y2;
    z3 = z1+z2;
    A1 = mc2cart(paz, ppl);
    A2 = mc2cart(taz, tpl);
  x1=A1$x;
  y1 = A1$y;
  z1 = A1$z;

  x2 = A2$x;
  y2 = A2$y;
  z2 = A2$z;



    dprodp =  (x1*x3+y1*y3+z1*z3)/(sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x3*x3+y3*y3+z3*z3));
    dprodt =  (x2*x3+y2*y3+z2*z3)/(sqrt(x2*x2+y2*y2+z2*z2)*sqrt(x3*x3+y3*y3+z3*z3));

    if(is.null(dprodt)) { print("NULL   dprodt"); return(1) }
    if(is.null(dprodp)) { print("NULL   dprodp"); return(1) }

    if(length(dprodt)<1) { print("No   dprodt"); return(1) }
    if(length(dprodp)<1) { print("N0   dprodp"); return(1) }

  ####   print(paste(sep=' ', "NULL   dprodt",dprodt, dprodp ));
      if(dprodt>dprodp)	{  ang2 = 1.0; }   else { ang2 = -1.0;}

      return(ang2);

}

