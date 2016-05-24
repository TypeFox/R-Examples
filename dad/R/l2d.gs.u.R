l2d.gs.u <-
function(x1,x2,check=FALSE)
{
	var1<-var(x1);
  var2<-var(x2);
  if(check)
    {if(var1<.Machine$double.eps |var2<.Machine$double.eps) 
      {stop("At least one variance is zero") 
      }
    }
  m1<-mean(x1);
	m2<-mean(x2);
	return(l2d.gp.u(m1,var1,m2,var2))    
}

