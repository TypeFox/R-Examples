`shade.col` <-
function(n, acol=c(1,0,0), bcol=c(1,1,1)  )
{

  if(missing(acol)) { acol=c(1,0,0) }
  if(missing(bcol)) { bcol=c(1,1,1) }
  
  if ((n <- as.integer(n[1])) > 0) 
    {
      r1 = acol[1];
      g1 = acol[2];
      b1 = acol[3];
      
      
      wr1 = bcol[1];
      wg1 = bcol[2];
      wb1 = bcol[3];
      
      dr1 = wr1-r1
      dg1 = wg1-g1
      db1 = wb1-b1
      
      hr1 = (wr1-r1)/n
      hg1 = (wg1-g1)/n
      hb1 = (wb1-b1)/n
      
      
      
      nr1 = seq(from=r1, length.out=n, by=hr1)
      nb1 = seq(from=b1, length.out=n, by=hb1)
      ng1 = seq(from=g1, length.out=n, by=hg1)
        
COL = rgb(red=nr1, green=ng1, blue=nb1)


      


	}	
	
    return(COL)

}

