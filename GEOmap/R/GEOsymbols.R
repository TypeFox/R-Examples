GEOsymbols<-function()
{

geotypes = c("contact", "anticline",
  "syncline",
  "OverTurned-ant","OverTurned-syn",
  "perp",
  "thrust",
  "normal",
  "dextral",
  "sinestral",
  "detachment", "bcars"
  )



  
  N = length(geotypes)
  ncol = 5
   nrow = round((N/ncol)+.5)



   dx = 1/ncol
    dy =  1/nrow

bx = dx*.2
by = dy*.2
plot(c(0,1), c(0,1), type='n',asp=1, axes=FALSE, xlab='', ylab='')

fin = par("fin")
pin = par("pin")
u = par("usr")
umm =   (u[4]-u[3])/pin[2]/25.4

G=list()
G$x=c(-1.0960,-0.9942,-0.8909,-0.7846,-0.6738,-0.5570,-0.4657,-0.3709,
-0.2734,-0.1740,-0.0734, 0.0246, 0.1218, 0.2169, 0.3086, 0.3956, 0.4641, 
0.5293, 0.5919, 0.6530, 0.7131)
G$y=c(-0.72392,-0.62145,-0.52135,-0.42599,-0.33774,-0.25896,-0.20759,
-0.16160,-0.11981,-0.08105,-0.04414,-0.00885, 0.02774, 0.06759, 0.11262, 
0.16480, 0.21487, 0.27001, 0.32895, 0.39044, 0.45319)



MYh = pin[1]/(u[2]-u[1])


   for(i in 1:length(geotypes) )
      {
        B =  RPMG::itoxyz(i, ncol, nrow, 1)
        x = (B$ix-1)*dx
        y = (B$iy-1)*dy
        rect(x , y , x+dx, y+dy, lty=1, col='white' )

        g1x = c(x+bx, x+dx-bx)
        g1y = c(y+by, y+dy-by)
        
     ###   segments(g1x[1], g1y[1],  g1x[2], g1y[2], col='blue')

       wig1x = RPMG::RESCALE(G$x, x+bx, x+dx-bx, -1, 1)
       wig1y = RPMG::RESCALE(G$y, y+by, y+dy-by, -1, 1)

        g = PointsAlong(wig1x, wig1y, N=3)


       text(x+bx*.1, y+by*.1, labels=geotypes[i], adj=c(0,0) )

        if(identical(geotypes[i], "contact"))
           {
             ###  
            ### segments(g1x[1],g1y[1],g1x[2],g1y[2], col='blue')
             lines(wig1x,wig1y,col='blue')
             
           }

        
        if(identical(geotypes[i], "anticline"))
           {
             ###  anticline(g1x,g1y,  N=1)
           ###  lines(wig1x,wig1y,col='blue')
             sk = 2
             SynAnticline(wig1x, wig1y, r1 = sk, r2 = sk, h1= 0, h2= 0,  N=2, syn=FALSE, endtol=.4)
           }

        if(identical(geotypes[i], "syncline"))
           {
             ###  anticline(g1x,g1y,  N=1)
                sk = 2
                SynAnticline(wig1x, wig1y, r1 = sk, r2 = sk, h1= 0, h2= 0,  N=2, syn=TRUE, endtol=.4)
            
           }
       if(identical(geotypes[i],"OverTurned"))
         {
           sk = 3

            OverTurned(wig1x,wig1y,  r1=0.8*sk, r2=sk, h1=.8*sk, h2=.8*sk, N=1, syn=FALSE, endtol=.2)

          
         }
        if(identical(geotypes[i],"OverTurned-ant"))
          {
           sk = 3

            OverTurned(wig1x,wig1y,  r1=0.8*sk, r2=sk, h1=.8*sk, h2=.8*sk, N=1, syn=FALSE, endtol=.2)

          }
        if(identical(geotypes[i],"OverTurned-syn"))
          {


            sk = 3

            OverTurned(wig1x,wig1y,  r1=0.8*sk, r2=sk, h1=.8*sk, h2=.8*sk, N=1, syn=TRUE, endtol=.2)

            
          }
 
      if(identical(geotypes[i],"perp"))
        {
          sk = 2
           ## segments(g1x[1], g1y[1],  g1x[2], g1y[2], col='blue')
          lines(wig1x,wig1y,col='blue')
          perpen(g$x,g$y, h=sk, rot=g$rot, lwd=1, col='blue' )

        }

        if(identical(geotypes[i],"thrust"))
          {
             sk = 2
              ## lines(wig1x,wig1y,col='blue')
            thrust(wig1x,wig1y, N=5, h=sk, REV=FALSE, col='blue', lty=1)
             

          }

        if(identical(geotypes[i],"normal"))
          {
             sk = 2
             lines(wig1x,wig1y,col='blue')
             normalfault(g$x,g$y,h=sk,hoff=sk, rot=g$rot , col='blue')


          }


             
        if(identical(geotypes[i],"dextral"))
          {
            sk = 1
           ###  g = PointsAlong(g1x,g1y, N=1)
            lines(wig1x,wig1y,col='blue')
            SSfault(g$x,g$y,h=sk,hoff=sk, rot=g$rot , col='blue')


          }
        if(identical(geotypes[i],"sinestral"))
          {
            sk = 1
          ###  g = PointsAlong(g1x,g1y, N=1)
            lines(wig1x,wig1y,col='blue')
            SSfault(g$x,g$y,h=sk,hoff=sk, rot=g$rot , col='blue', dextral=FALSE)

            
          }

        if(identical(geotypes[i],"detachment"))
          {
            
            sk = 2
    
           ### g = PointsAlong(g1x,g1y, N=1)
            lines(wig1x,wig1y,col='blue')
            
            horseshoe(g$x  , g$y , r1=sk, r2=sk-sk*.1, h2=0, h1=0, rot=g$rot , col='blue', fill=TRUE)
       
          }

      if(identical(geotypes[i],"bcars"))
          {
            
            sk = 3
    
           ### g = PointsAlong(g1x,g1y, N=1)
            lines(wig1x,wig1y,col='blue')
            
            bcars(g$x,g$y, h1=sk, h2=sk*.5, g$rot, col='black', border='black' )
       
          }
        

      }


}

### source("/home/lees/Progs/R_PAX/GEOmap/R/GEOsymbols.R")  ;  GEOsymbols()
### source("/home/lees/Progs/R_PAX/GEOmap/R/OverTurned.R")  ;  GEOsymbols()
