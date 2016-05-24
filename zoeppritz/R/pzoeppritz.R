`pzoeppritz` <-
function( chtype= "Amplitude" , alpha1, alpha2, beta1, beta2, rho1 ,rho2, chincw="P", choutkind="ALL")
  {
    ##X##    # %   plot the output of zoep
    ##X##   #   plotzpr( "Amplitude" , alpha1, alpha2, beta1, beta2, rho1 ,rho2, "SV", "S");

##  chtype = type of output: "Amplitude", "Potential", "Energy"
    if(missing(chtype)) { chtype= "Amplitude"}
    if(missing(chincw)) { chincw="P" }
    if(missing(choutkind)) { choutkind="ALL" }

    

    ###########  protect against upper/lower case and use only first letter
    lowtype = tolower(chtype)
    type1 = substr(lowtype,1,1)

    if(type1=="a" )
      {
        type = 1;
        chtype= "Amplitude"
      } else
    if( type1=="p")
      {

        type = 2;
        chtype= "Potential"
      }else
    if(type1=="e" )
      {

        type = 3;
        chtype="Energy"
      }

    chincw = toupper(chincw)
    if(chincw == "P")
      {
        incw = 1;

      }else
    if(chincw == "S")
      {
        incw = 2;
      }else
    if(chincw == "SV")
      {
        incw = 2;
      }

    
 choutkind = toupper(choutkind)
    if(choutkind=="SV")
      {
        outkind = 2;
      }else
    if(choutkind=="S")
      {
        outkind = 2;
      }else
    if(choutkind=="P")
      {
        outkind = 1;
      }
    if(choutkind=="ALL")
      {
        outkind = 5;
      }
    if(choutkind=="NONE")
      {
        outkind = 0;
      }
     
    


    A = zoeppritz( type , alpha1, alpha2, beta1, beta2, rho1 ,rho2, incw);
    A$chincw = chincw
    A$chtype  =  chtype

    
    A$alphacrit = NULL

    if(chincw=="P")
      {
        if(alpha1<alpha2)
          {
            A$alphacrit = 180*asin(alpha1/ alpha2)/pi
          }
        else
          {
            A$alphacrit = NULL

          }
       
      }

   if(chincw=="S")
      {
        if(beta1<beta2)
          {
            A$alphacrit = 180*asin(beta1/ beta2)/pi
          }
        else
          {
            A$alphacrit = NULL

          }
       
      }
##########   plotting decisions:

    if(outkind>4)
      {
        plotzoeppritz(A)
        
        u = par("usr")
        LL = list(x=c(u[1]+0.76*(u[2]-u[1]) , u[1]+0.98*(u[2]-u[1])   ) , y = c(u[3]+0.5*(u[4]-u[3]) , u[3]+0.7*(u[4]-u[3]) ) ) 
       piczoeppritz(LL, chincw)
      }
    else
      {

        if(outkind>0) 
           {
             plot(A$angle, A$rmat[,outkind], col=4, type="l", xlab="Angle",
                  ylab=chtype, ylim=c(min(c(-.1, min(A$rmat[, outkind]))) ,1.2 ) ) ;
             
             
             title(paste(sep="", chincw, " Incident/", choutkind, " Out" ) );
             
             text(0 , 1.15, paste('Layer 1: Vp=',alpha1,' Vs=',beta1, ' Rho=', rho1), pos=4)
             text(0, 1.05, paste('Layer 2: Vp=',alpha2,' Vs=',beta2, ' Rho=', rho2), pos=4)
             
           }

  }

    
    invisible(A)
  }

