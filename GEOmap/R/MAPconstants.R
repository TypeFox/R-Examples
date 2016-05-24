`MAPconstants` <-
function()
  {

  ###  MAPconstants()$A.MAPK
  ###  MAPconstants()$R.MAPK

###  most projection information is taken from the snyder book
    DEG2RAD=pi/180
    RAD2DEG=180/pi
    
    A.MAPK=6378206.4

###   page 13 in snyder
###  datum: clark 1866:   
    E2.MAPK=0.006768658

###   datum: GRS 80:
    E2.GRS80=0.0066943800

    
    E.MAPK=0.0822719
    E1.MAPK=0.993231340
    TwoE.MAPK=0.164543800
    R.MAPK=6378.2064 
  
    FEET2M= 0.3048
    M2FEET=1/FEET2M

    

rlist = list(DEG2RAD=DEG2RAD,
  RAD2DEG=RAD2DEG, A.MAPK=A.MAPK,
  E2.MAPK=E2.MAPK,
  E2.GRS80=E2.GRS80,
  E.MAPK=E.MAPK,
  E1.MAPK=E1.MAPK,
  TwoE.MAPK=TwoE.MAPK,
  R.MAPK=R.MAPK,
  FEET2M=FEET2M,
  M2FEET=M2FEET)

    
return(rlist)
    
  }

