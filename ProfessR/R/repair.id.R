repair.id<-function(sisroster, scrfn)
{

####   sisroster = list(ID=number,  lastname='last name of student',  fullname='full name of student')
 ####   scrfn = list(ID=number, nam="name on scantron")

  
  ###########   fix ID = repair.id
  
#########   fix ID/namews that do not match the roster list
######### sisroster -> roster from BLACK Board
#########  scrfn    -> list of scantron

###  scrfn = list(ID, nam, 
###  sisroster  list( lastname  , fullname


    scan.names = as.vector(scrfn$nam)
    scan.ID  = as.vector(scrfn$ID)

    if(length(sisroster$fullname)<1)
      {
        sisroster$fullname = paste(sep=' ',sisroster$firstname, sisroster$lastname)

      }

    BBlast = tolower(sisroster$lastname)
    BBfull = tolower(sisroster$fullname)
    BBID = sisroster$ID
    
  
  newid = as.vector(scan.ID)
    
  m1 = match(BBID, scan.ID)

  m2 = match(scan.ID, BBID )

  w1 = which(is.na( BBID[m2]))


## cbind(BBID[m2], scan.ID, scan.names)

 

 ##   BBlast

    
  if(length(w1)<1) { print("no problems"); return() }


    
  print("list of mismatches:")
  print( cbind( scan.ID[w1] , scan.names[w1]) )

    

  for(i in 1:length(w1))
    {

      thename = as.character( scan.names[w1[i]] )
      
      print(paste(sep=" ", i,  w1[i],  scan.names[w1[i]], "=SCANNED NAME"))

      
      Lname = tolower(unlist(strsplit(scan.names[w1[i]], split=" ")))


       
      if(length(Lname) == 1 )
        {
          print(paste("This is a One-Name ID, searching for repair:",Lname ))
          M = 0
          for(m in 1:length(BBlast))
            {
              tname = tolower(BBlast[m])
              gnum = grep(tname , Lname)
              ##print(paste(sep=' ', m, gnum, tname, Lname))
              
              if(length(gnum)>=1)
                {
                  
                  M = m
                  break


                }


            }

          Lname=c(BBlast[M], sisroster$firstname[M])

        }

      
      gN1 = grep(Lname[1], BBlast)
      gN2 = grep(Lname[2], BBlast)

      g1 = c(gN1, gN2) 

      if(length(g1)<1)
        {
          lameo = Lname[1]
          name15 = substr(lameo, 1,5)
          namend = substr(lameo, nchar(lameo)-5, nchar(lameo))
          
          g1 = grep(name15, BBlast)
          g2 = grep(namend, BBlast)
          g1 = c(g1, g2)
          if(length(g1)<1)
            {
              name15 = substr(lameo, 1,1)
              g1 = grep(name15, BBlast)
            }
        }


      
      if(length(g1)==1)
	{
          cat("---FOUND MATCH - FIXING---------", sep="\n")
          print(data.frame(cbind(scan.names[w1[i]], scan.ID[w1[i]], BBlast[g1], BBfull[g1]	, BBID[g1])))
          newid[w1[i]]=BBID[g1]
        }
      else
        {
          cat("---------------------------", sep="\n")

          cat("working on this dude:", sep="\n")

          cat( cbind( scan.ID[w1[i]] , scan.names[w1[i]]) , sep="\n" )
          cat("--------LIST--------", sep="\n")
          cat("possible matches:", sep="\n")

          pname = paste(sep=":", 1:length(g1), BBfull[g1])

          cat(pname, sep="\n")

          cho = readline("possible matches: ")

          if(cho=="") { next }
          else
            {
              icho = as.numeric(cho)
              
              newid[w1[i]]=BBID[g1[icho]]
            }

        }
    }

  return(newid)

}
