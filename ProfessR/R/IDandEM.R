IDandEM<-function(scrfn,sisroster, sel=1:2,  hnote="Exam Results", SEND=TRUE)
{
####   sisroster = list(ID=number,  lastname='last name of student',  fullname='full name of student')
 ####   scrfn = list(ID=number, nam="name on scantron")

 ####  sel=1:2,    select only a few to send (specify for a specific student)
  #### hnote="Exam Results",  subject line on E-mail
  #### SEND=TRUE         logical, if false, do not send
  
if(missing(sel)) { sel = 1:length(scrfn$ID) }
if(missing(hnote)) {  hnote="Exam Results"  }
if(missing(SEND)) { SEND=TRUE   }

  for(i in sel)
    {

      m1 = match(scrfn$ID[i],sisroster$ID)
     

      if(length(m1)>0 & !is.na(m1) )
        {

          em =  sisroster$email[m1]
          
          cat(paste(sep=" ", i,   scrfn$tfile[i], em), sep="   " )
          pfile = autoemail(em,scrfn$tfile[i], hnote=hnote )
          cat(paste(sep=" ",  pfile), sep="\n" )
          system(paste(sep=" ", "chmod +x",  pfile))

          ##
          if(SEND)
            {
              system(pfile)
            }

        }

    }


}


