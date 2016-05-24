`NEW.getUWSTAS` <-
function(PICS)
  {
    
    STAS = list(tag="", name="", comp="", c3="", phase="" , sec=0, err=0,  pol="", flg=0,  res=0 )
    k =  0 
    for( i in 1:length(PICS))
      {
        card = PICS[i]
        AS1  = unlist(strsplit(card," "))
        AS1  = unlist(strsplit(card,split="\\("))
        
####  the first section is the station/comp identifier, parse this:
        
        AS2 = AS1[1]
        AS3 = unlist(strsplit(AS2,"\\."))
        
        name = AS3[2]
        
        c3 = AS3[3]
        c3  = unlist(strsplit(c3,split=""))
        c3  = paste(collapse="", c3[c3!=" "])
        
        
#############   here convert all Z,U or whatever to V=vertical
        comp = fixcomps(c3)
        
        
### ppacs = which(AS1=="(p")
        
####  loop through the tokens on the DOT card
        for(m in 2:length(AS1))
          {
            tok1 = unlist(strsplit(AS1[m],split="\\)"))
            tok = unlist(strsplit(tok1,split=" "))

            tok = tok[tok!=""]
            tok[tok=="_"] = NA
            
            tok.kind = tok[1]

           
            ######### print(tok)
            
            if(identical(tok.kind, "P") | identical(tok.kind, "p"))
              {
                k = k +1
                
                STAS$name[k] = name
                STAS$c3[k] =c3
                STAS$comp[k] =  comp
                STAS$tag[k] = paste(sep=".", STAS$name[k], STAS$c3[k])
                STAS$phase[k] = tok[2]
                STAS$pol[k] = tok[3]
                STAS$sec[k] = as.numeric(tok[4])
                STAS$flg[k] = as.numeric(tok[5])
                STAS$err[k] = as.numeric(tok[6])
                STAS$res[k] = as.numeric(tok[7])

                if( identical(STAS$pol[k], "_") | identical(STAS$pol[k], ".")) {   STAS$pol[k] =NA  }
                
              }
            ###########   if there are other kinds of tokens, these are not saved here
          }
        
      }
        return(STAS)

  }

