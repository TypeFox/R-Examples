ConvP2Vdepthinv <-
function(rpd,v,r)
  
  
{
res=NaN
if((!is.na(rpd))&(rpd!=0)){
    
    
    ident=which(r==rpd) 
                        
                        
                        
    aboves=which(r>=rpd) 
    if(length(aboves)>0){
        
        anz=length(aboves)
        
           if(length(ident)==0){
              top=aboves[anz]
              bottom=aboves[anz]+1
           }else if(length(ident)==1){
              top=aboves[anz-1]
              bottom=aboves[anz]
           }else if(length(ident)==2){
              top=aboves[anz-2]
              bottom=aboves[anz-1]
           }else{
              stop('ConvP2Vdepthinv: unexpected value of ident')
        } 
        v1=approx(c(r[top],r[bottom]),c(v[top],v[bottom]),rpd)$y
        
        if( v1==0){
           res=NaN
        }else{
           res=rpd/v1
        } 
    }else{
        
        res=NaN
    } 
} 
p=res * pi/180  
return(p)
}

