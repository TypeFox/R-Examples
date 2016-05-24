InterpModel <-
function(model,newz=NULL,preserve=NULL)

{

eps = 10^-6 
  
if(is.null(newz)){
   
	
	dz=5 
	rpf=TransformS2Fz(model$rp,model$z[length(model$z)]-1,model$rp)[[2]]
	zf=seq(from=0,by=dz,to=rpf)
	zf=TransformF2Sz(zf,zf,model$rp)[[2]] 
	newz=zf
}
if(is.null(preserve)){preserve='preserve'}
preserve=tolower(preserve)



newz=sort(newz) 





if(preserve == 'preserve' | preserve == TRUE){
         
         
         
         discindy=which(diff(model$z)==0)
         discz=model$z[discindy] 
         
         discz=c(discz,
		model$conr,
                model$moho,
                model$d410,
                model$d520,
                model$d660,
                model$cmb,
                model$icb,
                model$dz)
         
         discz=discz[which(!is.na(discz))]
         
         discz=unique(sort(discz))
         
         
         
         

         for(indy in 1:length(discz)){
             indies=which(newz==discz[indy]); 
             newz[indies]=newz[indies]+NaN;  
         }
	 newz=newz[!is.na(newz)]


         
         newz=sort(c(newz, discz))


}


newmod=list()
newmod$z=numeric(0) 
newmod$vp=numeric(0)
newmod$vs=numeric(0)
newmod$rho=numeric(0)
newmod$qp=numeric(0) 
newmod$qs=numeric(0) 
newmod$conr=model$conr
newmod$moho=model$moho
newmod$d410=model$d410
newmod$d520=model$d520
newmod$d660=model$d660
newmod$cmb=model$cmb
newmod$icb=model$icb
newmod$dz=model$dz
newmod$dname=model$dname
newmod$year=model$year
newmod$rp=model$rp
newmod$name=model$name


indies=which(diff(model$z)==0) 
modstartindy=c(1,indies+1)
modendindy=c(indies,length(model$z))
pieceanz=length(modstartindy) 
indies=which(diff(newmod$z)==0)
newstartindy=c(1,indies+1)
newendindy=c(indies,length(newmod$z))
newpieceanz=length(newstartindy) 


method='linear'

for(piececnt in 1:pieceanz){
    
    modrange=modstartindy[piececnt]:modendindy[piececnt]
    

    if(newpieceanz>1){
       
       
       
       
       
       
       
       
       
       
       indies1=which(newmod$z[newstartindy]==model$z[modstartindy[piececnt]])  


       if(0!=length(indies1)){           
           
           
           
           firstnewsample=newstartindy[indies1[1]]

           
           
           subspace=(firstnewsample+1):length(newmod$z) 
           indies2=which(newmod$z[subspace]<=model$z[modendindy[piececnt]])

           if(0!=length(indies2)){
               
               
               lastnewsample=subspace[indies2[length(indies2)]]
               }else
               
               
               {lastnewsample=firstnewsample}


           
           
           newrange=firstnewsample:lastnewsample
                }else
	   
           
           
           {newrange=NULL}
	   }else{	          
       
       
             newrange=which((newz>=min(model$z[modrange] - eps))&(newz<=max(model$z[modrange] + eps)))
           }


    
    
    
    if(!is.null(newrange)){
       
       
       
       
       
       
       
       newmod$z=c(newmod$z,newz[newrange])

       
       if(sum(!is.na(model$vp[modrange])) > 0){
          
          interpresult=approx(model$z[modrange],model$vp[modrange],newz[newrange],method, rule = 2, ties = 'ordered')$y
          newmod$vp=c(newmod$vp,interpresult)
       }else
          
          {newmod$vp=c(newmod$vp,newz[newrange]*NaN)}


       
       if(sum(!is.na(model$vs[modrange])) > 0){
          
          interpresult=approx(model$z[modrange],model$vs[modrange],newz[newrange],method, rule = 2, ties = 'ordered')$y
          newmod$vs=c(newmod$vs,interpresult)
       }else{
          
          newmod$vs=c(newmod$vs,newz[newrange]*NaN)}

	
       if(sum(!is.na(model$rho[modrange])) > 0){
          
          interpresult=approx(model$z[modrange],model$rho[modrange],newz[newrange],method, rule = 2, ties = 'ordered')$y
          newmod$rho=c(newmod$rho,interpresult)
       }else{
          
          newmod$rho=c(newmod$rho,newz[newrange]*NaN)} 


       
       if(sum(!is.na(model$qp[modrange])) > 0){
          
          interpresult=approx(model$z[modrange],model$qp[modrange],newz[newrange],method, rule = 2, ties = 'ordered')$y
          newmod$qp=c(newmod$qp,interpresult)
       }else{
          
          newmod$qp=c(newmod$qp,newz[newrange]*NaN)}

       
       if(sum(!is.na(model$qs[modrange])) > 0){
          
          interpresult=approx(model$z[modrange],model$qs[modrange],newz[newrange],method, rule = 2, ties = 'ordered')$y
          newmod$qs=c(newmod$qs,interpresult)
       }else{
          
          newmod$qs=c(newmod$qs,newz[newrange]*NaN)}
}
}



return(newmod)}

