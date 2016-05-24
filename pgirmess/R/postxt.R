postxt<-function(cd="ul"){
    if(!is.numeric(cd)){ 
    cd<-switch(cd,
             ul = c(0.025,0.985),
             bl = c(0.025,0.025),
             ur = c(0.985,0.985),
             br = c(0.985,0.025),
             uc = c(0.5,0.985),
             bc = c(0.5,0.025),
             ml = c(0.025,0.5),
             mr = c(0.985,0.5),
             mc = c(0.5,0.5),
             stop("cd must be a numerical vector (length 2, 0<=values<=1, or choosen among 'ul','bl','ur','br','uc','bc','ml','mr,'mc'")
             )
    }
    else{
        if (any(length(cd)!=2,cd[1]>1,cd[1]<0,cd[2]>1,cd[2]<0)) {
        stop("if cd is a numerical vector, its length must be 2 and values comprised between 0 and 1")
        }
    }
    bxy<-par()$usr  
    xy<-list(x=bxy[1]+cd[1]*diff(bxy[1:2]),y=bxy[3]+cd[2]*diff(bxy[3:4]))
    xy
}
