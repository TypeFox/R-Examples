StripRepetitions <-
function(phase)
{

phase=strsplit(phase,'')[[1]]
phaseout=phase
repetitions=NULL

w=getOption('warn')
options(warn=-1)
if(!is.na(as.numeric(phase[length(phase)]))){
    
    repstring=NULL
    done=0
    len=length(phase)
    indy=len
    while( done==0){
    	if(!is.na(as.numeric(phase[indy]))){
            repstring=c(phase[indy],repstring)
            indy=indy-1
            if( indy==0 ){
                done=1
            } 
        }else{
            done=1
        } 
    } 
    repetitions=as.numeric(paste(repstring,collapse=''))
    phaseout=paste(phaseout[1:indy],collapse='')
}else{
    repetitions=0
} 
options(warn=w)
return(list(paste(phaseout,collapse=''),repetitions))
}

