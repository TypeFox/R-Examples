maglab <-
function(lims, n, log=FALSE, exptext=TRUE, crunch=TRUE, logpretty=TRUE, usemultloc=FALSE, multloc=c(1,2,5), prettybase=10, hersh=FALSE, trim=FALSE){
if(usemultloc & log==F){stop('If using multloc then log must be TRUE!')}
lims=lims/(prettybase/10)
if(log & usemultloc==F){lims=log10(lims)}
if(usemultloc==F){if(missing(n)){labloc=pretty(lims)}else{labloc=pretty(lims,n)}}
if(log){
    if(usemultloc==F){
	    labloc=labloc+log10(prettybase/10)
	    labloc=labloc[round(labloc -log10(prettybase/10),10) %% 1==0]
      if(min(labloc)>lims[1]){labloc=c(min(labloc)-1,labloc)}
	    if(max(labloc)<lims[2]){labloc=c(labloc,max(labloc)+1)}
      labloc=round(labloc,10)
      labloc=10^labloc
      tickloc=labloc
	}
    if(usemultloc){
        labloc={}
        for(i in 1:length(multloc)){labloc=c(labloc,multloc[i]*10^seq(ceiling(log10(lims[1]))-1,floor(log10(lims[2]))+1))}
        labloc=sort(labloc)
        tickloc={}
        for(i in 1:9){tickloc=c(tickloc,i*10^seq(ceiling(log10(lims[1]))-1,floor(log10(lims[2]))+1))}
        tickloc=sort(tickloc)
    }
    #annoyingly I get weird issues for some numbers (e.g 0.00035) if they are in an otherwise scientific format list, and this behaves differently to the formatting on the actual plots. Only way round this is to format each number individually.
    char={}
    if(exptext){for(i in 1:length(labloc)){char=c(char,format(labloc[i]))}}
    if(! exptext){for(i in 1:length(labloc)){char=c(char,format(log10(labloc[i])))}}
}else{
labloc=labloc*(prettybase/10)
tickloc=labloc
char={}
for(i in 1:length(labloc)){char=c(char,format(labloc[i]))}
}

if(log & usemultloc==F){lims=10^(lims)}
if(trim){
char=char[labloc>=lims[1] & labloc<=lims[2]]
labloc=labloc[labloc>=lims[1] & labloc<=lims[2]]
tickloc=tickloc[tickloc>=lims[1] & tickloc<=lims[2]]
}

check=grep('e',char)
if(length(check)>0){
    char=format(labloc,scientific=T)
    check=grep("0e+00",char,fixed=T)
    char[check]="0"
    if(hersh){
    check=grep("e+0",char,fixed=T)
    char[check]=sub('e+0','e+',char[check],fixed=T)
    check=grep("e-0",char,fixed=T)
    char[check]=sub('e-0','e-',char[check],fixed=T)
    check=grep('e+',char,fixed=T)
    char[check]=paste(sub('e+','\\mu10\\sp',char[check],fixed=T),'\\ep',sep='')
    check=grep('e-',char,fixed=T)
    char[check]=paste(sub('e-','\\mu10\\sp-',char[check],fixed=T),'\\ep',sep='')
    }else{
    check=grep('e+',char,fixed=T)
    char[check]=paste(sub('e+','*x*10^{',char[check],fixed=T),'}',sep='')
    check=grep('e-',char,fixed=T)
    char[check]=paste(sub('e-','*x*10^{-',char[check],fixed=T),'}',sep='')
    }
}
if(crunch){
check=grep('1*x*',char)
    if(length(check)>0){
        if(hersh){
            char[check]=sub('1\\mu','',char[check],fixed=T)
        }else{
            char[check]=sub('1*x*','',char[check],fixed=T)
        }
    }
}
if(hersh){exp=char}else{exp=parse(text=char)}
return=list(tickat=tickloc,labat=labloc,exp=exp)
}
