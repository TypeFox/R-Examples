
check.state<-function(state,statefips=FALSE){
check.state.aux<-function(state,statefips){
data("countyfips",envir = parent.frame())
assign("temp",countyfips)
assign("countyfips",temp)

state<-tolower(state)
if(!statefips & nchar(state)==2){
		state<-countyfips$statename[countyfips$acronym==tolower(state)][1]

}else if(statefips){
		state<-countyfips$statename[substr(countyfips$fips,1,2)==state][1]
	}else{
		}

if(state%in%unique(countyfips$statename)==FALSE & statefips==FALSE){
	return(NULL)
	}
state
}
out<-check.state.aux(state,statefips)
}