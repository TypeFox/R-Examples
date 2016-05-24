inf <-
function(inpr, ..., oppf=5, opdf=1, oprho=0, dfwarning=0,plot=0,others) {

#  ============== begin by defining dealwithrestd ==============================

dealwithrest<-function(inpr, therest) {

# This does the items given to inf that I don't format myself. So dealwithrest prints some things, and checks on existence of inpr$others.

# The reason some output parameters have "op" in front is to avoid confusion with internal variables of same name --
#       consider at some stage renaming the internals --- although rho is external too, so more complicated.

printstrcomp<-"if (is.null(inpr$item)) cat(paste('\nTo save this information, specify item when you invoke phyreg\n')) else
print(inpr$item)"

savestrcomp<-"if (is.null(inpr$item)) cat(paste('To save this information, specify item=TRUE or 1 when you invoke phyreg\n')) else
cat(paste('This information has been saved. Access as (phyreg-result)$item.\n')) "

nwrong<-0
for (ss in therest) {
 cat(paste("\n",ss,"\n"))
 if (ss %in% saveitemlist) eval(with(inpr,parse(text=gsub("item",ss,fixed=TRUE,x=savestrcomp))))  else 
 if (ss %in% printitemlist) eval(with(inpr,parse(text=gsub("item",ss,fixed=TRUE,x=printstrcomp)))) else
      {nwrong<-nwrong+1; cat(paste( "\n", ss, " is not recognised\n"))}
                }
   
    if (nwrong>0L) cat(paste("\nThe possible values are (each within quotation marks, though the first five can instead be of the form oppf=7):\n", paste(c(formatlist,  printitemlist, saveitemlist),collapse=", "),"\n"))

 }  # finish dealwithrest

# ========== done defining dealwithrest, and moving on  to define inputstoredinf : ================

#                NOTE only one defn now, at top level ==============


#  ============== done defining opsio, moving on to body of inf  ==============================


# first define some useful lists -- these might or might not be saved. Tell the user and explain how to
saveitemlist<-c( "lmshortx", "lmshortxz", "lmlongx", "lmlongxz", "paper", "linputs", "sinputs", "hinput")

# print these on request
printitemlist<-c("parmx", "parmxz","opfunccall","means", "fullphy", "usedphy", "originalIDs")

# these I format -- included here only to say to the user they are permitted values
 formatlist<-list("oppf","opdf","oprho","dfwarning","plot")

if (missing(inpr)) cat(paste("You need an argument in position 1 or called as inpr= that has been assigned as in x<-phyreg(... relevant parameters...)\n")) else if (class(inpr)!="phyreglm") {cat(paste("The first (or inpr=) argument, say x, should have been assigned as in x<-phyreg(... relevant parameters...)\n")) }

else {

# Copy the call, and lose the head, redefining as a list

mycall<-match.call()
mycall<-tail(as.list(mycall),-2)

# Check there is something in the list

if (length(mycall)==0) {
  cat(paste("You haven't said what information you'd like about the stored phyreg output:\n"))
  cat(paste("\nThe possible values are (each within quotation marks, though the first four can instead be of the form oppf=7):\n", paste(c(formatlist,  printitemlist, saveitemlist),collapse=", "),"\n"))}
  
 else {
 
# First, the things I have to format

req<-list()
if (!missing(oppf)) {if (oppf>0L) req$oppf<-oppf} else if ("oppf" %in% mycall) req$oppf<-5
if (!missing(opdf)) {if (opdf>0L) req$opdf<-opdf} else if ("opdf" %in% mycall) req$opdf<-1
if (!missing(oprho)) {if (oprho>0L) req$oprho<-oprho} else if ("oprho" %in% mycall) req$oprho<-5
if (!missing(dfwarning)) { if (dfwarning>0L) req$dfwarning<-dfwarning} else if ("dfwarning" %in% mycall) req$dfwarning<-1
if (!missing(plot)) { if (plot>0L) req$plot<-plot} else if ("plot" %in% mycall) req$plot<-1

if (sum(as.numeric(req))>0L) opsi(inpr,requests=req,originalIDs=inpr$originalIDs)

# Then, the rest.  

 formatlist<-list("oppf","opdf","oprho","dfwarning","plot")
 mycallrest<-mycall[!(names(mycall) %in% formatlist) & !(mycall %in% formatlist)]
 
 do.call(dealwithrest, list(inpr=inpr,therest=mycallrest))
 
 } # end of else that checks for something being asked
 
 }  # end of the else that checks the first argument is of class phyreglm

 }
