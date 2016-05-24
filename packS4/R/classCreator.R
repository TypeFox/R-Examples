classCreator <- function(name,field,type){

    lengthField <- length(field)
    nameLowerCase <- paste(tolower(substr(name,1,1)),substring(name,2),sep="")

    d1 <- ' # # # # # # # # # # # # # # # # # #\n# # #     Class  definition     # # #\n'
    d2 <- paste(' # # #',paste(rep(' ',12-floor(nchar(name)/2)),collapse=""),name,
      paste(rep(' ',13-ceiling(nchar(name)/2)),collapse=""),'# # #\n',sep="",collapse="")
    d3 <- "  # # # # # # # # # # # # # # # # # \n\n"
    d4 <- paste('.',name,'.validity <- function(object){
    cat("**** validity ',name,' <empty> ****\\n")
    return(TRUE)\n}\n\n\n',sep='')
    cat(d1,d2,d3,d4,sep="")

    t0 <- 'cat("### Class definition ###\\n"'
    t2 <- paste('setClass(\n    Class="',name,'",\n    representation=representation(\n',sep='')
    t3 <- paste('        ',field[-lengthField],'="',type[-lengthField],'",\n',sep='')
    t4 <- paste('        ',field[lengthField],'="',type[lengthField],'"\n',sep='')
    t5 <- '   ),\n    prototype=prototype(\n'
    t6 <- paste('        ',field[-lengthField],'=',type[-lengthField],'(),\n',sep='')
    t7 <- paste('        ',field[lengthField],'=',type[lengthField],'()\n',sep='')
    t8 <- paste('   ),\n    validity=.',name,'.validity\n)\n\n\n',sep='')
    cat(t2,t3,t4,t5,t6,t7,t8,sep="")

    e0 <- 'cat("### Constructor ###\\n")\n'
    e1 <- paste(nameLowerCase,' <- function(',sep="")
    e2 <- paste(field[-lengthField],',',sep='')
    e3 <- paste(field[lengthField],'){\n',sep='')
    e4 <- paste('    return(new("',name,'"',sep='')
    e5 <- paste(',',field,"=",field,sep='')
    e6 <- paste('))\n}\n\n\n')
    cat(e0,e1,e2,e3,e4,e5,e6,sep="")

    a0 <- 'cat("### Show ###\\n")\n'
    a1 <- paste('.',name,'.show <- function(object){\n     cat("   ~~~ Class:",class(object),"~~~ ")\n',sep='')
    a2 <- paste('    cat("\\n ~ ',field,' : ",object@',field,')\n',sep='')
    a3 <- paste('    cat("\\n")\n    return(invisible())\n}\nsetMethod(f="show",signature="',name,'",definition=.',name,'.show)\n\n\n',sep='')
    cat(a0,a1,a2,a3,sep="")

    b0 <- 'cat("### Getteur ###\\n")\n'
    b1 <- paste('.',name,'.get <- function(x,i,j,drop){\n    switch(EXPR=i,\n',sep='')
    b2 <- paste('       "',field,'"={return(x@',field,')},\n',sep='')
    b3 <- paste('       stop("[',name,':get] ",i," is not a \'',name,'\' slot")\n    )\n    return(invisible())\n}\n',sep='')
    b4 <- paste('setMethod(f="[",signature="',name,'",definition=.',name,'.get)\n\n\n',sep='')
    cat(b0,b1,b2,b3,b4,sep="")

    c0 <- 'cat("### Setteur ###\\n")\n'
    c1 <- paste('.',name,'.set <- function(x,i,j,value){\n    switch(EXPR=i,\n',sep='')
    c2 <- paste('       "',field,'"={x@',field,'<-value},\n',sep='')
    c3 <- paste('       stop("[',name,':set] ",i," is not a \'',name,'\' slot")\n',sep='')
    c4 <- paste('    )\n    validObject(x)\n    return(x)\n}\nsetMethod(f="[<-",signature="',name,'",definition=.',name,'.set)\n',sep='')
    cat(c0,c1,c2,c3,c4,sep="")
}


