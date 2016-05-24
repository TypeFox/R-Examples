WriteMatrix = 
function (x, File = "", sep = "|",row.names = TRUE , row.name=NULL, ID=TRUE, ID.name=NULL, col.names=TRUE, col.width=6, na=NA)
{
       x[is.na(x)]<-""
    x <- as.matrix(x)
if (sep!=" "){
x[,ncol(x)]=paste(x[,ncol(x)], sep,sep="")
colnames(x)[ncol(x)]=paste(colnames(x)[ncol(x)], sep,sep="")
}


if (is.null(rownames(x))) row.names=FALSE;      ##############################

if (row.names)  {   
                        rn<-format(rownames(x), justify="l");
                        x <- cbind(NA, x);
                        x[,1]<-rn;
                        if(!is.null(row.name)) colnames(x)[1]<-format(row.name, width=nchar(x[1,1]),justify="l");
                        x[,1]<-format(x[,1],justify="l");
                        }



if (is.null(colnames(x))) col.names=FALSE;      ##############################

if (ID)                 {
                        id<-format(1:(nrow(x)), justify="r");
                        x <- cbind(NA, x);
                        x[,1]<-id;                      
                        if(!is.null(ID.name)) colnames(x)[1]=ID.name;###############################
                        }

if(col.names)   x <- rbind(colnames(x), x);
                     
  x[is.na(x)]<-na
                                                                                   
if (ID) {apply(x[,-c(1,2)],2,format, width=col.width,justify = "r")->x[,-c(1,2)]
x[,1]<-format(x[,1], justify="left");
if(row.names) x[,2]<-format(x[,2], justify="left");
if(col.names||row.names)    x[,2]<-format(x[,2], justify="left");
}
else {
apply(x,2,format, width=col.width,justify = "r")->x
}
                                                                                                        
cat(c( t(x)), file = File, sep = c(rep(sep, ncol(x) - 1), "\n"))
}
