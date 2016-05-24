Calc.Expr <-
function(filename, con_group_name, tr_group_name, 
                    con_gene_name= "18s",tr_gene_name, 
                    repeatnum = 3, ... )
{
#load data
  rawdata<-read.csv(filename)
#subset data
  unknown<-   subset(rawdata, Type != "NTC")
  num    <-   length(unknown$Type)
  groupnum <- num/repeatnum
  con_group<- subset(unknown, Identifier == con_group_name)
  tr_group <- subset(unknown, Identifier == tr_group_name)
  con_group_con<- subset(con_group, Fluor == con_gene_name)
  tr_group_con <- subset(tr_group, Fluor == con_gene_name)
  con_group_tr <- subset(con_group, Fluor == tr_gene_name)
  tr_group_tr  <- subset(tr_group, Fluor == tr_gene_name)
#read CT value, convert the factors into CT number
  con.con<-con_group_con[6]
  con.con<-as.numeric(paste(con.con[[1]]))
  tr.con<- tr_group_con[6]
  tr.con<- as.numeric(paste(tr.con[[1]]))
  con.tr<- con_group_tr[6]
  con.tr<- as.numeric(paste(con.tr[[1]]))
  tr.tr <- tr_group_tr[6]
  tr.tr <- as.numeric(paste(tr.tr[[1]]))
#calculate 
  cal.CT <- cal_ct(con.con,tr.con,con.tr,tr.tr)
#output data
  names(cal.CT) <- c(con_group_name , tr_group_name)
#return data
  return(cal.CT)

}

