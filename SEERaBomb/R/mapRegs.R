mapRegs<-function(code=NA ){
  desc=c("1501"="San Francisco-Oakland SMSA (1973)",  "1502"="Connecticut (1973)",  "1520"="Metropolitan Detroit (1973)",
         "1521"="Hawaii (1973)",    "1522"="Iowa (1973)",    "1523"="New Mexico (1973)",
         "1525"="Seattle (Puget Sound) (1974)",    "1526"="Utah (1973)",    "1527"="Metropolitan Atlanta (1975)",
         "1529"="Alaska*",   "1531"="San Jose-Monterey",   "1535"="Los Angeles",    "1537"="Rural Georgia",
         "1541"="Greater California (excluding SF, LA & SJ)",     "1542"= "Kentucky",
         "1543"="Louisiana",      "1544"=  "New Jersey",       "1547"= "Greater Georgia (excluding AT and RG)")
  sym=c("1501"="sf","1502"="CT","1520"="dM","1521"="HI","1522"="IA","1523"="NM","1525"="sW","1526"="UT","1527"="aG",
        "1529"="AK", "1531"="sj", "1535"="la", "1537"="rG",
        "1541"="CA", "1542"="KY", "1543"="LA", "1544"="NJ", "1547"="GA")
  regs=data.frame(sym,desc)
  if (is.na(code[1])) return(regs) else return(regs[as.character(code),"sym"])
}  
