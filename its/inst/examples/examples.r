#a cautionary tale
newIts(start="2003-10-28",end="2004-10-05",by="month") #OK
newIts(start="2003-10-30",end="2004-10-05",by="month") #unexpected
newIts(start="2003-10-31",end="2004-10-05",by="month") #unexpected
newIts(start="2003-09-31",end="2004-10-05",by="month") #unexpected
newIts(start="2003-09-31",end="2005-05-05",by="day",period="month",find="last",extract=TRUE) #unexpected
newIts(start="2003-09-31",end="2005-05-05",by="DSTday",period="month",find="last",extract=TRUE) #OK
