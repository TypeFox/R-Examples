DateConversion<-function(DateVar,DateIn,DateOut) {
  #=======================================================
  #
  #  TITLE:     Conversion between Different Date Format
  #  FUNCTION:  DateConversion()
  #  AUTHOR:    Zhijie Zhang
  #  DATE:      07 JANUARY 2010
  #  CALLS:     
  #  NEEDS:
  #  NOTES: #%Y--yyyy;#%m--number;#%d--number
  #  DateVar-Specify your date variable
  #  DateIn- Specify the date format for your original variable,such as "%m/%d/%Y" 
  #  DateOut-Specify the date format that you hope to output,such as "%d/%m/%Y" 
  #  CHANGE HISTORY:
  #=======================================================
 DateNew<-format(as.Date(as.character(DateVar),format=DateIn),DateOut)
 return(DateNew)
}