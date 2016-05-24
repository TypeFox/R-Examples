fel.repeated <-
function(x,y=NULL,subjects=NULL,repeated=NULL,subjects.in="all",repeated.in="all" ,...) {

  if (subjects.in[1]=="all" & !is.null(subjects)) {
    subjects <- factor(subjects)
subjects.in <- levels(subjects)}
if (repeated.in[1]=="all" & !is.null(repeated)) {
  repeated <- factor(repeated)
  repeated.in <- levels(repeated)}

  if (is.null(subjects) & is.null(repeated)) {
    subjects2<-NULL
    subset<-NULL}
   else if (is.null(subjects)) {
      subjects2<-repeated
      subset<-repeated %in% repeated.in}
  else if (is.null(repeated)) {
    subjects2<-subjects
    subset<-subjects %in% subjects.in}
  else  {
    subjects2<-list("subjects"=subjects,"repeated"=repeated)
    subset<-((subjects %in% subjects.in) & (repeated %in% repeated.in))}
  fel(x,y,subjects=subjects2,subset=subset,...)
        }
