if(is.element("X",ls())) warning("X was replaced",call. = FALSE)
if(is.element("pdata",ls())) warning("pdata was replaced",call. = FALSE)
if(is.element("pos1.cell.counter",ls())) warning("pos1.cell.counter was replaced",call. = FALSE)

load("ACL394data.rda")
X$images$path<-factor(system.file('img', package='RcellData'))

X$data$QC<-rep(TRUE,times=dim(X$data)[1])
X$data<-subset(X$data,select=-c(time.min,f.total.y,n.tot,AF.nM))
X$QC.history<-list()
X$subset.history<-list()
X$transform<-list()
X$variables$transformed<-NULL
X$variables$merged<-NULL
X$variables$merged.by<-NULL
X$variables$as.factor<-c("pos","cellID","ucid")
X$variables$all<-setdiff(X$variables$all,c("time.min","f.total.y","n.tot","AF.nM"))