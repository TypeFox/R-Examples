

inlineCxxPlugin <- Rcpp:::Rcpp.plugin.maker(
    include.before = paste("#include<RcppEigen.h>",
                           "#include<covafill/Tree>",
                           "",
                           "typedef Eigen::Array<double,Eigen::Dynamic,1>  cVector;",
                           "typedef Eigen::MatrixXd cMatrix;",
                           "",
                           sep = "\n"),
    LinkingTo=c('covafillr','RcppEigen','Rcpp'),
    package = "covafillr"
)

