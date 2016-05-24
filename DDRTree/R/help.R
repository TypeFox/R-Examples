#'DDRTree - An algorithm to reduce dimensionality and learning principal graphs simultaneously
#'
#' This is an R and C code implementation of the DDRTree algorithm from Qi Mao, Li Wang et al. \cr \cr
#' Qi Mao, Li Wang, Steve Goodison, and Yijun Sun. Dimensionality Reduction via Graph Structure Learning.
#' The 21st ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD'15), 2015 \cr
#' \url{http://dl.acm.org/citation.cfm?id=2783309} \cr \cr
#' to perform dimension reduction and principal graph
#' learning simultaneously. Please cite this package and KDD'15 paper if you found DDRTree is useful for your research.
#'
#'@section Introduction:
#'The unprecedented increase in big-data causes a huge difficulty in data visualization and downstream analysis.
#'Conventional dimension reduction approaches (for example, PCA, ICA, Isomap, LLE, etc.) are limited in their ability
#'to explictly recover the intrinisic structure from the data as well as the discriminative feature representation,
#'both are important for scientific discovery. The DDRTree algorithm is a new algorithm to perform the following
#'three tasks in one setting: \cr
#'\cr
#'1. Reduce high dimension data into a low dimension space \cr
#'\cr
#'2. Recover an explicit smooth graph structure with local geometry only captured by distances of data
#'points in the low dimension space. \cr
#'\cr
#'3. Obtain clustering structures of data points in reduced dimension \cr
#'
#'@section Dimensionality reduction via graph structure learning:
#'
#'Reverse graph embedding is previously applied to learn the intrinisic graph structure in the original dimension.
#'The optimization of graph inference can be represented as:
#'\deqn{\mathop{min}_{f_g \in \mathcal{F}} \mathop{min}_{\{\mathbf{z}_1, ..., \mathbf{z}_M\}} \sum_{(V_i, V_j) \in
#'\mathcal{E}} b_{i,j}||f_g(\mathbf{z}_i) - f_g(\mathbf{z}_j)||^2} where
#'\eqn{f_g} is a function to map the instrinsic data space \eqn{\mathcal{Z} = \{\mathbf{z}_1, ..., \mathbf{z}_M\}} back to the input data space (reverse embedding)
#'\eqn{\mathcal{X} = \{ \mathbf{x}_1, ..., \mathbf{x}_N\}}. \eqn{V_i} is the the vertex of the instrinsic undirected graph
#'\eqn{\mathcal{G} = (\mathcal{V}, \mathcal{E})}. \eqn{b_{ij}} is
#'the edge weight associates with the edge set \eqn{\mathcal{E}}.
#'In order to learn the intrinsic structure from a reduced dimension, we need also to consider a term which includes
#'the error during the learning of the instrinsic structure. This strategy is incorporated as the following:
#'\deqn{\mathop{min}_{\mathcal{G} \in \hat{\mathcal{G}}_b}\mathop{min}_{f_g \in \mathcal{F}} \mathop{min}_{\{\mathbf{z}_1, ...,
#'\mathbf{z}_M\}} \sum_{i = 1}^N ||\mathbf{x}_i - f_g (\mathbf{z}_i)||^2 + \frac{\lambda}{2} \sum_{(V_i, V_j) \in \mathcal{E}}
#'b_{i,j}||f_g(\mathbf{z}_i) - f_g(\mathbf{z}_j)||^2} where \eqn{\lambda} is a non-negative parameter which controls the tradeoff between the data
#'reconstruction error and the reverse graph embedding.
#'
#'@section Dimensionality reduction via learning a tree:
#'The general framework for reducing dimension by learning an intrinsic structure in a low dimension requires a feasible set
#'\eqn{\hat{\mathcal{G}}_b} of graph and a mapping function \eqn{f_\mathcal{G}}. The algorithm uses minimum spanning tree as the feasible tree
#'graph structure, which can be solved by Kruskal' algoritm. A linear projection model \eqn{f_g (\mathbf{z}) = \mathbf{Wz}} is used as the mapping function.
#'Those setting results in the following specific form for the previous framework:
#'\deqn{\mathop{min}_{\mathbf{W}, \mathbf{Z}, \mathbf{B}} \sum_{i = 1}^N ||\mathbf{x}_i  - \mathbf{W}\mathbf{z}_i||^2 + \frac{\lambda}{2} \sum_{i,j}b_{i,j}||\mathbf{W} \mathbf{z}_i - \mathbf{W} \mathbf{z}_j||^2}
#'where \eqn{\mathbf{W} = [\mathbf{w}_1, ..., \mathbf{w}_d] \in
#'\mathcal{R}^{D \times d}} is an orthogonal set of \eqn{d} linear basis vectors. We can group tree graph \eqn{\mathbf{B}}, the orthogonal set of linear basis vectors
#'and projected points in reduced dimension \eqn{\mathbf{W}, \mathbf{Z}} as two groups and apply alternative structure optimization to optimize the tree graph.
#'This method is defined as DRtree (Dimension Reduction tree) as discussed by the authors.
#'
#'@section Discriminative dimensionality reduction via learning a tree:
#'In order to avoid the issues where data points scattered into different branches (which leads to lose of cluster information) and to
#' incorporate the discriminative information,another set of points \eqn{\{\mathbf{y}_k\}_{k = 1}^K} as the centers of \eqn{\{\mathbf{z}_i\}^N_{i = 1}} can be also introduced.
#' By so doing, the objective functions of K-means and the DRtree can be simulatenously minimized. The author further proposed a soft partition method
#' to account for the limits from K-means and proposed the following objective function:
#' \deqn{\mathop{min}_{\mathbf{W}, \mathbf{Z}, \mathbf{B}, \mathbf{Y}, \mathbf{R}} \sum_{i = 1}^N ||\mathbf{x}_i - \mathbf{W} \mathbf{z}_i||^2 +
#' \frac{\lambda}{2} \sum_{k, k'}b_{k, k'}||\mathbf{W} \mathbf{y}_k - \mathbf{W} \mathbf{y}_k'||^2 +
#' \gamma\Big[\sum_{k = 1}^K \sum_{i = 1}^N r_{i, k} ||\mathbf{z}_i - \mathbf{y}_k||^2 + \sigma \Omega (\mathbf{R})\Big]}
#' \deqn{s.t.\ \mathbf{W}^T \mathbf{W} = \mathbf{I}, \mathbf{B} \in \mathcal{B}, \sum_{k = 1}^K r_{i, k} = 1,
#' r_{i, k} \leq 0, \forall i, \forall k} where \eqn{\mathbf{R} \in \mathcal{R}^{N \times N}, \Omega(\mathbf{R}) = \sum_{i = 1}^N \sum_{k = 1}^k r_{i, k} log\ r_{i, k}} is the negative
#' entropy regularization which transforms the hard assignments used in K-means into soft assignments and \eqn{\sigma > 0} is the regulization parameter.
#'Alternative structure optimization is again used to solve the above problem by separately optimize each group \eqn{{\mathbf{W}, \mathbf{Z}, \mathbf{Y}}, {\mathbf{B}, \mathbf{R}}} until convergence.
#'
#'@section The actual algorithm of DDRTree:
#'\eqn{1.} \eqn{\mathbf{Input}}: Data matrix \eqn{\mathbf{X}}, parameters \eqn{\lambda, \sigma, \gamma} \cr
#'\eqn{2.} Initialize \eqn{\mathbf{Z}} by PCA \cr
#'\eqn{3.} \eqn{K = N, \mathbf{Y} = \mathbf{Z}} \cr
#'\eqn{4.} \eqn{\mathbf{repeat}}: \cr
#'      \eqn{\    5.}  \eqn{d_{k,k'} = ||\mathbf{y}_k - \mathbf{y}_{k'}||^2, \forall k, \forall k'} \cr
#'      \eqn{\     6.}   Obtain \eqn{\mathbf{B}} via Kruskal's algorithm \cr
#'      \eqn{\     7.}   \eqn{\mathbf{L} = diag(\mathbf{B1}) - \mathbf{B}} \cr
#'      \eqn{\     8.}   Compute \eqn{\mathbf{R}} with each element \cr
#'      \eqn{\     9.}   \eqn{\tau = diag(\mathbf{1}^T \mathbf{R})} \cr
#'      \eqn{\     10.}  \eqn{\mathbf{Q} = \frac{1}{\mathbf{1} + \gamma} \Big[\mathbf{I} + \mathbf{R} (\frac{1 + \gamma}{\gamma}(\frac{\lambda}{\gamma} \mathbf{L} +
#'     \tau) - \mathbf{R}^T \mathbf{R})^{-1} \mathbf{R}^T\Big]} \cr
#'      \eqn{\     11.}  \eqn{\mathbf{C} = \mathbf{X Q X}^T} \cr
#'      \eqn{\     12.}  Perform eigen-decomposition on \eqn{\mathbf{C}} such that \eqn{\mathbf{C} = \mathbf{U} \wedge \mathbf{U}^T} and \eqn{diag(\wedge)} is sorted in a descending order \cr
#'      \eqn{\     13.} \eqn{\mathbf{W} = \mathbf{U}(:, 1:d)} \cr
#'      \eqn{\     14.} \eqn{\mathbf{Z} = \mathbf{W}^T \mathbf{X Q}} \cr
#'      \eqn{\     15.} \eqn{\mathbf{Y} = \mathbf{Z R} (\frac{\lambda}{\gamma} \mathbf{L} + \tau)^{-1}} \cr
#'\eqn{16.} \eqn{\mathbf{Until}} Convergence \cr
#'
#'@section Implementation of DDRTree algorithm:
#'We implemented the algorithm mostly in Rcpp for the purpose of efficiency. It also has extensive optimization
#'for sparse input data. This implementation is originally based on the matlab code provided from the author of DDRTree paper.
#'
#'@docType package
#'@name DDRTree
#'@import irlba
#'@importFrom Rcpp evalCpp
#'@useDynLib DDRTree
#'@aliases DDRTree DDRTree-package
NULL