# http://stackoverflow.com/questions/12391195/include-data-examples-in-developing-r-packages

#' Blood Transfusion Service Center Data Set 
#'
#' Data taken from the Blood Transfusion Service Center in Hsin-Chu City in Taiwan.
#' To demonstrate the RFMTC marketing model (a modified version of RFM), this study adopted the donor database
#' of Blood Transfusion Service Center in Hsin-Chu City in Taiwan.
#' The center passes their blood transfusion service bus to one university in Hsin-Chu City to gather blood
#' donated about every three months.
#' To build a FRMTC model, we selected 748 donors at random from the donor database.
#' These 748 donor data, each one included R (Recency - months since last donation),
#' F (Frequency - total number of donation), M (Monetary - total blood donated in c.c.),
#' T (Time - months since first donation), and a binary variable representing whether he/she donated blood in March 2007
#' (1 stand for donating blood; 0 stands for not donating blood).
#' 
#' The variables are as follows:
#'
#' \itemize{
#' \item R. Recency - months since last donation
#' \item F. Frequency - total number of donations
#' \item M. Monetary - total blood donated in c.c. (mL)
#' \item T. Time - months since first donation
#' \item y. a binary variable representing whether he/she donated blood in March 2007 (1=yes; 0 =no)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name UCI.transfusion
#' @usage data(UCI.transfusion)
#' @format A data frame with 748 rows and 5 variables
#' @source Dataset downloaded from the UCI Machine Learning Repository.
#' \url{http://archive.ics.uci.edu/ml/datasets/Blood+Transfusion+Service+Center}
#' 
#' Original Owner and Donor:
#' Prof. I-Cheng Yeh
#' Department of Information Management
#' Chung-Hua University
#' Hsin Chu, Taiwan 30067, R.O.C.
#' e-mail: icyeh 'at' chu.edu.tw
#' @references Yeh, I-Cheng, Yang, King-Jang, and Ting, Tao-Ming, "Knowledge discovery on RFM model using Bernoulli sequence",
#' Expert Systems with Applications, 2008. DOI: 10.1016/j.eswa.2008.07.018
#' @references Lichman, M. (2013). UCI Machine Learning Repository [\url{http://archive.ics.uci.edu/ml}].
#' Irvine, CA: University of California, School of Information and Computer Science.
NULL


#' Breast Cancer Wisconsin (Diagnostic) Data Set
#'
#' Features are computed from a digitized image of a fine needle aspirate (FNA) of a breast mass.
#' They describe characteristics of the cell nuclei present in the image.
#' 
#' Separating plane described above was obtained using Multisurface Method-Tree (MSM-T)
#' [K. P. Bennett, "Decision Tree Construction Via Linear Programming." Proceedings of the 4th Midwest Artificial
#' Intelligence and Cognitive Science Society, pp. 97-101, 1992], a classification method which uses linear programming to
#' construct a decision tree. Relevant features were selected using an exhaustive search in the space of 1-4 features
#' and 1-3 separating planes.
#' 
#' The actual linear program used to obtain the separating plane in the 3-dimensional space is that described in:
#' [K. P. Bennett and O. L. Mangasarian: "Robust Linear Programming Discrimination of Two Linearly Inseparable Sets",
#' Optimization Methods and Software 1, 1992, 23-34].
#' 
#' The variables are as follows:
#'
#' \itemize{
#' \item ID number
#' \item Diagnosis (1 = malignant, 0 = benign) 
#' \item Ten real-valued features are computed for each cell nucleus
#' }
#'
#' @docType data
#' @keywords datasets
#' @name UCI.BCD.Wisconsin
#' @usage data(UCI.BCD.Wisconsin)
#' @format A data frame with 569 rows and 32 variables
#' @source Dataset downloaded from the UCI Machine Learning Repository.
#' \url{http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)}
#' 
#' Creators:
#' 
#' 1. Dr. William H. Wolberg, General Surgery Dept.
#' University of Wisconsin, Clinical Sciences Center
#' Madison, WI 53792
#' wolberg 'at' eagle.surgery.wisc.edu
#' 
#' 2. W. Nick Street, Computer Sciences Dept.
#' University of Wisconsin, 1210 West Dayton St., Madison, WI 53706
#' street 'at' cs.wisc.edu 608-262-6619
#' 
#' 3. Olvi L. Mangasarian, Computer Sciences Dept.
#' University of Wisconsin, 1210 West Dayton St., Madison, WI 53706
#' olvi 'at' cs.wisc.edu
#' 
#' Donor: Nick Street
#' @references W.N. Street, W.H. Wolberg and O.L. Mangasarian. Nuclear feature extraction for breast tumor diagnosis.
#' IS&T/SPIE 1993 International Symposium on Electronic Imaging: Science and Technology, volume 1905, pages 861-870, San Jose, CA, 1993.
#' @references Lichman, M. (2013). UCI Machine Learning Repository [\url{http://archive.ics.uci.edu/ml}].
#' Irvine, CA: University of California, School of Information and Computer Science. 
NULL


#' ISOLET Data Set (ABC)
#'
#' This data set was generated as follows. 150 subjects spoke the name of each letter of the alphabet twice.
#' Hence, we have 52 training examples from each speaker.
#' 
#' To reduce package size, only the 3 first letters are included here. The full dataset can be obtained
#' from \url{http://archive.ics.uci.edu/ml/datasets/ISOLET}.
#' 
#' The features are described in the paper by Cole and Fanty cited below.
#' The features include spectral coefficients; contour features, sonorant features, pre-sonorant features,
#' and post-sonorant features. Exact order of appearance of the features is not known.
#'
#' @docType data
#' @keywords datasets
#' @name UCI.ISOLET.ABC
#' @usage data(UCI.ISOLET.ABC)
#' @format A data frame with 900 rows and 618 variables
#' @source Dataset downloaded from the UCI Machine Learning Repository.
#' \url{http://archive.ics.uci.edu/ml/datasets/ISOLET}
#' 
#' Creators:
#' 
#' Ron Cole and Mark Fanty
#' Department of Computer Science and Engineering,
#' Oregon Graduate Institute, Beaverton, OR 97006.
#' cole 'at' cse.ogi.edu, fanty 'at' cse.ogi.edu
#' 
#' Donor:
#' 
#' Tom Dietterich
#' Department of Computer Science
#' Oregon State University, Corvallis, OR 97331
#' tgd 'at' cs.orst.edu
#' @references Fanty, M., Cole, R. (1991). Spoken letter recognition. In Lippman, R. P., Moody, J.,
#' and Touretzky, D. S. (Eds). Advances in Neural Information Processing Systems 3. San Mateo, CA: Morgan Kaufmann.
#' [\url{http://rexa.info/paper/bee6de062d0d168c5c5b5290b11cd6b12ca8472e}]
#' @examples
#' # NB: 50 iterations isn't enough in this case,
#' # it was chosen so that the example runs fast enough on CRAN check farm
#' data(UCI.ISOLET.ABC);
#' X=as.matrix(sN.normalizeDF(as.data.frame(UCI.ISOLET.ABC[,1:617])));
#' y=as.matrix(UCI.ISOLET.ABC[,618]-1);
#' myMLP=sN.MLPtrain(X=X,y=y,hidden_layer_size=20,it=50,lambda=0.5,alpha=0.5);
#' myPrediction=sN.MLPpredict(nnModel=myMLP,X=X,raw=FALSE);
#' table(y,myPrediction);
NULL
