#' @name thomas01a
#' @title 2x2 face recognition confusion matrix for Observer A
#' @description This data set contains the results of a full-report face recognition experiment reported in Thomas (2001). For Observer A, the two channels are degree of eye separation and nose length. 
#' @docType data
#' @usage data(thomas01a)
#' @format a \code{matrix} instance, containing counts for all stimulus-response combinations. Each row corresponds to a different stimulus presentation (in the order aa, ab, ba, bb) and each column in that row represents the frequency of each response (in the order aa, ab, ba, bb).
#' @source Thomas, R. D. (2001). Characterizing perceptual interactions in face identification using multidimensional signal detection theory. In M.Wenger & J.T. Townsend (Eds.) Computational, geometric, and process perspectives on facial cognition: Contexts and challenges. Hillsdale, NJ: Erlbaum.
#' @author Robin D. Thomas 
NULL

#' @name thomas01b
#' @title 2x2 face recognition confusion matrix for Observer B
#' @description This data set contains the results of a full-report face recognition experiment reported in Thomas (2001). For Observer B, the two channels are degree of eye separation and mouth width. 
#' @docType data
#' @usage data(thomas01b)
#' @format a \code{matrix} instance, containing counts for all stimulus-response combinations. Each row corresponds to a different stimulus presentation (in the order aa, ab, ba, bb) and each column in that row represents the frequency of each response (in the order aa, ab, ba, bb).
#' @source Thomas, R. D. (2001). Characterizing perceptual interactions in face identification using multidimensional signal detection theory. In M.Wenger & J.T. Townsend (Eds.) Computational, geometric, and process perspectives on facial cognition: Contexts and challenges. Hillsdale, NJ: Erlbaum.
#' @author Robin D. Thomas 
NULL

#' @name wo89xt
#' @title Cross-tabulated concurrent detection data
#' @description This data set contains a slightly coarse-grained version of Table 1 from Wickens and Olzak (1989). 
#' For each of four possible combinations of stimuli, participants gave a graded confidence judgement (collapsed here to 1-4) on both dimensions concurrently.
#' A rating of 1 corresponded to "definitely absent" and a rating of 4 corresponded to "definitely present".
#' @usage data(wo89xt)
#' @format an \code{xtabs} instance, containing counts for all stimulus-response combinations. 
#' For each of 4 Stim levels (where NN = absent+absent, LN = low-frequency signal+absent, NH = absent+high-frequency signal, LH = low-frequency signal+high-frequency signal), 
#' there is a 4x4 table giving the frequency of each rating.
#' @source Wickens, T. D., & Olzak, L. A. (1989). The statistical analysis of concurrent detection ratings. Perception & psychophysics, 45(6), 514-528.
#' @author Thomas D. Wickens and Lynn A. Olzak
NULL

#' @name silbert12
#' @title 2x2 phoneme confusion matrix 
#' @description Confusion matrix from speech perception experiment probing confusions between noise-masked tokens of English [p],[b],[f], and [v] (observer 3 in Ref.)
#' @usage data(silbert12)
#' @format A \code{matrix} instance, containing counts for all stimulus-response combinations. Rows correspond to stimuli, columns to responses
#' @source Silbert, N. H. (2012). Syllable structure and integration of voicing and manner of articulation information in labial consonant identification. Journal of the Acoustical Society of America, 131(5), 4076-4086.
#' @author Noah H. Silbert
NULL

#' @name silbert09a
#' @title 2x2 Frequency vs. Duration confusion matrix
#' @description Confusion matrix from auditory perception experiment, in which listeners identified noise stimuli varying across frequency range and duration (Experiment 1, Observer 3 in Ref.)
#' @usage data(silbert09a)
#' @format A \code{matrix} instance, containing counts for all stimulus-response combinations. Rows correspond to stimuli, columns to responses
#' @source Silbert, N. H., Townsend, J. T., & Lentz, J. J. (2009). Independence and separability in the perception of complex nonspeech sounds. Attention, Perception, & Psychophysics, 71(8), 1900-1915.
#' @author Noah H. Silbert
NULL

#' @name silbert09b
#' @title 2x2 Pitch vs. Timbre confusion matrix
#' @description Confusion matrix from auditory perception experiment, in which listeners identified 13-component harmonic stimuli varying across fundamental frequency and location of spectral prominence (Experiment 2, Observer 7 in Ref..
#' @usage data(silbert09b)
#' @format A \code{matrix} instance, containing counts for all stimulus-response combinations. Rows correspond to stimuli, columns to responses
#' @source Silbert, N. H., Townsend, J. T., & Lentz, J. J. (2009). Independence and separability in the perception of complex nonspeech sounds. Attention, Perception, & Psychophysics, 71(8), 1900-1915.
#' @author Noah H. Silbert
NULL

#' @name thomas15a
#' @title 3x3 face recognition confusion matrix for Observer A
#' @description This data set contains the results of a 3x3 full-report face recognition experiment reported in Thomas et al (2015). The two channels are degree of eye separation and nose width, with three levels on each dimension.
#' @docType data
#' @usage data(thomas15a)
#' @format an \code{xtabs} instance, containing counts for all stimulus-response combinations. 
#' The first two dimensions consist of the response counts on each level of nose width and eye separation, respectively, and the third dimension indexes the stimulus.
#' @source Thomas, R. D., Altieri, N. A., Silbert, N. H., Wenger, M. J., & Wessels, P. M. (2015). Multidimensional signal detection decision models of the uncertainty task: Application to face perception. Journal of Mathematical Psychology, 66, 16-33.
#' @author Robin D. Thomas 
NULL

#' @name thomas15b
#' @title 3x3 face recognition confusion matrix for Observer B
#' @description This data set contains the results of a 3x3 full-report face recognition experiment reported in Thomas et al (2015). The two channels are degree of eye separation and nose width, with three levels on each dimension.
#' @docType data
#' @usage data(thomas15b)
#' @format an \code{xtabs} instance, containing counts for all stimulus-response combinations. 
#' The first two dimensions consist of the response counts on each level of nose width and eye separation, respectively, and the third dimension indexes the stimulus.
#' @source Thomas, R. D., Altieri, N. A., Silbert, N. H., Wenger, M. J., & Wessels, P. M. (2015). Multidimensional signal detection decision models of the uncertainty task: Application to face perception. Journal of Mathematical Psychology, 66, 16-33.
#' @author Robin D. Thomas 
NULL
