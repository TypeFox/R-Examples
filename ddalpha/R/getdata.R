#all_datasets = c("baby","banknoten","biomed","bloodtransfusion","breast_cancer_wisconsin","bupa","chemdiab_1vs2","chemdiab_1vs3","chemdiab_2vs3","cloud","crabB_MvsF","crabF_BvsO","crabM_BvsO","crabO_MvsF","crab_BvsO","crab_MvsF","cricket_CvsP","diabetes","ecoli_cpvsim","ecoli_cpvspp","ecoli_imvspp","gemsen_MvsF","glass","groessen_MvsF","haberman","heart","hemophilia","indian_liver_patient_1vs2","indian_liver_patient_FvsM","iris_setosavsversicolor","iris_setosavsvirginica","iris_versicolorvsvirginica","irish_ed_MvsF","kidney","pima","plasma_retinol_MvsF","segmentation","socmob_IvsNI","socmob_WvsB","tae","tennis_MvsF","tips_DvsN","tips_MvsF","uscrime_SvsN","vertebral_column","veteran_lung_cancer","vowel_MvsF","wine_1vs2","wine_1vs3","wine_2vs3")

getdata <- function (name){
  data(list = name, envir = environment())
  return(get(name))
}
