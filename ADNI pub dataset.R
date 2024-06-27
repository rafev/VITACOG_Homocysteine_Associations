ADNI_MRIs = read.csv("ADNI1_Screening_1.5T_4_18_2023.csv")
ADNI_MRIs = ADNI_MRIs[,2:11]
ADNI_MRIs$RID = gsub("(.*?_.*?_)", "", ADNI_MRIs$Subject)



ADNI_Meth = as.data.frame(readxl::read_xlsx("ADNI_DNA_Methylation_SampleAnnotation_20170530.xlsx"))
ADNI_Meth = ADNI_Meth[order(ADNI_Meth$RID, ADNI_Meth$Edate),]
ADNI_Meth$RID = stringr::str_pad(ADNI_Meth$RID, 4, side = "left", pad = 0)
ADNI_Meth$CollectYr = NA

patient_IDs = unique(ADNI_Meth$RID)

ADNI_Meth = lapply(patient_IDs, function(x) {
  dat = ADNI_Meth[ADNI_Meth$RID == x,]
  
  rep_num = nrow(dat)
  
  for (i in 1:rep_num) {
    j = i-1
    
    # conditions to assign year
    if (i == 1) {
      dat[i,9] = 0
    }
    if (i > 1) {
      if (dat[i,4] == dat[j,4]) {
        dat[i,9] = dat[j,9]
      } else {
        dat[i,9] = dat[j,9]+1
      }
    }
    
  }
  
  data.table(dat)
})
ADNI_Meth = rbindlist(ADNI_Meth)

rm(patient_IDs)


ADNI_Meth_final = ADNI_Meth[ADNI_Meth$RID %in% ADNI_MRIs$RID,]

Patient_count_BL = length(unique(ADNI_Meth_final$RID))

Patient_count_2yr = length(unique(ADNI_Meth_final$RID %in% ADNI_Meth_final[ADNI_Meth_final$CollectYr == 2,2]))
