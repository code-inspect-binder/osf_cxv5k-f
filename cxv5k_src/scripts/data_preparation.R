
####PREPARATION VARIABLES METALEXA#####
library(questionr)
library(tidyverse)
library(countrycode)

DATA<- read.csv2("data/raw_data_Kimmoun_et_al.csv", dec = c(".", ","))

#unique ID
DATA$IDU<-paste(DATA$PMID,sep="_",DATA$author_name)

#Year of recruitement
DATA$TIME<- ifelse(!is.na(DATA$start) & !is.na(DATA$end),  
                          round((DATA$start+DATA$end)/2), 
                   ifelse(!is.na(DATA$start) & is.na(DATA$end), 
                            DATA$start, 
                    ifelse(is.na(DATA$start) & !is.na(DATA$end), 
                            DATA$end, NA)))

DATA$TIME[is.na(DATA$TIME)]<-DATA$published[is.na(DATA$TIME)]-5

#####################################
##### Patients characteristics  #####
#####################################

## Ejection fraction ##

DATA$EF <-  ifelse(is.na(DATA$"all_mean_EF")==T, 
                   DATA$all_median_EF, DATA$"all_mean_EF")            
# taking median value if no mean value

# weighted mean EF
Mp<-subset(DATA, select=c(mean_EF_1, mean_EF_2,mean_EF_3,mean_EF_4,mean_EF_5,mean_EF_6, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt,number_5_ttt,number_6_ttt))

meanp<-function(monvec){
  return(wtd.mean(x=monvec[1:6], w=monvec[7:12], na.rm=T))
}

DATA$mean_p<-apply(Mp,FUN=meanp,MARGIN=1)

#Weighted mean from median considered as surrogate of mean
Mep<-subset(DATA, select=c(median_EF_1, median_EF_2,median_EF_3,median_EF_4, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt))

meanmep<-function(monvec){
  return(wtd.mean(x=monvec[1:4], w=monvec[5:8], na.rm=T))
}
DATA$meanme_p<-apply(data.frame(Mep),FUN=meanp,MARGIN=1)

#EF_all
DATA$EF_all<-DATA$EF
DATA$EF_all[is.na(DATA$EF_all)==T] <- DATA$mean_p[is.na(DATA$EF_all)==T]
DATA$EF_all[is.na(DATA$EF_all)==T] <- DATA$meanme_p[is.na(DATA$EF_all)==T]


## Heart failure ##

DATA$HF_1_n<- round(DATA$HF_1*DATA$number_1_ttt)/100
DATA$HF_2_n<- round(DATA$HF_2*DATA$number_2_ttt)/100
DATA$HF_3_n<- round(DATA$HF_3*DATA$number_3_ttt)/100
DATA$HF_4_n<- round(DATA$HF_4*DATA$number_4_ttt)/100
DATA$HF_5_n<- round(DATA$HF_5*DATA$number_5_ttt)/100
DATA$HF_6_n<- round(DATA$HF_6*DATA$number_6_ttt)/100
DATA$HF_all_n<- round(DATA$HF_all*DATA$number_ttt)/100

DATA$HF_all_n_2 <- apply(DATA[,c("HF_1_n","HF_2_n","HF_3_n","HF_4_n","HF_5_n","HF_6_n")],1,sum,na.rm=T)

DATA$HF_all_n_2[DATA$HF_all_n_2==0] <- NA 

DATA$HF_all_n[is.na(DATA$HF_all_n)==T] <- DATA$HF_all_n_2[is.na(DATA$HF_all_n)==T]

DATA$HF_all_rate <- DATA$HF_all_n/DATA$number_ttt


## Hypertension ##

DATA$HT_1_n<- round(DATA$HT_1*DATA$number_1_ttt)/100
DATA$HT_2_n<- round(DATA$HT_2*DATA$number_2_ttt)/100
DATA$HT_3_n<- round(DATA$HT_3*DATA$number_2_ttt)/100
DATA$HT_4_n<- round(DATA$HT_4*DATA$number_2_ttt)/100
DATA$HT_5_n<- round(DATA$HT_5*DATA$number_2_ttt)/100
DATA$HT_6_n<- round(DATA$HT_6*DATA$number_2_ttt)/100
DATA$all_HT_n<- round(DATA$all_HT*DATA$number_ttt)/100

DATA$all_HT_n_2 <- apply(DATA[,c("HT_1_n","HT_2_n","HT_3_n","HT_4_n","HT_5_n","HT_6_n")],1,sum,na.rm=T)

DATA$all_HT_n_2[DATA$all_HT_n_2==0] <- NA 

DATA$all_HT_n[is.na(DATA$all_HT_n)==T] <- DATA$all_HT_n_2[is.na(DATA$all_HT_n)==T]

DATA$all_HT_rate <- DATA$all_HT_n/DATA$number_ttt


## Diabetes mellitus ##

DATA$DM_1_n<- round(DATA$DM_1*DATA$number_1_ttt)/100
DATA$DM_2_n<- round(DATA$DM_2*DATA$number_2_ttt)/100
DATA$DM_3_n<- round(DATA$DM_3*DATA$number_3_ttt)/100
DATA$DM_4_n<- round(DATA$DM_4*DATA$number_4_ttt)/100
DATA$DM_5_n<- round(DATA$DM_5*DATA$number_5_ttt)/100
DATA$DM_6_n<- round(DATA$DM_6*DATA$number_6_ttt)/100
DATA$all_DM_n<- round(DATA$all_DM*DATA$number_ttt)/100

DATA$all_DM_n_2 <- apply(DATA[,c("DM_1_n","DM_2_n","DM_3_n","DM_4_n","DM_5_n","DM_6_n")],1,sum,na.rm=T)

DATA$all_DM_n_2[DATA$all_DM_n_2==0] <- NA 

DATA$all_DM_n[is.na(DATA$all_DM_n)==T] <- DATA$all_DM_n_2[is.na(DATA$all_DM_n)==T]

DATA$all_DM_rate <- DATA$all_DM_n/DATA$number_ttt


## Coronary heart disease ##

DATA$CAD.IHD_1_n<- round(DATA$CAD.IHD_1*DATA$number_1_ttt)/100
DATA$CAD.IHD_2_n<- round(DATA$CAD.IHD_2*DATA$number_2_ttt)/100
DATA$CAD.IHD_3_n<- round(DATA$CAD.IHD_3*DATA$number_3_ttt)/100
DATA$CAD.IHD_4_n<- round(DATA$CAD.IHD_4*DATA$number_4_ttt)/100
DATA$CAD.IHD_5_n<- round(DATA$CAD.IHD_5*DATA$number_5_ttt)/100
DATA$CAD.IHD_6_n<- round(DATA$CAD.IHD_6*DATA$number_6_ttt)/100
DATA$all_CAD.IHD_n<- round(DATA$all_CAD.IHD*DATA$number_ttt)/100

DATA$all_CAD.IHD_n_2 <- apply(DATA[,c("CAD.IHD_1_n","CAD.IHD_2_n","CAD.IHD_3_n","CAD.IHD_4_n","CAD.IHD_5_n","CAD.IHD_6_n")],1,sum,na.rm=T)

DATA$all_CAD.IHD_n_2[DATA$all_CAD.IHD_n_2==0] <- NA 

DATA$all_CAD.IHD_n[is.na(DATA$all_CAD.IHD_n)==T] <- DATA$all_CAD.IHD_n_2[is.na(DATA$all_CAD.IHD_n)==T]

DATA$all_CAD.IHD_rate <- DATA$all_CAD.IHD_n/DATA$number_ttt


## Chronic obstructive respiratory disease ##

DATA$COPD_1_n<- round(DATA$COPD_1*DATA$number_1_ttt)/100
DATA$COPD_2_n<- round(DATA$COPD_2*DATA$number_2_ttt)/100
DATA$COPD_3_n<- round(DATA$COPD_3*DATA$number_3_ttt)/100
DATA$COPD_4_n<- round(DATA$COPD_4*DATA$number_4_ttt)/100
DATA$COPD_5_n<- round(DATA$COPD_5*DATA$number_5_ttt)/100
DATA$COPD_6_n<- round(DATA$COPD_6*DATA$number_6_ttt)/100
DATA$all_COPD_n<- round(DATA$all_COPD*DATA$number_ttt)/100

DATA$all_COPD_n_2 <- apply(DATA[,c("COPD_1_n","COPD_2_n","COPD_3_n","COPD_4_n","COPD_5_n","COPD_6_n")],1,sum,na.rm=T)

DATA$all_COPD_n_2[DATA$all_COPD_n_2==0] <- NA 

DATA$all_COPD_n[is.na(DATA$all_COPD_n)==T] <- DATA$all_COPD_n_2[is.na(DATA$all_COPD_n)==T]

DATA$all_COPD_rate <- DATA$all_COPD_n/DATA$number_ttt


## Atrial fibrillation ##

DATA$Afib_1_n<- round(DATA$Afib_1*DATA$number_1_ttt)/100
DATA$Afib_2_n<- round(DATA$Afib_2*DATA$number_2_ttt)/100
DATA$Afib_3_n<- round(DATA$Afib_3*DATA$number_3_ttt)/100
DATA$Afib_4_n<- round(DATA$Afib_4*DATA$number_4_ttt)/100
DATA$Afib_5_n<- round(DATA$Afib_5*DATA$number_5_ttt)/100
DATA$Afib_6_n<- round(DATA$Afib_6*DATA$number_6_ttt)/100
DATA$all_Afib_n<- round(DATA$all_Afib*DATA$number_ttt)/100

DATA$all_Afib_n_2 <- apply(DATA[,c("Afib_1_n","Afib_2_n","Afib_3_n","Afib_4_n","Afib_5_n","Afib_6_n")],1,sum,na.rm=T)

DATA$all_Afib_n_2[DATA$all_Afib_n_2==0] <- NA 

DATA$all_Afib_n[is.na(DATA$all_Afib_n)==T] <- DATA$all_Afib_n_2[is.na(DATA$all_Afib_n)==T]

DATA$all_Afib_rate <- DATA$all_Afib_n/DATA$number_ttt


## ACEi and/or ARB ##

DATA$all_ACE_ARB_adm <- apply(cbind(DATA$"A_all_ACE",DATA$"A_all_ARB",DATA$"A_all_ACE_ARB"),1,max,na.rm=T)
DATA$all_ACE_ARB_adm[DATA$all_ACE_ARB_adm=="-Inf"] <- NA

DATA$ACE_ARB_adm_1 <- apply(cbind(DATA$"A_ACE_1", DATA$"A_ARB_1", DATA$"A_ACE_ARB_1"), 1 , max,na.rm=T)
DATA$ACE_ARB_adm_1[DATA$ACE_ARB_adm_1=="-Inf"] <- NA
DATA$ACE_ARB_adm_2 <- apply(cbind(DATA$"A_ACE_2", DATA$"A_ARB_2", DATA$"A_ACE_ARB_2"), 1 , max,na.rm=T)
DATA$ACE_ARB_adm_2[DATA$ACE_ARB_adm_2=="-Inf"] <- NA
DATA$ACE_ARB_adm_3 <- apply(cbind(DATA$"A_ACE_3", DATA$"A_ARB_3", DATA$"A_ACE_ARB_3"), 1 , max,na.rm=T)
DATA$ACE_ARB_adm_3[DATA$ACE_ARB_adm_3=="-Inf"] <- NA
DATA$ACE_ARB_adm_4 <- apply(cbind(DATA$"A_ACE_4", DATA$"A_ARB_4", DATA$"A_ACE_ARB_4"), 1 , max,na.rm=T)
DATA$ACE_ARB_adm_4[DATA$ACE_ARB_adm_4=="-Inf"] <- NA
DATA$ACE_ARB_adm_5 <- apply(cbind(DATA$"A_ACE_5", DATA$"A_ARB_5", DATA$"A_ACE_ARB_5"), 1 , max,na.rm=T)
DATA$ACE_ARB_adm_5[DATA$ACE_ARB_adm_5=="-Inf"] <- NA
DATA$ACE_ARB_adm_6 <- DATA$"A_ACE_ARB_6"

DATA$ACE_ARB_adm_1_n<- round(DATA$ACE_ARB_adm_1*DATA$number_1_ttt)/100

DATA$ACE_ARB_adm_2_n<- round(DATA$ACE_ARB_adm_2*DATA$number_2_ttt)/100

DATA$ACE_ARB_adm_3_n<- round(DATA$ACE_ARB_adm_3*DATA$number_3_ttt)/100

DATA$ACE_ARB_adm_4_n<- round(DATA$ACE_ARB_adm_4*DATA$number_4_ttt)/100

DATA$ACE_ARB_adm_5_n<- round(DATA$ACE_ARB_adm_5*DATA$number_5_ttt)/100

DATA$ACE_ARB_adm_6_n<- round(DATA$ACE_ARB_adm_5*DATA$number_6_ttt)/100

DATA$all_ACE_ARB_adm_n<- round(DATA$all_ACE_ARB_adm*DATA$number_ttt)/100

DATA$all_ACE_ARB_adm_n_2 <- apply(DATA[,c("ACE_ARB_adm_1_n",
                                          "ACE_ARB_adm_2_n",
                                          "ACE_ARB_adm_3_n",
                                          "ACE_ARB_adm_4_n",
                                          "ACE_ARB_adm_5_n",
                                          "ACE_ARB_adm_6_n")],1,sum,na.rm=T)

DATA$all_ACE_ARB_adm_n_2[DATA$all_ACE_ARB_adm_n_2==0] <- NA 

DATA$all_ACE_ARB_adm_n[is.na(DATA$all_ACE_ARB_adm_n)==T] <- DATA$all_ACE_ARB_adm_n_2[is.na(DATA$all_ACE_ARB_adm_n)==T]

DATA$all_ACE_ARB_adm_rate <- DATA$all_ACE_ARB_adm_n/DATA$number_ttt


## Loop diuretics ##

DATA$A_diuretic_1_n<- round(DATA$A_diuretic_1*DATA$number_1_ttt)/100
DATA$A_diuretic_2_n<- round(DATA$A_diuretic_2*DATA$number_2_ttt)/100
DATA$A_diuretic_3_n<- round(DATA$A_diuretic_3*DATA$number_3_ttt)/100
DATA$A_diuretic_4_n<- round(DATA$A_diuretic_4*DATA$number_4_ttt)/100
DATA$A_diuretic_5_n<- round(DATA$A_diuretic_5*DATA$number_5_ttt)/100
DATA$A_all_diuretic_n<- round(DATA$A_all_diuretic*DATA$number_ttt)/100

DATA$A_all_diuretic_n_2 <- apply(DATA[,c("A_diuretic_1_n","A_diuretic_2_n","A_diuretic_3_n","A_diuretic_4_n","A_diuretic_5_n")],1,sum,na.rm=T)

DATA$A_all_diuretic_n_2[DATA$A_all_diuretic_n_2==0] <- NA 

DATA$A_all_diuretic_n[is.na(DATA$A_all_diuretic_n)==T] <- DATA$A_all_diuretic_n_2[is.na(DATA$A_all_diuretic_n)==T]

DATA$A_all_diuretic_rate <- DATA$A_all_diuretic_n/DATA$number_ttt


## Digoxin ##

DATA$A_digoxin_1_n<- round(DATA$A_digoxin_1*DATA$number_1_ttt)/100
DATA$A_digoxin_2_n<- round(DATA$A_digoxin_2*DATA$number_2_ttt)/100
DATA$A_digoxin_3_n<- round(DATA$A_digoxin_3*DATA$number_3_ttt)/100
DATA$A_digoxin_4_n<- round(DATA$A_digoxin_4*DATA$number_4_ttt)/100
DATA$A_all_digoxin_n<- round(DATA$A_all_digoxin*DATA$number_ttt)/100

DATA$A_all_digoxin_n_2 <- apply(DATA[,c("A_digoxin_1_n","A_digoxin_2_n","A_digoxin_3_n","A_digoxin_4_n")],1,sum,na.rm=T)

DATA$A_all_digoxin_n_2[DATA$A_all_digoxin_n_2==0] <- NA 

DATA$A_all_digoxin_n[is.na(DATA$A_all_digoxin_n)==T] <- DATA$A_all_digoxin_n_2[is.na(DATA$A_all_digoxin_n)==T]

DATA$A_all_digoxin_rate <- DATA$A_all_digoxin_n/DATA$number_ttt


## Beta-blockers ##

DATA$A_B_1_n<- round(DATA$A_B_1*DATA$number_1_ttt)/100
DATA$A_B_2_n<- round(DATA$A_B_2*DATA$number_2_ttt)/100
DATA$A_B_3_n<- round(DATA$A_B_3*DATA$number_3_ttt)/100
DATA$A_B_4_n<- round(DATA$A_B_4*DATA$number_4_ttt)/100
DATA$A_B_5_n<- round(DATA$A_B_5*DATA$number_5_ttt)/100
DATA$A_B_6_n<- round(DATA$A_B_6*DATA$number_6_ttt)/100
DATA$A_all_B_n<- round(DATA$A_all_B*DATA$number_ttt)/100

DATA$A_all_B_n_2 <- apply(DATA[,c("A_B_1_n","A_B_2_n","A_B_3_n","A_B_4_n","A_B_5_n","A_B_6_n")],1,sum,na.rm=T)

DATA$A_all_B_n_2[DATA$A_all_B_n_2==0] <- NA 

DATA$A_all_B_n[is.na(DATA$A_all_B_n)==T] <- DATA$A_all_B_n_2[is.na(DATA$A_all_B_n)==T]

DATA$A_all_B_rate <- DATA$A_all_B_n/DATA$number_ttt


###############
#MAIN OUTCOMES#
###############

####Data management : Major variables####

#### 30-day death

DATA$death_30_1[is.na(DATA$death_30_1)==T] <- round((DATA$death_30_1_.[is.na(DATA$death_30_1)==T]*DATA$number_1_follow_up[is.na(DATA$death_30_1)==T])/100)

DATA$death_30_2[is.na(DATA$death_30_2)==T] <- round((DATA$death_30_2_.[is.na(DATA$death_30_2)==T]*DATA$number_2_follow_up[is.na(DATA$death_30_2)==T])/100)

DATA$death_30_3[is.na(DATA$death_30_3)==T] <- round((DATA$death_30_3_.[is.na(DATA$death_30_3)==T]*DATA$number_3_follow_up[is.na(DATA$death_30_3)==T])/100)

DATA$death_30_4[is.na(DATA$death_30_4)==T] <- round((DATA$death_30_4_.[is.na(DATA$death_30_4)==T]*DATA$number_4_follow_up[is.na(DATA$death_30_4)==T])/100)

DATA$death_30_5[is.na(DATA$death_30_5)==T] <- round((DATA$death_30_5_.[is.na(DATA$death_30_5)==T]*DATA$number_5_follow_up[is.na(DATA$death_30_5)==T])/100)

DATA$death_30_6[is.na(DATA$death_30_6)==T] <- round((DATA$death_30_6_.[is.na(DATA$death_30_6)==T]*DATA$number_6_follow_up[is.na(DATA$death_30_6)==T])/100)

DATA$death_30_all[is.na(DATA$death_30_all)==T] <- round((DATA$death_30_all_.[is.na(DATA$death_30_all)==T]*DATA$number_follow_up[is.na(DATA$death_30_all)==T])/100)

DATA$death_30_all_2 <- apply(DATA[,c("death_30_1","death_30_2", "death_30_3" ,"death_30_4","death_30_5","death_30_6")],1,sum,na.rm=T)

DATA$death_30_all[is.na(DATA$death_30_all)==T] <- DATA$death_30_all_2[is.na(DATA$death_30_all)==T]

DATA$death_30_all[DATA$death_30_all==0] <- NA 

DATA$death_30_rate <- DATA$death_30_all/DATA$number_follow_up


### One-year death

DATA$death_365_1[is.na(DATA$death_365_1)==T] <- round((DATA$death_365_1_.[is.na(DATA$death_365_1)==T]*DATA$number_1_follow_up[is.na(DATA$death_365_1)==T])/100)

DATA$death_365_2[is.na(DATA$death_365_2)==T] <- round((DATA$death_365_2_.[is.na(DATA$death_365_2)==T]*DATA$number_2_follow_up[is.na(DATA$death_365_2)==T])/100)

DATA$death_365_3[is.na(DATA$death_365_3)==T] <- round((DATA$death_365_3_.[is.na(DATA$death_365_3)==T]*DATA$number_3_follow_up[is.na(DATA$death_365_3)==T])/100)

DATA$death_365_4[is.na(DATA$death_365_4)==T] <- round((DATA$death_365_4_.[is.na(DATA$death_365_4)==T]*DATA$number_4_follow_up[is.na(DATA$death_365_4)==T])/100)

DATA$death_365[is.na(DATA$death_365)==T] <- round((DATA$death_365_.[is.na(DATA$death_365)==T]*DATA$number_follow_up[is.na(DATA$death_365)==T])/100)

DATA$death_365_all_2 <- apply(DATA[,c("death_365_1","death_365_2", "death_365_3", "death_365_4")],1,sum,na.rm=T)

DATA$death_365_all_2[DATA$death_365_all_2==0] <- NA 

DATA$death_365[is.na(DATA$death_365)==T] <- DATA$death_365_all_2[is.na(DATA$death_365)==T]

DATA$death_365_rate <- DATA$death_365/DATA$number_follow_up


### 30-day readmission

DATA$rehosp_30_all.cause_1[is.na(DATA$rehosp_30_all.cause_1)==T] <- round((DATA$rehosp_30_all.cause_1_.[is.na(DATA$rehosp_30_all.cause_1)==T]*DATA$number_1_follow_up[is.na(DATA$rehosp_30_all.cause_1)==T])/100)

DATA$rehosp_30_all.cause_2[is.na(DATA$rehosp_30_all.cause_2)==T] <- round((DATA$rehosp_30_all.cause_2_.[is.na(DATA$rehosp_30_all.cause_2)==T]*DATA$number_2_follow_up[is.na(DATA$rehosp_30_all.cause_2)==T])/100)

DATA$rehosp_30_all.cause_3[is.na(DATA$rehosp_30_all.cause_3)==T] <- round((DATA$rehosp_30_all.cause_3_.[is.na(DATA$rehosp_30_all.cause_3)==T]*DATA$number_3_follow_up[is.na(DATA$rehosp_30_all.cause_3)==T])/100)

DATA$rehosp_30_all.cause_4[is.na(DATA$rehosp_30_all.cause_4)==T] <- round((DATA$rehosp_30_all.cause_4_.[is.na(DATA$rehosp_30_all.cause_4)==T]*DATA$number_4_follow_up[is.na(DATA$rehosp_30_all.cause_4)==T])/100)

DATA$rehosp_30_all.cause_5[is.na(DATA$rehosp_30_all.cause_5)==T] <- round((DATA$rehosp_30_all.cause_5_.[is.na(DATA$rehosp_30_all.cause_5)==T]*DATA$number_5_follow_up[is.na(DATA$rehosp_30_all.cause_5)==T])/100)

DATA$rehosp_30_all.cause_all[is.na(DATA$rehosp_30_all.cause_all)==T] <- round((DATA$rehosp_30_all.cause_all_.[is.na(DATA$rehosp_30_all.cause_all)==T]*DATA$number_follow_up[is.na(DATA$rehosp_30_all.cause_all)==T])/100)

DATA$rehosp_30_all.cause_all_2 <- apply(DATA[,c("rehosp_30_all.cause_1","rehosp_30_all.cause_2","rehosp_30_all.cause_3","rehosp_30_all.cause_4", "rehosp_30_all.cause_5")],1,sum,na.rm=T)

DATA$rehosp_30_all.cause_all_2[DATA$rehosp_30_all.cause_all_2==0] <- NA 

DATA$rehosp_30_all.cause_all[is.na(DATA$rehosp_30_all.cause_all)==T] <- DATA$rehosp_30_all.cause_all_2[is.na(DATA$rehosp_30_all.cause_all)==T]

DATA$rehosp_30_rate <- DATA$rehosp_30_all.cause_all/DATA$number_follow_up


#### One-year readmission

DATA$rehosp_365_1[is.na(DATA$rehosp_365_1)==T] <- round((DATA$rehosp_365_1_.[is.na(DATA$rehosp_365_1)==T]*DATA$number_1_follow_up[is.na(DATA$rehosp_365_1)==T])/100)

DATA$rehosp_365_2[is.na(DATA$rehosp_365_2)==T] <- round((DATA$rehosp_365_2_.[is.na(DATA$rehosp_365_2)==T]*DATA$number_2_follow_up[is.na(DATA$rehosp_365_2)==T])/100)

DATA$rehosp_365_3[is.na(DATA$rehosp_365_3)==T] <- round((DATA$rehosp_365_3_.[is.na(DATA$rehosp_365_3)==T]*DATA$number_3_follow_up[is.na(DATA$rehosp_365_3)==T])/100)

DATA$rehosp_365_4[is.na(DATA$rehosp_365_4)==T] <- round((DATA$rehosp_365_4_.[is.na(DATA$rehosp_365_4)==T]*DATA$number_4_follow_up[is.na(DATA$rehosp_365_4)==T])/100)

DATA$rehosp_365[is.na(DATA$rehosp_365)==T] <- round((DATA$rehosp_365_.[is.na(DATA$rehosp_365)==T]*DATA$number_follow_up[is.na(DATA$rehosp_365)==T])/100)

DATA$rehosp_365_all_2 <- apply(DATA[,c("rehosp_365_1","rehosp_365_2","rehosp_365_3","rehosp_365_4")],1,sum,na.rm=T)

DATA$rehosp_365_all_2[DATA$rehosp_365_all_2==0] <- NA 

DATA$rehosp_365[is.na(DATA$rehosp_365)==T] <- DATA$rehosp_365_all_2[is.na(DATA$rehosp_365)==T]

DATA$rehosp_365_rate <- DATA$rehosp_365/DATA$number_follow_up


######################
#sensitivity analysis#
######################

#for continuous variables, in order to minimize exclusion of studies, medians were considered as mean.
#When median for overall groups were not available, weighted mean were calculated from subgroups with medians considered as mean. 

### Systolic blood pressure ###

#median considered as mean when mean not available
DATA$SBP <-  ifelse(is.na(DATA$"all_mean_SBP")==T, DATA$all_median_SBP, DATA$"all_mean_SBP")     

#weighted mean from mean when overall mean was not available
Mp1<-subset(DATA, select=c(mean_SBP_1, mean_SBP_2,mean_SBP_3,mean_SBP_4,mean_SBP_5, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt, number_5_ttt))

meanp<-function(monvec){
  return(wtd.mean(x=monvec[1:5], w=monvec[6:10], na.rm=T))
}
DATA$mean_p1<-apply(Mp1,FUN=meanp,MARGIN=1)


#weighted mean from mean when overall median considered as mean was not available
Mep1<-subset(DATA, select=c(median_SBP_1, median_SBP_2,median_SBP_3,median_SBP_4,median_SBP_5,median_SBP_6, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt,number_5_ttt,number_6_ttt))

meanmep<-function(monvec){
  return(wtd.mean(x=monvec[1:6], w=monvec[7:12], na.rm=T))
}
DATA$meanme_p1<-apply(Mep1,FUN=meanp,MARGIN=1)

#Creation of a unique variable

DATA$SBP_all<-DATA$SBP
DATA$SBP_all[is.na(DATA$SBP_all)==T] <- DATA$mean_p1[is.na(DATA$SBP_all)==T] 
DATA$SBP_all[is.na(DATA$SBP_all)==T] <- DATA$meanme_p1[is.na(DATA$SBP_all)==T] 


### Heart rate ###

#median considered as mean when mean not "available"
DATA$HR <-  ifelse(is.na(DATA$"all_mean_HR")==T, DATA$all_median_HR, DATA$"all_mean_HR")  

#weighted mean from mean when overall mean was not available
Mp2<-subset(DATA, select=c(mean_HR_1, mean_HR_2,mean_HR_3,mean_HR_4,mean_HR_5,mean_HR_6, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt, number_5_ttt,number_6_ttt))

meanp<-function(monvec){
  return(wtd.mean(x=monvec[1:6], w=monvec[7:12], na.rm=T))
}
DATA$mean_p2<-apply(Mp2,FUN=meanp,MARGIN=1)

#weighted mean from mean when overall median considered as mean was not available
Mep2<-subset(DATA, select=c(median_HR_1, median_HR_2,median_HR_3,median_HR_4,median_HR_5,median_HR_6, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt,number_5_ttt,number_6_ttt))
meanmep<-function(monvec){
  return(wtd.mean(x=monvec[1:6], w=monvec[7:12], na.rm=T))
}
DATA$meanme_p2<-apply(Mep2,FUN=meanmep,MARGIN=1)

#Creation of a unique variable
DATA$HR_all<-DATA$HR
DATA$HR_all[is.na(DATA$HR_all)==T] <- DATA$mean_p2[is.na(DATA$HR_all)==T]
DATA$HR_all[is.na(DATA$HR_all)==T] <- DATA$meanme_p2[is.na(DATA$HR_all)==T]


### AGE ###

#median considered as mean when mean not available
DATA$AGE <-  ifelse(is.na(DATA$all_mean_age_years)==T, DATA$all_median_age, DATA$all_mean_age_years)  

#weighted mean from mean when overall mean was not available
Mp3<-subset(DATA, select=c(mean_age_1, mean_age_2,mean_age_3,mean_age_4,mean_age_5, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt,number_5_ttt))

meanp<-function(monvec){
  return(wtd.mean(x=monvec[1:5], w=monvec[6:10], na.rm=T))
}
DATA$mean_p3 <- apply(Mp3,FUN=meanp,MARGIN=1)

#weighted mean from mean when overall median considered as mean was not available
Mep3<-subset(DATA, select=c(median_age_1, median_age_2,median_age_3,median_age_4, number_1_ttt,number_2_ttt,number_3_ttt,number_4_ttt))

meanmep<-function(monvec){
  return(wtd.mean(x=monvec[1:4], w=monvec[5:8], na.rm=T))
}
DATA$meanme_p3<-apply(Mep3,FUN=meanmep,MARGIN=1)

#Creation of a unique variable
DATA$AGE_all<-DATA$AGE
DATA$AGE_all[is.na(DATA$AGE_all)==T] <- DATA$mean_p3[is.na(DATA$AGE_all)==T]
DATA$AGE_all[is.na(DATA$AGE_all)==T] <- DATA$meanme_p3[is.na(DATA$AGE_all)==T]

### SEX ###

DATA$gender_male_1[is.na(DATA$gender_male_1)==T] <- round((DATA$gender_male_1_.[is.na(DATA$gender_male_1)==T]*DATA$number_1_ttt[is.na(DATA$gender_male_1)==T])/100)

DATA$gender_male_2[is.na(DATA$gender_male_2)==T] <- round((DATA$gender_male_2_.[is.na(DATA$gender_male_2)==T]*DATA$number_2_ttt[is.na(DATA$gender_male_2)==T])/100)

DATA$gender_male_3[is.na(DATA$gender_male_3)==T] <- round((DATA$gender_male_3_.[is.na(DATA$gender_male_3)==T]*DATA$number_3_ttt[is.na(DATA$gender_male_3)==T])/100)

DATA$gender_male_4[is.na(DATA$gender_male_4)==T] <- round((DATA$gender_male_4_.[is.na(DATA$gender_male_4)==T]*DATA$number_4_ttt[is.na(DATA$gender_male_4)==T])/100)

DATA$gender_male_5[is.na(DATA$gender_male_5)==T] <- round((DATA$gender_male_5_.[is.na(DATA$gender_male_5)==T]*DATA$number_5_ttt[is.na(DATA$gender_male_5)==T])/100)

DATA$gender_male_6[is.na(DATA$gender_male_6)==T] <- round((DATA$gender_male_6_.[is.na(DATA$gender_male_6)==T]*DATA$number_6_ttt[is.na(DATA$gender_male_6)==T])/100)

DATA$all_gender_male[is.na(DATA$all_gender_male)==T] <- round((DATA$all_gender_male_.[is.na(DATA$all_gender_male)==T]*DATA$number_ttt[is.na(DATA$all_gender_male)==T])/100)

DATA$all_gender_male_all_2 <- apply(DATA[,c("gender_male_1","gender_male_2","gender_male_3","gender_male_4","gender_male_5","gender_male_6")],1,sum,na.rm=T)

DATA$all_gender_male[is.na(DATA$all_gender_male)==T] <- DATA$all_gender_male_all_2[is.na(DATA$all_gender_male)==T]

DATA$all_gender_male[DATA$all_gender_male==0] <- NA 

DATA$all_gender_male_rate <- DATA$all_gender_male/DATA$number_ttt

DATA_clean <- subset(DATA, select = c("IDU","PMID", "study.name", "author_name", "type_Obsv_0_trial_1", "mono0_multi1", "published", "number_follow_up", "d_adm_dis", "rehosp_30_all.cause_all", "death_30_all", "rehosp_365", "death_365", "number_ttt", "all_gender_male", "TIME", "EF_all", "HF_all_n", "all_HT_n", "all_DM_n", "all_CAD.IHD_n", "all_COPD_n", "all_Afib_n", "all_ACE_ARB_adm_n", "all_ACE_ARB_adm_rate", "A_all_diuretic_n", "A_all_diuretic_rate", "A_all_digoxin_n", "A_all_digoxin_rate", "A_all_B_n", "A_all_B_rate", "SBP_all", "AGE_all", "HR_all", "all_gender_male_rate"))

################
## Adding information about countries involved in each study
################

PMID <- unique(DATA_clean$PMID)

countries <- read.csv2("data/countries.csv", na.strings =c(""," ","NA")) 
# for each study it has 1 in the countries column where it occured and NA otherwise

# creates a list of vector with one vector per study                      
countries_list <- countries[countries$PMID %in% PMID,] %>% # using only studies included in data
  dplyr::select(-c(Author_name:verif)) %>%
  group_by(PMID) %>% 
  nest() %>%
  mutate(list = map(data, function(df) {colnames(df)[which(df == 1)]}))


max(sapply(countries_list$list, function(x) length(x))) # 31 countries max per study

# creates a dataset with country names for each study                    
countries_data <- as.data.frame(matrix(NA, ncol = 32, nrow = nrow(countries_list), 
                                       dimnames = list(NULL, c("PMID",paste0("country", 1:31)))))
countries_data$PMID <- PMID
for(i in 1:nrow(countries_data)) {
  if(length(countries_list$list[[i]]) > 0){
    countries_data[i,2:(length(countries_list$list[[i]])+1)] <- countries_list$list[[i]]
  }
}

# joins with dataset and export
new_DATA <- full_join(DATA_clean, countries_data, by = "PMID")

# creates continent variables
new_DATA$continent <- countrycode(new_DATA$country1, origin = "country.name",
                                  destination = "continent")
new_DATA <- new_DATA %>% mutate(continent = ifelse(country1 == "USA" | country1 == "Canada", "North_America", 
                                                   ifelse(continent == "Americas", "South_America",
                                                          ifelse(country1 == "Armenia" | country1 == "Turkey", "Europe",
                                                                 continent))))

# checking if multi-countries studies are mono-continental
multicount <- new_DATA %>% filter(!is.na(country2)) %>% 
  dplyr::select(c(PMID, contains("country"), continent)) %>% 
  mutate(multicontinent = c("Yes", "No", "No", "No", "No", "Yes", "Yes", "No",
                            "Yes", "Yes", "Yes", "Yes","Yes", "Yes", "No", "Yes",
                            "Yes", "Yes", "Yes", "No", "No", "No", "Yes", "Yes")) # manually creating multicontinent variable

new_DATA <- full_join(new_DATA, multicount %>% dplyr::select(PMID, multicontinent), by = "PMID")

write.csv2(new_DATA, "data/clean_data_Kimmoun_et_al.csv")
