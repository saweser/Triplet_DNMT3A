

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/AML_triplets/RNA_Seq/TimeToRelapse")


library(survminer)
library(survival)
library(dplyr)
library(readxl)

# CN patients
anno<-readRDS("/data/htp/A07/AML_triplets/RNA_Seq/AnnoAndCounts/anno_CN.rds")%>%
  mutate(stat=ifelse(is.na(TimeToRelapse)==TRUE, 0, 1)) %>%
  filter(SampleType=="Diagnosis")

# all patients from RNA Seq
anno<-readRDS("/data/htp/A07/AML_triplets/RNA_Seq/AnnoAndCounts/anno50patients.rds")%>%
  mutate(stat=ifelse(is.na(TimeToRelapse)==TRUE, 0, 1)) %>%
  filter(SampleType=="Diagnosis")

rownames(anno)<-anno$PatientID
# add column stat: 1 means  patients stayed in study , 0 means patients were censored meaning no data for the time to relapse


fit <- survfit(Surv(TimeToRelapse, stat) ~ DNMT3A, data = anno)
surv<-ggsurvplot(fit, data = anno, pval=TRUE)+
  labs(y= "Relapse probability")
surv

# plot will split by DNMT3A

##################################################################################################
# plot the cox model
cox <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A,  data = anno)
ggadjustedcurves(cox, data=anno, variable = "DNMT3A")
# without variable does not split into DNMT3A Status


cox <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A +strata(DNMT3A),  data = anno)
ggadjustedcurves(cox, data=anno, variable = "DNMT3A")

#cox <- coxph(Surv(TimeToRelapse, stat) + strata(DNMT3A),  data = anno)
#ggadjustedcurves(cox, data=anno, variable = "DNMT3A")

# what exactely does strata do? 
# this does look exactely the same as the kaplan meyer

ggforest(cox, data = anno)



ggforest(fit.coxph, data = anno)



ggsurvplot(surv_fit(cox, data=anno), conf.int=TRUE, legend.labs=c("DNMT3A_WT", "DNMT3A_MUT"),  palette= "#2E9FDF",
           ggtheme = theme_minimal())


#### stratify for the treatment approach

therapy<-read_excel("/data/htp/A07/AML_triplets/auxilliary_data/Supplement_Exome_Greif/SupS2.xlsx") %>%
  mutate(UPN=gsub("-", "\\.", UPN)) %>%
  filter(!is.na(UPN)) %>%
  mutate(therapy_cat=c(rep("AD",11), rep("M",21), rep("T",18)))

colnames(therapy)[9]<-"Induction_Therapy"
colnames(therapy)[8]<-"ELN_Risk_Group"
colnames(therapy)[7]<-"AML_Type"

anno<-readRDS("/data/htp/A07/AML_triplets/RNA_Seq/AnnoAndCounts/anno_CN.rds") %>%
  left_join(therapy[,c(1,7,8,9,12)], by=c("PatientID"="UPN")) %>%
  mutate(stat=ifelse(is.na(TimeToRelapse)==TRUE, 0, 1)) %>%
  filter(SampleType=="Diagnosis") %>%
  mutate(age_cat=ifelse(AgeAtDiagnosis>60, "old", "young" ))

rownames(anno)<-anno$SampleName
# add column stat: 1 means  patients stayed in study , 0 means patients were censored meaning no data for the time to relapse

fit <- survfit(Surv(TimeToRelapse, stat) ~ DNMT3A + Induction_Therapy, data = anno)
surv<-ggsurvplot(fit, data = anno, pval=TRUE)+
  labs(y= "Relapse probability")
surv

# reason why we need cox is i can include covariates in the model!
# kaplan meyer will only add strata!

fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A, data = anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "DNMT3A")

fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A + strata (Induction_Therapy), data = anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "DNMT3A")

fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A + Induction_Therapy, data = anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "DNMT3A")

fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A + Induction_Therapy + age_cat, data = anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "DNMT3A")

fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A + Induction_Therapy + ELN_Risk_Group, data = anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "DNMT3A")


fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A + Induction_Therapy+ age_cat + FLT3.Status, data = anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "DNMT3A")


# Gender has a significant effect on Time to Relapse
fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A + Gender, data = anno)
fit.coxph <- coxph(Surv(TimeToRelapse, stat) ~ Gender, data = anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "Gender")


ggforest(fit.coxph, data=anno)
ggadjustedcurves(fit.coxph, data=anno, variable = "Induction_Therapy")

# Gender makes it worse
# Age does not seem to make a difference
# ELN and AML_Type  a bit worse
# 

summary(fit.coxph)
cox.zph(fit.coxph)








#### on the cohort of 664 patients?
#P values were calculated from a univariate Cox proportional hazards model
#stratified for trial arm, using the Wald test.




# find out how they filered in the 664 patients !
# here the DNMT3A mut is significant


# 50 patients paper: time to relapse diff in DNMT3A persistent vs non persistent ones (19:8)
anno<-read_excel("/data/htp/A07/AML_triplets/auxilliary_data/Supplement_Exome_Greif/SupS2.xlsx")%>%




#### save the session Info #####################################################################
writeLines(capture.output(sessionInfo()), "sessionInfoAnno.txt")





####!!! try with the whole cohort once all 50 patients!!!
#### 6leave out the step the patients in counts only


#### cox regression instead of kaplan meyer







https://www.datacamp.com/community/tutorials/survival-analysis-R