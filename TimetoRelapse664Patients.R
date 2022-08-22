########## Survial Analysis ##################################################################################
# Patients from Metzeler et al. 664 AML patients

##### 664 patients from Metzeler et al. ! ###########################################################################
rm(list=ls());  # empty workspace
setwd("/data/htp/A07/AML_triplets/RNA_Seq/TimeToRelapse")


library(survminer)
library(survival)

library(dplyr)
library(readxl)



# DNMT3A mutated patients
mut<-read_excel("/data/htp/A07/AML_triplets/auxilliary_data/Supplement_Exome_Greif/SupS3.xlsx", sheet = 1) %>%
  mutate(`% VAF* Rem` = as.numeric(`% VAF* Rem`)) %>%
  filter(Gene=="DNMT3A") %>%
  mutate(Persistent=ifelse(`% VAF* Rem` >5, "Persistent","Non")) %>%
  slice(-24)

mut$Persistent[17]<-"Non"
# patient CN-036, we dont have data for remission --> no data for persitant
# CN-046 has 2 DNMT3A mutations, one persistent on non persistent. this patient counts as persistent. 


var<-read_excel("/data/htp/A07/AML_triplets/auxilliary_data/Supplement_Exome_Greif/SupS2.xlsx", sheet = 1) %>%
  slice(1:50) %>%
  mutate(DNMT3A=ifelse(UPN %in% mut$UPN, "MUT", "WT")) %>%
  left_join(mut[,c(1,16)], by="UPN") %>%
  mutate(stat=ifelse(is.na(`Time to Relapse`)==TRUE, 0, 1)) %>%
  mutate(`Time to Relapse`=as.numeric(`Time to Relapse`)) %>%
  mutate(Persistent = replace_na(Persistent, "WT")) %>%
  mutate(Age= ifelse(`Age at Diagnosis`>60, "Old", "Young"))

rownames(var)<-var$UPN
colnames(var)[2]<-"Age_at_Diagnosis"
colnames(var)[5]<-"Karyotype_at_Relapse"
colnames(var)[7]<-"AML_Type"
colnames(var)[8]<-"ELN_Risk_Group"
colnames(var)[9]<-"Induction_Therapy"
colnames(var)[10]<-"TimeToRelapse"
colnames(var)[11]<-"Evolutionary_Pattern"

var<-var %>% mutate(Age_at_Diagnosis=as.numeric(Age_at_Diagnosis))

var<-var %>% as.data.frame(var)
# add column stat: 1 means  patients stayed in study , 0 means patients were censored meaning no data for the time to relapse
############################################################################################################

##### Kaplan Meyer ##########################################################################################
#### DNMT3A MUT vs WT
fit <- survfit(Surv(TimeToRelapse, stat) ~ DNMT3A, data = var)
surv<-ggsurvplot(fit, data = var, pval=TRUE)
surv
# plot will split by DNMT3A

#### DNMT3A MUT Persistent vs DNMT3A MUT non persistent
var2<-var %>% filter(DNMT3A=="MUT")

fit <- survfit(Surv(TimeToRelapse, stat) ~ Persistent, data = var2)
surv<-ggsurvplot(fit, data = var2, pval=TRUE)
surv
# there is a difference between the DNMT3A persistent and non persistent ones
# same as in the Exome Paper Greif et al.

##### DNMT3A MUT persistent against WT?
var3<-var %>% filter(Persistent %in% c("Persistent", NA))

fit <- survfit(Surv(TimeToRelapse, stat) ~ DNMT3A, data = var3)
surv<-ggsurvplot(fit, data = var3, pval=TRUE)+
  labs(y= "Relapse probability")
surv
# also persistnt DNMT3A only against the WT patients not significant!
##################################################################################################

##### Cox proportional Hazards Model #############################################################

#### which covariates to include #################################################################
covariates <- c("Age_at_Diagnosis", "Gender",  "Karyotype_at_Relapse", "AML_Type", "ELN_Risk_Group", "Induction_Therapy", "Evolutionary_Pattern",
                "DNMT3A","Persistent")
# Age at diagnosis (numeric)
cox <- coxph(Surv(TimeToRelapse, stat) ~ Age_at_Diagnosis,  data = var)
summary(cox)
# not significant

# Age (categorical: Old, Young)
cox <- coxph(Surv(TimeToRelapse, stat) ~ Age_at_Diagnosis,  data = var)
summary(cox)
# not significant

# Gender
cox <- coxph(Surv(TimeToRelapse, stat) ~ Gender,  data = var)
summary(cox)
# not significant

# Karyotype_at_Relapse
cox <- coxph(Surv(TimeToRelapse, stat) ~ Karyotype_at_Relapse,  data = var)
summary(cox)
# too many different karyotypes 

# AML Type
cox <- coxph(Surv(TimeToRelapse, stat) ~ AML_Type,  data = var)
summary(cox)
# not significant

# ELN Risk Group
cox <- coxph(Surv(TimeToRelapse, stat) ~ ELN_Risk_Group,  data = var)
summary(cox)
# SIGNIFICANT

# Induction_Therapy
cox <- coxph(Surv(TimeToRelapse, stat) ~ Induction_Therapy,  data = var)
summary(cox)
# not significant

# Evolutionary Pattern
cox <- coxph(Surv(TimeToRelapse, stat) ~ Evolutionary_Pattern,  data = var)
summary(cox)
# SIGNIFICANT

# DNMT3A
cox <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A,  data = var)
summary(cox)
# not significant

# Persistent
cox <- coxph(Surv(TimeToRelapse, stat) ~ Persistent, data = var)
summary(cox)
curve<-ggadjustedcurves(cox, data=var, variable = "Persistent")
curve
ggsave("DNMT3A_Persistent_ExomePatients.pdf", plot = curve, width= 150, height= 100, units= "mm", device = "pdf",
       path = "/data/htp/A07/AML_triplets/RNA_Seq/TimeToRelapse")

# SIGNIFICANT
# but small samples size, only 8 for non, this can drive pvalues up
# persistent vs non is significant 
# non vs wt is not!
# there is a visible difference but not quantified.

# Multivariate model
cox <- coxph(Surv(TimeToRelapse, stat) ~ DNMT3A+ ELN_Risk_Group+ Evolutionary_Pattern+ Persistent,  data = var)
summary(cox)
# in the context of all covariates only the Evolutionary Pattern shows significance
# this means that the contribution of the Evolutionary pattern is higher than the others, 
# this causes the others to miss the 10 % mark
# also the samples size for non peristent is only 8
# still the persistent has an influence and can be shown in an univariate model!
# why it throws NA for WT now i dont understand

ggadjustedcurves(cox, data=var, variable = "Evolutionary_Pattern")
##############################################################################################################

#### test for proportional hazards assumption ###############################################################
test.ph <- cox.zph(cox)
test.ph
# non of them is significant, therefore all covariates have proprtional hazards.
#############################################################################################################



# add column stat: 1 means  patients stayed in study , 0 means patients were censored meaning no data for the time to relapse

# reason why we need cox is i can include covariates in the model!
# kaplan meyer will only add strata!

#### on the cohort of 664 patients?
#P values were calculated from a univariate Cox proportional hazards model
#stratified for trial arm, using the Wald test.

#### save the session Info #####################################################################
writeLines(capture.output(sessionInfo()), "SessionInfo_TimeToRelapeExomePatients.txt")
