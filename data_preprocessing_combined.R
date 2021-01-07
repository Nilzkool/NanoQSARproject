#########################################################################
# Script to  pre-process dataset combined from all the case studies.
#########################################################################
### clear workspace 

## Clear Variables
rm(list = ls())

## Clear console
cat("/014") 

## Load libraries
library(readxl)
library(glmulti)
library(glmnet)
library(xlsx)

#########################################################################
folder<- 'C:/Folder name/'
df_DNA<-read_excel(paste0(folder,'all case studies_1.xlsx'), sheet='DNA damage')
df_via<-read_excel(paste0(folder,'all case studies_1.xlsx'), sheet='Viability')

#########################################################################
# Subset req_cols
req2_DNA<-c('C_type','Conc','P_size', 'Z_average','Z_average_CCM','Cryst_str', 'Coating','Time','DNA_damage')
req2_via<-c('C_type','Conc','P_size', 'Z_average','Z_average_CCM','Cryst_str', 'Coating', 'Viability')

df2_DNA<-df_DNA[,req2_DNA]
df2_via<-df_via[,req2_via]
#########################################################################
# Process data type
df2_DNA$Cryst_str<-as.factor(df2_DNA$Cryst_str)
df2_DNA$Coating<-as.factor(df2_DNA$Coating)
df2_DNA$C_type<-as.factor(df2_DNA$C_type)
#df2_DNA$Conc<-as.factor(df2_DNA$Conc)

df2_via$Cryst_str<-as.factor(df2_via$Cryst_str)
df2_via$Coating<-as.factor(df2_via$Coating)
df2_via$C_type<-as.factor(df2_via$C_type)
#df2_via$Conc<-as.factor(df2_via$Conc)

sup_data<-list('DNA' = df2_DNA,
               'via' = df2_via)

#########################################################################
# Save processed datasets
saveRDS(sup_data,file=paste0(folder,"processed_data.Rda"))

#####################################################"""
# Train validation split (70:30)
# Select df2_via or df2_DNA
dfReq<-df2_via
smp_size <- floor(0.70 * nrow(dfReq))

## set the seed to make your partition reproducible
set.seed(2)
train_ind <- sample(seq_len(nrow(dfReq)), size = smp_size)


train <- dfReq[train_ind, ]
test <- dfReq[-train_ind, ]

## save accordingly
dev_data<-list('train' = train, 'test' = test)
saveRDS(dev_data,file=paste0(folder,"Via_dev_data.Rda")) #or "DNA_dev_data.Rda"

