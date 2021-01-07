#########################################################################
# Script to  pre-process supplementary files 1,2 and 3 before analysis
#########################################################################
### clear workspace 

## Clear Variables
rm(list = ls())

## Clear console
cat("/014") 

## Load libraries
library(readxl)
library(xlsx)
#########################################################################
## Load the case study files and appropriate sheets
folder<- 'C:/Folder name/'
df1<-read_excel(paste0(folder,'Supplementary data file 1.xlsx'))
df2_DNA<-read_excel(paste0(folder,'Supplementary data file 2.xlsx'), sheet='DNA damage')
df2_via<-read_excel(paste0(folder,'Supplementary data file 2.xlsx'), sheet='Viability')

df3_DNA<-read_excel(paste0(folder,'Supplementary data file 3.xlsx'), sheet='DNA damage')
df3_via<-read_excel(paste0(folder,'Supplementary data file 3.xlsx'), sheet='Viability')

#########################################################################
# Select required fields 

req1<-c('C_type', 'Conc','P_size', 'Zeta_po','Feret_min','Aspect_Ratio','Z_average','PTA_size','Z_average_CCM', 'Viability','DNA_damage')

req2_DNA<-c('P_size', 'Z_average','Z_average_CCM','Cryst_str', 'DNA_damage', 'Viability')
req2_via<-c('P_size', 'Z_average','Z_average_CCM','Cryst_str', 'Viability')

req3_DNA<-c('P_size', 'Zeta_po', 'Z_average','Z_average_CCM','Cryst_str','Coating','DNA_damage')
req3_via<-c('P_size', 'Zeta_po', 'Z_average','Z_average_CCM','Cryst_str','Coating','Viability')


df1<-df1[,req1]
df2_DNA<-df2_DNA[,req2_DNA]
df2_via<-df2_via[,req2_via]
df3_DNA<-df3_DNA[,req3_DNA]
df3_via<-df3_via[,req3_via]

sup_data<-list('sup1' = df1, 
               'sup2_DNA' = df2_DNA,
               'sup2_via' = df2_via,
               'sup3_DNA' = df3_DNA,
               'sup3_via' = df3_via)

#########################################################################
# Modify column types
df1$C_type<-as.factor(df1$C_type)
df1$S_type<-as.factor(df1$S_type)

#########################################################################
# Save processed datasets
folder<- 'C:/Users/u0113548/Desktop/Extra PhD projects/Nanotoxicity project/Analysis/Analysis 6 May 2020/'
saveRDS(sup_data,file=paste0(folder,"processed_data.Rda"))

#####################################################"""
# Train validation split (70:30), only for df1 or supplementary file 1
dfReq<-df1 
smp_size <- floor(0.70 * nrow(dfReq))

## set the seed to make your partition reproducible
set.seed(2)
train_ind <- sample(seq_len(nrow(dfReq)), size = smp_size)
train <- dfReq[train_ind, ]
test <- dfReq[-train_ind, ]

## save
dev_data<-list('train' = train, 'test' = test)
saveRDS(dev_data,file=paste0(folder,"sup1_dev_data.Rda"))