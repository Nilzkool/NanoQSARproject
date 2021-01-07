#########################################################################
# Script for ML analysis on supplementary files 1,2 and 3
#########################################################################
### clear workspace 

## Clear Variables
rm(list = ls())

## Clear console
cat("\014") 

## Load libraries
library(readxl)
library(ggplot2)
library(grid)
library(gridExtra)
library(psych)
require(reshape2)
library(glmulti)
library(glmnet)
library(caret)
library(xlsx)
library(praznik)
library(factoextra)
library(Rtsne)
library(cluster)
library(gtools)
library(ggpubr)
library(randomForest)
library(Metrics)


#########################################################################
# Load processed dataset
folder<-'C:/Folder name/'
dev_data<-readRDS(file = paste0(folder,"sup1_dev_data.Rda"))
dev_train<- dev_data$train
dev_train$C_type<-NULL
dev_train$Conc<-NULL
dev_train<-na.omit(dev_train)

# NOTE: for supplementary files 2 and 3
# dev_data<-readRDS(file = paste0(folder,"processed_data.Rda"))
# dev_train<- dev_data['sup2_DNA'] # Note for viability use 'sup2_via'
########################################################################
# Generate MRMR rankings

# Convert continuous variables to quantiles
for (col in colnames(dev_train))
{   if (col=='DNA_damage' | col=='Viability')
{
  dev_train[[col]]<-quantcut(dev_train[[col]], q=4,na.rm=TRUE)
  dev_train[[col]]<- factor(dev_train[[col]])
} 
  else if (col=='Z_average_CCM')
  {
    dev_train[[col]]<-quantcut(dev_train[[col]], q=2,na.rm=TRUE)
    dev_train[[col]]<- factor(dev_train[[col]])
  }
  else {
    dev_train[[col]]<- factor(dev_train[[col]])
    
  }
  
}

#Select desired response: Viability or DNA_damage
resp = 'Viability'

if (resp=='DNA_damage')
{ dev_train$Viability<-NULL} else
{dev_train$DNA_damage<-NULL }

# Compute MRMR rankings
x<-dev_train[, -which(names(dev_train) == resp)]
y<-dev_train[,c(resp)]

mrmr<- MRMR(x, y[[resp]], k = ncol(dev_train)-1, threads = 0)
param_order_list<-names(mrmr$selection)

#Prepare MRMR plot
dfPlot<-data.frame(Parameter=param_order_list, Score=mrmr$score)
rownames(dfPlot) <- c()
dfPlot$Parameter <- factor(dfPlot$Parameter)

# MRMR plot
N=nrow(x)+1
title<-paste0(paste0('Feature importance for ', resp,'(N=',as.character(N),')'))
p<-ggplot(dfPlot, aes(x = Parameter, y = Score, group=1)) + 
  geom_line()+geom_point()+theme_bw()+ coord_flip()+
  scale_x_discrete(limits= rev(param_order_list))+
  xlab("")+
  theme_bw()+
  theme(axis.text.y =element_text(face=c('bold'),size=11,angle = 0))+
  ylab("MRMR score")+
  geom_hline(yintercept = 0,  colour="black",linetype = "dashed")+
  scale_y_continuous(breaks = round(seq(min(dfPlot$Score), max(dfPlot$Score)+0.1, by = 0.1),1)) +
  ggtitle(title)
p

#########################################################################
# ML modelling

# Select desired response: 'DNA_damage' or 'Viability'
resp<-'DNA_damage' 
dev_train<- dev_data$train

# NOTE: test set was only created for supplemtary file 1
# for files 2 and 3, test set is ignored as results are reported on 
# the entire datset. We use one case here for code generalizability

if("test" %in% colnames(dev_data)) 
  {dev_valid<-dev_data$valid} else
  {dev_valid<- dev_train[1,]} 


# Stepwise modelling, define sequence of parameter insertion 
# starting with C_type and Conc as default parameters
param_order_list_orig<-param_order_list
param_order_list<-append(c('C_type','Conc'), param_order_list)

rsq <- function (x, y) cor(x, y) ^ 2

# Pre-porcess response data by dividing by 100
train<- dev_data$train
valid<-dev_data$test
train[[resp]]<-train[[resp]]/100
valid[[resp]]<-valid[[resp]]/100

# Define output table
dfModResults<-data.frame(matrix(ncol = 9))
colnames(dfModResults)<-c('Response', 'Model' ,'Inputs', 'Train_R2','Train_RMSE' ,'Train_MAE', 'Valid_R2', 'Valid_MAE','Valid_RMSE')

# Define type of modelling: Linear Regression (LR) and Random Forests (RF)
models<-c('LR', 'RF')

# Define max no features to be added in stepwise modelling
n_feats=5

# Loop and generate results
for (modeltyp in models)
{ 
  for (i in 2:n_feats)
  { 
    params<-param_order_list[1:i]
    
    # prep train and valid sets
    dfTrain<-na.omit(train[,append(params, resp)])
    dfValid<-na.omit(valid[,append(params, resp)])
    
    # preprocess and fit model
    if (modeltyp=='LR')
    { nums <- unlist(lapply(dfTrain, is.numeric))  
    nums[length(nums)]<-FALSE
    pp<-preProcess(dfTrain[,nums], method = c("center", "scale"))
    dfTrain[,nums]<-predict(pp,dfTrain[,nums])
    dfValid[,nums]<-predict(pp,dfValid[,nums])
    model<-lm(paste(resp,"~ ."), data=dfTrain)
    
    } else {
      x_train<-dfTrain[,-which(names(dfTrain)==resp)]
      y_train<-dfTrain[[resp]]
      model<-randomForest(x=x_train, y=y_train,ntree=100)
    }
    
    #Make  model preds
    pred_train<-predict(model, newdata =dfTrain[,-which(names(dfTrain)==resp)])
    pred_valid<-predict(model, newdata =dfValid[,-which(names(dfValid)==resp)])
    
    # Evaluate preds
    train_r2<-rsq(pred_train,dfTrain[[resp]] )
    train_MAE<-mae(pred_train,dfTrain[[resp]])
    train_rmse<-rmse(pred_train,dfTrain[[resp]])
    
    valid_r2<-rsq(pred_valid,dfValid[[resp]] )
    valid_MAE<-mae(pred_valid,dfValid[[resp]])
    valid_rmse<-rmse(pred_valid,dfValid[[resp]])
    #valid_se <- sd(abs(pred_valid - dfValid[[resp]]))/sqrt(nrow(dfValid))
    
    #Update results
    newRow<- data.frame(Response = resp, 
                        Model = modeltyp,
                        Inputs = paste(params,collapse ='+'),
                        Train_R2 = train_r2, 
                        Train_MAE = train_MAE,
                        Train_RMSE=train_rmse,
                        Valid_R2 = valid_r2,
                        Valid_MAE = valid_MAE,
                        Valid_RMSE = valid_rmse)
    
    dfModResults<-rbind(dfModResults, newRow)
  }
}

# Final results table
dfModResults<-dfModResults[2:nrow(dfModResults),]



#########################################################################
# Y randomnization test of final model in supplementary file 1 
# Select desired response: 'DNA_damage' or 'Viability'
resp<-'DNA_damage'

if (resp=='DNA_damage')
{
  req_cols<-c('C_type', 'Conc','Aspect_ratio')
} else
{req_cols<-c('C_type', 'Conc','Zeta_po')}

# Reload dataset
folder<-'C:/folder name/'
dev_data<-readRDS(file = paste0(folder,"sup1_dev_data.Rda"))
train<- dev_data$train
valid<-dev_data$test
train[[resp]]<-train[[resp]]/100
valid[[resp]]<-valid[[resp]]/100

rsq <- function (x, y) cor(x, y) ^ 2

# Orininal train
dfTrain<-na.omit(train[,append(params, resp)])
dfValid<-na.omit(valid[,append(params, resp)])


x_train<-dfTrain[,-which(names(dfTrain)==resp)]
y_train<-dfTrain[[resp]]
model<-randomForest(x=x_train, y=y_train,ntree=100)

#Make  model predictions
pred_train<-predict(model, newdata =dfTrain[,-which(names(dfTrain)==resp)])
pred_valid<-predict(model, newdata =dfValid[,-which(names(dfValid)==resp)])

train_r2_orig<-rsq(pred_train,dfTrain[[resp]] )
valid_r2_orig<-rsq(pred_valid,dfValid[[resp]] )

# define k randomnizations
k=100
train_r2_list<-c()
valid_r2_list<-c()
for (i in 1:k)
{ # randomnize response
  dfTrain<-na.omit(train[,append(params, resp)])
  dfValid<-na.omit(valid[,append(params, resp)])
  dfTrain[[resp]]<-sample(dfTrain[[resp]])
  dfValid[[resp]]<-sample(dfValid[[resp]])
  
  x_train_orig<-dfTrain[,-which(names(dfTrain)==resp)]
  y_train<-dfTrain[[resp]]
  model<-randomForest(x=x_train, y=y_train,ntree=100)
  
  #Make  model predictions randomnized response
  pred_train<-predict(model, newdata =dfTrain[,-which(names(dfTrain)==resp)])
  pred_valid<-predict(model, newdata =dfValid[,-which(names(dfValid)==resp)])
  
  train_r2<-rsq(pred_train,dfTrain[[resp]] )
  valid_r2<-rsq(pred_valid,dfValid[[resp]] )
  
  train_r2_list<-append(train_r2_list,train_r2)
  valid_r2_list<-append(valid_r2_list,valid_r2)
}

# Analyze results

describe(train_r2_list)
describe(valid_r2_list)