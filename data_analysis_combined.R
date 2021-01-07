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
# Select 'DNA_dev_data.Rda' or 'Via_dev_data.Rda'
req_file='DNA_dev_data.Rda'
# MRMR rankings for sup 1 train data
folder<- 'C:/Folder name/'
dev_data<-readRDS(file = paste0(folder,req_file))
dev_train<- dev_data$train
dev_train$C_type<-NULL
dev_train$Conc<-NULL

if (req_file=='DNA_dev_data.Rda')
{ dev_train$Time<-NULL
}
dev_train<-na.omit(dev_train)

########################################################################
# Generate MRMR rankings

# Convert continuous variables to quantiles
for (col in colnames(dev_train))
{   if (col=='DNA_damage' | col=='Viability')
   {
    dev_train[[col]]<-quantcut(dev_train[[col]], q=4,na.rm=TRUE)
    dev_train[[col]]<- factor(dev_train[[col]])
    } 
    else if (col=='Z_average_CCM' | col=='P_size')
    {
      dev_train[[col]]<-quantcut(dev_train[[col]], q=4,na.rm=TRUE)
      dev_train[[col]]<- factor(dev_train[[col]])
    }

}
  
# Desired response (Viability or DNA_damage) is selected based on file name 
if (req_file=='Via_dev_data.Rda')
{ resp = 'Viability' } else
{ resp= 'DNA_damage'}

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
title<-paste0(paste0('Feature importance for ', resp,' (train N=',as.character(N),')'))
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

# Desired response (Viability or DNA_damage) is selected based on file name 
if (req_file=='Via_dev_data.Rda')
{ resp='Viability'
  dev_file = 'Via_dev_data.Rda' } else
  { dev_file= 'DNA_dev_data.Rda'
    resp='DNA_damage'
    time_req=0}

# Reload dataset
dev_data<-readRDS(file = paste0(folder,dev_file))
dev_train<- dev_data$train
dev_valid<-dev_data$valid


# Stepwise modelling, define sequence of parameter insertion 
# starting with C_type and Conc as default parameters
param_order_list_orig<-param_order_list
param_order_list<-append(c('C_type','Conc','Time'), param_order_list)

rsq <- function (x, y) cor(x, y) ^ 2

#Pre-porcess response data by dividing by 100train<- dev_data$train
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
for (modeltyp in models)
{ 
  for (i in 3:n_feats)
  { #print(i)
    #modeltyp='LR'
    #i=2
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
# Y randomnization test of final model
resp<-'DNA_damage'

if (resp=='DNA_damage')
{
  req_cols<-c('C_type', 'Conc','Aspect_ratio') # Features based on final model for 'DNA_damage'
} else
{req_cols<-c('C_type', 'Conc','Zeta_po')} # Features based on final model for Viability

# Reload dataset
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

#Make  model preds
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