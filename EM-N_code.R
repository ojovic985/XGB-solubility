
### You have to do this input loading part for yourself ####
...## Load your input training set LogS vector and label that vector as "B1"
...## Load your input test set LogS vector and label that vector as "Btest"
...## Load your input training set scaled molecular descriptor matrix and label it as "Fe2"
...## Load your input test set scaled molecular descriptor matrix according to training set and label it as "Xb"
...## Load your input database scaled molecular descriptor matrix according to training set and label it as "Ge1"
...## Matrices Fe2, Xb and Ge1 must have the same number of variables and exactly the same ordered variables 
###

B.pred_test.error<- array(rep(NA,length(Btest)*kfold*kfold),dim=c(length(Btest),kfold,kfold))
B.pred_test.error2<- array(rep(NA,length(Btest)*kfold*kfold),dim=c(nrow(Ge1),kfold,kfold))
colnames(Fe2)<-seq(1,ncol(Fe2),by=1)
colnames(Xb)<-seq(1,ncol(Xb),by=1)
colnames(Ge1)<-seq(1,ncol(Ge1),by=1)
xgb_train = xgb.DMatrix(data = Fe2, label = B1)
xgb_test = xgb.DMatrix(data = Xb, label = Btest)
Ge1<-as.matrix(Ge1)

max_depth_opt<- ### input your optimized max_depth parameter
eta_opt<- ### input your optimized eta parameter
ser_tv1<- ### input your optimized number of rounds parameter


set.seed(71) #set your seed in order to be able to exactly reproduce your results in the future

ABS_endpoint<-matrix(rep(NA, length(B1)*kfold),nrow=kfold)
Error_endpoint <-matrix(rep(NA, length(B1)*kfold),nrow=kfold)
Tr<-seq(1,length(B1),by=1)
Sys.time()

kfold<-20

for (ii in 1:kfold) {
cv = xgb.cv(data = xgb_train, nrounds = ser_tv1, nthread = 2, nfold = kfold, max_depth = max_depth_opt, eta = eta_opt, objective = "reg:squarederror",prediction=TRUE)

telta2<-cv[8]$pred
ABS_endpoint[ii,]<-abs(telta2-B1)

Error_Btrain<- ABS_endpoint[ii,]
Error_xgb_train = xgb.DMatrix(data = Fe2, label = Error_Btrain)
Error_model_xgboost = xgb.cv(data = Error_xgb_train, nthread = 2, nfold = kfold, max.depth = max_depth_opt, eta=eta_opt, nrounds = ser_tv1, objective = "reg:squarederror",prediction=TRUE)
Error_endpoint[ii,]<- abs(Error_model_xgboost[8]$pred)
AA<- Error_model_xgboost[7]
for (j in 1:kfold) {
xgb_train1 = xgb.DMatrix(data = Fe2[-AA$folds[[j]],], label = Error_Btrain[Tr[-AA$folds[[j]]]])
model_XGB<-xgb.train(data = xgb_train1, nthread = 2, max.depth = max_depth_opt, eta=eta_opt, nrounds = ser_tv1, objective = "reg:squarederror")

B.pred_test.error[,ii,j]<-predict(model_XGB,Xb)
B.pred_test.error2[,ii,j]<-predict(model_XGB,Ge1)
}
}

error_te<-array(nrow(Xb))
error_te_database<-array(nrow(Ge1))

beta<-0 # Set sensitivity factor (beta)

for (k in 1:nrow(Ge1)) {
error_te_database[k]<-beta+(mean(abs(B.pred_test.error2[k,,])))} # produces sigma-test for database

for (k in 1:nrow(Xb)) {
error_te[k]<-beta+(mean(abs(B.pred_test.error[k,,])))} # produces sigma-test for test set

alfa_value<-array(ncol(Error_endpoint))
for (i in 1:ncol(Error_endpoint)) {

alfa_value[i]<-mean(ABS_endpoint[,i])/ (beta+ (mean(Error_endpoint[,i])) )}  # produces alfa values (NC-s)

### Not necessary part, but for training set prediction intervals if one is interested too
ABS_mean<- array(ncol(Error_endpoint))
Error_function<-array(ncol(Error_endpoint))
for (i in 1:ncol(Error_endpoint)) {
Error_function[i]<-beta + (mean(Error_endpoint[,i]))
ABS_mean[i]<-mean(ABS_endpoint[,i])}


#### If your training data set does NOT have the number of training samples as 99 + whole number x 100  then:

consta<-100 # 99% conf.
y2<-(length(B1)- (((length(B1)+1) %/% consta)-1))/(length(B1)+1)
y1<-(length(B1)- (((length(B1)+1) %/% consta)))/(length(B1)+1)
x2<-(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta])
x1<-(sort(alfa_value,decreasing=TRUE)[((length(B1)+1) %/% consta)+1])
fact99<-((((1-(1/consta))-y1)*(x2-x1)/(y2-y1))+x1)/sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta]

consta<-20 # 95% conf.
y2<-(length(B1)- (((length(B1)+1) %/% consta)-1))/(length(B1)+1)
y1<-(length(B1)- (((length(B1)+1) %/% consta)))/(length(B1)+1)
x2<-(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta])
x1<-(sort(alfa_value,decreasing=TRUE)[((length(B1)+1) %/% consta)+1])
fact95<-((((1-(1/consta))-y1)*(x2-x1)/(y2-y1))+x1)/sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta]


consta<-10 # 90% conf.
y2<-(length(B1)- (((length(B1)+1) %/% consta)-1))/(length(B1)+1)
y1<-(length(B1)- (((length(B1)+1) %/% consta)))/(length(B1)+1)
x2<-(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta])
x1<-(sort(alfa_value,decreasing=TRUE)[((length(B1)+1) %/% consta)+1])
fact90<-((((1-(1/consta))-y1)*(x2-x1)/(y2-y1))+x1)/sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta]

consta<-5 # 80% conf.
y2<-(length(B1)- (((length(B1)+1) %/% consta)-1))/(length(B1)+1)
y1<-(length(B1)- (((length(B1)+1) %/% consta)))/(length(B1)+1)
x2<-(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta])
x1<-(sort(alfa_value,decreasing=TRUE)[((length(B1)+1) %/% consta)+1])
fact80<-((((1-(1/consta))-y1)*(x2-x1)/(y2-y1))+x1)/sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% consta]

fact99
fact95
fact90
fact80


#### Otherwise, If your training data set has the number of training samples as 99 + whole number x 100  then:
fact99<-1
fact95<-1
fact90<-1
fact80<-1


### Not necessary part, but for training set prediction intervals if one is interested too
## training set prediction intervals
telta_half_width_99perc<- Error_function*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 100])*fact99 
telta_half_width_95perc<- Error_function*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 20])*fact95
telta_half_width_90perc<- Error_function*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 10])*fact90
telta_half_width_80perc<- Error_function*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 5])*fact80
## training set prediction intervals



## test set prediction intervals
pred_y_half_width_99perc<- error_te*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 100])*fact99
pred_y_half_width_95perc<- error_te*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 20])*fact95
pred_y_half_width_90perc<- error_te*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 10])*fact90
pred_y_half_width_80perc<- error_te*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 5])*fact80
## test set prediction intervals


## database prediction intervals
pred_y_half_width_99perc_database<- error_te_database*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 100])*fact99
pred_y_half_width_95perc_database<- error_te_database*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 20])*fact95
pred_y_half_width_90perc_database<- error_te_database*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 10])*fact90
pred_y_half_width_80perc_database<- error_te_database*(sort(alfa_value,decreasing=TRUE)[(length(B1)+1) %/% 5])*fact80
## database prediction intervals


# means and medians prediction intervals of your test set
mean(pred_y_half_width_99perc) 
mean(pred_y_half_width_95perc)
mean(pred_y_half_width_90perc)
mean(pred_y_half_width_80perc)
median(pred_y_half_width_99perc)
median(pred_y_half_width_95perc)
median(pred_y_half_width_90perc)
median(pred_y_half_width_80perc)
# means and medians prediction intervals of your test set


# Error rate analysis of your test set
stat_pred_y_99<-0
vect_pred_y_99<-array(rep(0),length(Btest))
vect_pred_y_95<-array(rep(0),length(Btest))
vect_pred_y_90<-array(rep(0),length(Btest))
vect_pred_y_80<-array(rep(0),length(Btest))
stat_pred_y_95<-0
stat_pred_y_90<-0
stat_pred_y_80<-0
for (i in 1:length(Btest)) {
if (abs(pred_y[i]-Btest[i])>pred_y_half_width_99perc[i]) {
stat_pred_y_99<-stat_pred_y_99+1
vect_pred_y_99[i]<-1}
if (abs(pred_y[i]-Btest[i])>pred_y_half_width_95perc[i]) {
stat_pred_y_95<-stat_pred_y_95+1
vect_pred_y_95[i]<-1}
if (abs(pred_y[i]-Btest[i])>pred_y_half_width_90perc[i]) {
stat_pred_y_90<-stat_pred_y_90+1
vect_pred_y_90[i]<-1}
if (abs(pred_y[i]-Btest[i])>pred_y_half_width_80perc[i]) {
stat_pred_y_80<-stat_pred_y_80+1
vect_pred_y_80[i]<-1}
}

stat_pred_y_99
stat_pred_y_95
stat_pred_y_90
stat_pred_y_80

stat_pred_y_99/length(Btest)
stat_pred_y_95/length(Btest)
stat_pred_y_90/length(Btest)
stat_pred_y_80/length(Btest)
# Error rate analysis of your test set



# means and medians prediction intervals of your database
mean(pred_y_half_width_99perc_database) 
mean(pred_y_half_width_95perc_database)
mean(pred_y_half_width_90perc_database)
mean(pred_y_half_width_80perc_database)
median(pred_y_half_width_99perc_database)
median(pred_y_half_width_95perc_database)
median(pred_y_half_width_90perc_database)
median(pred_y_half_width_80perc_database)
# means and medians prediction intervals of your database


