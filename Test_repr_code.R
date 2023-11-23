#AFTER you complete Info section 1 and obtain two files:
#one file with 1444 Padel descriptors as "file_with_1444_Padel_descriptors.csv"
#and the other file with 17 rdkit descriptors as "file_with_17_rdkit_descriptors.csv"

# THEN: Open R in the same directory of that files:

# V<-read.csv("file_with_17_rdkit_descriptors.csv",sep=',',header=TRUE) # file with 17 rdkit #descriptors, check the separation symbol and whether the file has a header, also check whether the #first column is the first descriptors, if it isnt, but if it is a simple column with row numbers, that #column must be removed, so that the first column is the first molecular descriptor

# V1<-read.csv('file_with_1444_Padel_descriptors.csv',sep=',',header=TRUE) # file with 1444 #Padel descriptors, check the separation symbol and whether the file has a header, also check #whether the first column is the first descriptors, if it isn't, but if it is a simple column with row #numbers, that column must be removed, so that the first column is the first molecular descriptor

#Fe1<-cbind(V,V1) 

#dim(Fe1) ### the second dimension must be 1461, as 17+1444 = 1461 mol. descriptors

#Just before we continue, These V and V1 are YOUR OWN FILES IF YOU CREATED THEM. #But we already have our own test file with 1461 such descriptors, so you can use our test file to #test our results. 

#In that case you can also load our file as:

Fe1<-read.csv("AqSolDB-n_test_set_mol_descr.csv",sep=",",header=TRUE)
Fe1<-Fe1[,-1] #Removing non-descriptor column!
dim(Fe1) # checking dimensions
# [1]  220 1461 dimension of 1461 descriptors is OK, we can continue, otherwise we would NOT #be able to continue

Fe1<-as.matrix(Fe1)

F<-Fe1[,c("MolLogP","XLogP","ATS0p","ZMIC1","piPC2","MDEC.33","MPC7","piPC1","piPC3","ZMIC2","piPC6","AATS5p","MolMR","TpiPC","AATS6v","piPC4","AATS1i","piPC10","GATS2c","AATS1e","AATS4v","TWC","nH","ATS1m","Mi","MPC8")] #Selecting 26 AqSolDB-n descrip.


### Loading mean vector and standard deviation vector when scaling the descriptor data: ####

meanvect<-read.csv("meanvect.csv",header=FALSE)
meanvect<-as.vector(t(meanvect))
sdvect<-read.csv("sdvect.csv",header=FALSE)
sdvect<-as.vector(t(sdvect))

### 

for (i in 1:ncol(F)) {
F[,i]<- (F[,i]- meanvect[i])/sdvect[i] }

Y<-read.csv("AqSolDB-n_test_set_LogS.csv",header=TRUE) # Loading
B2<-as.vector(t(Y))

library(xgboost) #Loading XGB library
bst<-xgb.load("FSTI-XGB_1619") #Loading our AqSolDB-N model

pred_y<-predict(bst,F)

library(caret) # Loading caret library
caret::RMSE(B2,pred_y)

# [1] 0.5943423 

#Bravo, you correctly reproduced our result of our internal independent test set (Table 4), last row #RMSEV value.

