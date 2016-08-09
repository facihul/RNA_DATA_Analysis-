##  This module explains noncoding RNA_Seq data clasification. 
# It classifies multiple classes (4 stages of breast cancer data) based on patient pathological stages. 
# data were prepare before for this module. 


library(e1071)

#data <- read.delim(file.choose(), sep = '\t') 
load("/Users/mdfacihulazam/Desktop/ref/final_data/Non_Coding_RNA.RData")
data1<- t(nonCodingRNA)

labels <- c(rep(0,106),rep(1,393),rep(2,290),rep(3,171)) # Class Labels
#labels <- c(rep(0,107),rep(1,395),rep(2,291),rep(3,173))
fold <- 5 # 5 fold cross validation 
ste <- 960/fold  # each fold size of data 
indx <- c(1:length(labels))
#indx <- row.names(data1)
j<-1
acc <- c()
#sample_test <- c()
count <- 1
for(i in 1:fold){
  start <- count
  count <- i*ste
  sample_test <- indx[start:count]   # test index
  sample_train <- indx[-c(start:count)]  # training index
  
  train_in <- data1[sample_train,] # training sample 
  train_out <- labels[sample_train] # training  sample label
  
  test_in <- data1[sample_test,] # test sample
  test_out <- labels[sample_test] # test sample
  
  model <- e1071::svm(train_in, train_out,kernel="radial") # SVM model generate 
  predsvm1 <- as.character(predict(model,test_in)) # class prediction 
  tab <- table(predsvm1,test_out) # confusion matrix 
  acc[j] <- (tab[1,1]/(ste+1))*100 
  j<-j+1
} 
acc_mean <- mean(acc) # mean of accuracy
Sd_acc = sd(acc)/sqrt(length(acc)) # mean Standered error




