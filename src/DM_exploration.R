###1. read data and some cleaning
df_otu <- read.csv("otutable.txt", sep="\t", row.names=1)
df_meta <- read.csv("task.txt", sep="\t")
# filter unmatched samples (some baseline control)
keep_otu <- intersect(colnames(df_otu), df_meta$X.SampleID)
df_otu <- df_otu[ ,keep_otu]
df_meta <- df_meta[df_meta$X.SampleID %in% keep_otu, c(1,2)]
rm(keep_otu)

###2. clean data (prepare into analysis-ready format)
row.names(df_meta) <- df_meta[,1]
df_meta <- df_meta[, -1, drop=F] 
colnames(df_meta) <- "label"
original_otu <- row.names(df_otu) 
fuso_otu <- original_otu[grep("Fusobacterium", original_otu)]
fuso_otu_id <- unname(sapply(fuso_otu, function(x){paste(unlist(strsplit(x, split="_"))[1:2], collapse = "_")}))
row.names(df_otu) <- sapply(row.names(df_otu), function(x){paste(unlist(strsplit(x, split="_"))[1:2], collapse = "_")})

otu_count <- t(df_otu)
otu_count <- merge(df_meta, otu_count, by=0)
row.names(otu_count) <- otu_count$Row.names
otu_count <- otu_count[, -1]


###3. feature selection
# 3.1 sample size (rowSum): mean 4501, median 4297, max 11339, min(237)
# indicating some samples are extremely sparse, can't be used
# arbitrary cutoff: >2000 (129 / 156 remaining)
par(mfrow=c(1,2))
plot(density(rowSums(otu_count[, -1])), main="density of library size", col="orange")
abline(v=2000, col="red")
boxplot(rowSums(otu_count[, -1]), main="Boxplot of library size", col="grey")
abline(h=2000, col="red")
f.otu.count <- otu_count[rowSums(otu_count[, -1])>2000, ]

# 3.2 remove low abundance species (biological reasoning)
plot(density(colSums(f.otu.count[,-1])))
quantile(colSums(f.otu.count[,-1]))
### so many variables appear very few times
# 0~1: not reliable, itu is not sensitive at low-abundance level
# 1~20: the appearance is too low to be biological informative (0.5% total popu)
#f.otu.count <- f.otu.count[, c(TRUE, colSums(f.otu.count[,-1])>20)]
f.otu.count <- f.otu.count[, c(TRUE, colSums(f.otu.count[,-1]>5)>5)]

#3.3 remove rerundant feature: cutoff 0.9
library(caret)
library(mlbench)
correlationMatrix <- cor(f.otu.count[,-1])
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9)
sum(colnames(f.otu.count[,-1]) == colnames(correlationMatrix)) #make sure they are in same order
f.otu.count <- f.otu.count[, -c(highlyCorrelated+1)] # remove those with >=0.9 correlation
### Important: though those features are not helpful in prediction model. 
### They are very useful in functional analysis because co-existing (+0.9) or compitition (-0.9) are important features to know in biology 
### May even contribute to the treatment
length(fuso_otu)
fuso202 <- intersect(colnames(f.otu.count), fuso_otu_id)
# after filtering 202 variables remaining
# 7 of them are "fuso" genus, consistent with the 

# confirm the paper finding that they are enriched:
aggregate(.~label, f.otu.count[,c("label",fuso202)], mean)
all_dif <- aggregate(.~label, f.otu.count, mean)
exp_dif <- abs(all_dif[1, -1] - all_dif[2, -1])
exp_dif <- sort(exp_dif, decreasing = T)

# corrplot
library(corrplot)
corrplot(cor(f.otu.count[,-1]), method="circle", title="Corrplot of features", tl.pos = "n", mar=c(0,0,1,0))

# pca plot
library(ggfortify)
df <- f.otu.count[, -1]
pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res, data = f.otu.count, colour = 'label')



############ DM analysis
#4.1 SVM by LIBSVM
library(e1071) 
run_svm <- function(input_df, input_kernel="sigmoid", input_seed=100) {
  set.seed(input_seed)
  index <- sample(2,nrow(input_df),replace = TRUE,prob=c(0.7,0.3))
  traindata <- input_df[index==1,]
  testdata <- input_df[index==2,]
  temp_svm <- svm(label~., data=traindata, kernel=input_kernel, probability=T)
  
  pred1 <- predict(temp_svm, traindata[,-1])
  table1 <- table(pred=pred1, true=traindata[,1])
  accuracy1 <- sum(pred1==traindata[,1])/length(pred1)
  pred2 <- predict(temp_svm, testdata[,-1])
  table2 <- table(pred=pred2, true=testdata[,1])
  accuracy2 <- sum(pred2==testdata[,1])/length(pred2)
  print("Prediction for training data")
  print(accuracy1)
  print(table1)
  print("Prediction for testing data")
  print(accuracy2)
  print(table2)
  
  # data for ROC
  x.svm.prob <- predict(temp_svm, type="prob", newdata=testdata, probability = TRUE)
  x.svm.prob.rocr <- prediction(attr(x.svm.prob, "probabilities")[,1], testdata$label)
  x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
  return(x.svm.perf)
}
#plot(temp_svm, testdata, NR_041364.1~NR_028961.1) #only can plot 2D

#4.1.1 compare kernel
a1<-run_svm(f.otu.count, input_kernel="linear")
a2<-run_svm(f.otu.count, input_kernel="polynomial")
a3<-run_svm(f.otu.count, input_kernel="radial")
a4<-run_svm(f.otu.count, input_kernel="sigmoid")
rm(a1,a2,a3,a4)

#4.1.2 var-selection and AUC
### selection by RF
### ref: http://r-statistics.co/Variable-Selection-and-Importance-With-R.html
library(party)
cf1 <- cforest(label ~ . , data= f.otu.count, control=cforest_unbiased(mtry=2,ntree=50)) 
fss <- sort(varimp(cf1), decreasing = T)

m5_otu <- f.otu.count[,c("label",names(fss[1:5]))]
m10_otu <- f.otu.count[,c("label",names(fss[1:10]))]
m20_otu <- f.otu.count[,c("label",names(fss[1:20]))]
m50_otu <- f.otu.count[,c("label",names(fss[1:50]))]

a5_rad <- run_svm(m5_otu, input_kernel="radial")
a10_rad <- run_svm(m10_otu, input_kernel="radial")
a20_rad <- run_svm(m20_otu, input_kernel="radial")
a50_rad <- run_svm(m50_otu, input_kernel="radial")
afull_rad <- run_svm(f.otu.count, input_kernel="radial")
a5_sig <- run_svm(m5_otu, input_kernel="sigmoid")
a10_sig <- run_svm(m10_otu, input_kernel="sigmoid")
a20_sig <- run_svm(m20_otu, input_kernel="sigmoid")
a50_sig <- run_svm(m50_otu, input_kernel="sigmoid")
afull_sig <- run_svm(f.otu.count, input_kernel="sigmoid")

# ROC
plot(a5_rad, col="orange", main="ROC curves of SVM")
plot(a10_rad, col="orange", add=TRUE)
plot(a20_rad, col="orange", add=TRUE)
plot(a50_rad, col="orange", add=TRUE)
plot(afull_rad, col="orange", add=TRUE)
plot(a5_sig, col="red", add=TRUE)
plot(a10_sig, col="red", add=TRUE)
plot(a20_sig, col="red", add=TRUE)
plot(a50_sig, col="red", add=TRUE)
plot(afull_sig, col="red", add=TRUE)
legend(0.5, 0.4, c("radial", "sigmoid"), c("orange","red"))

### 4.1.3 5-folder cross validation by SVM with "sigmoid" kernel
k_cross <- function(input_data,k=5, kernel, seed=123) {
  set.seed(seed)
  all_index <- sample(seq(1:dim(f.otu.count)[1]))
  index_cutoff <- seq(1, dim(input_data)[1], length.out = k+1)
  accu_train=c()
  accu_test=c()
  for (i in 1:k) {
    testdata <- input_data[all_index[floor(index_cutoff[i]):floor(index_cutoff[i+1])], ]
    traindata <- input_data[-all_index[floor(index_cutoff[i]):floor(index_cutoff[i+1])], ]
    temp_svm <- svm(label~., data=traindata, kernel=kernel)
    pred1 <- predict(temp_svm, traindata[,-1])
    accuracy1 <- sum(pred1==traindata[,1])/length(pred1)
    pred2 <- predict(temp_svm, testdata[,-1])
    accuracy2 <- sum(pred2==testdata[,1])/length(pred2)
    accu_train=c(accu_train, accuracy1)
    accu_test=c(accu_test, accuracy2)
  }
  print("Mean accuracy in training data:")
  print(mean(accu_train))
  print("mean accuracy in test data:")
  print(mean(accu_test))
}

k_cross(m5_otu, kernel="sigmoid")
k_cross(m5_otu, kernel="radial")

k_cross(m10_otu, kernel="sigmoid")
k_cross(m10_otu, kernel="radial", seed=999)

k_cross(m20_otu, kernel="sigmoid")
k_cross(m20_otu, kernel="radial", seed=999)

k_cross(m50_otu, kernel="sigmoid")
k_cross(m50_otu, kernel="radial")

k_cross(f.otu.count, kernel="sigmoid")
k_cross(f.otu.count, kernel="radial")

# plot
n_var=c(5,10,20,50,202)
rad_train = c(0.76,0.82,0.86,0.82,0.92)
rad_test = c(0.70,0.78,0.70,0.70,0.63)
sig_train = c(0.75,0.75,0.77,0.88,0.84)
sig_test = c(0.71, 0.74, 0.73, 0.71, 0.65)
plot(n_var, rad_train, type = "n", ylim = c(0.5, 1), xlim = c(0, 210), main="5-CV mean accuracy", ylab="accuracy", xlab="number of variables")
lines(n_var, rad_train, col = "red")
lines(n_var, rad_test, col="red")
lines(n_var, sig_test, col="black")
lines(n_var, sig_train, col="black")
text(locator(), labels = c("rad_train", "rad_test", "sig_train", "sig_test"), col=c("red","red","black", "black"))


### 4.2 naive bayesian
k_cross_bayes <- function(input_data,k=5, kernel, seed=123) {
  set.seed(seed)
  all_index <- sample(seq(1:dim(f.otu.count)[1]))
  index_cutoff <- seq(1, dim(input_data)[1], length.out = k+1)
  accu_train=c()
  accu_test=c()
  for (i in 1:k) {
    testdata <- input_data[all_index[floor(index_cutoff[i]):floor(index_cutoff[i+1])], ]
    traindata <- input_data[-all_index[floor(index_cutoff[i]):floor(index_cutoff[i+1])], ]
    temp_bayes <- naiveBayes(label~., data=traindata, laplace = 3)
    pred1 <- predict(temp_bayes, traindata[,-1])
    accuracy1 <- sum(pred1==traindata[,1])/length(pred1)
    pred2 <- predict(temp_bayes, testdata[,-1])
    accuracy2 <- sum(pred2==testdata[,1])/length(pred2)
    accu_train=c(accu_train, accuracy1)
    accu_test=c(accu_test, accuracy2)
  }
  print("Mean accuracy in training data:")
  print(mean(accu_train))
  print("mean accuracy in test data:")
  print(mean(accu_test))
}

m30_otu <- f.otu.count[,c("label",names(fss[1:30]))]
m40_otu <- f.otu.count[,c("label",names(fss[1:40]))]
m60_otu <- f.otu.count[,c("label",names(fss[1:60]))]
m70_otu <- f.otu.count[,c("label",names(fss[1:70]))]
m80_otu <- f.otu.count[,c("label",names(fss[1:80]))]
m80_otu <- f.otu.count[,c("label",names(fss[1:80]))]
m100_otu <- f.otu.count[,c("label",names(fss[1:100]))]

k_cross_bayes(m10_otu)
k_cross_bayes(m20_otu)
k_cross_bayes(m30_otu)
k_cross_bayes(m40_otu)
k_cross_bayes(m50_otu)
k_cross_bayes(m60_otu)
k_cross_bayes(m70_otu)
k_cross_bayes(m80_otu)
k_cross_bayes(m90_otu)
k_cross_bayes(m100_otu)
k_cross_bayes(f.otu.count)

n_var = c(10,20,30,40,50,60,70,100,202)
accu_train = c(0.71,0.75,0.74,0.75,0.78,0.77,0.79,0.81,0.87)
accu_test = c(0.67,0.68,0.68,0.67,0.67,0.65,0.64,0.68,0.64)

plot(n_var, accu_train,type="l", col="red", ylim=c(0.5,1), ylab="average accuracy", main="5-VC of Naive Bayesian")
lines(n_var, accu_test, type="l", col="black")
text(locator(), labels = c("Bayes_train", "Bayes_test"), col=c("red","black"))

### 4.3 Neural network
library("neuralnet")
k_cross_nn <- function(input_data,k=5, hidden=3, seed=123) {
  set.seed(seed)
  all_index <- sample(seq(1:dim(f.otu.count)[1]))
  index_cutoff <- seq(1, dim(input_data)[1], length.out = k+1)
  accu_train=c()
  accu_test=c()
  for (i in 1:k) {
    testdata <- input_data[all_index[floor(index_cutoff[i]):floor(index_cutoff[i+1])], ]
    traindata <- input_data[-all_index[floor(index_cutoff[i]):floor(index_cutoff[i+1])], ]
    temp_nn <- neuralnet(label~.,data=traindata, hidden=hidden,act.fct = "logistic", linear.output = FALSE)
    pred1 <- predict(temp_nn, traindata[,-1])
    pred1 <- ifelse(pred1[,2]>=0.5, "Tumor", "Healthy")
    accuracy1 <- sum(pred1==traindata[,1])/length(pred1)
    pred2 <- predict(temp_nn, testdata[,-1])
    pred2 <- ifelse(pred2[,2]>=0.5, "Tumor", "Healthy")
    accuracy2 <- sum(pred2==testdata[,1])/length(pred2)
    accu_train=c(accu_train, accuracy1)
    accu_test=c(accu_test, accuracy2)
  }
  print("Mean accuracy in training data:")
  print(mean(accu_train))
  print("mean accuracy in test data:")
  print(mean(accu_test))
}


### random forest
library(randomForest)
test_rf <- function(maxnode){
  temp_rf <- randomForest(label~., data=traindata, maxnode=maxnode)
  pred1 <- predict(temp_rf, traindata[,-1])
  accuracy1 <- sum(pred1==traindata[,1])/length(pred1)
  pred2 <- predict(temp_rf, testdata[,-1])
  accuracy2 <- sum(pred2==testdata[,1])/length(pred2)
  print(accuracy1)
  print(accuracy2)
}



########### Archive: code backup
## Bias-reduction logistic: to handle linear separatility
# fit bias-reduction logit
# https://stats.stackexchange.com/questions/11109/how-to-deal-with-perfect-separation-in-logistic-regression
library(arm)
new_logistic <- bayesglm(label~., data=traindata, family="binomial")
pred1 <- predict(new_logistic, traindata[, -1])
optCutOff=optimalCutoff(traindata$label, pred1)

