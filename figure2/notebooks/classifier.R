library(e1071)
library(reshape2)
library(MASS)
a<-read.table("../raw/data", check.names = F, sep = "\t", header = T)
a<-read.table("../raw/data", sep = "\t", header = T)
a1<-a[,c(4:dim(a)[2])] / a[,3]
b1<-factor(a$Pathology)
a2<-cbind(a1, b1)
colnames(a2)[dim(a2)[2]]<-"pathol"
###############################################################

# a2 is the input dataframe
# a1 (the first part of the a2): gene expression value
# b1 (the last column of the a2): pathology annotation
# I tried both whole genes and top1000 highly variable gene (from Seurat)
# We used the top1000 since the performance was the best

###############################################################
total_res<-c()
###svm
res<-svm(pathol~., a2)
pred<-predict(res, a1)
a3<-a$Pathology == pred
res1 <- length(which(a3)) / length(pred)
total_res<-c(total_res, res1)
names(total_res)[1] <- "svm"
###lda
res<-lda(pathol~., a2, tol=1e-20)
pred<-predict(res, a1)
a3<-a$Pathology == pred$class
res1 <- length(which(a3)) / length(pred$class)
total_res<-c(total_res, res1)
names(total_res)[2] <- "lda"
###
#multinomial logistic regression
library(nnet)
res<-multinom(pathol~., a2)
pred<-predict(res, a1)
a3<-a$Pathology == pred
res1 <- length(which(a3)) / length(pred)
total_res<-c(total_res, res1)
names(total_res)[3] <- "multi_logistic_reg"
#################
## decision tree
library(rpart)
res<-rpart(pathol~., a2)
pred<-predict(res, a1)
m<-apply(pred, 1, function(x) which(x == max(x)))
m1<-colnames(pred)[m]
a3<-a$Pathology == m1
res1 <- length(which(a3)) / length(pred)
total_res<-c(total_res, res1)
names(total_res)[4] <- "decision_tree"
###naive bayesian
res<-naiveBayes(pathol~., a2)
pred<-predict(res, a1)
a3<-a$Pathology == pred
res1 <- length(which(a3)) / length(pred)
total_res<-c(total_res, res1)
names(total_res)[5] <- "naive_bayesian"
###random forest
res<-train(pathol~., a2, method = 'rf')
pred<-predict(res, a1)
a3<-a$Pathology == pred
res1 <- length(which(a3)) / length(pred)
total_res<-c(total_res, res1)
names(total_res)[6] <- "random_forest"
write.table(total_res, "self_validation_pathol_celltype.csv", sep = ",", quote = F) (edited) 
