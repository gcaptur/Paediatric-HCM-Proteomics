# Paediatric-HCM-Proteomics

# This code was used to analyse the plasma proteomics data reported in the paper: 

# "A plasma proteomic biomarker panel can identify hypertrophic cardiomyopathy and predict risk of sudden cardiac death in children"

# Normality tests (exemplar analyte code)

a5<- read.csv("~/Directory/a5.csv")				

plot(density(a5$ALDOA));
shapiro.test(a5$ALDOA);
qqnorm(a5$ALDOA);
qqline(a5$ALDOA, col = 2);

# Parametric and non-parametric difference  (exemplar analyte code)


a1<- read.csv("~/Directory/a1.csv")
a1$GROUP <- as.factor(a1$status)
GROUP <-a1$GROUP 
table(GROUP)

t.test(a1$ALDOA~a1$GROUP, data=a1, p.adjust.method = "holm")
wilcox.test(a1$ALDOA~a1$GROUP, data=a1, p.adjust.method = "holm")
describeBy(a1$A2AP, a1$GROUP)

res.aovIGHG1 <- aov(a1$IGHG1 ~a1$GROUP, data = a1)
summary(res.aovIGHG1)
TukeyHSD(res.aovIGHG1)

# Diagnostic panel by SVM

# Training set 
data <- read.csv("~/Directory/data.csv")
head(data,5)
data$status <- as.factor(data$status)
status <- data$status
table(data$status)
prop.table(table(data$status))
rose <- ROSE(status~., data = data, seed = 1)$data
table(rose$status)
prop.table(table(rose$status))

# Oversample 
rose <- ovun.sample(status~., data = data, method = "over", N = 300)$data 
table(rose$status)
prop.table(table(rose$status))

# Over-undersample
rose <- ovun.sample(status~., data = data, method = "both", p=0.5, N=1000, seed = 1)$data
table(rose$status)
prop.table(table(rose$status))

index <- c(1:nrow(rose))
test.index <- sample(index, size = (length(index)/3))
training <- rose[-test.index ,]
test <- rose[test.index ,]
prop.table(table(training$status))
prop.table(table(test$status))

x <- subset(training, select=-status)
y <- training$status
xx <- subset(test, select=-status)
yy <- test$status
table(training$status)
prop.table(table(training$status))

svm_model1 <- svm(x,y)

summary(svm_model1)
pred <- predict(svm_model1,x, decision.values = TRUE)
table(pred,y)
svm_tune <- tune(svm, train.x=x, train.y=y, 
              kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))
print(svm_tune)

svm_model_after_tune <- svm(x,y, kernel="radial", cost=10, gamma=2, decision.values = TRUE)
summary(svm_model_after_tune)

pred <- predict(svm_model_after_tune,x, decision.values = TRUE) 
table(pred,y)
summary(pred)
pred
y

pred2<-prediction(as.numeric(pred), as.numeric(training$status))
curve<-performance(pred2,"auc")
rates <- performance(pred2, "tpr", "fpr") 
summary(pred2)
summary(rates)
pred2
plot(rates)
print(str(curve))
print(svm_model_after_tune)

rates2 <- performance(pred2, measure = "tpr", x.measure = "fpr")
plot(rates2)
abline(a=0, b= 1) 

perf1 <- performance(pred2, "prec", "rec")
plot(perf1)

perf2 <- performance(pred2, "sens", "spec")
plot(perf2)

# Validation set 

best.svm.pred <- predict(svm_model_after_tune,xx, decision.values = TRUE)
table(best.svm.pred,yy)
table(Prediction = best.svm.pred, Truth = test$status)

pred3<-prediction(as.numeric(best.svm.pred), as.numeric(test$status))
curve<-performance(pred3,"auc")
rates3 <- performance(pred3, "tpr", "fpr") 
summary(pred3)
best.svm.pred

yy 

plot(rates3)
print(str(curve))
print(svm_model_after_tune)

rates3 <- performance(pred3, measure = "tpr", x.measure = "fpr")
plot(rates3)
abline(a=0, b= 1) 

perf3 <- performance(pred3, "prec", "rec")
plot(perf3)
print(perf3)
summary(perf3)

# Boxplot of prediction biomarker scores for AUC

a6 <- read.csv("~/Directory/training.csv")

par(mar=c(8,8,4,2)) 
boxplot(a6$value ~ a6$group, data = a6, names = c("Control","HCM"), col = c("#0000ff22","royalblue2"),
ylab = "Biomarker Score", cex.axis=1.5, cex.lab=1.5, outline = FALSE)

a6 <- read.csv("~/Directory/validation.csv")

par(mar=c(8,8,4,2)) 
boxplot(a6$value ~ a6$group, data = a6, names = c("Control","HCM"), col = c("#0000ff22","royalblue2"),
ylab = "Biomarker Score", cex.axis=1.5, cex.lab=1.5, outline = FALSE)

# Boxplot of analyte differences HCM vs control

a6 <- read.csv("~/Directory/a6.csv")

par(mar=c(8,8,4,2)) 
boxplot(a6$TLN1 ~ a6$status, data = a6, names = c("Control","HCM"), col = c("#0000ff22","royalblue2"),
ylab = "GST01 (ppm)", cex.axis=1.5, cex.lab=1.5, outline = FALSE)


