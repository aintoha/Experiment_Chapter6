source("functions_file1.R")
library(MASS)
library(cvTools)
attach(Boston)
mydata <- as.data.frame(Boston[,13])
#############################################
##############################
#grid search
#bw: bandwidth
#loss_in: empircal loss in the in-sample
#est_density_in: estimated density in the in-sample
#loss: empirical loss in out-of-sample
#est_density_out: estimate density in the out-of-sample
#se_loss_in: 
set.seed(5202)
bw  <- 2^seq(from=-4, to= 4, by=1)
length_bandwidth <- length(bw)
numouterfolds <- 3
folds <- cvFolds(NROW(mydata), K=numouterfolds)
loss_in <- matrix(NA, numouterfolds, length_bandwidth)
loss_out <- rep(NA, numouterfolds)
se_loss_out <- rep(NA,numouterfolds)

# others_data <- rep(NA, numouterfolds)
# test_data <- rep(NA, numouterfolds)
# 
# for(outerfolds in 1:numouterfolds){
#   
#   others_data[outerfolds] <- list(splitfolder(mydata, folds, numberoffolds = outerfolds)[[1]]) #select the other
#   test_data[outerfolds] <- list(splitfolder(mydata, folds, numberoffolds = outerfolds)[[2]]) # select the test set
#   
# }  
# others_data  
# test_data

set.seed(5202)
for(outerfolds in 1:numouterfolds){
  others_data <- splitfolder(mydata, folds, numberoffolds = outerfolds)[[1]] #select the other
  test_data <- splitfolder(mydata, folds, numberoffolds = outerfolds)[[2]] # select the test set
  for(h in 1:length(bw)){
    loss_in[outerfolds,h] <- unlist(V.logloss(logh=log(bw[h]), obsdata=others_data, testdata= others_data))
  }
}

loss_in
#mean loss for each inner folds
mean_loss_in <- apply(loss_in, 2,mean)
mean_loss_in
#the parameter is performed 5 times and the models are selected
minh <- which(mean_loss_in == min(mean_loss_in), arr.ind = T)

#model evaluation
for(outerfolds in 1:numouterfolds){
  
  others_data <- splitfolder(mydata, folds, numberoffolds = outerfolds)[[1]]
  test_data <- splitfolder(mydata, folds, numberoffolds = outerfolds)[[2]]
  
  loss_out[outerfolds] <- unlist(V.logloss(logh=log(bw[minh]), obsdata=others_data, testdata=test_data))
  se_loss_out[outerfolds] <- se.logloss(data=others_data, testdata=test_data, logbw=log(bw[minh]))
  
}
mean(loss_out)
mean(se_loss_out)

plot(log(bw), mean_loss_in, type="l", main="", 
     xlab= "log bandwidth", ylab="Loss", ylim=c(2.5, 4))
legend("bottomright", legend=c("Lower bound mean", "Loss", "Upper bound mean"),
       col=c("red","black", "blue"), lty=2:1, cex=0.7)
dev.off()
