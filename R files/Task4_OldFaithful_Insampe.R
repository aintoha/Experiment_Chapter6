#In this script, we will compare two methods. 
#One method is using grid search
#the other method is a combination of grid search + gradient descent
#One -  which is fastest?
#Second - which give the smallest error 


source("functions_file1.R")
library(MASS)
library(cvTools)
library(mlr3)
library(mlr3proba)
data <- read.table("old_faithful.txt")
mydata <- as.data.frame(data$V2)
#############################################
##############################
#grid search,
#bw: bandwidth
#loss_in: empircal loss in the in-sample
#est_density_in: estimated density in the in-sample
#loss: empirical loss in out-of-sample
#est_density_out: estimate density in the out-of-sample
#se_loss_in: 
set.seed(5202)
bw  <- c(0.0001, 0.0100000, 0.3252632, 0.6405263, 0.9557895, 1.2710526, 1.5863158, 1.9015789, 2.2168421, 2.5321053, 2.8473684, 3.1626316,
         3.4778947, 3.7931579, 4.1084211, 4.4236842, 4.7389474, 5.0542105, 5.3694737, 5.6847368, 6.0000000)
length_bandwidth <- length(bw)
numouterfolds <- 3
folds <- cvFolds(NROW(mydata), K=numouterfolds)
loss_in <- matrix(NA, numouterfolds, length_bandwidth)
loss_out <- rep(NA, numouterfolds)
se_loss_out <- rep(NA,numouterfolds)
minh <- rep(NA, numouterfolds)

set.seed(5202)
for(outerfolds in 1:numouterfolds){
  others_data <- splitfolder(mydata, folds, numberoffolds = outerfolds)[[1]] #select the other
  test_data <- splitfolder(mydata, folds, numberoffolds = outerfolds)[[2]] # select the test set
  for(h in 1:length(bw)){
    loss_in[outerfolds,h] <- unlist(V.logloss(logh=log(bw[h]), obsdata=others_data, testdata= others_data))
  }
  minh[outerfolds] = bw[which.min(loss_in)]
  loss_out[outerfolds] =  unlist(V.logloss(logh=log(minh[outerfolds]), obsdata=others_data, testdata=test_data))
}

minh
loss_out
loss_in


png(filename= "task2_insample.png", width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(bw), loss_in[1,], type = "b", lwd = 2, pch = 16, 
     xlab  = "log bandwidth", ylab = "Out-of-sample empirical loss",main = "Task 1")
lines(log(bw), loss_in[1,], type = "b", col = "red", lwd = 2, lty = 4, pch = 18)
lines(log(bw), loss_in[1,], type = "b", col = "blue", lwd = 2, lty = 3, pch = 17)
legend("topright", 
       legend = c("Fold 1", "Fold 2", "Fold 3"), 
       col = c("black", "red", "blue"), 
       pch = c(16,18,17), 
       lty = c(1, 4, 3),
       lwd = 2, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
dev.off()

