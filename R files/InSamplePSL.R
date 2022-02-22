#Simulation experiment for out-of-sample


library(mlr3proba)
library(distr6)
library(mlr3tuning)
library(paradox)
library(mlr3)
library(future)
library(MASS)
library(data.table)

#load data
#----------
set.seed(5202)
Data1 <- data.frame("A" = rnorm(200, mean = 2, sd = 1))
Data2 <- data.frame("B" = c(rnorm(100, mean = 1, sd = 1), rnorm(100, mean = 10, sd = 2)))
oldfaithful <- read.table("old_faithful.txt")
attach(Boston)
Boston <- Boston
auto <- read.table("auto.Data", header = T)
head(auto)
energy <- read.table("energyefficiency.Data", header = T)

#set Task 
#------------
Task1 = TaskDens$new(id = "A", backend = Data1$A)
Task2 = TaskDens$new(id = "B", backend = Data2$B)
Task3 = TaskDens$new(id = "boston", backend = Boston$medv)
Task4 = TaskDens$new(id = "old", backend = oldfaithful$V2)
Task5 = TaskDens$new(id = "energy", backend = energy$Y2)
Task6 = TaskDens$new(id = "auto", backend = auto[,1] )

#set learrner 
#--------------
lrn = lrn("dens.kde", kernel = "Norm")

#set tuning parameter
#------------------------
ps = ParamSet$new(list(
  ParamDbl$new("bandwidth", lower = 0.01, upper = 6)
))
tuner = tnr("grid_search", resolution = 20)

#define the resampling 
#----------------------------
cv1 = rsmp("insample") #inner
cv2 = rsmp("cv", folds = 3) #outer

set.seed(5202)
#set tuning learner 
NormLrn_task1 = AutoTuner$new(
  learner = lrn,
  resampling = cv1,
  measure = msr("dens.squared"),
  search_space = ps,
  terminator = trm("evals", n_evals = 20),
  tuner = tuner,
  store_tuning_instance = TRUE
)
NormLrn_task1$train(Task1)$tuning_result
NormLrn_task1$tuning_result

#d. set up benchmark
#----------------------
set.seed(5202)
grid_task1 = resample(
  task = Task1,
  learner = NormLrn_task1,
  resampling = cv2, 
  store_models=TRUE
)
grid_task1$aggregate()
lapply(grid_task1$learners, function(x) x$tuning_result)
rep_task1 = lapply(grid_task1$learners, function(x) x$model)

#extract from first fold
dat1_task1 = data.table(rep_task1[[1]]$tuning_instance$archive$data())
dat1_task1 = dat1_task1[order(dat1_task1$bandwidth)]
dat2_task1 = data.table(rep_task1[[2]]$tuning_instance$archive$data())
dat2_task1 = dat2_task1[order(dat2_task1$bandwidth)]
dat3_task1 = data.table(rep_task1[[3]]$tuning_instance$archive$data())
dat3_task1 = dat3_task1[order(dat3_task1$bandwidth)]

png(filename= "task1_insample_psl.png", width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(dat1_task1$bandwidth), dat1_task1$dens.squared, type = "b", lwd = 2, pch = 16, 
     xlab  = "log bandwidth", ylab = "Out-of-sample empirical PSL",main = "Task 1")
lines(log(dat1_task1$bandwidth), dat2_task1$dens.squared, type = "b", col = "red", lwd = 2, lty = 4, pch = 18)
lines(log(dat1_task1$bandwidth), dat3_task1$dens.squared, type = "b", col = "blue", lwd = 2, lty = 3, pch = 17)
legend("bottomright", 
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

mean_bandwith_task1 = apply(matrix(c(dat1_task1$bandwidth, dat2_task1$bandwidth, dat3_task1$bandwidth), ncol = 3), 1, mean)
mean_loss_task1 = apply(matrix(c(dat1_task1$dens.squared, dat2_task1$dens.squared, dat3_task1$dens.squared), ncol = 3), 1, mean)



###################################################################
#Task 2
#------

#set tuning learner 
NormLrn_task2 = AutoTuner$new(
  learner = lrn,
  resampling = cv1,
  measure = msr("dens.squared"),
  search_space = ps,
  terminator = trm("evals", n_evals = 20),
  tuner = tuner,
  store_tuning_instance = TRUE
)

#d. set up benchmark
#----------------------
set.seed(5202)
grid_task2 = resample(
  task = Task2,
  learner = NormLrn_task2,
  resampling = cv2, 
  store_models=TRUE
)
grid_task2$aggregate()
rep_task2 = lapply(grid_task2$learners, function(x) x$model)
rep_task2

lapply(grid_task2$learners, function(x) x$tuning_result)
as.data.table(grid_task2)$learner[[1]]$tuning_result


grid = resample(
  task = Task2,
  learner = lrn("dens.kde", kernel = "Norm", bandwidth = 0.9557895 ),
  resampling = cv2
)

#extract from first fold
dat1_task2 = data.table(rep_task2[[1]]$tuning_instance$archive$data())
dat1_task2 = dat1_task2[order(dat1_task2$bandwidth)]
dat2_task2 = data.table(rep_task2[[2]]$tuning_instance$archive$data())
dat2_task2 = dat2_task2[order(dat2_task2$bandwidth)]
dat3_task2 = data.table(rep_task2[[3]]$tuning_instance$archive$data())
dat3_task2 = dat3_task2[order(dat3_task2$bandwidth)]

png(filename= "task2_insample_psl.png",  width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(dat1_task2$bandwidth), dat1_task2$dens.squared, type = "b", lwd = 2, pch = 16, 
     xlab  = "log bandwidth", ylab = "Out-of-sample empirical PSL", main = "Task 2")
lines(log(dat1_task2$bandwidth), dat2_task2$dens.squared, type = "b", col = "red", lwd = 2, lty = 4, pch = 18)
lines(log(dat1_task2$bandwidth), dat3_task2$dens.squared, type = "b", col = "blue", lwd = 2, lty = 3, pch = 17)
legend("bottomright", 
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

mean_bandwith_task2 = apply(matrix(c(dat1_task2$bandwidth, dat2_task2$bandwidth, dat3_task2$bandwidth), ncol = 3), 1, mean)
mean_loss_task2 = apply(matrix(c(dat1_task2$dens.squared, dat2_task2$dens.squared, dat3_task2$dens.squared), ncol = 3), 1, mean)


#####################################################################
#Task 3 

#set tuning learner 
set.seed(5202)
NormLrn_task3 = AutoTuner$new(
  learner = lrn,
  resampling = cv1,
  measure = msr("dens.squared"),
  search_space = ps,
  terminator = trm("evals", n_evals = 20),
  tuner = tuner,
  store_tuning_instance = TRUE
)

#d. set up benchmark
#----------------------
set.seed(5202)
grid_task3 = resample(
  task = Task3,
  learner = NormLrn_task3,
  resampling = cv2, 
  store_models=TRUE
)
grid_task3$aggregate()
rep_task3 = lapply(grid_task3$learners, function(x) x$model)
rep_task3

lapply(grid_task3$learners, function(x) x$tuning_result)

set.seed(5202)
grid3 = resample(
  task = Task3,
  learner = lrn("dens.kde", kernel = "Norm", bandwidth = 1.586316),
  resampling = cv2
)
grid3$aggregate()

as.data.table(grid_task3)$learner[[1]]$tuning_result


#extract from first fold
dat1_task3 = data.table(rep_task3[[1]]$tuning_instance$archive$data())
dat1_task3 = dat1_task3[order(dat1_task3$bandwidth)]
dat2_task3 = data.table(rep_task3[[2]]$tuning_instance$archive$data())
dat2_task3 = dat2_task3[order(dat2_task3$bandwidth)]
dat3_task3 = data.table(rep_task3[[3]]$tuning_instance$archive$data())
dat3_task3 = dat3_task3[order(dat3_task3$bandwidth)]

png(filename= "task3_insample_psl.png", width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(dat1_task3$bandwidth), dat1_task3$dens.squared, type = "b", lwd = 2, pch = 16, 
     xlab  = "log bandwidth", ylab = "Out-of-sample empirical PSL", main = "Task 3")
lines(log(dat1_task3$bandwidth), dat2_task3$dens.squared, type = "b", col = "red", lwd = 2, lty = 4, pch = 18)
lines(log(dat1_task3$bandwidth), dat3_task3$dens.squared, type = "b", col = "blue", lwd = 2, lty = 3, pch = 17)
legend("bottomright", 
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

mean_bandwith_task3 = apply(matrix(c(dat1_task3$bandwidth, dat3_task3$bandwidth, dat3_task3$bandwidth), ncol = 3), 1, mean)
mean_loss_task3 = apply(matrix(c(dat1_task3$dens.squared, dat3_task3$dens.squared, dat3_task3$dens.squared), ncol = 3), 1, mean)

##################################################################
#task 4
#------

#set tuning learner 
NormLrn_task4 = AutoTuner$new(
  learner = lrn,
  resampling = cv1,
  measure = msr("dens.squared"),
  search_space = ps,
  terminator = trm("evals", n_evals = 20),
  tuner = tuner,
  store_tuning_instance = TRUE
)

#d. set up benchmark
#----------------------
grid_task4 = resample(
  task = Task4,
  learner = NormLrn_task4,
  resampling = cv2, 
  store_models=TRUE
)
grid_task4$aggregate()
rep_task4 = lapply(grid_task4$learners, function(x) x$model)
rep_task4
lapply(grid_task4$learners, function(x) x$tuning_result)

as.data.table(grid_task4)$learner[[1]]$tuning_result


gridk4 = resample(
  task = Task4,
  learner = NormLrn_task4,
  resampling = cv2, 
  store_models=TRUE
)

#extract from first fold
dat1_task4 = data.table(rep_task4[[1]]$tuning_instance$archive$data())
dat1_task4 = dat1_task4[order(dat1_task4$bandwidth)]
dat2_task4 = data.table(rep_task4[[2]]$tuning_instance$archive$data())
dat2_task4 = dat2_task4[order(dat2_task4$bandwidth)]
dat3_task4 = data.table(rep_task4[[3]]$tuning_instance$archive$data())
dat3_task4 = dat3_task4[order(dat3_task4$bandwidth)]

png(filename= "task4_insample_psl.png", width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(dat1_task4$bandwidth), dat1_task4$dens.squared, type = "b", lwd = 2, pch = 16, 
     xlab  = "log bandwidth", ylab = "Out-of-sample empirical PSL", main = "Task 4")
lines(log(dat1_task4$bandwidth), dat2_task4$dens.squared, type = "b", col = "red", lwd = 2, lty = 4, pch = 18)
lines(log(dat1_task4$bandwidth), dat3_task4$dens.squared, type = "b", col = "blue", lwd = 2, lty = 3, pch = 17)
legend("bottomright", 
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

mean_bandwith_task4 = apply(matrix(c(dat1_task4$bandwidth, dat2_task4$bandwidth, dat3_task4$bandwidth), ncol = 3), 1, mean)
mean_loss_task4 = apply(matrix(c(dat1_task4$dens.squared, dat2_task4$dens.squared, dat3_task4$dens.squared), ncol = 3), 1, mean)


#########################################################
# Task 5
#-----------

#set tuning learner 
NormLrn_task5 = AutoTuner$new(
  learner = lrn,
  resampling = cv1,
  measure = msr("dens.squared"),
  search_space = ps,
  terminator = trm("evals", n_evals = 20),
  tuner = tuner,
  store_tuning_instance = TRUE
)

#d. set up benchmark
#----------------------
grid_task5 = resample(
  task = Task5,
  learner = NormLrn_task5,
  resampling = cv2, 
  store_models=TRUE
)
grid_task5$aggregate()
rep_task5 = lapply(grid_task5$learners, function(x) x$model)
rep_task5
lapply(grid_task5$learners, function(x) x$tuning_result)
as.data.table(grid_task5)$learner[[1]]$tuning_result


#extract from first fold
dat1_task5 = data.table(rep_task5[[1]]$tuning_instance$archive$data())
dat1_task5 = dat1_task5[order(dat1_task5$bandwidth)]
dat2_task5 = data.table(rep_task5[[2]]$tuning_instance$archive$data())
dat2_task5 = dat2_task5[order(dat2_task5$bandwidth)]
dat3_task5 = data.table(rep_task5[[3]]$tuning_instance$archive$data())
dat3_task5 = dat3_task5[order(dat3_task5$bandwidth)]

png(filename= "task5_insample_psl.png",  width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(dat1_task5$bandwidth), dat1_task5$dens.squared, type = "b", lwd = 2, pch = 16, 
     xlab  = "log bandwidth", ylab = "Out-of-sample empirical PSL", main = "Task 5")
lines(log(dat1_task5$bandwidth), dat2_task5$dens.squared, type = "b", col = "red", lwd = 2, lty = 4, pch = 18)
lines(log(dat1_task5$bandwidth), dat3_task5$dens.squared, type = "b", col = "blue", lwd = 2, lty = 3, pch = 17)
legend("bottomright", 
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

mean_bandwith_task5 = apply(matrix(c(dat1_task5$bandwidth, dat2_task5$bandwidth, dat3_task5$bandwidth), ncol = 3), 1, mean)
mean_loss_task5 = apply(matrix(c(dat1_task5$dens.squared, dat2_task5$dens.squared, dat3_task5$dens.squared), ncol = 3), 1, mean)


#####################################################################
#task 6
#------

#set tuning learner 
NormLrn_task6 = AutoTuner$new(
  learner = lrn,
  resampling = cv1,
  measure = msr("dens.squared"),
  search_space = ps,
  terminator = trm("evals", n_evals = 20),
  tuner = tuner,
  store_tuning_instance = TRUE
)

#d. set up benchmark
#----------------------
set.seed(5202)
grid_task6 = resample(
  task = Task6,
  learner = NormLrn_task6,
  resampling = cv2, 
  store_models=TRUE
)
grid_task6$aggregate(msr("dens.squared"))
grid_task6$score()
rep_task6 = lapply(grid_task6$learners, function(x) x$model)
rep_task6
lapply(grid_task6$learners, function(x) x$tuning_result)

as.data.table(grid_task6)$learner[[1]]$tuning_result

#extract from first fold
dat1_task6 = data.table(rep_task6[[1]]$tuning_instance$archive$data())
dat1_task6 = dat1_task6[order(dat1_task6$bandwidth)]
dat2_task6 = data.table(rep_task6[[2]]$tuning_instance$archive$data())
dat2_task6 = dat2_task6[order(dat2_task6$bandwidth)]
dat3_task6 = data.table(rep_task6[[3]]$tuning_instance$archive$data())
dat3_task6 = dat3_task6[order(dat3_task6$bandwidth)]

png(filename= "task6_insample_psl.png",  width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(dat1_task6$bandwidth), dat1_task6$dens.squared, type = "b", lwd = 2, pch = 16, 
     xlab  = "log bandwidth", ylab = "Out-of-sample empirical PSL", main = "Task 6")
lines(log(dat1_task6$bandwidth), dat2_task6$dens.squared, type = "b", col = "red", lwd = 2, lty = 4, pch = 18)
lines(log(dat1_task6$bandwidth), dat3_task6$dens.squared, type = "b", col = "blue", lwd = 2, lty = 3, pch = 17)
legend("bottomright", 
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

mean_bandwith_task6 = apply(matrix(c(dat1_task6$bandwidth, dat2_task6$bandwidth, dat3_task6$bandwidth), ncol = 3), 1, mean)
mean_loss_task6 = apply(matrix(c(dat1_task6$dens.squared, dat2_task6$dens.squared, dat3_task6$dens.squared), ncol = 3), 1, mean)
# #plot
# set.seed(5202)
# grid6 = resample(
#   task = Task6,
#   learner = lrn("dens.kde", kernel = "Norm", bandwidth = 1.901579 ),
#   resampling = cv2
# )
# grid6$score()

##########################################################################################
#Tabulate 
#--------

result_table_tune_psl = data.table(mean_bandwith_task1,log(mean_bandwith_task1), mean_loss_task1, mean_loss_task2, mean_loss_task3, mean_loss_task4, 
                               mean_loss_task5, mean_loss_task6)
write.csv(result_table_tune_psl, "C:\\Users\\Ain Toha\\Dropbox\\Simulation\\resultInSamplePSL.csv")

png("InSample_psl.png",  width = 5, height = 5, units = "in", res= 150, pointsize=11)
plot(log(mean_bandwith_task1), mean_loss_task1, type = "l")
lines(log(mean_bandwith_task1), mean_loss_task2, type = "l", col = "red")
lines(log(mean_bandwith_task1), mean_loss_task3, type = "l", col = "blue")
lines(log(mean_bandwith_task1), mean_loss_task4, type = "l", col = "green")
lines(log(mean_bandwith_task1), mean_loss_task5, type = "l", col = "orange")
lines(log(mean_bandwith_task1), mean_loss_task6, type = "l", col = "purple")
dev.off()
