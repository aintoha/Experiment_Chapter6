#This file contain various functions
#-----------------------------------


#1. gaussiandensity
#==================
#Description
#Usage: gaussiandensity(x, xtest)
#Arguments: 
##  x: sample observations which the estimate is computed
##  xtest: a vector of of test values which we want to estimate their probability density 
##  bw: bandwidth

#Details:
#The vectors of x and xtest are replicated. x is replicated by the 
#length of xtest and then transposed, while xtest is replicated by 
#the length of vector x. D is the difference of each x and each xtest 
#D = (x_i - x_j)
#The kernel function is found and it should be in a matrix
#of a dimension nx*nxtest where nx is the length of vector x and ntest
#is the length of xtest. The mean of the kernel (fx) is the estimated pdf 
#of a gaussian density. fx is a vector of length xtest
#xtest must be more than 1

#Value: estimated density for the point xtest

#Example1:
##  x <- c(1:20)
##  xtest <- c(5.5)
##  bw <- 1
## gaussiandensity(x, xtest, bw)

#Example3:
##  x <- c(1:5)
##  xtest <- c(1.5, 2.5)
##  bw <- 1
##  gaussiandensity(x, xtest, logbw)


#Example2:
## ptm <- proc.time()
##  x <- runif(20, min =1, max= 50)
##  xtest <- c(5, 6.7, 63)
##  logbw <- log(0.4)
##  ptm <- proc.time()
##  
##  proc.time() - ptm
#Author: Ain Toha

gaussiandensity <- function(x, xtest, logbw){
  
  t.bw <- exp(logbw)
  nx <- length(x)
  nxtest <- length(xtest)
  norm <- 1/(sqrt(2*pi)* t.bw)
  xrep <- replicate(nxtest, x)
  txrep <- t(xrep)
  xtestrep <- replicate(nx, xtest)
  D <- txrep - xtestrep
  #D <- sapply(x, function(x, y) (x -y), y = xtest)
  kern <- exp(-(D^2)/(2*t.bw^2))*norm
  fx <- rowMeans(kern)
  #fx <- apply(kern, 1, mean)
  return(fx)
  
}

#--------------------------------------------------------------------------------------------------------

#2 Vectorize Gaussian density
#=============================
#function: v.gaussiandensity
#Description
#Usage: gaussiandensity(x, xtest)
#Arguments: 
##  x: sample observations which the estimate is computed
##  xtest: a vector of of test values which we want to estimate their probability density 
##  bw: bandwidth

#Details:
#The vectors of x and xtest are replicated. x is replicated by the 
#length of xtest and then transposed, while xtest is replicated by 
#the length of vector x. D is the difference of each x and each xtest 
#D = (x_i - x_j)
#The kernel function is found and it should be in a matrix
#of a dimension nx*nxtest where nx is the length of vector x and ntest
#is the length of xtest. The mean of the kernel (fx) is the estimated pdf 
#of a gaussian density. fx is a vector of length xtest
#xtest must be more than 1

#Value: estimated density for the point xtest

v.gaussiandensity <- function(x, xtest, logbw){
  
  nx <- length(x)
  nxtest <- length(xtest)
  xrep <- replicate(nxtest, x)
  txrep <- t(xrep)
  xtestrep <- replicate(nx, xtest)
  D <- txrep - xtestrep
  t.bw <- exp(logbw)
  kern <- lapply(t.bw, function(t.bw, D) exp(-(D^2)/(2*t.bw^2))*(1/(sqrt(2*pi)* t.bw)),
                 D=D)
  fx <- lapply(kern, rowMeans)
  return(fx) 
  
  
}

#----------------------------------------------------------------------------------------------------------

#3. Empirical Log loss of a Gaussian density 
#=====================================

#Function:logloss
#Descrition: This function find the -log loss of nonparametric density for 
#            univariate case. 
#            L = 1/M sum_{j=1}^{M}[-log[1/N sum_{i=1}^{N} 1/(h*sqrt(2pi))
#                exp{-1/2 * (t_{ij}/h)^2}]]  
#            where training set is X = {X_1, ..., X_M}
#            where test set is X* = {X*_1, ..., X*_M}
#            h = bandwidth
#            t_ij = X_i - X*_j

#Usage: logloss(bw, x, xtest)

#Arguments: 
##bw: bandwidth
##x: training set
##xtest: test set 

#Values: the loss value

#Example:
#logloss(bw=log(2), x=c(1:10), xtest =c(0.1, 5.6))

logloss <- function(logh, obsdata, testdata){
  
  dens <- gaussiandensity(logbw=logh, x=obsdata, xtest=testdata)
  loss <- mean(-log(dens))
  
attr(loss, "gradient") <- t.gradient.logloss(logbw=logh, x=obsdata, xtest=testdata)
attr(loss, "hessian") <- t.hessian.logloss(logbw = logh, x=obsdata, xtest=testdata)
  
  return(loss)
  
}

#--------------------------------------------------------------------------------------------------------

#4. Vectorization of log loss
#function: V.logloss

#Description: Vectorized transformed log loss

#Usage: V.logloss(logbw, x, xtest)

#Arguments:
##logbw: a vector of log bandwidth 
##x: A vector of training set
##xtest: A vector of test set

#Value: 
#The output will be a list 

#Example: 
#x <- c(1:10)
#xtest <- c(0.1, 1, 2.3, 4.5,12)
#logbw <- log(c(0.1, 1.1, 2))
#V.logloss(logbw, x, xtest)

V.logloss <- function(logh, obsdata, testdata){
  
  dens <- v.gaussiandensity(x=obsdata, xtest =testdata, logbw=logh)
  logdens <- lapply(dens,  function(dens) -log(dens))
  loss <- lapply(logdens, mean) 
  return(loss)
  
}




#---------------------------------------------------------------------------------------------------------

#5. standard error of Gaussian density log-loss
#================================================
#Function: se.logloss

#Description: Find the standard error of the logloss. The bandwidth need to be tranformed into log before 
#             using this function              

#Usage: se.logloss(logbw, x, xtest)

#Arguments:
##  logbw: log bandwidth
##  x: training set
##  xtest: test set 

#Example1:
##  x <- c(1:20)
##  xtest <- c(5.5, 4.5)
##  bw <- 2
## se.logloss (data=x, testdata=xtest, log

se.logloss <- function(data, testdata, logbw){
  
  M <- length(testdata)
  loss <- -log(gaussiandensity(x=data, xtest=testdata, logbw= logbw))
  empirical_loss <- logloss(obsdata=data, testdata=testdata, logh=logbw)
  diff <- (loss - empirical_loss)^2
  standard_error <- sqrt(sum(diff)/(M*(M-1)))
  
  return(standard_error)
  
}


#------------------------------------------------------------------------------------------------------------

#6.Gradient of Gaussian density log-loss
#==============================================

#Function: gradient.logloss

#Description: This function calculate the gradient of the -log loss density
#             -log(f) w.r.t the bandwidth numerically using the training 
#             data and test data by specifying the bandwidth

#Usage: gradient_logloss(x, xtest,bw)

#Arguments: 
##x: training set
##xtest: test set 
##bw: bandwidth

#value: gradient of -log(fx) w.r.t to the bandwidth

#Example: 
#x <- c(1, 3, 4, 6, 8, 9, 10, 11, 12, 12.3)
#xtest <- c(4.2, 5.9)
#bw <- 0.4
# gradient_logloss(bw, x, xtest)

gradient.logloss <- function(bw, x, xtest){  
  
  nx <- length(x)
  nxtest <- length(xtest)
  xrep <- replicate(nxtest, x)
  txrep <- t(xrep)
  xtestrep <- replicate(nx, xtest)
  t_ij <- txrep - xtestrep
  a <- rowSums(exp(-(t_ij ^2)/(2*bw^2)))
  b <- rowSums((t_ij ^2) * exp(-(1/2)*(t_ij/bw)^2))
  ratio_ab <- b/a
  sum_ratio_ab <- sum(ratio_ab)
  c <- 1/(nxtest * bw^3)
  d <- 1/bw
  gradient <- d - c*sum_ratio_ab
  return(gradient)
  
}

#------------------------------------------------------------------------------------------------------

#function: t.gradient.logloss

#Description: This function calculate the gradient of the -log loss density
#             L-log(f) w.r.t the log(bandwidth) numerically using the training 
#             data and test data.  A log(bandwidth) must be supplied such that
#             logbw = log(bandwidth)
#             dL/d(logbw) = dL/d(log(bw)) = dL/d(bw) * 
#                             d(bw)/d(logbw)

#Usage: t.gradient.logloss(bw, x, xtest)

#Arguments: 
##bw: bandwidth
##x: training set
##xtest: test set 

#value: gradient of -log(fx) w.r.t to the log(bandwidth)

#Example: 
#x <- c(1, 3, 4, 6, 8, 9, 10, 11, 12, 12.3)
#xtest <- c(4.2, 5.9)
#logbw <- log(1)
#t.gradient.logloss(logbw, x, xtest)

t.gradient.logloss <- function(logbw, x, xtest){  
  
  t.bw <- exp(logbw) 
  nx <- length(x)
  nxtest <- length(xtest)
  xrep <- replicate(nxtest, x)
  txrep <- t(xrep)
  xtestrep <- replicate(nx, xtest)
  D <- txrep - xtestrep
  
  a <- rowSums(exp(-(D^2)/(2*t.bw^2)))
  b <- rowSums((D^2) * exp(-(D^2)/(2*t.bw^2)))
  ratio_ab <- b/a
  sum_ratio_ab <- sum(ratio_ab)
  c <- 1/(nxtest * t.bw^2)
  gradient <- 1 - (c*sum_ratio_ab)
  return(gradient)
  
}


#--------------------------------------------------------------------------------------------------------
#7. 2nd derivative of the log loss
#=================================
#function: hessian_logloss

#Description: This function calculate the second derivative of the log 
#             loss function 

#Usage: hessian_logloss(bw, x, xtest)

#Arguments:
##  bw: bandwidth
##  x: training set
##  xtest: test set

#Value: The value of the second the derivative

#Example:
#hessian.logloss(bw=1, x=c(1:10), xtest=0.4)
#0.02600986

hessian.logloss <- function(bw, x, xtest){
  
  nx <- length(x)
  nxtest <- length(xtest)
  bw <- bw
  xrep <- replicate(nxtest, x)
  txrep <- t(xrep)
  xtestrep <- replicate(nx, xtest)
  D <- txrep - xtestrep
  a <- rowSums(exp(-(D^2)/(2*bw^2)))
  b <- rowSums((D^2) * exp(-(D^2)/(2*bw^2)))
  c <-rowSums((D^4) * exp(-(D^2)/(2*bw^2)))
  b2 <- b^2
  a2 <- a^2
  ratio_ba <- b/a
  sum_ratio_ba <- sum(ratio_ba)
  ratio_ba2 <- b2/a2
  first <- -1/(bw^2)
  second <- (3/((bw^4)*nxtest)) * sum_ratio_ba
  ratio_ca <- c/a
  sum_ratio_ca <- sum(ratio_ca)
  third <- (1/((bw^6)*nxtest)) * sum(ratio_ba2)
  four <- -(1/((bw^6)*nxtest))*sum_ratio_ca
  result <- first + second + third+four
  return(result)
}


#---------------------------------------------------------------------------------------------------------


#function: t.hessian_logloss

#Description: This function calculate the second derivative of the log loss
#             function. 
#             h = log(gamma)
#             gamma = exp(h)
#             d(gamma)/dh = exp(h)

#Usage: t.hessian_logloss(logbw, x, xtest)

#Arguments: 
##    logbw: The input the bandwidth must be log of the bandwidth
##    x: training set
##    xtest: test set

#Value: 

#Example: 
#   t.hessian_logloss(logbw=log(1), x=c(1:10), xtest=0.4)
# 0.0183706

t.hessian.logloss<- function(logbw, x, xtest){
  
  h <- exp(logbw)
  nx <- length(x)
  nxtest <- length(xtest)
  xrep <- replicate(nxtest, x)
  txrep <- t(xrep)
  xtestrep <- replicate(nx, xtest)
  D <- txrep - xtestrep
  a <- rowSums(exp(-(D^2)/(2*h^2)))
  b <- rowSums((D^2) * exp(-(D^2)/(2*h^2)))
  c <- rowSums((D^4) * exp(-(D^2)/(2*h^2)))
  a2 <- a^2
  b2 <- b^2
  ratio_ba<- b/a
  ratio_ba2 <- b2/a2
  ratio_ca <- c/a
  sum_ratio_ba <- sum(ratio_ba)
  first <- (2/(nxtest * (h^2)))*sum_ratio_ba
  second <- (1/(nxtest*(h^4)))*sum(ratio_ba2 - ratio_ca)
  result <- first + second
  return(result)
}

#-----------------------------------------------------------------------
#MISE GAUSSIAN 

#----------------------------------------------------------------------
#ise

#Function: gaussian.ise

#Description: The function computes the integrated squared error (ISE)by
#             Rudemo (182) and Bowman (1984). Unlike the original method that uses
#             LOOCV, this method uses k-fold cross validation. required the
#             function fxsquared 

#Usage: gaussia.ise(x, xtest, bw)

#Arguments: 
#### x: training data
#### xtest : test data 
#### bw: bandwidth

#Value: a numeric value 

#Examples: 
##  set.seed(1001099)
##  data <- c(1:10, 1, 2, 0.14)
##  testdata <- c(0.14, 0.65)
##  bw <- 0.001
##  gaussian.ise(x = data, xtest = testdata, bw=bw)

fxsquared <- function(x, bw){  
  
  h <- bw
  N <- length(x)
  D <- sapply(x, function(x, y) ((x - y)/h), y=x)
  expD <- exp(-(D^2)/4)
  totalx <- colSums(expD)
  totalrepx <- sum(totalx)
  fx.squared <- totalrepx/(2*(N^2)*h*sqrt(pi))
  
  return(fx.squared)
}

gaussian.density <- function(x, xtest, bw){
  
  nx <- length(x)
  nxtest <- length(xtest)
  norm <- 1/(sqrt(2*pi)*bw)
  xrep <- replicate(nxtest, x)
  txrep <- t(xrep)
  xtestrep <- replicate(nx, xtest)
  D <- txrep - xtestrep
  #D <- sapply(x, function(x, y) ((x -y)/bw), y = xtest)
  expD <- exp(-(D^2)/(2*(bw^2)))
  kern <- expD*norm
  fx <- rowMeans(kern)
  return(fx)
  
}


gaussian.ise <- function(x, xtest, bw){
  
  fx2 <- fxsquared(x=x, bw=bw)
  dens <- gaussian.density(x =x, xtest = xtest, bw=bw)
  loss <- -2*dens + fx2
  emp.loss <- mean(loss)
  
# attr(emp.loss, "gradient") <- gradient.gaussian.ise(x=x , xtest = xtest,bw =bw)
  return(emp.loss)
}


#--------------------------------------------------------------------
#Function: se.gaussian.ise

#Description: Return the standard error of the empirical loss using the ise as the loss

#Arguments:
##x: training data 
##xtest: test data
##bw: bandwidth

#Value: value of the standard error 

#Example: 
## set.seed(10010)
## x <- runif(10)
## xtest <- c(0.44, 0.22)
## bw <- 2
## se.gaussian.ise(x, xtest, bw)

se.gaussian.ise <- function(x, xtest, bw){
  
  M <- length(xtest)
  dens <- gaussian.density(x=x, xtest=xtest, bw=bw)
  fx2 <- fxsquared(x=x, bw=bw)
  loss <- -2*dens + fx2
  
  emp.loss <- -2*mean(dens) + fx2
  dif <- sum((loss - emp.loss)^2)
  se <- sqrt(dif/(M*(M-1)))
  return(se)
  
}


#---------------------------------------------------------------------------
#Gradient empirical out-of-sample squared loss 
#Gaussian kernel

#Description: This is the empirical out-of-sample squared loss 
#             using the Gaussian kernel
#Arguments:
#x: a vector of training data (must be more than one. if only 
#           one then make 0 as another point)
#xtest : a vector of test data (must be more than one. if only 
#           one then make 0 as another point)
#bw: the bandwidth

#Value: a numerical value that is the gradient 

#Example: 
# x = c(1,2)
# xtest = c(1, 3)
# bw = 1
# gradient.gaussian.ise(x, xtest, bw=h)

gradient.gaussian.ise <- function(x, xtest, bw){
  
  h <- bw
  norm <- 1/(sqrt(2*pi))
  N <- length(x)
  M <- length(xtest)
  
  #functions for the int f(x)^2
  repx <- rep(x)
  Nrepx <- length(repx)
  D1 <- sapply(x, function(x,y) (x-y), y=repx)
  expD1 <- exp(-(1/4)*((D1/h)^2))
  totalx1 <- sum(expD1)
  A <- -(totalx1)/(2*(h^2)*(N^2)*sqrt(pi))
  
  expD2 <- (D1^2)*(expD1)
  totalx2 <- sum(expD2)
  B <- totalx2/(4*(h^4)*(N^2)*sqrt(pi))
  lh <- A + B 
  
  #functions for the train and test
  xrep <- replicate(M, x)
  txrep <- t(xrep)
  xtestrep <- replicate(N, xtest)
  Dtest <- txrep - xtestrep
  expDtest <- exp(-(Dtest^2)/(2*(h^2)))
  #row= xtest; column= x
  totalxtest <- sum(expDtest)
  C <- (2*totalxtest*norm)/(N*M*(h^2))
  
  #expDtest1 <- sapply(Dtest, function(x,y) (x^2)*y, y=expDtest)
  expDtest1 <- (Dtest^2)*expDtest
  rowtotalxtest1 <- rowSums(expDtest1)
  totalxtest1 <- sum(rowtotalxtest1)
  E <- (-2*totalxtest1*norm)/((h^4) * N * M)
  
  gh <- C+ E
  gradient <- lh + gh
  return(gradient)
  
}

#-----------------------------------------------------------------------
#Gradient empirical out-of-sample squared loss 
#Gaussian kernel

#Description: This is the empirical out-of-sample squared loss 
#             using the Gaussian kernel
#Arguments:
#traindata: a vector of training data (must be more than one. if only 
#           one then make 0 as another point)
#testdata : a vector of test data (must be more than one. if only 
#           one then make 0 as another point)
#h: the bandwidth

#Value: a numerical value that is the gradient 

#Example: 
#traindata = c(1,0)
#testdata = c(2.5, 0)
#h = 0.1
#gradient.squared.loss(traindata = traindata, testdata = testdata, h=1)


gradient.squared.loss <- function(traindata, testdata, h){
  
  #this is for the term that include training and test 
  N <- length(traindata)
  M <- length(testdata)
  diff1 <- sapply(traindata, function(x,y) (x-y), y=testdata)
  exp_diff1 <- exp(-(1/2)*(diff1/h)^2)
  sum_exp_diff1 <- (2/(sqrt(2*pi)^(h^2)))*sum(exp_diff1)
  total_first <- sum_exp_diff1*(1/(N*M) - mean(diff1^2)/(h^2))
  
  #this is for the training part int f(x)^2

  rep_train <- rep(traindata)
  diff2 <- sapply(traindata, function(x,y) (x-y), y=rep_train)
  exp_diff2 <- exp(-(1/4)*(diff1/h)^2)
  sum_exp_diff2 <- (1/(2*sqrt(pi)^(h^2)))*sum(exp_diff2)
  total_second <- sum_exp_diff2*(1/(N^2) + mean(diff2^2)/2)
  
  gradient_squared_loss <- total_first - total_second
  return(gradient_squared_loss)
  
}



#----------------------------------------------------------------------------------------------------------

#Function: t.area
#==================
#Finding the area under the curve using the trapezium rule (can be used 
#for uniformed and non-uniformed grid)

#Usage: t.area(xtrain, ytrain)

#arguments: 
##xtrain: vector of x-axis values
##ytrain: vector of y-axis values

#Details:
#trapezoid rule:
#the x-axis values will be sorted first (if has not been sorted before)
#the y-axis values will be sorted according to the x-axis
#let the sorted x-axis be xtrain={x_1, ..., x_n}
#let [a,b] in x={x_1,..., x_n} where a=x_1(first) and b=x_n(last)
#the equation is difx/2 * (y_1 + 2y_2 + ... + 2y_(n-1) + y(n)) where
#difx = (x_n - x_1)/(n-1) and n is the length of x-axis
#Note that the length of x-axis and y-axis must be the same

#values: 
##trapezoidarea: area under the Y curve

#Example1: 
# set.seed(100)
# n <- 512
# d <- 0.3
# xtrain <- seq.int(1, by = d, length.out = n)
# ytrain <- sapply(xtrain, function(x) x +2.5)
# t.area(xtrain, ytrain)

#Example2:
# xtrain <- c(1:5)
# ytrain <- c(3:7)
# t.area(xtrain, ytrain)

#Author: Ain Toha

t.area <- function(x, y){
  
  xs <- sort(x, decreasing = FALSE)
  ys <- y[order(x)]
  lengthx <- length(xs)  
  dx <- 2:lengthx
  area <- (xs[dx] - xs[dx-1]) %*% (ys[dx] + ys[dx-1])/2
  return(area)
}


#-------------------------------------------------------------------------------------------------------

#8. Function: splitfolder
#==================


#Description: Split the folder into training and validation set.

#Usage: spliting(dataset, folder, numberoffolds)

#Arguments: 
###dataset: dataset(observation) can either be a dataframe or matrix
###folders: folder that as been named before
###numberoffolds: an integer giving the number of groups which the data (observation)
#               should be split

#Value:As a list
###[[1]]: group for training set
###[[2]]: group for validation test

#Example: 
#   set.seed(100)
#   x <- runif(100, min = 0, max=1000)
#   x <- as.data.frame(x)
#   library(cvTools) #install.packages("cvTools")
#   fold1 <- cvFolds(NROW(x), K= 5)
#   splitfolder(x, fold1, 2)

#Author: Ain Toha

splitfolder <- function(dataset, folders, numberoffolds){
  
  trainset <- dataset[folders$subsets[folders$which != numberoffolds], ] 
  #Set the training set
  validationset <- dataset[folders$subsets[folders$which == numberoffolds], ] 
  #Set the validation set
  return(list(trainset, validationset))
  
}

#-------------------------------------------------------------------------------------------------


#9. function: splitdata
#========================


#Description: this function split data into training vector and test vector
#             based on ratio. The data must be in vector format

#Usage: splitdata(data, ratio)

#Arguments: 
##  data: The data set used must be in matrix format
##  ratio: Ratio for the training set (numeric)


#Values: 
## The output is a list
## splitdata(data, ratio)[[1]] is the training set
## splitdata(data, ratio)[[2]] is the test set



##Example
# set.seed(100)
# data <- rnorm(100, mean=1, sd = 2)
# ratio <- 0.7

splitdata <- function(data, ratio){
  
  smp_size <- floor(ratio * length(data))
  train_ind <- sample(seq_len(length(data)), size = smp_size, replace = F)
  train <- data[train_ind]
  test <- data[-train_ind]
  
  return(list(train, test))
  
}    

#10. Bootstrap
#===========

#Function: bootstrap

#Description: This function is used in replicating the original data. This
#             is done by bootstrap. In bootstrap, the 

#Arguments:


bootstrap <- function(data,numbags){
  
  set.seed(100)
  data <- as.matrix(data)
  training_positions <-replicate(numbags,sample(nrow(data), 
                                                size=nrow(data), replace=T))
  resample <- matrix(NA, nrow(data), numbags)
  for(i in 1:numbags){
    resample[,i] <- as.matrix(data[training_positions[,i]])
  }
  return(resample)
  
}  
