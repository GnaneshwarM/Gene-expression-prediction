install.packages("ggplot2")
install.packages("dplyr")
install.packages("corrplot")
install.packages("caTools")

library(caTools)
library(corrplot)
library("tidyverse")
library(ggplot2)
library(dplyr)

#Set the directory
setwd("/Users/gnani/Downloads/Candy Statics Assignment") 

# get the CSV file
genedata <- read.csv("genedata.csv")

max(genedata$x3)

head(genedata)

# Line Plot for all the  genes

# Line plot for gene x1
ggplot(data = genedata, aes(x = genedata[,c(1)], y = x1))+
  geom_line(color = "#00AFBB", size = 1)

# Line plot for gene x2
ggplot(data = genedata, aes(x = genedata[,c(1)], y = x2))+
  geom_line(color = "#FC4E07", size = 1)

# Line plot for gene x3
ggplot(data = genedata, aes(x = genedata[,c(1)], y = x3))+
  geom_line(color = "#612CCB", size = 1)

# Line plot for gene x4
ggplot(data = genedata, aes(x = genedata[,c(1)], y = x4))+
  geom_line(color = "#830B87", size = 1)

# Line plot for gene x5
ggplot(data = genedata, aes(x = genedata[,c(1)], y = x5))+
  geom_line(color = "#6F6F66", size = 1)




#basic histogram
#Histogram for x1 in genedata
hist(genedata$x1 , col="lightblue")
abline(v = mean(genedata$x1), lwd= 2)#adds vertical line at the mean value
text(1.6 , 30 , "mean = 1.40" , adj= c(0,0)) #adds text

#Histogram for x2 in genedata
hist(genedata$x2 , col="lightblue")
abline(v = mean(genedata$x2), lwd= 2)  #adds vertical line at the mean value
text(1.5 , 40 , "mean = 1.122" , adj= c(0,0)) #adds text

#Histogram for x3 in genedata
hist(genedata$x3 , col="lightblue")
abline(v = mean(genedata$x3), lwd= 2)  #adds vertical line at the mean value
text(1.7 , 40 , "mean = 1.121" , adj= c(0,0)) #adds text

#Histogram for x4 in genedata
hist(genedata$x4 , col="lightblue")
abline(v = mean(genedata$x4), lwd= 2)  #adds vertical line at the mean value
text(1.5 , 40 , "mean = 1.415" , adj= c(0,0)) #adds text

#Histogram for x5 in genedata
hist(genedata$x5 , col="lightblue")
abline(v = mean(genedata$x5), lwd= 2)  #adds vertical line at the mean value
text(0.5 , 50 , "mean = 1.599" , adj= c(0,0)) #adds text

genes <- genedata[, c(2,3,4,5,6)]
## Scatter Matrix
head(genes)
lapply(genedata[, 2:6], sd) # Standard Deviation
lapply(genedata[, 2:6], var) # variance

genes <- genedata[, c(2,3,4,5,6)]
pairs(genes[,1:5], pch = 19)

sapply(genedata[, 2:6], sd) # Standard Deviation
sapply(genedata[, 2:6], var) # variance

# Compute the mean of each column
sapply(genes, mean)
# Compute the median of each column
sapply(genes, median)
# Compute the median absolute deviation
sapply(genes, mad)


# Compute the minimum value
sapply(genes, min)

# Compute the maximum value
sapply(genes, max)

# Range
sapply(genes, range)

#Interquartile range
sapply(genes, quantile)
sapply(genes, IQR)


# correlation Plot

library(corrplot)
M <-cor(genes)
corrplot(M)
 
# heat Map for the genes
palatte = colorRampPalette(c("green","white","red")) ((10))
heatmap(x=genes.cor, col=palatte, symm = TRUE)

#boxplot
boxplot(genedata$x1,genedata$x2,genedata$x3,genedata$x4,genedata$x5)# 



#correlation plot
M <-cor(genes)
corrplot(M,method = c("square") , type = c("full"))



# Missing Value detection

missmap(genedata,col=c('yellow','black'),y.at=1,y.labels="",legend=TRUE)

#corrplot

corrplot(cor(genes),type="upper",method="color",addCoef.col = "black",tl.col = "black", number.cex = 0.6,mar=c(0,0,0,0))

# Scatter Matrix
X = rbind( matrix(genedata$`x1`,1,301) , matrix(genedata$`x2`,1,301) , matrix(genedata$`x3`,1,301) , 
           matrix(genedata$`x4`,1,301) , matrix(genedata$`x4`,1,301) )

par(mfcol=c(2,5))

# Plots for all the Variables
plot( X[2,] , X[1,] , main="x1,x2")
plot( X[3,] , X[1,] , main="x1,x3")
plot( X[4,] , X[1,] , main="x1,x4")
plot( X[5,] , X[1,] , main="x1,x5")

plot( X[3,] , X[2,] , main="x2,x3")
plot( X[4,] , X[2,] , main="x2,x4")
plot( X[5,] , X[2,] , main="x2,x5")

plot( X[4,] , X[3,] , main="x3,x4")
plot( X[5,] , X[3,] , main="x3,x5")

plot( X[5,] , X[4,] , main="x4,x5")


##PCA (Principal component analysis) ##

## X contains all the variable except the dependent term
mean_X_rows = matrix( rowMeans( X[1:5,] ) , 5 , 301 ) # Mean matrix of rowise means

X_features_de_mean = X[1:5,] - mean_X_rows # remove mean of the features from X. Note: features are from rows 1 to 4

Cx = (1/301) * ( X_features_de_mean %*% t(X_features_de_mean) )

#eigen Values and vectors
Eig = eigen(Cx)
Eig_values = Eig$values
Eig_vectors = Eig$vectors

print(Eig_values)

print(Eig_vectors)



# defining the dimention by k = 2
k=2
P = Eig_vectors[1:k,]

Y_k2 = P %*% X_features_de_mean # Y for k = 2

plot( Y_k2[2,] , Y_k2[1,] , main="1,2", col = c("red", "blue")) # for k = 2

# Task 3 

genedata <- read.csv("genedata.csv")

gene_model <- genedata[, c(1,4,5,6)] # copies the column time(min), X3,x4 and x5


# Test train split using sample.split(function)
set.seed(123)
split = sample.split(gene_model$x3, SplitRatio = 0.8)

train_set= as.matrix(subset(gene_model, split == TRUE))
test_set= as.matrix(subset(gene_model,split == FALSE ))

time_train = cbind(train_set[,c(1)])

y = cbind(train_set[,c(2)])

x4_train =cbind(train_set[,c(3)]) 

x5_train =cbind(train_set[,c(4)]) 



train_set_for_prediction = cbind(x4_train,x5_train)

# training dataset which contains all the 8 variable
X1 = cbind(x4_train,x4_train^2,x4_train^3,x4_train^4,x5_train,x5_train^2,x5_train^3,x5_train^4)

time_test = cbind(test_set[,c(1)])

y_test = cbind(test_set[,c(2)])

x4_test =cbind(test_set[,c(3)]) 

x5_test =cbind(test_set[,c(4)]) 

# testing dataset which contains all the 8 variable
X2 = cbind(x4_test,x4_test^2,x4_test^3,x4_test^4,x5_test,x5_test^2,x5_test^3,x5_test^4)



# creating the list to store the MSE, SSE, error values of the first iteratiom
iteration <- 1
l_MSE <- vector("list", iteration)
l_SSE <- vector("list", iteration)
l_error <- vector("list", iteration)

# While loop for the first iteration
while (iteration < 9) {
  X <- cbind(X1[,c(iteration)])
  X_test<- cbind(X2[,c(iteration)])
  thetaHat = solve(t(X) %*% X) %*% t(X) %*% y  
  y_Hat = X_test %*% thetaHat # compute y hat
  error = y_test - y_Hat # computes error
  l_error[[iteration]] <- error
  SSE = norm(error, type = "2")^2
  l_SSE[[iteration]] <- SSE
  MSE = SSE/61
  l_MSE[[iteration]] <- MSE
  iteration  <-  iteration + 1
}

#MSE of the first iteration
MSE_of_iteration1 = min(unlist(list(l_MSE)))      

# index postion of the MSE
index_of_mse1 = which.min(unlist(list(l_MSE)))  

value_of_first_x_train = X1[,c(index_of_mse1)]

value_of_first_x_test = X2[,c(index_of_mse1)]

# removing the first iteration values from the second iteration  
X_train_iteration2 = subset(X1, select = -c(index_of_mse1) )

X_test_iteration2 = subset(X2, select = -c(index_of_mse1) )

# creating the listd for the second iteration 
iteration2 <- 1
l_MSE2 <- vector("list", iteration2)
l_SSE2 <- vector("list", iteration2)
l_thetahat_2 <- vector("list", iteration2)
l_yhat_2 <- vector("list", iteration2)

#The second iteration for the second column
while (iteration2 < 8) {
  X_2 <- cbind(value_of_first_x_train,X_train_iteration2[,c(iteration2)])
  X_2test<- cbind(value_of_first_x_test,X_test_iteration2[,c(iteration2)])
  thetaHat_2 = solve(t(X_2) %*% X_2) %*% t(X_2) %*% y  
  l_thetahat_2[[iteration2]] <- thetaHat_2
  y_Hat_2 = X_2test %*% thetaHat_2 # compute y hat
  l_yhat_2[[iteration2]]  <- y_Hat_2
  error_2 = y_test - y_Hat_2 # computes error
  SSE_2 = norm(error_2, type = "2")^2
  l_SSE2[[iteration2]] <- SSE_2
  MSE_2 = SSE_2/61
  l_MSE2[[iteration2]] <- MSE_2
  iteration2  <-  iteration2 + 1
}

# smallest of the MSE in the list MSE,SSE
MSE_of_iteration2= min(unlist(list(l_MSE2)))      

SSE = min(unlist(list(l_SSE2)))

index_of_mse2 = which.min(unlist(list(l_MSE2)))  


value_of_second_x_train = X_train_iteration2[,c(index_of_mse2)]

value_of_second_x_test = X_test_iteration2[,c(index_of_mse2)]


final_X = cbind(value_of_first_x_train,value_of_second_x_train)

final_x_test = cbind(value_of_first_x_test,value_of_second_x_test)

thetaHat_final = solve(t(final_X) %*% final_X) %*% t(final_X) %*% y  


y_hat_final = final_x_test %*% thetaHat_final

plot(time_test,y_test, col="blue") #time versus the given data 'y' which is x3 testing data
lines(time_test,y_hat_final, col="red")#  prediction line for time versus predicted values

## parameter Covariance Matrix

n = length(y_test)

sigma_2 = SSE/( n - 1 ) # sample variance

cov_thetaHat = 1 * (solve(t(final_X) %*% final_X))# covariance Matrix
cov_thetaHat_inv = (t(final_X) %*% final_X) * (1/sigma_2)  # inverse of parameter covariance matrix
det_cov_thetaHat = det(cov_thetaHat) # determinent of cov_thetaHat

no_points  = 20
number_of_parameters = 2
theta_1_range = seq(thetaHat_final[1]-0.005 , thetaHat_final[1]+0.005 , length=no_points) # range of parameter 1
theta_2_range = seq(thetaHat_final[2]-0.005 , thetaHat_final[2]+0.005 , length=no_points) # range of parameter 2

p_thetaHat_D = matrix(0 , no_points , no_points) # initializing an emoty p.d.f matrix

for(r in 1:20){
  for(c in 1:20){
    
    theta_12 = matrix( c( theta_1_range[r] , theta_2_range[c] ) , number_of_parameters , 1)
    thetaHat_theta = theta_12 - thetaHat_final
    
    p_thetaHat_D[r,c] = ( 1/sqrt( ( (2*pi)^number_of_parameters ) * det_cov_thetaHat) ) * exp( -0.5 * t(-thetaHat_theta) %*% cov_thetaHat_inv %*% -thetaHat_theta )
  }
}

contour(theta_1_range, theta_2_range, p_thetaHat_D) # contour plot
persp(theta_1_range, theta_2_range, p_thetaHat_D , theta = 20 , phi = 50) # 3d plot p.d.f

## Prediction on training data
y_hat_pred_train = final_X %*% thetaHat_final

## Computing 95% CI with error bars
number_of_parameters =2 

var_y_hat = matrix(0 , n , 1)

for( i in 1:n){
  X_i = matrix( X[i,] , 1 , number_of_parameters ) # X[i,] creates a vector. Convert it to matrix
  var_y_hat[i,1] = X_i %*% cov_thetaHat %*% t(X_i) # same as sigma_2 * ( X_i %*% ( solve(t(X) %*% X)  ) %*% t(X_i) )
}

CI = 2 * sqrt(var_y_hat) # Confidance interval

dim(y_hat_final)


plot(time_test, y_hat_final , type = "l")
segments(time_test, y_hat_final-CI, time_test, y_hat_final+CI) # Adds error bars to the indivigual data points
segments(time_test, mean(y_hat_final)-CI, time_test, mean(y_hat_final)+CI) # Adds error bars to the indivigual data points





# Validating using new test , train splitting approach
split = sample.split(gene_model$x3, SplitRatio = 2/3)
train_set_for_validation= as.matrix(subset(gene_model, split == TRUE))
test_set_for_validation= as.matrix(subset(gene_model,split == FALSE ))

time_train_val = cbind(train_set_for_validation[,c(1)])
y_train_val = cbind(train_set_for_validation[,c(2)])

x4_train_val =cbind(train_set_for_validation[,c(3)]) 

x5_train_val =cbind(train_set_for_validation[,c(4)]) 



# training dataset which contains all the 8 variable
X1_val = cbind(x4_train_val,x4_train_val^2,x4_train_val^3,x4_train_val^4,x5_train_val,x5_train_val^2,x5_train_val^3,x5_train_val^4)

time_test_val = cbind(test_set_for_validation[,c(1)])

y_test_val = cbind(test_set_for_validation[,c(2)])

x4_test_val =cbind(test_set_for_validation[,c(3)]) 

x5_test_val =cbind(test_set_for_validation[,c(4)]) 

# testing dataset which contains all the 8 variable
X2_val = cbind(x4_test_val,x4_test_val^2,x4_test_val^3,x4_test_val^4,x5_test_val,x5_test_val^2,x5_test_val^3,x5_test_val^4)


X_val <- cbind(X1_val[,c(index_of_mse1)],X1_val[,c(index_of_mse2)])
X_test_val<- cbind(X2_val[,c(index_of_mse1)],X2_val[,c(index_of_mse2)])
thetaHat_val = solve(t(X_val) %*% X_val) %*% t(X_val) %*% y_train_val 
dim(thetaHat_val)
y_Hat_val = X_test_val %*% thetaHat_val


error_val = y_test_val - y_Hat_val # computes error
SSE_val = norm(error_val, type = "2")^2
MSE_val = SSE_val/61
MSE_val

plot(time_test_val,y_test_val, col="blue")
lines(time_test_val,y_Hat_val, col="red")





## Approximate Bayesian Computation
install.packages("EasyABC")
library("EasyABC")
observedData =  final_X
observedSummary = c(mean(final_X), sd(final_X))


#geene model for ABC
gene_model_bayes <- function(x){
  + c( final_x_test[1] + final_x_test[2] + runif(61,0,0.1), final_x_test[1] * final_x_test[2] + runif(61,0,0.1) )
}

#prior list to create estimated list
gene_prior = list(c("unif",0,1),c("normal",mean(final_x_test),sd(final_x_test)))

n = 61

# applying rejection ABC
ABC_rej<-ABC_rejection(model=gene_model_bayes, prior=gene_prior, nb_simul=2)




  
  
  
