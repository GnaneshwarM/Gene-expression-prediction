genedata <- read.csv("genedata.csv")

gene_model <- genedata[, c(4,5,6)] # copies the column of X3,x4 and x5

str(gene_model)


split = sample.split(gene_model$x3, SplitRatio = 0.8)
train_set= as.matrix(subset(gene_model, split == TRUE))
test_set= as.matrix(subset(gene_model,split == FALSE ))


y = cbind(train_set[,c(1)])

x4_train =cbind(train_set[,c(2)]) 

x5_train =cbind(train_set[,c(3)]) 

train_set_for_prediction = cbind(x4_train,x5_train)

# training dataset which contains all the 8 variable
X1 = cbind(x4_train,x4_train^2,x4_train^3,x4_train^4,x5_train,x5_train^2,x5_train^3,x5_train^4)

y_test = cbind(test_set[,c(1)])

x4_test =cbind(test_set[,c(2)]) 

x5_test =cbind(test_set[,c(3)]) 

X2 = cbind(x4_test,x4_test^2,x4_test^3,x4_test^4,x5_test,x5_test^2,x5_test^3,x5_test^4)



iteration <- 1
l_MSE <- vector("list", iteration)
l_SSE <- vector("list", iteration)
l_error <- vector("list", iteration)
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

MSE_of_Param_1 = min(unlist(list(l_MSE)))      


SSE_of_Param_1 = min(unlist(list(l_SSE)))

index_of_mse1 = which.min(unlist(list(l_MSE)))  





value_of_first_x_train = X1[,c(index_of_mse1)]

value_of_first_x_test = X2[,c(index_of_mse1)]



X_train_iteration2 = subset(X1, select = -c(index_of_mse1) )

X_test_iteration2 = subset(X2, select = -c(index_of_mse1) )


iteration2 <- 1
l_MSE2 <- vector("list", iteration2)
l_SSE2 <- vector("list", iteration2)
l_thetahat_2 <- vector("list", iteration2)
l_yhat_2 <- vector("list", iteration2)
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


MSE_of_Param_2= min(unlist(list(l_MSE2)))      


SSE_of_Param_2 = min(unlist(list(l_SSE2)))

index_of_mse2 = which.min(unlist(list(l_MSE2)))  




value_of_second_x_train = X_train_iteration2[,c(index_of_mse2)]

value_of_second_x_test = X_test_iteration2[,c(index_of_mse2)]



final_X = cbind(value_of_first_x_train,value_of_second_x_train)

final_x_test = cbind(value_of_first_x_test,value_of_second_x_test)

thetaHat_final = solve(t(final_X) %*% final_X) %*% t(final_X) %*% y  

y_hat_final = final_x_test %*% thetaHat_final

## parameter Covariance Matrix

n = length(y_test)
sigma_2 = SSE_of_Param_2/( n - 1 )

cov_thetaHat = 1 * (solve(t(final_X) %*% final_X))
cov_thetaHat_inv = (t(final_X) %*% final_X) * (1/sigma_2)
det_cov_thetaHat = det(cov_thetaHat) # determinent of cov_thetaHat
cov_thetaHat
cov_thetaHat_inv
det_cov_thetaHat

no_points  = 20
number_of_parameters = 2
theta_1_range = seq(0.009 , 0.150 , length=no_points)
theta_2_range = seq(0.001 , 0.006 , length=no_points)

p_thetaHat_D = matrix(0 , no_points , no_points)

for(r in 1:20){
  for(c in 1:20){
    
    theta_12 = matrix( c( theta_1_range[r] , theta_2_range[c] ) , number_of_parameters , 1)
    thetaHat_theta = theta_12 - thetaHat_final
    p_thetaHat_D[r,c] = ( 1/sqrt( ( (2*pi)^number_of_parameters ) * det_cov_thetaHat) ) * exp( -0.5 * t(-thetaHat_theta) %*% cov_thetaHat_inv %*% -thetaHat_theta )
    
    
  }
}





## Computing 95% CI with error bars
number_of_parameters =2 

var_y_hat = matrix(0 , n , 1)

for( i in 1:n){
  X_i = matrix( X[i,] , 1 , number_of_parameters ) # X[i,] creates a vector. Convert it to matrix
  var_y_hat[i,1] = X_i %*% cov_thetaHat %*% t(X_i) # same as sigma_2 * ( X_i %*% ( solve(t(X) %*% X)  ) %*% t(X_i) )
}

CI = 2 * sqrt(var_y_hat) # Confidance interval

plot(final_x_test, y_hat_final , type = "l")
segments(final_x_test, y_hat_final-CI, final_x_test, y_hat_final+CI) # Adds error bars to the indivigual data points




## Prediction 

y_hat_pred = train_set_for_prediction %*% thetaHat_final






## Approximate Bayesian Computation
n = 20000
fit = data.frame(shape = runif(n, 1, 6), scale = runif(n, 1,10), summary1 = rep(NA, n), summary2 = rep(NA, n), distance = rep(NA, n))

for (i in 1:n){
  prediction <- model(fit[i,1:2])
  deviation = sqrt(sum((prediction- observedSummaryStatistics)^2))
  fit[i,3:5] = c(prediction, deviation)
}

plot(fit[fit[,5] < 1.5, 1:2], xlim = c(1,6), ylim = c(1,10), col = "lightgrey", main = "Accepted parameters for different values of epsilon")
points(fit[fit[,5] < 1, 1:2],  pch = 18, col = "gray")
points(fit[fit[,5] < 0.5, 1:2],  pch = 8, col = "red")

legend("topright", c("< 1.5", "< 1", "< 0.5"), pch = c(1,18,8), col = c("lightgrey", "gray", "red"))

abline(v = 2)
abline(h = 5) 