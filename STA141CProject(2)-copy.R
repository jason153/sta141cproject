library(dplyr)
library(tidyr)
library(Matrix)
data <- read.csv("train.csv", header=TRUE)
data1 =cbind(1, data$LotArea, data$TotalBsmtSF, data$X1stFlrSF,data$X2ndFlrSF, data$GarageArea, data$SalePrice)
data1

X <-  data1[, 1:6] # X and y
X
y <- as.matrix(data1[,7])
y

start_time <- Sys.time()
XX <- crossprod(X)
XX.lu <- expand(lu(XX))
L <- XX.lu$L
U <- XX.lu$U
P <- XX.lu$P
bs <- forwardsolve(L, crossprod(P, crossprod(X, y)))
betahat <- backsolve(U, bs)
betahat

sigma2hat <- crossprod(y - X %*% betahat) / (dim(X)[1] - 6)
sigma2hat

#betacov <- sigma2hat * solve(XX)  

betacov <- as.numeric(sigma2hat) * chol2inv(XX)
betase <- sqrt(diag(betacov))
betase
end_time <- Sys.time()
end_time - start_time

#Cholesky
start_time <- Sys.time()
yXXy <- crossprod(data1)
yXXy
# Cholesky decomposition of Grammian matrix
cholR <- chol(yXXy)
# regression coefficients
betahat <- backsolve(cholR[1:6, 1:6], cholR[1:6, 7])
betahat
# s.e. of regression coefficients
sigma2hat <- cholR[7, 7]^2 / (dim(data1)[1] - 6)
betacov <- sigma2hat * chol2inv(cholR[1:6, 1:6])
betase <- sqrt(diag(betacov))
betase
end_time <- Sys.time()
end_time - start_time

#QR
start_time <- Sys.time()
x.qr <- qr(data1[, 1:6])
x.qr$rank
x.qr$pivot
# regression coefficients
betahat <- forwardsolve(x.qr$qr, yXXy[order(x.qr$pivot), 7],
                        upper.tri = TRUE, transpose = TRUE)
betahat <- backsolve(x.qr$qr, betahat)[x.qr$pivot]
betahat
# s.e. of regression coefficients
sigma2hat <- crossprod(Xy[, 7] - Xy[, 1:6] %*% betahat) / (dim(Xy)[1] - 6)
betacov <- as.numeric(sigma2hat) * chol2inv(x.qr$qr)
betase <- sqrt(diag(betacov))
betase <- betase[x.qr$pivot]
betase

end_time <- Sys.time()
end_time - start_time

#LU Time difference of 0.006774902 secs
#Cho Time difference of 0.004496098 secs
#QR Time difference of 0.005033016 secs

testdata <- read.csv("Desktop/STA141Cproject/test.csv",sep = ",")
data2 =cbind(1, testdata$LotArea, testdata$TotalBsmtSF, testdata$X1stFlrSF,testdata$X2ndFlrSF, testdata$GarageArea)
data2
testX <-  data2[, 1:6] # X and y
testX <- data.frame(testX)

model = lm(y~X)

#predicts the future values
predict(model,newdata = testX)

#AIC
nrow(data1)*(log(2*pi)+1+log((sum(model$residuals^2)/nrow(data1))))+((length(model$coefficients)+1)*2)
AIC(model)

