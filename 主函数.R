error <- data.frame(iter=c(0:10))

# library(foreach)
# library(parallel)
# library(doParallel)
#system.time({
  # 注册核心集群
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  # 启动并行计算
  res <- foreach(k=1:1000,.combine = 'cbind')%dopar%{
    # Suppose X1,X2,X3 follows multivariate uniform distribution with correlation rho_jk = 0.5^|j-k|, epsilon follows N(0,1)
    # Generate data, first generate multivariate normal, then use F() transform to Uniform(0,1)
    # Define correlation matrix
    rho <- 0.5
    corr_matrix <- matrix(0, nrow=3, ncol=3)
    for(i in 1:3){
      for(j in 1:3){
        corr_matrix[i,j] <- rho^abs(i-j)
      }
    }
    # Generate multivariate sample
    library(MASS)
    n <- 10000 # internal sample size
    N <- 10000 # external sample size
    normal_samples_ex <- mvrnorm(n, mu=rep(0,3), Sigma=corr_matrix) # external data
    normal_samples_in <- mvrnorm(n, mu=rep(0,3), Sigma=corr_matrix) # internal data
    # Transform to uniform(0,1)
    uniform_samples_ex <- pnorm(normal_samples_ex)
    uniform_samples_in <- pnorm(normal_samples_in)
    # Regression coefficient
    library(quantreg)
    tau <- 0.25
    beta0 <- 1
    beta1 <- 1
    beta2 <- 1
    beta3 <- 1
    p <- 3
    beta.real <- c(beta0+qnorm(tau),beta1,beta2,beta3)
    # Get summary information from external data
    x1.ex <- uniform_samples_ex[,1]
    x2.ex <- uniform_samples_ex[,2]
    x3.ex <- uniform_samples_ex[,3]
    epsilon.ex <- rnorm(n)
    y.ex <- cbind(rep(1,n),x1.ex,x2.ex,x3.ex) %*% rbind(beta0,beta1,beta2,beta3) + epsilon.ex
    data.ex <- data.frame(x1=x1.ex,x2=x2.ex,x3=x3.ex,y=y.ex)
    summary.info <- as.matrix(coef(rq(y~.,tau=tau,data=data.ex)))
    # Inference based on internal data
    x1.in <- uniform_samples_in[,1]
    x2.in <- uniform_samples_in[,2]
    x3.in <- uniform_samples_in[,3]
    epsilon.in <- rnorm(n)
    y.in <- cbind(rep(1,n),x1.in,x2.in,x3.in) %*% rbind(beta0,beta1,beta2,beta3) + epsilon.in
    data.in <- data.frame(x1=x1.in,x2=x2.in,x3=x3.in,y=y.in)
    beta.single <- as.matrix(coef(rq(y~.,tau=tau,data=data.in)))
    # Internal data
    X <- cbind(rep(1,n),as.matrix(data.in[-dim(data.in)[2]]))
    Y <- as.matrix(data.in[dim(data.in)[2]])
    
    # Define smooth function
    H <- function(x){
      if (x<=0) {return(0)}
      else if (x>=1) {return(1)}
      else {return(1/2+15/16*(x-2/3*x^3+1/5*x^5))}
    }
    # Define randomness of summary info
    V <- diag(1,p+1)
    # Define bandwidth
    h <- 1/5*n^(-1/5)
    # Define iteration function
    iter <- function(beta_0,h){
      # Initialize U_k and V_k
      U_k <- numeric(length = ncol(X))
      V_k <- matrix(0, ncol = ncol(X), nrow = ncol(X))
      library(MASS)
      library(numDeriv)
      # Calculate U_k
      for (i in 1:n) {
        temp <- c((Y[i] - X[i,] %*% beta_0) / h)
        U_k <- U_k + X[i,] * (H(temp) + tau-1 + Y[i] / h * grad(H,temp))
      }
      # Calculate V_k
      for (i in 1:n) {
        temp <- c((Y[i] - X[i,] %*% beta_0) / h)
        V_k <- V_k + X[i,] %*% t(X[i,]) * (1/h * grad(H,temp))
      }
      return(ginv(V_k)%*%U_k)
    }
    # Iteration
    beta.hat <- beta.single
    beta.process <- as.matrix(beta.hat,nrow=4) # 记录迭代的beta值
    error.process <- data.frame(iter=0,error=sum((beta.hat-beta.real)^2)) # 记录误差
    max_iter <- 10
    tol <- 1e-5
    for (i in 1:max_iter) {
      if (sum((beta.hat-beta.real)^2)<tol){
        beta.process <- cbind(beta.process,beta.hat)
        error.process <- rbind(error.process,c(i,sum((beta.hat-beta.real)^2)))
        break
      }
      else{
        beta.hat <- iter(beta_0=beta.hat,h=h)
        beta.process <- cbind(beta.process,beta.hat)
        error.process <- rbind(error.process,c(i,sum((beta.hat-beta.real)^2)))
      }
    }
    error[paste0('error',k)] <- error.process$error
  }
  stopCluster(cl)
})
#res <- as.data.frame(res)
#res$mean <- apply(res[,2:dim(error)[2]],1,mean)







error <- data.frame(iter=c(0:10))
for (k in 1:10) {
  # Suppose X1,X2,X3 follows multivariate uniform distribution with correlation rho_jk = 0.5^|j-k|, epsilon follows N(0,1)
  # Generate data, first generate multivariate normal, then use F() transform to Uniform(0,1)
  # Define correlation matrix
  rho <- 0.5
  corr_matrix <- matrix(0, nrow=3, ncol=3)
  for(i in 1:3){
    for(j in 1:3){
      corr_matrix[i,j] <- rho^abs(i-j)
    }
  }
  # Generate multivariate sample
  library(MASS)
  n <- 10000 # internal sample size
  N <- 10000 # external sample size
  normal_samples_ex <- mvrnorm(n, mu=rep(0,3), Sigma=corr_matrix) # external data
  normal_samples_in <- mvrnorm(n, mu=rep(0,3), Sigma=corr_matrix) # internal data
  # Transform to uniform(0,1)
  uniform_samples_ex <- pnorm(normal_samples_ex)
  uniform_samples_in <- pnorm(normal_samples_in)
  # Regression coefficient
  library(quantreg)
  tau <- 0.25
  beta0 <- 1
  beta1 <- 1
  beta2 <- 1
  beta3 <- 1
  p <- 3
  beta.real <- c(beta0+qnorm(tau),beta1,beta2,beta3)
  # Get summary information from external data
  x1.ex <- uniform_samples_ex[,1]
  x2.ex <- uniform_samples_ex[,2]
  x3.ex <- uniform_samples_ex[,3]
  epsilon.ex <- rnorm(n)
  y.ex <- cbind(rep(1,n),x1.ex,x2.ex,x3.ex) %*% rbind(beta0,beta1,beta2,beta3) + epsilon.ex
  data.ex <- data.frame(x1=x1.ex,x2=x2.ex,x3=x3.ex,y=y.ex)
  summary.info <- as.matrix(coef(rq(y~.,tau=tau,data=data.ex)))
  # Inference based on internal data
  x1.in <- uniform_samples_in[,1]
  x2.in <- uniform_samples_in[,2]
  x3.in <- uniform_samples_in[,3]
  epsilon.in <- rnorm(n)
  y.in <- cbind(rep(1,n),x1.in,x2.in,x3.in) %*% rbind(beta0,beta1,beta2,beta3) + epsilon.in
  data.in <- data.frame(x1=x1.in,x2=x2.in,x3=x3.in,y=y.in)
  beta.single <- as.matrix(coef(rq(y~.,tau=tau,data=data.in)))
  # Internal data
  X <- cbind(rep(1,n),as.matrix(data.in[-dim(data.in)[2]]))
  Y <- as.matrix(data.in[dim(data.in)[2]])
  
  # Define smooth function
  H <- function(x){
    if (x<=0) {return(0)}
    else if (x>=1) {return(1)}
    else {return(1/2+15/16*(x-2/3*x^3+1/5*x^5))}
  }
  # Define randomness of summary info
  V <- diag(1,p+1)
  # Define bandwidth
  h <- 1/5*n^(-1/5)
  # Define iteration function
  iter <- function(beta_0,h){
    C <- 0.5*n*N*(ginv(V)+t(ginv(V)))
    # Initialize U_k and V_k
    U_k <- numeric(length = ncol(X))
    V_k <- matrix(0, ncol = ncol(X), nrow = ncol(X))
    library(MASS)
    library(numDeriv)
    # Calculate U_k
    for (i in 1:n) {
      temp <- c((Y[i] - X[i,] %*% beta_0) / h)
      U_k <- U_k + X[i,] * (H(temp) + tau-1 + Y[i] / h * grad(H,temp))
    }
    # Calculate V_k
    for (i in 1:n) {
      temp <- c((Y[i] - X[i,] %*% beta_0) / h)
      V_k <- V_k + X[i,] %*% t(X[i,]) * (1/h * grad(H,temp))
    }
    return(ginv(V_k)%*%(U_k+C%*%(beta_0-summary.info)))
  }
  # Iteration
  beta.hat <- beta.single
  beta.process <- as.matrix(beta.hat,nrow=4) # 记录迭代的beta值
  error.process <- data.frame(iter=0,error=sum((beta.hat-beta.real)^2)) # 记录误差
  max_iter <- 10
  tol <- 1e-5
  for (i in 1:max_iter) {
    if (sum((beta.hat-beta.real)^2)<tol){
      beta.process <- cbind(beta.process,beta.hat)
      error.process <- rbind(error.process,c(i,sum((beta.hat-beta.real)^2)))
      break
    }
    else{
      beta.hat <- iter(beta_0=beta.hat,h=h)
      beta.process <- cbind(beta.process,beta.hat)
      error.process <- rbind(error.process,c(i,sum((beta.hat-beta.real)^2)))
    }
  }
  print(k)
  error[paste0('error',k)] <- error.process$error
}
error$mean <- apply(error[,2:dim(error)[2]],1,mean)


