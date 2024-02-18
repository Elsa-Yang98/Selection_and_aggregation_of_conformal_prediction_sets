library("quantregForest")
library(MASS)
library(mvtnorm)
#########################################################################################
####################EFCP method##########################################################
#########################################################################################
efcp.fun=function(X,Y,X0,lambda=seq(0,100,length=100),alpha=0.1){
  n=dim(X)[1]
  d=dim(X)[2]
  n0=dim(X0)[1]
  D1=1:floor(n/2)
  D2=(floor(n/2)+1):n
X1=X[D1,]
X2=X[D2,]
Y1=Y[D1]
Y2=Y[D2]
n2=length(Y2)

out.gnet = glmnet(X1,Y1,alpha=0,lambda=lambda) 
lambda = out.gnet$lambda
RR2=matrix(0,n2,length(lambda))
YY2=matrix(rep(Y2,length(lambda)),nrow=length(Y2))
RR2=abs(YY2-cbind(rep(1,n2),X2)%*%coef(out.gnet))
q2=apply(RR2,2,quantile,probs=(1-alpha)*(1+1/n2))
index=which(q2==min(q2))[1]
#chosen.lambda=lambda[index]

R2=abs(Y2-cbind(rep(1,n2),X2)%*%coef(out.gnet)[,index])
q3=quantile(R2,probs=(1-alpha)*(1+1/n2))
up=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]+q3
lo=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]-q3

return(list(up=up,lo=lo))   
}

#########################################################################################
####################VFCP method##########################################################
#########################################################################################
vfcp.fun=function(X,Y,X0,lambda=seq(0,100,length=100),alpha=0.1){
  n=dim(X)[1]
  d=dim(X)[2]
  n0=dim(X0)[1]
  D1=1:floor(n/3)
  D2=(floor(n/3)+1):floor(2*n/3)
  D3=(floor(2*n/3)+1):n
  X1=X[D1,]
  X2=X[D2,]
  X3=X[D3,]
  Y1=Y[D1]
  Y2=Y[D2]
  Y3=Y[D3]
  n2=length(Y2)
  n3=length(Y3)

  out.gnet = glmnet(X1,Y1,alpha=0,lambda=lambda) 
  lambda = out.gnet$lambda
  RR2=matrix(0,n2,length(lambda))
  YY2=matrix(rep(Y2,length(lambda)),nrow=length(Y2))
  RR2=abs(YY2-cbind(rep(1,n2),X2)%*%coef(out.gnet))
  q2=apply(RR2,2,quantile,probs=(1-alpha)*(1+1/n2))
  index=which(q2==min(q2))[1]
  #chosen.lambda=lambda[index]
  
  R3=abs(Y3-cbind(rep(1,n3),X3)%*%coef(out.gnet)[,index])
  q3=quantile(R3,probs=(1-alpha)*(1+1/n3))
  up=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]+q3
  lo=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]-q3
  
  return(list(up=up,lo=lo))   
}



my.param.fun = function(x, y, x0) {
  n = nrow(x); n0 = nrow(x0)
  out = my.lm.funs$train(x,y)
  fit = matrix(my.lm.funs$predict(out,x),nrow=n)
  pred = matrix(my.lm.funs$predict(out,x0),nrow=n0)
  m = ncol(pred)
  
  x1 = cbind(rep(1,n0),x0)
  q = qt(1-alpha/2, n-d-1)
  lo = up = matrix(0,n0,m)
  
  for (j in 1:m) {
    sig.hat = sqrt(sum((y - fit[,j])^2)/(n-ncol(x1)))
    g = diag(x1 %*% chol.solve(out$chol.R[[j]], t(x1)))
    lo[,j] = pred[,j] - sqrt(1+g)*sig.hat*q
    up[,j] = pred[,j] + sqrt(1+g)*sig.hat*q
  }
  
  # Return proper outputs in proper formatting
  return(list(pred=pred,lo=lo,up=up,fit=fit))
}

#########################################################################################
####################Naive method#########################################################
#########################################################################################
naive.fun=function(X,Y,X0,alpha=0.1){
  x=X
  y=Y
  x0=X0
  n = nrow(x); n0 = nrow(x0)
  beta=ginv(t(x)%*%x)%*%t(x)%*%y
  R=abs(y-x%*%beta)
  q=quantile(R,probs=(1-alpha)*(1+1/n))
  
  up=X0%*%beta+q
  lo=X0%*%beta-q
  return(list(up=up,lo=lo))  
}


#########################################################################################
####################Linear method########################################################
#########################################################################################
ginverselm.funs=function (intercept = TRUE, lambda = 0) 
{
  
  m = length(lambda)
  for (j in 1:m) 
    train.fun = function(x, y, out = NULL) {
      n = nrow(x)
      p = ncol(x)
      v = rep(1, p)
      if (!is.null(out)) {
        chol.R = out$chol.R
      }
      else {
        chol.R = vector(mode = "list", length = m)
        for (j in 1:m) {
          chol.R[[j]] = crossprod(x) 
        }
      }
      beta = matrix(0, p, m)
      for (j in 1:m) {
        beta[, j] = ginv(chol.R[[j]])%*%t(x) %*% y
      }
      return(list(beta = beta, chol.R = chol.R))
    }
  predict.fun = function(out, newx) {
    
    return(newx %*% out$beta)
  }
  special.fun = function(x, y, out) {
    n = nrow(x)
    p = ncol(x)
    
    res = y - x %*% out$beta
    for (j in 1:m) {
      s = diag(x %*% ginv(out$chol.R[[j]])%*% t(x))
      res[, j] = res[, j]/(1 - s)
    }
    return(res)
  }
  active.fun = function(out) {
    p =  nrow(out$beta)
    m = ncol(out$beta)
    act.list = vector(length = m, mode = "list")
    for (i in 1:m) act.list[[i]] = 1:p
    return(act.list)
  }
  return(list(train.fun = train.fun, predict.fun = predict.fun, 
              special.fun = special.fun, active.fun = active.fun))
}


my.ginverselm.funs = ginverselm.funs(lambda=0)

ginverse.fun = function(x, y, x0,alpha=0.1) {
  n = nrow(x); n0 = nrow(x0)
  out = my.ginverselm.funs$train(x,y)
  fit = matrix(my.ginverselm.funs$predict(out,x),nrow=n)
  pred = matrix(my.ginverselm.funs$predict(out,x0),nrow=n0)
  m = ncol(pred)
  
  x1 = x0
  #q = qt(1-alpha/2, n-d)
  q = qnorm(1-alpha/2)
  lo = up = matrix(0,n0,m)
  
  for (j in 1:m) {
    #sig.hat = sqrt(sum((y - fit[,j])^2)/(n-ncol(x1)))
    sig.hat = sqrt(sum((y - fit[,j])^2)/(n))
    g = diag(x1 %*% ginv(out$chol.R[[j]])%*%t(x1))
    lo[,j] = pred[,j] - sqrt(1+g)*sig.hat*q
    up[,j] = pred[,j] + sqrt(1+g)*sig.hat*q
  }
  
  # Return proper outputs in proper formatting
  return(list(pred=pred,lo=lo,up=up,fit=fit))
}
#########################################################################################
####################All methods wrapper example##########################################
df=3 # degrees of freedom
l=60    #number of dimensions 
lambda_seq=seq(0,200,l=100)
dim=round(seq(5,300,l=l))
alpha=0.1
n=200   #number of training samples
n0=100  #number of prediction points
nrep=100 #number of independent trials
rho=0.5
#######################################################################
four_methods_comparison = function(nrep=100, l=60, lambda_seq=seq(0,200,l=100), dim=round(seq(5,300,l=l)), n=200, n0=100, rho=0.5, alpha=0.1,setting = "linear", df=3){
  cov.nested=len.nested=cov.our=len.our=cov.naive=len.naive=cov.param=len.param=matrix(0,nrep,l)
  out.nested.up=out.nested.lo=out.our.up=out.our.lo=out.naive.up=out.our.lo=out.naive.lo=out.param.up=out.param.lo=matrix(0,n0,l)
  for(i in 1:nrep){
    cat(i,"\n")
    for (r in 1:l){
      d=dim[r]
      set.seed(i)
      
      Sigma=matrix(rho,d,d)
      diag(Sigma)=rep(1,d)
      X=rmvt(n,Sigma,df)	#multivariate t distribution
      X0=rmvt(n0,Sigma,df)
      if (setting=="nonlinear"){
      eps1=rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
      eps2=rt(n,df)*(1+sqrt(X[,1]^4+X[,2]^4))
      Y=rpois(n,sin(X[,1])^2 + cos(X[,2])^4+0.01 )+0.03*X[,1]*eps1+25*(runif(n,0,1)<0.01*eps2)
      eps01=rt(n0,df)*(1+sqrt(X0[,1]^2+X0[,2]^2))
      eps02=rt(n0,df)*(1+sqrt(X0[,1]^4+X0[,2]^4))
      Y0=rpois(n0,sin(X0[,1])^2 + cos(X0[,2])^4+0.01 )+0.03*X0[,1]*eps01+25*(runif(n0,0,1)<0.01)*eps02
      } else if(setting=="linear"){
        beta=rep(1:5,d/5)
        eps=rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
        Y=X%*%beta+eps
        eps0=rt(n0,df)*(1+sqrt(X0[,1]^2+X0[,2]^2))
        Y0=X0%*%beta+eps0
      } else{
        print("error in data setting")
        break
      }
      ##Linear##
      out.param=ginverse.fun(X,Y,X0,alpha=alpha)
      out.param.lo[,r]=out.param$lo
      out.param.up[,r]=out.param$up
      cov.param[i,r]= mean(out.param.lo[,r] <= Y0 & Y0 <= out.param.up[,r]) 
      len.param[i,r]=mean(out.param.up[,r]-out.param.lo[,r]) 
      ##VFCP##
      out.our=vfcp.fun(X,Y,X0,lambda=lambda_seq,alpha=alpha)
      out.vfcp.up[,r]=out.vfcp$up
      out.vfcp.lo[,r]=out.vfcp$lo
      cov.vfcp[i,r]= mean(out.vfcp.lo[,r] <= Y0 & Y0 <= out.vfcp.up[,r]) 
      len.vfcp[i,r]=mean(out.vfcp.up[,r]-out.vfcp.lo[,r]) 
      ##naive##
      out.naive=naive.fun(X,Y,X0,alpha=alpha)
      out.naive.up[,r]=out.naive$up
      out.naive.lo[,r]=out.naive$lo
      cov.naive[i,r]= mean(out.naive.lo[,r] <= Y0 & Y0 <= out.naive.up[,r]) 
      len.naive[i,r]=mean(out.naive.up[,r]-out.naive.lo[,r]) 
      ##EFCP
      out.nested=nested.fun(X,Y,X0,lambda=lambda_seq,alpha=alpha)
      out.nested.up[,r]=out.nested$up
      out.nested.lo[,r]=out.nested$lo
      cov.nested[i,r]= mean(out.nested.lo[,r] <= Y0 & Y0 <= out.nested.up[,r]) 
      len.nested[i,r]=mean(out.nested.up[,r]-out.nested.lo[,r]) 
    }
  }
  return(list(dim,cov.param, cov.naive, cov.vfcp,cov.nested,len.param, len.naive, len.vfcp,len.nested))
}


conf_CQR_prelim <- function(X1, Y1, X2, Y2, beta_grid, mtry, ntree, alpha = 0.1){
  n2 <- length(Y2)
  quant_reg_fit <- quantregForest(X1, Y1, ntree = ntree, mtry = mtry)
  qhat_beta_grid <- predict(quant_reg_fit, newdata = X2, what = c(beta_grid, 1/2, 1 - beta_grid))
  colnames(qhat_beta_grid) <- paste0("beta",c(beta_grid, 1/2, 1-beta_grid))
  width_beta <- colMeans(abs(qhat_beta_grid[,paste0("beta",1-beta_grid)] - qhat_beta_grid[,paste0("beta", beta_grid)]),na.rm=TRUE)
  width_beta <- unname(width_beta)
  #### CQR score computation
  cqr_left_score <- qhat_beta_grid[,paste0("beta",beta_grid)] - Y2
  cqr_right_score <- Y2 - qhat_beta_grid[,paste0("beta",1 - beta_grid)]
  cqr_score <- pmax(cqr_left_score, cqr_right_score)
  ### CQR score computation complete
  ### CQR-m score computation
  cqr_m_left_denom <- qhat_beta_grid[,paste0("beta",1/2)] - qhat_beta_grid[,paste0("beta",beta_grid)]
  cqr_m_left_denom <- abs(cqr_m_left_denom) + 1e-08
  cqr_m_left <- (qhat_beta_grid[,paste0("beta",beta_grid)] - Y2)/cqr_m_left_denom
  cqr_m_right_denom <- qhat_beta_grid[,paste0("beta",1-beta_grid)] - qhat_beta_grid[,paste0("beta",1/2)]
  cqr_m_right_denom <- abs(cqr_m_right_denom) + 1e-08
  cqr_m_right <- (Y2 - qhat_beta_grid[,paste0("beta",1-beta_grid)] )/cqr_m_right_denom
  cqr_m_score <- pmax(cqr_m_left, cqr_m_right)
  ### CQR-m score computation complete
  ### CQR-r score computation
  cqr_r_denom <- qhat_beta_grid[,paste0("beta",1-beta_grid)] - qhat_beta_grid[,paste0("beta",beta_grid)]
  cqr_r_denom <- abs(cqr_r_denom) + 1e-08
  cqr_r_left <- (qhat_beta_grid[,paste0("beta",beta_grid)] - Y2)/cqr_r_denom
  cqr_r_right <- (Y2 - qhat_beta_grid[,paste0("beta",1-beta_grid)] )/cqr_r_denom
  cqr_r_score <- pmax(cqr_r_left, cqr_r_right)
  ### CQR-r score computation complete
  ### Computation of quantiles for each method
  quant <- width <- matrix(0, nrow = 3, ncol = length(beta_grid))
  rownames(quant) <- rownames(width) <- c("CQR", "CQR-m", "CQR-r")
  colnames(quant) <- colnames(width) <- paste0("beta", beta_grid)
  quant["CQR",] <- apply(cqr_score, 2, quantile, probs = (1 - alpha)*(1 + 1/n2), na.rm = TRUE)
  quant["CQR-m",] <- apply(cqr_m_score, 2, quantile, probs = (1-alpha)*(1 + 1/n2),  na.rm = TRUE)
  quant["CQR-r",] <- apply(cqr_r_score, 2, quantile, probs = (1-alpha)*(1 + 1/n2), na.rm = TRUE)
  ### Computing the width for each method
  width["CQR",] <- 2*quant["CQR",] + width_beta
  width["CQR-m",] <- (1 + quant["CQR-m",])*width_beta
  width["CQR-r",] <- (1 + 2*quant["CQR-r",])*width_beta
  opt_ind <- arrayInd(which.min(width), dim(width))
  ret <- list(call= match.call())
  ret$opt_threshold <- unname(quant[opt_ind])
  ret$cqr_method <- rownames(width)[opt_ind[,1]]
  ret$beta <- unname(beta_grid[opt_ind[,2]])  
  ret$width <- unname(min(width))
  ret$width_beta <- unname(width_beta[opt_ind[,2]])
  ret$mtry <- mtry
  ret$ntree <- ntree
  ret$forest <- quant_reg_fit
  return(ret)
}

conf_CQR_reg <- function(x, y, split, beta_grid, mtry_grid, ntree_grid, method = "efficient", alpha = 0.1){
  y_data <- as.vector(y)
  x_data <- as.matrix(x)
  N <- length(y_data)
  if(method == "efficient"){
    ## If method is efficient there is only one split required.
    ## Check is split is a vector of length 1.
    stopifnot(length(split) == 1)
    
    ## Divide the sample into two parts randomly
    n1 <- ceiling(N*split[1])
    n2 <- N - n1
    I1 <- sample.int(N, n1, replace = FALSE)
    X1 <- x_data[I1,]
    X1 = as.matrix(X1)
    Y1 <- y_data[I1]
    X2 <- x_data[-I1,]
    X2 = as.matrix(X2)
    Y2 <- y_data[-I1]
  } else if(method == "valid"){
    stopifnot(length(split) == 2)
    n1 <- ceiling(N*split[1])
    n2 <- ceiling({N - n1}*split[2])
    n3 <- N - n1 - n2
    I1 <- sample.int(N, n1, replace = FALSE)
    I2 <- sample.int(N-n1, n2, replace = FALSE)
    
    X1 <- x_data[I1,]
    X1 = as.matrix(X1)
    Y1 <- y_data[I1]
    Xtmp <- x_data[-I1,]
    Xtmp <- as.matrix(Xtmp)
    Ytmp <- y_data[-I1]
    X2 <- Xtmp[I2,]
    X2 = as.matrix(X2)
    Y2 <- Ytmp[I2]
    X3 <- Xtmp[-I2,]
    X3 = as.matrix(X3)
    Y3 <- Ytmp[-I2]
  }
  mtry_ntree_grid <- expand.grid(mtry_grid, ntree_grid)
  num_params <- nrow(mtry_ntree_grid)
  opt_width <- diff(range(y_data))
  for(idx in 1:num_params){
    mtry <- mtry_ntree_grid[idx, 1]
    ntree <- mtry_ntree_grid[idx, 2]
    ret_mtry_ntree <- conf_CQR_prelim(X1, Y1, X2, Y2, beta_grid = beta_grid, mtry = mtry, ntree = ntree, alpha = alpha)
    if(ret_mtry_ntree$width <= opt_width){
      ret <- ret_mtry_ntree
      opt_width <- ret_mtry_ntree$width
    }
  }
  if(method == "valid"){
    forest <- ret$forest
    beta <- ret$beta
    cqr_method <- ret$cqr_method
    predictions_D3 <- predict(forest, newdata = X3, what = c(beta, 1/2, 1 - beta))
    if(cqr_method == "CQR"){
      ret$opt_threshold <- quantile(pmax(predictions_D3[,1] - Y3, Y3 - predictions_D3[,3]), probs = (1-alpha)*(1+1/n3))
      ret$width <- 2*ret$opt_threshold + ret$width_beta    
    } else if(cqr_method == "CQR-m"){
      low <- (predictions_D3[,1] - Y3)/(abs(predictions_D3[,2] - predictions_D3[,1]) + 1e-08)
      up <- (Y3 - predictions_D3[,3])/(abs(predictions_D3[,3] - predictions_D3[,2]) + 1e-08)
      ret$opt_threshold <- quantile(pmax(low, up), probs = (1-alpha)*(1 + 1/n3))
      ret$width <- (1 + ret$opt_threshold)*ret$width_beta
    } else if(cqr_method == "CQR-r"){
      denom <- abs(predictions_D3[,3] - predictions_D3[,1]) + 1e-08
      low <- (predictions_D3[,1] - Y3)/denom
      up <- (Y3 - predictions_D3[,3])/denom
      ret$opt_threshold <- quantile(pmax(low, up), probs = (1-alpha)*(1 + 1/n3))
      ret$width <- (1 + 2*ret$opt_threshold)*ret$width_beta
    }
  }
  if(ret$cqr_method == "CQR"){
    forest <- ret$forest
    beta <- ret$beta
    quant <- ret$opt_threshold
    pred_set_verify <- function(xx, yy){
      qhat_final <- predict(forest, newdata = xx, what = c(beta, 1 - beta))
      return(pmax(qhat_final[,1] - yy, yy - qhat_final[,2]) <= quant)
    }      
    ret$pred_set <- pred_set_verify
  } else if(ret$cqr_method == "CQR-m"){
    forest <- ret$forest
    beta <- ret$beta
    quant <- ret$opt_threshold
    pred_set_verify <- function(xx, yy){
      qhat_final <- predict(forest, newdata = xx, what = c(beta, 1/2, 1 - beta))
      low <- (qhat_final[,1] - yy)/(abs(qhat_final[,2] - qhat_final[,1]) + 1e-08)
      up <- (yy - qhat_final[,3])/(abs(qhat_final[,3] - qhat_final[,2]) + 1e-08)
      return(pmax(low, up) <= quant)
    }      
    ret$pred_set <- pred_set_verify    
  } else if(ret$cqr_method == "CQR-r"){
    forest <- ret$forest
    beta <- ret$beta
    quant <- ret$opt_threshold
    pred_set_verify <- function(xx, yy){
      qhat_final <- predict(forest, newdata = xx, what = c(beta, 1 - beta))
      low <- (qhat_final[,1] - yy)/(abs(qhat_final[,2] - qhat_final[,1]) + 1e-08)
      up <- (yy - qhat_final[,2])/(abs(qhat_final[,2] - qhat_final[,1]) + 1e-08)
      return(pmax(low, up) <= quant)
    }      
    ret$pred_set <- pred_set_verify    
  }
  ret$method <- method
  return(ret)
}







