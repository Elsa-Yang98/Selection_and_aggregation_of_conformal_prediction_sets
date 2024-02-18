source('methods_functions.r')
df <- 3
d <- 3
l.lambda <- 100
lambda_seq <- seq(0,200,l=l.lambda)
nset <- c(50,100,500,1000,5000)
alpha <- 0.1 #miscoverage level
n0 <- 100  #number of prediction points
nrep <- 100 #number of independent trials
rho <- 0.5

evaluations <- expand.grid(1:nrep, nset, c("efficient", "valid"))
no_eval <- nrow(evaluations)
width_mat <- cov_mat <- data.frame(number = rep(0, no_eval), 
                                   rep = evaluations[,1], 
                                   nset = evaluations[,2],
                                   method = evaluations[,3])
colnames(width_mat) <- colnames(cov_mat) <- c("number", "rep", "sample size", "method")

Sigma <- matrix(rho,d,d)
diag(Sigma) <- rep(1,d) #covariance matrix for X
  

for(idx in 1:nrow(evaluations)){
  set.seed(evaluations[idx, 1])
  if(idx%%1 == 0){
    print(idx)
  }
  n <- evaluations[idx, 2]  
  
X <- rmvt(n,Sigma,df)	#multivariate t distribution
  
  eps1 <- rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
  eps2 <- rt(n,df)*(1+sqrt(X[,1]^4+X[,2]^4))
  Y <- rpois(n,sin(X[,1])^2 + cos(X[,2])^4+0.01 )+0.03*X[,1]*eps1+25*(runif(n,0,1)<0.01)*eps2
  
  X0 <- rmvt(n0,Sigma,df)
  eps01 <- rt(n0,df)*(1+sqrt(X0[,1]^2+X0[,2]^2))
  eps02 <- rt(n0,df)*(1+sqrt(X0[,1]^4+X0[,2]^4))
  Y0 <- rpois(n0,sin(X0[,1])^2 + cos(X0[,2])^4+0.01 )+0.03*X0[,1]*eps01+25*(runif(n0,0,1)<0.01)*eps02
  
  width_mat[idx,3] <- cov_mat[idx, 3] <- n
  method <- evaluations[idx, 3]
  width_mat[idx,4] <- cov_mat[idx, 4] <- method
  width_mat[idx, 2] <- cov_mat[idx, 2] <- evaluations[idx, 1]
  
    if(method == "valid"){
    split <- c(1/2, 1/2)
  } else {
    split <- 1/2
  }  
  
  beta_grid <- seq(1e-03, 4, length = 20)*alpha
  mtry_grid <- unique(ceiling(seq(1/10, 1, length = 20)*d))
  ntree_grid <- seq(50, 400, by = 50)
  
  tmp <- try(conf_CQR_reg(X, Y, split = split, beta_grid, mtry_grid, ntree_grid, method = method, alpha = alpha))
  
  while (class(tmp)=="try-error"){
    
    tmp <- try(conf_CQR_reg(X, Y, split = split, beta_grid, mtry_grid, ntree_grid, method = method, alpha = alpha),silent=TRUE)
    
  }
  width_mat[idx, 1] <- tmp$width
  cov_mat[idx, 1] <- mean(tmp$pred_set(X0, Y0))
}




width_efcp <- width_vfcp <- sd_width_efcp <- sd_width_vfcp <- NULL
for(n in nset){
  TMP <- width_mat[evaluations[,3] == "efficient", ]
  TMP_prime <- TMP[TMP[,3] == n,]
  
  TMP <- width_mat[evaluations[,3] == "valid", ]
  TMP_prime_vfcp <- TMP[TMP[,3] == n,]
  TMP_prime_vfcp_clean =TMP_prime_vfcp[ TMP_prime_vfcp[,1]<=10^5,1]
  
  width_efcp <- c(width_efcp, mean(TMP_prime[,1] / TMP_prime_vfcp[,1]))
  sd_width_efcp <- c(sd_width_efcp, sd(TMP_prime[,1]/ TMP_prime_vfcp[,1])/sqrt(nrep))

  width_vfcp <- c(width_vfcp, mean(TMP_prime_vfcp[,1] / TMP_prime_vfcp[,1]))
  sd_width_vfcp <- c(sd_width_vfcp, sd(TMP_prime_vfcp[,1]/ TMP_prime_vfcp[,1])/sqrt(nrep))
}


df.cov=data.frame(factor(nset),cov_efcp,cov_vfcp)
df.cov_sd=data.frame(factor(nset),sd_cov_efcp,sd_cov_vfcp)

names=c("EFCP","VFCP")
col_names = c("V1","V2","sd","Method")
`EFCP` <- cbind(df.cov[,c(1,2)],df.cov_sd[,2])
`EFCP`$"Method"=names[1]


`VFCP` <- cbind(df.cov[,c(1,3)],df.cov_sd[,3])
`VFCP`$"Method"=names[2]

force_bind = function(df1, df2) {
  colnames(df2) = colnames(df1)
  rbind(df1, df2)
}
cov = force_bind(`EFCP`, `VFCP`)
cov$Var="Coverage"


df.width=data.frame(factor(nset),width_efcp,width_vfcp)
df.width_sd=data.frame(factor(nset),sd_width_efcp,sd_width_vfcp)


`EFCP` <- cbind(df.width[,c(1,2)],df.width_sd[,2])
`EFCP`$"Method"=names[1]
`VFCP` <- cbind(df.width[,c(1,3)],df.width_sd[,3])
`VFCP`$"Method"=names[2]

width = force_bind(`EFCP`, `VFCP`)
width$Var="Width"
synthetic = force_bind(cov,width)
colnames(synthetic) <- col_names
synthetic$Name="Synthetic Data"

save(synthetic, file = 'synthetic_to_plot.RData')


