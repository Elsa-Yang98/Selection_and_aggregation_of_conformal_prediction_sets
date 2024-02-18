source("methods_functions.R")

alpha=0.1
nrep=100
dim = c("blog","protein", "concrete","news", "kernel","superconduct")

evaluations <- expand.grid(1:nrep, dim, c("efficient", "valid"))
no_eval <- nrow(evaluations)
width_mat <- cov_mat <- data.frame(number = rep(0, no_eval), 
                                   rep = evaluations[,1], 
                                   dim = evaluations[,2],
                                   method = evaluations[,3])
colnames(width_mat) <- colnames(cov_mat) <- c("number", "rep", "dim", "method")
for(idx in 1:nrow(evaluations)){
  set.seed(evaluations[idx, 1])
  if(idx%%1 == 0){
    print(idx)
  }
  d <- evaluations[idx, 2]  
  if (d =="protein"){
    protein=as.matrix(read.csv("real/CASP.csv",header=TRUE))
    X_full=protein[,2:9]
    Y_full=protein[,1]
  }else if (d=="blog"){
    blog=as.matrix(read.csv("real/blogData_train.csv",header=FALSE))
    X_full=blog[,1:280]
    Y_full=blog[,281]
  }else if (d=="concrete"){
    concrete=as.matrix(read.csv("real/concrete.csv",header=FALSE,sep=","))
    X_full=concrete[,1:8]
    Y_full=concrete[,9]
  }else if (d=="news"){
    news=as.matrix(read.csv("real/OnlineNewsPopularity.csv",header=FALSE))
    X_full=news[,1:59]
    Y_full=news[,60]
  }else if (d=="kernel"){
    kernel=as.matrix(read.csv("real/sgemm_product.csv",header=FALSE))
    X_full=kernel[,1:14]
    Y_full=apply(kernel[,15:18],1,mean)
  }else if (d=="superconduct"){
    superconduct=as.matrix(read.csv("real/train.csv",header=FALSE))
    X_full=superconduct[,1:81]
    Y_full=superconduct[,82]
  }
  
  subSampleSize = 1000; 
  rows=sample(nrow(X_full))
  X_sample = X_full[rows[1:subSampleSize],]
  Y_sample = Y_full[rows[1:subSampleSize]]
  X=X_sample[1:768,]
  X0=X_sample[769:1000,]
  Y=Y_sample[1:768]
  Y0=Y_sample[769:1000]
  
  width_mat[idx,3] <- cov_mat[idx, 3] <- d
  method <- evaluations[idx, 3]
  width_mat[idx,4] <- cov_mat[idx, 4] <- method
  width_mat[idx, 2] <- cov_mat[idx, 2] <- evaluations[idx, 1]
  
  if(method == "valid"){
    split <- c(1/2, 1/2)
  } else {
    split <- 1/2
  }  
  
  p = dim(X0)[2]
  beta_grid <- seq(1e-03, 4, length = 10)*alpha
  mtry_grid <- unique(ceiling(seq(1/10, 1, length = 10)*p))
  ntree_grid <- seq(100, 400, by = 100)
  
  tmp =   try(conf_CQR_reg(X, Y, split = split, beta_grid, mtry_grid, ntree_grid, method = method, alpha = alpha))
  
  while (class(tmp)=="try-error"){
    
    tmp =   try(conf_CQR_reg(X, Y, split = split, beta_grid, mtry_grid, ntree_grid, method = method, alpha = alpha),silent=TRUE)
    
  }
  width_mat[idx, 1] <- tmp$width
  cov_mat[idx, 1] <- mean(tmp$pred_set(X0, Y0))
}

width_efcp <- width_vfcp <- sd_width_efcp <- sd_width_vfcp <- NULL
for(d in dim){
  TMP <- width_mat[evaluations[,3] == "efficient", ]
  TMP_prime <- TMP[TMP[,3] == d,]
  
  TMP <- width_mat[evaluations[,3] == "valid", ]
  TMP_prime_vfcp <- TMP[TMP[,3] == d,]
  
  
  width_efcp <- c(width_efcp, mean(TMP_prime[,1] / TMP_prime_vfcp[,1]))
  sd_width_efcp <- c(sd_width_efcp, sd(TMP_prime[,1]/ TMP_prime_vfcp[,1])/sqrt(nrep))
  

  width_vfcp <- c(width_vfcp, mean(TMP_prime_vfcp[,1] / TMP_prime_vfcp[,1]))
  sd_width_vfcp <- c(sd_width_vfcp, sd(TMP_prime_vfcp[,1]/ TMP_prime_vfcp[,1])/sqrt(nrep))
}

cov_efcp <- cov_vfcp <-sd_cov_efcp <- sd_cov_vfcp <- NULL
for(d in dim){
  TMP <- cov_mat[evaluations[,3] == "efficient", ]
  TMP_prime <- TMP[TMP[,3] == d,]
  cov_efcp <- c(cov_efcp, mean(TMP_prime[,1]))
  sd_cov_efcp <- c(sd_cov_efcp, sd(TMP_prime[,1])/sqrt(nrep))
  
  TMP <- cov_mat[evaluations[,3] == "valid", ]
  TMP_prime <- TMP[TMP[,3] == d,]
  cov_vfcp <- c(cov_vfcp, mean(TMP_prime[,1]))
  sd_cov_vfcp <- c(sd_cov_vfcp, sd(TMP_prime[,1])/sqrt(nrep))
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

real = rbind(cov,width)
colnames(real) <- col_names
real$Name = "Real Data"

save(real, file = 'real_to_plot.RData')

