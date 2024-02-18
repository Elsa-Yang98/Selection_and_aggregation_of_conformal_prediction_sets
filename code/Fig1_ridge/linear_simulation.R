library(glmnet,MASS,mvtnorm, gsubfn)
source("methods_functions.R")
df = 3
list(dim,cov.param, cov.naive, cov.vfcp,cov.nested,len.param, len.naive, len.vfcp,len.nested) = four_methods_comparison(nrep=100, l=60, lambda_seq=seq(0,200,l=100), dim=round(seq(5,300,l=l)), n=200, n0=100, rho=0.5, alpha=0.1, setting="linear", df=3)

n=nrow(len.param)
len.param = len.param/len.vfcp
len.naive = len.naive/len.vfcp
len.nested = len.nested/len.vfcp
len.vfcp = len.vfcp/len.vfcp
df.cov3=data.frame(dim,apply(cov.param,2,mean),apply(cov.naive,2,mean),apply(cov.vfcp,2,mean),apply(cov.nested,2,mean))
df.cov3_sd=data.frame(dim,apply(cov.param,2,sd),apply(cov.naive,2,sd),apply(cov.vfcp,2,sd),apply(cov.nested,2,sd))/sqrt(n)
df.len3=data.frame(dim,apply(len.param,2,mean),apply(len.naive,2,mean),apply(len.vfcp,2,mean),apply(len.star,2,mean),apply(len.cv5,2,mean), apply(len.nested,2,mean))
df.len3_sd=data.frame(dim,apply(len.param,2,sd),apply(len.naive,2,sd),apply(len.vfcp,2,sd),apply(len.star,2,sd),apply(len.cv5,2,sd), apply(len.nested,2,sd))/sqrt(n)


df = 5
list(dim,cov.param, cov.naive, cov.vfcp,cov.nested,len.param, len.naive, len.vfcp,len.nested) = four_methods_comparison(nrep=100, l=60, lambda_seq=seq(0,200,l=100), dim=round(seq(5,300,l=l)), n=200, n0=100, rho=0.5, alpha=0.1, setting="linear", df=5)
len.param = len.param/len.vfcp
len.naive = len.naive/len.vfcp
len.nested = len.nested/len.vfcp
len.vfcp = len.vfcp/len.vfcp
df.cov5=data.frame(dim,apply(cov.param,2,mean),apply(cov.naive,2,mean),apply(cov.vfcp,2,mean),apply(cov.nested,2,mean))
df.cov5_sd=data.frame(dim,apply(cov.param,2,sd),apply(cov.naive,2,sd),apply(cov.vfcp,2,sd),apply(cov.nested,2,sd))/sqrt(n)
df.len5=data.frame(dim,apply(len.param,2,mean),apply(len.naive,2,mean),apply(len.vfcp,2,mean),apply(len.nested,2,mean))
df.len5_sd=data.frame(dim,apply(len.param,2,sd),apply(len.naive,2,sd),apply(len.vfcp,2,sd),apply(len.nested,2,sd))/sqrt(n)


seq=seq(2,60,by=2)

df.cov3=df.cov3[seq,]
df.len3=df.len3[seq,]
df.cov5=df.cov5[seq,]
df.len5=df.len5[seq,]

df.cov3_sd=df.cov3_sd[seq,]
df.len3_sd=df.len3_sd[seq,]
df.cov5_sd=df.cov5_sd[seq,]
df.len5_sd=df.len5_sd[seq,]
#######################################################################
####################save data##########################################
#######################################################################
names=c("Linear","Naive","VFCP","EFCP")
force_bind = function(df1, df2) {
  colnames(df2) = colnames(df1)
  rbind(df1, df2)
}
format_data = function(data_df, data_sd, Moment_name, Var_name){
  `Linear` <- cbind(data_df[,c(1,2)],data_sd[,2])
  `Linear`$"Method"=names[1]
  `Naive` <- cbind(data_df[,c(1,3)],data_sd[,3])
  `Naive`$method=names[2]
  `Conformal` <- cbind(data_df[,c(1,4)],data_sd[,4])
  `Conformal`$Method=names[3]
  `Nested` <- cbind(data_df[,c(1,5)],data_sd[,5])
  `Nested`$Method=names[4]
  output=force_bind(`Linear`,`Naive`,`Conformal`,`Nested`)
  output$Moment=Moment_name
  output$Var=Var_name
  return(output)
}

a3_cov=format_data(df.cov3, df.cov3_df, "3rd moment", "Coverage")
a5_cov=format_data(df.cov5, df.cov5_df, "5th moment", "Coverage")
a3_len=format_data(df.len3, df.len3_sd, "3rd moment", "Width ratio")
a5_len=format_data(df.len5, df.len5_sd, "5th moment", "Width ratio")

data_linear_fm=force_bind(a3_cov,a5_cov,a3_len,a5_len)
colnames(data_linear_fm) <- c("V1","V2","sd","Method")
save(data_linear_fm,file="linear_data_to_plot.RData")