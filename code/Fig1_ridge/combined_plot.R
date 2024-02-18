rm(list=ls())
library("ggplot2")
load("non_linear_data_to_plot.RData")
load("linear_data_to_plot.RData")

data_nonlinear_fm$Linear= "Non Linear"
data_linear_fm$Linear= "Linear"
str(data_linear_fm)
str(data_nonlinear_fm)
data_plot <- rbind(data_nonlinear_fm, data_linear_fm)
str(data_plot)
data_plot$Method<- as.factor(data_plot$Method)
data_plot$Method <- factor(data_plot$Method, levels = c("EFCP", "VFCP","Linear","Naive"))
str(data_plot)
names=c("Linear","Naive","VFCP","CV*","CV-5-fold","EFCP")
colors_manual=c("dodgerblue4","slategrey","darkorchid3","red","turquoise3","grey23")
cbbPalette <- c("#000000","red",  "#009E73", "#0072B2", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
dummy_coverage <- data.frame(Var=c("Coverage"),Z= 0.9)
ggplot(data=data_plot, aes(x=V1, y=V2, group=Method,color=Method, linetype=Method))+geom_hline(data=dummy_coverage, aes(yintercept=Z), col="grey30")  +
  geom_ribbon(aes(ymin=V2-sd, ymax=V2+sd, fill=Method), alpha=.2,colour = NA)+
  geom_line( aes(color=Method))+ xlab("Dimension")+ ylab("")+theme(axis.text=element_text(size=11),axis.title=element_text(size=13,face="bold"),strip.text = element_text(size=13,face = "bold"))+
  facet_grid(Var~Linear+Moment,scales = "free") +scale_color_manual(values=cbbPalette)+ scale_linetype_manual(values=c(1, 5,4,2))+theme(legend.position="bottom", legend.text=element_text(size=13),legend.title=element_text(size=13))