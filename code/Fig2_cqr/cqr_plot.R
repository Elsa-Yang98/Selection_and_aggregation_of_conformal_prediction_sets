library(ggplot2)
load('real_to_plot.RData')
load('synthetic_to_plot.RData')

levels(real$V1)[6] <- "conduct"
cqr = rbind(synthetic, real)
cqr[,5][cqr[,5] == 'Width']='Width ratio'
dummy_coverage <- data.frame(Var=c("Coverage"),Z= 0.9)
cbbPalette <- c("blue","red",  "#000000", "#0072B2", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data=cqr, aes(x=V1, y=V2, group=Method,color=Method, linetype=Method))+geom_hline(data=dummy_coverage, aes(yintercept=Z), col="grey30")  +
  geom_ribbon(aes(ymin=V2-sd, ymax=V2+sd, fill=Method), alpha=.2,colour = NA)+
  geom_line( aes(color=Method))+ xlab("Dimension")+ ylab("")+theme(axis.text=element_text(size=11),axis.title=element_text(size=13,face="bold"),strip.text = element_text(size=13,face = "bold"))+
  facet_grid(Var~Name,scales = "free") + scale_linetype_manual(values=c(1, 5,4,2))+theme(legend.position="bottom", legend.text=element_text(size=13),legend.title=element_text(size=13))+ylim(0,1)+scale_color_manual(values=cbbPalette)
