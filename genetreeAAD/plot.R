require(ggplot2); require(scales); require(reshape2)

aai = read.csv("GORGv1_16SSAGs_aai_summary.csv.xz",sep="\t")
aai$Mean.AAI = (100 - aai$Mean.AAI)/100
head(aai)
aai2 = aai[,c(3:4,1:2,5:8)]
names(aai2)=names(aai)
aai = rbind(aai,aai2)
rm(aai2)

newd = read.csv('all-closest-dist.txt.xz',sep="\t",h=F)
names(newd)=c("gene","query","neighbor","treedist","neighbor","treedist","neighbor","treedist")

newd = rbind(newd[,c(1,2,3,4)],newd[,c(1,2,5,6)],newd[,c(1,2,7,8)])
nrow(newd)


ming = read.csv("minimum-AAD-per-gene.csv.xz")[,2:4]

newd = merge(newd,aai,by.x=c("query","neighbor"),by.y=c("Genome.A","Genome.B"))
nrow(newd)
newd = merge(newd,ming,by.x=c(3,1),by.y=1:2)
nrow(newd)
head(newd)


ggplot(aes(y=Mean.AAI,x=min,color=treedist+10^-5),
       data=newd[newd$gene %in% c("dnaK", "recA", "gyrA", "ychF" )
                 & newd$treedist < 1,])+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  #annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.5)+
  facet_wrap(~gene,nrow=1)+
  scale_color_gradient(high = "#4080A0",low="#F02000",name="Gene tree dist",
                       breaks=c(10^-5,10^-3,10^-2,10^-1,10^0),
                       labels=c(    0,10^-3,10^-2,10^-1,10^0),
                        trans="log10")+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="log10") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="log10")+
  theme(panel.spacing = unit(0,"pt"),
        legend.position=c(.914,.16),
        legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.5),xlim=c(0.001,0.5))+
  guides(colour = guide_colorbar(title.position="top",title.hjust = 0.5))
ggsave("dminA_vs_dminG_all_new_selected2_log_closest3.pdf",width = 9,height=3.4)


require(tidyverse)
newd[newd$gene %in% c("trxA" ) & newd$treedist < 1,] %>%
  mutate(delta = Mean.AAI- min) %>%
  ggplot(aes(x=delta))+geom_density()


dnak = newd[newd$gene %in% c("trxA" ) & newd$treedist < 1,] %>%
  mutate(delta = (Mean.AAI-min)) 

ggplot(data=dnak,aes(x=delta))+geom_density(adjust = 2)+#coord_cartesian(xlim=c(0,0.02)) +
  stat_function(fun=function(x) dgamma(x,shape=mean(dnak$delta)^2/var(dnak$delta),rate=mean(dnak$delta)/var(dnak$delta)),color="red") +
  stat_function(fun=function(x) dexp(x,log(2)/median(dnak$delta)),color="blue")

ggplot(aes(sample=delta),data=dnak)+
  geom_qq_line(distribution = stats::qgamma,dparams = c(shape=mean(dnak$delta)^2/var(dnak$delta),rate=mean(dnak$delta)/var(dnak$delta)))+ 
  stat_qq(distribution = stats::qgamma,dparams = c(shape=mean(dnak$delta)^2/var(dnak$delta),rate=mean(dnak$delta)/var(dnak$delta)))

ggplot(aes(sample=delta),data=dnak)+
  geom_qq_line(distribution = stats::qexp,dparams = c(rate=log(2)/median(dnak$delta)))+ 
  stat_qq(distribution = stats::qexp,dparams = c(rate=log(2)/median(dnak$delta)))


ggplot(aes(x=delta),data=dnak)+stat_ecdf()+
  geom_vline(xintercept = qexp(1-c(0.001,0.01,0.05),rate=log(2)/median(dnak$delta)),color="red")+
  stat_function(fun=function(x) pexp(x,log(2)/median(dnak$delta)),color="blue")
  

ggplot(aes(y=Mean.AAI,x=min,color=cut(delta,qexp(1-c(0,0.001,0.01,0.05,1),rate=log(2)/median(delta)))),data=dnak)+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  #annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.5)+
  facet_wrap(~gene,nrow=1)+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="log10") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="log10")+
  theme(panel.spacing = unit(0,"pt"),
        legend.position=c(.614,.16),
        legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.5),xlim=c(0.001,0.5))+
  scale_colour_viridis_d(direction = -1,name="p-value")


ggplot(aes(y=Mean.AAI,x=min,color=treedist+10^-5),
       data=newd[newd$gene %in% c("dnaK", "recA", "gyrA", "ychF" )
                 & newd$treedist < 1,])+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  geom_point(size=0.5)+
  annotate("rect", xmin = 0, xmax = 0.02, ymin = 0.06, ymax =0.5, alpha = .15)+
  facet_wrap(~gene,nrow=1)+
  scale_color_gradient(high = "#4080A0",low="#F02000",name="Gene tree dist",
                       breaks=c(10^-5,10^-3,10^-2,10^-1,10^0),
                       labels=c(    0,10^-3,10^-2,10^-1,10^0),
                       trans="log10")+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="identity") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="identity")+
  theme(panel.spacing = unit(0,"pt"),
        legend.position=c(.914,.16),
        legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.4),xlim=c(0.001,0.4))+
  guides(colour = guide_colorbar(title.position="top",title.hjust = 0.5))
ggsave("dminA_vs_dminG_all_new_selected2_closest3.pdf",width = 9,height=3.4)


ggplot(aes(y=Mean.AAI,x=min,color=treedist+10^-5),
       data=newd[# newd$gene %in% c("dnaK", "recA", "gyrA", "gyrB" )&
                   newd$treedist < 1,])+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.2)+
  facet_wrap(~gene,ncol=12)+
  scale_color_gradient(high = "#4080A0",low="#F02000",name="gene tree distance",
                       breaks=c(10^-5,10^-3,10^-2,10^-1,10^0),
                       labels=c(    0,10^-3,10^-2,10^-1,10^0),
                       trans="log10")+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="log10") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="log10")+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom",legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.5),xlim=c(0.001,0.5))
ggsave("dminA_vs_dminG_new_log_closest3.png",width = 14,height=18,dpi = "retina")
ggsave("dminA_vs_dminG_new.eps",width = 14,height=15)
ggsave("dminA_vs_dminG_new_log_closest3.pdf",width = 14,height=15)

#write.csv(
  with(gyrA,gyrA[Mean.AAI-min>.1 & min<0.05 & V3<0.01,c(1,2,3,11)])
  #,file="gyrA-examples.csv")
