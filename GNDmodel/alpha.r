require(ggplot2);require(scales); require(reshape2); library(ggpubr)



## Read the blastn results
bd = read.csv('220916_sag_pair_summary_shared_orfs_ani.csv.xz')
#bd = read.csv('220926_sag_pair_summary_shared_orfs_ani.csv.xz')
nrow(bd)
head(bd)

# Divide GND into bins of size 0.1%
bd$gndbin=cut(bd$mean_pident/100,breaks=(670:1000)/1000, labels =round((1-(671:1000)/1000+0.0005)*100,digits=5))
dcast(data=bd[,c("X99_pid_500.1500bp_orthologs","gndbin")],formula=gndbin~.)
nrow(bd)
head(bd)
nrow(bd[bd$X99_pid_500.1500bp_orthologs == 0,])

## Read completeness stats
co = read.csv('Table_S1_80pct-comp_contains-16S.csv')
co$Genome_completeness_. = co$Genome_completeness_./100
head(co)

# Merge blastn results with completeness
bdrm = merge(merge(bd,co[,c(1,12)],by.x = "qsag", by.y="SAG"),co[,c(1,12)],by.x = "ssag", by.y="SAG")
nrow(bdrm)
names(bdrm)[33:34] = c("CompA","CompB")
head(bdrm)
bdrm[is.na(bdrm$X99),]

# Remove odd pairs where the number of hits is too low for the GND level.  
p1 = ggplot(aes(y=  total_hits  , color =total_hits> 800-50*(100-mean_pident),
           x=(100-mean_pident)/100),data=bdrm)+ 
  #stat_function(fun = function(x) 100*(1-x*5),color="black")+
  #         color= (total_hits> 50*mean_pident-4200), x=as.numeric(as.character(gndbin))/100),data=bdrm)+
  stat_function(fun = function(x) 800-5000*x,color="black")+
  geom_point(alpha=0.2,size=0.3)+ 
  #facet_wrap(.~(is.na(ani_aln_coverage_ab) |  is.na(ani_aln_coverage_ba) |(ani_aln_coverage_ab<0.03 | ani_aln_coverage_ba<0.03)))+
  scale_x_continuous(name="GND",labels=percent)+
  scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
  scale_color_manual(name="Include?",values = c("#E855B0","#11A026","#000577"))+
  scale_y_continuous(name="Number of blastn hits")+
  theme_bw()+
  theme(legend.position = c(.8,.7),panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0,0.29), ylim=c(0,2000))

p1
ggsave("Filtering-low-hits.pdf",width = 6,height = 4)
ggsave("Filtering-low-hits.png",width = 6,height = 4)


bdrmi=bdrm[bdrm$total_hits> 800-50*(100-bdrm$mean_pident),]

### This is what we used
alphas=quantile(with(bdrmi[bdrmi$mean_pident!= 100 & bdrmi$stdev_pident!=0 &  bdrmi$mean_pident <95 & bdrmi$mean_pident > 80  &
                          !is.na(bdrmi$total_hits)   & bdrmi$total_hits > 4,
                        ,], 
                     1/( (stdev_pident/100)/(1-mean_pident/100) )^2 ),c(0.1,0.5,0.9))

alphas
p2= ggplot(aes(x=1/( (stdev_pident/100)/(1-mean_pident/100) )^2 ),
       data= bdrmi[ bdrmi$stdev_pident!=0 & bdrmi$mean_pident != 100& bdrmi $ mean_pident > 0 & 
                      !is.na(bdrmi$total_hits) & !is.na(bdrmi$mean_pident) &
                      bdrmi$total_hits>4,] )+
  geom_histogram(binwidth =0.5)+
  scale_x_continuous(name=expression(alpha),trans = "identity")+
  geom_vline(xintercept = alphas,color="red",linetype=2)+
  theme_bw()+
  facet_grid( .~cut(100-mean_pident,100-c(70,80,95,100)))+
  xlab("alpha")+coord_cartesian(xlim = c(0.1,25))
p2
ggsave("alpha-estimate-nolog.pdf",width=6.5,height = 3.8)

ggplot(aes(x=1/( (stdev_pident/100)/(1-mean_pident/100) )^2 ),
       data= bdrmi[ bdrmi$stdev_pident!=0 & bdrmi$mean_pident != 100& bdrmi $ mean_pident > 0 & 
                      !is.na(bdrmi$total_hits) & !is.na(bdrmi$mean_pident) &
                      bdrmi$total_hits>4,] )+
  geom_histogram(binwidth =0.05)+
  scale_x_continuous(name=expression(alpha),trans = "log10")+
  geom_vline(xintercept = alphas,color="red",linetype=2)+
  theme_bw()+
  facet_grid( .~cut(100-mean_pident,100-c(70,80,95,100)))+
  xlab("alpha")#+coord_cartesian(xlim = c(0.1,100))
#ggsave("alpha-estimate.pdf",width=6.5,height = 3.8)


ggarrange(p2, p1, ncol = 1, labels = c("A", "B"),heights = c(2/3,1))
ggsave(paste("SupplementaryS",5,".pdf",sep=""),width=6.5,height = 8)
ggsave(paste("SupplementaryS",5,".png",sep=""),width=6.5,height = 8)


###
### NOW YOU HAVE TO RUN alpha.nb in Mathematica to produce the following files. 
### Make sure the following alphas are used
alphas

# Read the results of the model, as produced by the Mathematica 
model=rbind(
  data.frame(read.csv('alpha7.77-c1-gl904-orf1.csv'),alpha=7.77,c=1,gl=904,orf=1),
  data.frame(read.csv('alpha5.28-c1-gl904-orf1.csv'),alpha=5.28,c=1,gl=904,orf=1),
  data.frame(read.csv('alpha3.31-c1-gl904-orf1.csv'),alpha=3.31,c=1,gl=904,orf=1))
model$GND = round(model$GND,5)
tail(model)

bdrmim = melt(bdrmi,measure.vars = 10:23)
bdrmim = bdrmim[grepl("500" ,bdrmim$variable),]
bdrmim$adjusted = with(bdrmim,value*(1/sgene_count_500.1500bp/CompB+1/qgene_count_500.1500bp/CompA)/2 )
bdrmim$similarity=sub("_.*","",bdrmim$variable)
head(bdrmim)
levels(bdrmim$variable)

vapply(c("99.9","99.5","99","98.5","98","97.5","97"),function(x) {
  xx=paste("X",x,sep="");
  v=paste(xx,"_pid_500.1500bp_orthologs",sep="");
  fn=paste("model-emp-",x,".csv",sep="");
  print(v);
  write.csv(merge(dcast(GND~alpha,data=model[,c("GND",xx,"alpha")],value.var = xx),
                merge(dcast(gndbin~"real_mean",data=bdrmim[bdrmim$variable==v,c("gndbin","adjusted")],fun.aggregate = mean),
                      dcast(gndbin~"real_count",data=bdrmim[bdrmim$variable==v,c("gndbin","adjusted")])),
                by.x="GND",by.y="gndbin"),fn);
  fn;
},c("1"))


ggplot(aes(y=  adjusted ,
           color="Data", x=1-mean_pident/100),
       data=bdrmim[bdrmim$similarity=="X99",])+
  geom_point(alpha=0.5,size=0.5)+ #geom_smooth(se=F,aes(color="Data (fit)"),size=1)+
  geom_line(aes(y=`.`,color="Data (mean)",x=as.numeric(as.character(gndbin))/100),
            data=dcast(gndbin+variable~.,data=bdrmim[bdrmim$similarity=="X99",c("gndbin","variable", "adjusted")],fun.aggregate = mean),
            size=1)+
  geom_ribbon(aes(ymin=`3.31`,y=`5.28`,ymax=`7.77`,x=GND/100,color="Model (80% CI alpha)"),
              data=dcast(GND~alpha,data=model[,c(1,5,11)],value.var = "X99"),
              size=0.2,alpha=0.3,fill="#EE60BB",show.legend = F)+
  geom_line(aes(y=X99/orf/c,x=GND/100, color="Model (median alpha)"),
            data=model[model$alpha==5.28,],
            size=1,alpha=0.8)+
  scale_x_continuous(name="GND",labels=percent)+
  scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
  scale_color_manual(name="",values = c("gray50","#1055EE","#FF60AA","#BB3333","#770010"))+
  scale_y_continuous(name="Shared genes (adjusted for incompleteness)",labels=percent,breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme_bw()+
  theme(legend.position = c(.8,.7),panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0,0.27))
ggsave("shared-genes-percent-newformula-99-210819.pdf",width=6.5,height = 4.5)
ggsave("shared-genes-percent-newformula-99-210819.png",width=6.5,height = 4.5)



# !(is.na(bdrmani_aln_coverage_ab) |  is.na(ani_aln_coverage_ba) |(ani_aln_coverage_ab<0.03 | ani_aln_coverage_ba<0.03))
ggplot(aes(y=  adjusted ,
           color="Data", x=1-mean_pident/100),
       data=bdrmim)+
  facet_wrap(~similarity)+
  geom_point(alpha=0.5,size=0.5)+ #geom_smooth(se=F,aes(color="Data (fit)"),size=1)+
  geom_line(aes(y=`.`,color="Data (mean)",x=as.numeric(as.character(gndbin))/100),
            data=dcast(gndbin+similarity~.,data=bdrmim[,c("gndbin","similarity", "adjusted")],fun.aggregate = mean),
            size=.4)+
  geom_ribbon(aes(ymin=`3.31`,y=`5.28`,ymax=`7.77`,x=GND/100,color="Model (80% CI alpha)"),
              data=dcast(GND+similarity~alpha,data=melt(model[,c(1,3:9,11)],measure.vars = 2:8,variable.name = "similarity")[,c("GND","alpha","similarity","value")],value.var = "value"),
              size=0.2,alpha=0.3,fill="#EE60BB",show.legend = F)+
  geom_line(aes(y=`5.28`,x=GND/100, color="Model (median alpha)"),
               data=dcast(GND+similarity~alpha,data=melt(model[model$alpha==5.28,c(1,3:9,11)],measure.vars = 2:8,variable.name = "similarity")[,c("GND","alpha","similarity","value")],value.var = "value"),
            size=1,alpha=0.8)+
  scale_x_continuous(name="GND",labels=percent)+
  scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
  scale_color_manual(name="",values = c("gray50","#1055EE","#FF60AA","#BB3333","#770010"))+
  scale_y_continuous(name="Shared genes (adjusted for incompleteness)",labels=percent,breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme_bw()+
  theme(legend.position = c(.8,.2),panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0,0.27))
#ggsave("shared-genes-percent-newformula-multiple-210819.pdf",width=6.5,height = 6.5)
ggsave("shared-genes-percent-newformula-multiple-210819.png",width=6.5,height = 6.5)


ds = merge(dcast(variable+GND~alpha,data=melt(model[,c(1,3:9,11)],measure.vars = 2:8)[,c("GND","alpha","variable","value")],value.var = "value"),
      dcast(gndbin+similarity~"real",data=bdrmim[,c("gndbin","similarity", "adjusted")],fun.aggregate = mean),
      by.x=c("variable","GND"),by.y=c("similarity","gndbin"))
ds=melt(ds,measure.vars = 3:5,variable.name = "alpha",value.name = "model")

write.csv(x = ds,file="model-emp-all.csv",row.names = F)
head(ds)

ggplot(aes(y=real-model,x=GND/100,color=sub("X","",variable)),data=ds)+geom_line()+
  scale_color_brewer(palette = "Dark2",name="similarity")+
  scale_x_continuous(lim=c(0,0.13),labels = percent)+
  facet_wrap(~alpha,labeller = label_both,nrow=2)+
  theme_bw()+
  theme(legend.position = c(.8,.2))
ggsave("divergence-multiple-210819.pdf",width=6.5,height = 6.5)
ggsave("divergence-multiple-210819.png",width=6.5,height = 6.5)



ggplot(aes(y=real-model,x=GND/100,color=as.factor(100-as.numeric(sub("X","",variable)))),data=ds[ds$alpha == 5.28 & ds$variable %in% c("X99","X98.5","X98"),])+
  geom_line(size=0.8)+
  scale_color_manual(name="Gene ND threshold",values=c("red","#009900","blue"))+
  scale_x_continuous(lim=c(0,0.13),labels = function(x) x*100,breaks = (0:13)/100,"Genomic nucleotid difference, %")+
  #facet_wrap(~alpha,labeller = label_both,nrow=2)+
  theme_classic()+
  geom_hline(yintercept = 0,color="grey")+
  theme(legend.position = c(.8,.85))
ggsave("divergence-main-210819.pdf",width=6,height = 4)






