require(ggplot2);require(scales)

require(chngpt)

d = read.table("combined_data.txt",header=T)

# S16
fitd16=chngptm(formula.1=qdistAstral16-qdist16~1, formula.2=~pdist, data=d , type="stegmented",
               family="gaussian")



t.test(with(d[d$pdist<fitd16$chngpt,],qdistAstral16-qdist16))

t.test(with(d[d$pdist>=fitd16$chngpt,],qdistAstral16-qdist16))
     
ggplot(d,aes(x=pdist,y=qdistAstral16-qdist16)) +
  geom_point(alpha=1/3,color="blue") + 
  geom_hline(yintercept = 0,linetype=2) + 
  geom_smooth(aes(group = pdist<fitd16[[2]][5]),method="lm",se=F,color="red")+
  theme_classic() + xlab("mean pairwise AAD") + ylab("quartet distance (adjusted)") 
ggsave("astral_vs_S16.pdf",width=7,height = 4)


# S23
fitd23=chngptm(formula.1=qdistAstral23-qdist23~1, formula.2=~pdist, data=d , type="stegmented",
            family="gaussian")

ggplot(d,aes(x=pdist,y=qdistAstral23-qdist23)) +
  geom_point(alpha=1/3,color="blue") + 
  geom_hline(yintercept = 0,linetype=2) + 
  geom_smooth(aes(group = pdist<fitd23[[2]][5]),method="lm",se=F,color="red")+
  theme_classic() + xlab("mean pairwise AAD") + ylab("quartet distance (adjusted)") 
ggsave("astral_vs_S23.pdf",width=7,height = 4)
