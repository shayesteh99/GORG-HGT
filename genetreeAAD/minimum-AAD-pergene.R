require(reshape2)

newd = read.csv('all-closest-dist.txt',sep="\t",h=F)
names(newd)=c("gene","query","neighbor","treedist","neighbor","treedist","neighbor","treedist")

newd = rbind(newd[,c(1,2,3,4)],newd[,c(1,2,5,6)],newd[,c(1,2,7,8)])
nrow(newd)


ming = dcast(query+gene.x~., 
             data=merge(merge(newd[newd$gene=="accA",c(1,2)],
                              aai[,c(1,3,6)],by.x=2,by.y="Genome.B"),
                        newd[newd$gene=="accA",c(1,2)],by.y=2,by.x="Genome.A"), 
             value.var = "Mean.AAI",fun.aggregate = min)
head(ming)
ming=ming[,c(2,1,3)]
names(ming)[1] = c("gene")
head(ming)

for (x in unique(newd$gene)[-1]
     #c("dnaK", "recA", "gyrA", "gyrB" )
) {
  tmp = dcast(query+gene.x~., 
              data=merge(merge(newd[newd$gene==x,c(1,2)],aai[,c(1,3,6)],by.x=2,by.y="Genome.B"),
                         newd[newd$gene==x,c(1,2)],by.y=2,by.x="Genome.A"), value.var = "Mean.AAI",
              fun.aggregate = min)
  tmp=tmp[,c(2,1,3)]
  names(tmp)[1] = c("gene")
  ming = rbind(ming,tmp)
  rm(tmp)
  print(x)
}
head(ming)
nrow(ming)
names(ming)[3]="min"

nrow(ming)
write.csv(ming,"minimum-AAD-per-gene.csv")

aaim = recast(aai[,c("Genome.A","Genome.B","Mean.AAI")],Genome.A~Genome.B,
                        fun.aggre=function(x)(x[1]))
nn = as.matrix(aaim[,1])
aaim=as.matrix(aaim[,-1])
rownames(aaim) = nn
head(aaim)
aaim["AG-359-E10","AG-918-K09"]

#m=merge(dcast(aai,Genome.A~.,value.var = "Mean.AAI",fun.aggregate = min),
#        dcast(aai,Genome.B~.,value.var = "Mean.AAI",fun.aggregate = min),by=1)
#m$min=apply(m[,2:3],1,min)
#head(m)
#summary(m$min)
