## run Tcheb Funcitons

setwd("~/Documents/MathCareer/Stephan Helfrich/")
source("TchebFunctions.R")

##START: Read in data 
fname <- "RunningEx-NDF" #running example or
fname <- "m25.1-NDF" #other example instances: m25.1, m25.2, m25.3, m25.4, m25.5 
baseNDF <- read.csv(paste0("Instances/",fname,".txt"),header=TRUE,sep=",")

### PRIMARY PLOTTING TOOLS

## Run Box-Based Criteria-Space Search Algorithm for NDF instance
policy=1                                   #policy 1 recommended for simplicity, but replace for exploration
outlist <- SimulatedPrimal(baseNDF,policy) #outlist is a list of output containing:
LNPlist = outlist[[1]]                     #list of data frames containing local nadir points (LNPs) per iteration
NDF = outlist[[2]]                         #data frame for nondominated frontier (NDF) including which iteration it was found
outlist[[3]]                               #simply prints the number of simulated integer programs (IPs) which would have been solved

## Plot final decomposition using Algorithm 1
s = length(LNPlist) #recommended to use final step; use subsequent section for choosing an intermediate step
LNP <- LNPlist[[s]] 
NDF[,c("w1","w2","w3")] <- t(apply(NDF[,c("y1","y2","y3")], 1, kernelWeight)) #compute kernel weight per image 
LNP[,c("w1","w2","w3")] <- t(apply(LNP[,c("y1","y2","y3")], 1, kernelWeight)) #compute kernel weight per lnp 
subNDF <- subset(NDF,StepFound<=s)                  #only include images found by step s
perimset <- ComputeAllPerimeters(subNDF,LNP)        #compute perimeter sets 
trilist <- TriangulatePerimeterSet(subNDF,perimset) #triangulate -> list of data frames containing 
triangledf <- trilist[[1]]                          #data frame for individual triangles
areadf <- trilist[[2]]                              #data frame for area per triangle
linesegdf <- trilist[[3]]                           #data frame for line segments of the perimeters

# 3 plotting options:
# (a) Complete Plot without labels
g1 <- ggplot() +
  geom_polygon(data=triangledf, aes(x=lambda1,y=lambda2,fill=image,group=tri),alpha=0.5) +
  geom_segment(data = outlinedf, aes(x=startx,y=starty,xend=endx,yend=endy),color="black",size=1)+
  geom_segment(data = linesegdf, aes(x=startx,y=starty,xend=endx,yend=endy),size=0.5,color="black") +
  theme_bw() +
  scale_fill_gradientn(colours=rainbow(length(unique(triangledf$image))), #length(unique(triangledf$image))
                       limits=c(1,length(unique(triangledf$image))))+ #length(unique(triangledf$image))
  #ggtitle(fname) +
  ggtitle("Instance m25.1") +
  xlab(expression(lambda[1])) +
  ylab(expression(lambda[2])) +
  theme(legend.position="none",#right or none
        legend.title = element_text(size=10),legend.text = element_text(size=8),
        axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))
g1

# (b) Complete plot with labels
g2<-g1+geom_point(data=subset(subNDF,id>0), 
                  aes(x=w1, y=w2), colour="black",size=4.5,shape = 21,fill="white") + #size=3
  geom_text(data=subset(subNDF,id>0), 
            aes(x=w1, y=w2, label = id), size=2.4)  #size=2

g2

# (c) Interactive using plotly
g3 <- ggplot() +
  geom_polygon(data=triangledf,
               aes(x=lambda1,y=lambda2,fill=image,group=tri,
                   text=paste0("Image: ",image)),alpha=0.5,color="gray") +
  geom_segment(data = outlinedf, aes(x=startx,y=starty,xend=endx,yend=endy),color="black",size=1)+
  geom_point(data = perimset,size=1, aes(x=w1,y=w2,
                                         text=paste0("Image: (",y1,",",y2,",",y3,")\n",
                                                     Kway,"-way LNP\n",
                                                     "C1: ",C1,"\n",
                                                     "C2: ",C2,"\n",
                                                     "C3: ",C3,"\n"))) +
  theme_bw() +
  scale_fill_gradientn(colours=rainbow(length(unique(triangledf$image))), #length(unique(triangledf$image))
                       limits=c(1,length(unique(triangledf$image))))+ #length(unique(triangledf$image))
  #labs(title=fname,x=expression("lambda_1"),y=expression("\lambda_2")) + 
  #ggtitle(fname) +
  xlab(expression(lambda[1])) +
  ylab(expression(lambda[2])) +
  theme(legend.position="none",#right or none
        legend.title = element_text(size=10),legend.text = element_text(size=8),
        axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))

ggplotly(g3,tooltip = 'text')

### PLOTTING APPROXIMATIONS

# Plot Inner or Outer Approximations at a single intermediate step
s = 10 #choose an intermediate step in order to evaluate an any-time approximation
LNP <- LNPlist[[s]] 
NDF[,c("w1","w2","w3")] <- t(apply(NDF[,c("y1","y2","y3")], 1, kernelWeight)) #compute kernel weights per image
LNP[,c("w1","w2","w3")] <- t(apply(LNP[,c("y1","y2","y3")], 1, kernelWeight)) #compute kernel weights per lnp
subNDF <- subset(NDF,StepFound<=s)                  #only include images found by step s

outer <- SingleStepOuterApprox(subNDF,LNP,s)
outertriangledf <- outer[[1]]
inner <- SingleStepInnerApprox(subNDF,LNP,s)
innertriangledf <- inner[[1]]
#toggle comments below depending on which approximation is of interest
g4 <- ggplot() +
  # for outer:
  ggtitle(paste("Step",s,"Outer Approximation")) +
  geom_polygon(data=outertriangledf, 
               aes(x=lambda1,y=lambda2,fill=image,group=tri),alpha=0.5) + #0.3
  # for inner: 
  #ggtitle(paste("Step",s,"Inner Approximation")) +
  #geom_polygon(data=innertriangledf, 
  #             aes(x=lambda1,y=lambda2,fill=image,group=tri),alpha=0.5,color="black") +
  geom_segment(data = outlinedf, aes(x=startx,y=starty,xend=endx,yend=endy),color="black",size=1) +
  scale_fill_gradientn(colours=rainbow(length(unique(NDF$id))-3), #length(unique(triangledf$image))
                       limits=c(1,length(unique(NDF$id))-3)) + #length(unique(triangledf$image))
  xlab(expression(lambda[1])) +
  ylab(expression(lambda[2])) +
  theme(legend.position="none",#right or none
        legend.title = element_text(size=10),legend.text = element_text(size=8),
        axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))

g4

# Plot Inner or Outer Approximations over run time
Approx <- SummarizeInnerOuterApprox(NDF,LNPlist) #takes some time to run
InnerApprox <- Approx[[1]]
OuterApprox <- Approx[[2]]
AllApprox <- Approx[[3]]
# optional aesthetics commented out below
g5 <- ggplot() +
  geom_line(data=OuterApprox, aes(x=Step, y=TotalArea, group=image, color=image),linetype="dashed") +
  #geom_point(data=OuterApprox, aes(x=Step, y=TotalArea, color=image)) +
  geom_line(data=InnerApprox, aes(x=Step, y=TotalArea, group=image, color=image)) +
  #geom_point(data=InnerApprox, aes(x=Step, y=TotalArea, color=image)) +
  xlab("Step") + ylab("Approximated Area") + 
  #ggtitle("Instance m25.1") +
  theme(legend.position="none",#right or none
        legend.title = element_text(size=10),legend.text = element_text(size=8),
        #axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))
g5

### SIMULATED EXPERIMENTS

## Evaluate Policies
eval_policies <- data.frame(Step=0, Gap=0, Policy=0)
repby_policies <- data.frame(Step=0, Dmean=0, Dmax=0, Policy=0)
for(policy in 1:5) {
  print(paste("Evaluating policy",policy))
  outlist <- SimulatedPrimal(baseNDF,policy)
  LNPlist = outlist[[1]]
  NDF = outlist[[2]]
  Approx <- SummarizeInnerOuterApprox(NDF,LNPlist)
  AllApprox <- Approx[[3]]
  gap_sum <- AllApprox %>% group_by(Step) %>% summarise(Gap = sum(OuterMinusInner))
  gap_sum$Policy <- policy
  eval_policies <- rbind(eval_policies,gap_sum)
  
  repmetrics = representationMetrics(NDF)
  repmetrics$Policy <- policy
  repby_policies <- rbind(repby_policies,repmetrics)
}
eval_policies <- eval_policies[-1,]
eval_policies$Policy <- as.factor(eval_policies$Policy)

repby_policies <- repby_policies[-1,]
repby_policies$Policy <- as.factor(repby_policies$Policy)

#Plotting Outer-Inner
g1 <- ggplot() +
  geom_line(data=eval_policies, aes(x=Step, y=Gap, group=Policy, colour=Policy)) +
  xlab("Step") + ylab("Approximation Gap (Sum)") + 
  #ggtitle(paste("Policy",policy)) +
  theme(legend.position="right",#right or none
        legend.title = element_text(size=10),legend.text = element_text(size=8),
        #axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))
g1

#Plot representation metrics
g2 <- ggplot() + 
  geom_line(data=repby_policies, aes(x=Step, y=Dmean, group=Policy, colour=Policy)) +
  geom_line(data=repby_policies, aes(x=Step, y=Dmax, group=Policy, colour=Policy),linetype="dashed") +
  xlab("Step") + ylab("Representation Metrics") + 
  #ggtitle(paste("Policy",policy)) +
  theme(legend.position="right",#right or none
        legend.title = element_text(size=10),legend.text = element_text(size=8),
        #axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))
g2

#Summarize ranks across policies
rankdf <- data.frame(Step=0, Metric="Label",
                     Policy1 = NA, Policy2 = NA, Policy3=NA, Policy4=NA, Policy5=NA)
for(s in 2:max(eval_policies$Step)) {
  if(s %in% eval_policies$Step) {
    sub_eval = subset(eval_policies,Step==s)
  } else {
    sprev = max(eval_policies$Step[eval_policies$Step<s])
    sub_eval = subset(eval_policies,Step==sprev)
  }
  
  if(dim(sub_eval)[1]>0) {
    sub_eval$Rank=NA
    for(i in 1:dim(sub_eval)[1]) {
      val=sub_eval$Gap[i]
      sub_eval$Rank[i]=sum(sub_eval$Gap<val)+1
    }
    addrow = data.frame(Step=s, Metric="Gap",
                        Policy1 = sub_eval$Rank[which(sub_eval$Policy==1)], 
                        Policy2 = sub_eval$Rank[which(sub_eval$Policy==2)], 
                        Policy3 = sub_eval$Rank[which(sub_eval$Policy==3)], 
                        Policy4 = sub_eval$Rank[which(sub_eval$Policy==4)], 
                        Policy5 = sub_eval$Rank[which(sub_eval$Policy==5)])
    rankdf = rbind(rankdf,addrow)
  }
  sub_rep = subset(repby_policies,Step==s)
  for(p in 1:5) {
    if(!(p%in%sub_rep$Policy)) {
      sprev = max(repby_policies$Step[repby_policies$Step<s & repby_policies$Policy==p])
      sub_rep = rbind(sub_rep,subset(repby_policies,repby_policies$Step==sprev & repby_policies$Policy==p))
    }
  }
  if(dim(sub_rep)[1]>0) {
    sub_rep$DmeanRank=NA
    sub_rep$DmaxRank=NA
    for(i in 1:dim(sub_rep)[1]) {
      val=sub_rep$Dmean[i]
      sub_rep$DmeanRank[i]=sum(sub_rep$Dmean<val)+1
      val=sub_rep$Dmax[i]
      sub_rep$DmaxRank[i]=sum(sub_rep$Dmax<val)+1
    }
    addrow1 = data.frame(Step=s, Metric="Dmean",
                         Policy1 = sub_rep$DmeanRank[which(sub_rep$Policy==1)], 
                         Policy2 = sub_rep$DmeanRank[which(sub_rep$Policy==2)], 
                         Policy3 = sub_rep$DmeanRank[which(sub_rep$Policy==3)], 
                         Policy4 = sub_rep$DmeanRank[which(sub_rep$Policy==4)], 
                         Policy5 = sub_rep$DmeanRank[which(sub_rep$Policy==5)])
    addrow2 = data.frame(Step=s, Metric="Dmax",
                         Policy1 = sub_rep$DmaxRank[which(sub_rep$Policy==1)], 
                         Policy2 = sub_rep$DmaxRank[which(sub_rep$Policy==2)], 
                         Policy3 = sub_rep$DmaxRank[which(sub_rep$Policy==3)], 
                         Policy4 = sub_rep$DmaxRank[which(sub_rep$Policy==4)], 
                         Policy5 = sub_rep$DmaxRank[which(sub_rep$Policy==5)])
    rankdf = rbind(rankdf,addrow1,addrow2)
  }
}
rankdf = rankdf[-1,]

summary_rank = data.frame(Policy = 0, Metric = "Label", PM="Label", Rank = 0, Freq = 0)
avg_rank = data.frame(Policy = 0, Metric = "Label", PM="Label", AvgRank = 0)
for(p in 1:5) {
  for(m in unique(rankdf$Metric)) {
    Rsum = 0
    Freqsum = 0
    for(r in 1:5) {
      freq = sum(rankdf[which(rankdf$Metric==m),paste0("Policy",p)]==r)
      summary_rank=rbind(summary_rank,c(p,m,paste0("P",p,"-",m),r,freq))
      Rsum = Rsum+r*freq
      Freqsum = Freqsum+freq
    }
    avg_rank = rbind(avg_rank,c(p,m,paste0("P",p,"-",m),Rsum/Freqsum))
  }
}
summary_rank = summary_rank[-1,]
summary_rank$Policy = as.factor(summary_rank$Policy)
summary_rank$PM = as.factor(summary_rank$PM)
summary_rank$Rank = as.numeric(summary_rank$Rank)
summary_rank$Freq = as.numeric(summary_rank$Freq)
avg_rank = avg_rank[-1,]
avg_rank$Policy = as.factor(avg_rank$Policy)
avg_rank$PM = as.factor(avg_rank$PM)
avg_rank$AvgRank = as.numeric(avg_rank$AvgRank)


g3 <- ggplot(data=summary_rank) +
  geom_point(aes(x=PM, y=Rank, colour=Policy, size=Freq)) +
  geom_point(data=avg_rank, aes(x=PM, y=AvgRank), shape=18, size=3) +
  xlab("") + 
  ylab("Ranked Performance") +
  scale_y_continuous(trans = "reverse") +
  #scale_y_discrete(limits=rev)
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0))


g3 #save 600x400

## Run Compromise Region
maxlambda=0.5
initLambda=data.frame(w1=c(1-2*maxlambda,maxlambda,maxlambda),
                      w2=c(maxlambda,1-2*maxlambda,maxlambda),
                      w3=c(maxlambda,maxlambda,1-2*maxlambda))
outlist2 <- CompromisePrimal(baseNDF,maxlambda,initLambda)
LNPlist = outlist2[[1]]
NDF = outlist2[[2]]
outlist2[[3]]

LNP <- LNPlist[[length(LNPlist)]] #length(LNPlist)
#Compute kernel weights per image and LNP
NDF[,c("w1","w2","w3")] <- t(apply(NDF[,c("y1","y2","y3")], 1, kernelWeight))
LNP[,c("w1","w2","w3")] <- t(apply(LNP[,c("y1","y2","y3")], 1, kernelWeight))

subNDF <- subset(NDF,StepFound>0)
#Compute perimeter sets and triangulate
perimset <- ComputeAllPerimeters(subNDF,LNP)
trilist <- TriangulatePerimeterSet(subNDF,perimset)
triangledf <- trilist[[1]]
areadf <- trilist[[2]]
linesegdf <- trilist[[3]]

g2 <- ggplot() +
  geom_polygon(data=triangledf, aes(x=lambda1,y=lambda2,fill=image,group=tri),alpha=0.5) +
  geom_segment(data = outlinedf, aes(x=startx,y=starty,xend=endx,yend=endy),color="black",size=1)+
  geom_segment(data = linesegdf, aes(x=startx,y=starty,xend=endx,yend=endy),size=0.5,color="black") +
  theme_bw() +
  scale_fill_gradientn(colours=rainbow(max(subNDF$id)), #length(unique(triangledf$image))
                       limits=c(1,max(subNDF$id)))+ #length(unique(triangledf$image))
  #ggtitle(fname) +
  #ggtitle("Instance m25.1") +
  xlab(expression(lambda[1])) +
  ylab(expression(lambda[2])) +
  theme(legend.position="none",#right or none
        legend.title = element_text(size=10),legend.text = element_text(size=8),
        axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))
g2

g3<-g2+geom_polygon(data=initLambda, aes(x=w1,y=w2),fill=NA,color="black") +
  geom_point(data=subset(subNDF,id>0), 
             aes(x=w1, y=w2), colour="black",size=4.5,shape = 21,fill="white") + #size=3
  geom_text(data=subset(subNDF,id>0), 
            aes(x=w1, y=w2, label = id), size=2.4)  #size=2

g3

