## Tchebychev Functions

library(ggplot2)
library(plotly)

#Plot
outlinedf <- data.frame(startx=c(0,0,1),starty=c(0,1,0),endx=c(0,1,0),endy=c(1,0,0))

InitializeNDF <- function(NDF) {
  NDF <- NDF[,c("id","y1","y2","y3")]
  utopia <- sapply(NDF[,c("y1","y2","y3")],min)-1
  #shift to positive orthant
  NDF$y1 <- NDF$y1-utopia[1]
  NDF$y2 <- NDF$y2-utopia[2]
  NDF$y3 <- NDF$y3-utopia[3]
  #add dummy images
  initialNadir <- c(max(NDF$y1),
                    max(NDF$y2),
                    max(NDF$y3))
  bigM <- 100*10^(ceiling(log10(max(initialNadir))))
  dummy <- NDF[1,]
  for(i in 1:3) {
    dummy[,c("y1","y2","y3")] <- utopia
    dummy$id <- -i
    dummy[,paste0("y",i)] <- initialNadir[i]+bigM #BigM
    NDF<-rbind(dummy,NDF)
  }
  return(NDF)
}

InitializeLNP <- function(NDF) {
  LNP <- data.frame(y1=NDF$y1[NDF$id==-1],
                    y2=NDF$y2[NDF$id==-2],
                    y3=NDF$y3[NDF$id==-3],
                    C1=I(list(-1)),
                    C2=I(list(-2)),
                    C3=I(list(-3)),
                    Confirmed=FALSE,
                    StepConfirmed=0)
  return(LNP)  
}

ComputeKway <- function(LNP) {
  LNP$Kway <- 0
  for(n in 1:dim(LNP)[1]) {
    LNP$Kway[n]<-length(unlist(LNP$C1[n])) + length(unlist(LNP$C2[n])) + length(unlist(LNP$C3[n]))
  }
  return(LNP)
}

KlamrothUpdate <- function(LNP, NDF, yrow) {
  A <- which(yrow$y1<LNP$y1 & yrow$y2<LNP$y2 & yrow$y3<LNP$y3)
  P <- LNP[1,]
  for(n in 1:dim(LNP)[1]) {
    if(yrow$y1==LNP$y1[n] && yrow$y2<LNP$y2[n] && yrow$y3<LNP$y3[n]) 
      LNP$C1[[n]] = c(LNP$C1[[n]],yrow$id)
    if(yrow$y1<LNP$y1[n] && yrow$y2==LNP$y2[n] && yrow$y3<LNP$y3[n]) 
      LNP$C2[[n]] = c(LNP$C2[[n]],yrow$id)
    if(yrow$y1<LNP$y1[n] && yrow$y2<LNP$y2[n] && yrow$y3==LNP$y3[n]) 
      LNP$C3[[n]] = c(LNP$C3[[n]],yrow$id)
  }
  for(i in A) {
    for(j in 1:3) {
      m_j_n = 0
      for(k in setdiff(1:3,j)) {
        ids <- LNP[[i,paste0("C",k)]]
        min_y_j <- min(NDF[NDF$id%in%ids, paste0("y",j)])
        m_j_n = max(m_j_n,min_y_j)
      }
      if(yrow[,paste0("y",j)]>m_j_n) {
        newn <- LNP[i,]
        newn[,paste0("y",j)] <- yrow[,paste0("y",j)]
        newn[,paste0("C",j)] <- yrow$id
        for(k in setdiff(1:3,j)) {
          ids <- LNP[[i,paste0("C",k)]]
          ids <- ids[NDF[NDF$id%in%ids, paste0("y",j)]<yrow[,paste0("y",j)]]
          newn[[paste0("C",k)]] <- I(list(ids))
        }
        P<-rbind(P,newn)
        
      }
    }
  }
  P<-P[-1,]
  if(length(A)>0) LNP <- LNP[-A,]
  if(dim(P)[1]>0) LNP<-rbind(LNP,P)
  return(LNP)
}

KlamrothUpdateMaximal <- function(LNP, NDF, yrow) {
  A <- which(yrow$y1<LNP$y1 & yrow$y2<LNP$y2 & yrow$y3<LNP$y3)
  P <- LNP[1,]
  for(n in 1:dim(LNP)[1]) {
    if(yrow$y1==LNP$y1[n] && yrow$y2<=LNP$y2[n] && yrow$y3<=LNP$y3[n]) 
      LNP$C1[[n]] = c(LNP$C1[[n]],yrow$id)
    if(yrow$y1<=LNP$y1[n] && yrow$y2==LNP$y2[n] && yrow$y3<=LNP$y3[n]) 
      LNP$C2[[n]] = c(LNP$C2[[n]],yrow$id)
    if(yrow$y1<=LNP$y1[n] && yrow$y2<=LNP$y2[n] && yrow$y3==LNP$y3[n]) 
      LNP$C3[[n]] = c(LNP$C3[[n]],yrow$id)
  }
  for(i in A) {
    for(j in 1:3) {
      m_j_n = 0
      for(k in setdiff(1:3,j)) {
        ids <- LNP[[i,paste0("C",k)]]
        min_y_j <- min(NDF[NDF$id%in%ids, paste0("y",j)])
        m_j_n = max(m_j_n,min_y_j)
      }
      if(yrow[,paste0("y",j)]>m_j_n) {
        newn <- LNP[i,]
        fullC <- c(unlist(newn$C1),unlist(newn$C2),unlist(newn$C3))
        newn[,paste0("y",j)] <- yrow[,paste0("y",j)]
        newn[,paste0("C",j)] <- I(list(yrow$id))
        for(k in setdiff(1:3,j)) {
          ids <- LNP[[i,paste0("C",k)]]
          ids <- ids[NDF[NDF$id%in%ids, paste0("y",j)]<=yrow[,paste0("y",j)]]
          newn[[paste0("C",k)]] <- I(list(ids))
          #ids <- LNP[[i,paste0("C",k)]]
          #ids <- ids[
          #inds <- which(NDF$id %in% fullC)
          addToCj <- which(NDF$id %in% fullC &
                             NDF[, paste0("y",j)]==yrow[,paste0("y",j)] &
                             NDF$y1<=newn$y1 &
                             NDF$y2<=newn$y2 &
                             NDF$y3<=newn$y3 )
          newn[[paste0("C",j)]] <- I(list(unique(c(unlist(newn[[paste0("C",j)]]),NDF$id[addToCj]))))
        }
        P<-rbind(P,newn)
        
      }
    }
  }
  P<-P[-1,]
  if(length(A)>0) LNP <- LNP[-A,]
  if(dim(P)[1]>0) LNP<-rbind(LNP,P)
  return(LNP)
}

SimulatedPrimal <- function(baseNDF,policy) {
  NDF <- InitializeNDF(baseNDF)
  NDF$StepFound = 0
  NDF$StepFound[1:3] = 1
  
  LNP <- InitializeLNP(NDF)
  LNP <- ComputeKway(LNP)
  LNPlist <- list("step0"=LNP)
  IPsolved = 1
  while(sum(LNP$Confirmed==FALSE)>0) {
    #policy for choosing LNP
    if(policy==1) {
      #Policy 1: First In First Out
      test = which(LNP$Confirmed==FALSE)[1] 
    } else if(policy==2) {
      #Policy 2: Max Volume (big M for dummy images replaced with 1)
      lnp_copy = LNP
      bigM = max(lnp_copy$y1)
      lnp_copy$y1[lnp_copy$y1==bigM] = 1
      lnp_copy$y2[lnp_copy$y2==bigM] = 1
      lnp_copy$y3[lnp_copy$y3==bigM] = 1
      lnp_copy$vol = apply(lnp_copy[,c("y1","y2","y3")],1,prod)
      maxvol = max(lnp_copy$vol[lnp_copy$Confirmed==FALSE])
      test = which(lnp_copy$Confirmed==FALSE & lnp_copy$vol==maxvol)[1]
    } else if(policy==3) {
      #Policy 3: Max-Min Weight, Interior-First
      lnp_copy = LNP
      lnp_copy[,c("w1","w2","w3")] <- t(apply(lnp_copy[,c("y1","y2","y3")], 1, kernelWeight))
      lnp_copy$minW = apply(lnp_copy[,c("w1","w2","w3")],1,min)
      maxminW = max(lnp_copy$minW[lnp_copy$Confirmed==FALSE])
      test = which(lnp_copy$Confirmed==FALSE & lnp_copy$minW==maxminW)[1]
    } else if(policy==4) {
      #Policy 4: Min-Min Weight, Boundary-First
      lnp_copy = LNP
      lnp_copy[,c("w1","w2","w3")] <- t(apply(lnp_copy[,c("y1","y2","y3")], 1, kernelWeight))
      lnp_copy$minW = apply(lnp_copy[,c("w1","w2","w3")],1,min)
      minminW = min(lnp_copy$minW[lnp_copy$Confirmed==FALSE])
      test = which(lnp_copy$Confirmed==FALSE & lnp_copy$minW==minminW)[1]
    } else if(policy==5) {
      lnp_copy = LNP
      lnp_copy$SumArea = 0
      trilist <- SingleStepOuterApprox(NDF,LNP,IPsolved)
      outer <- trilist[[2]]
      if(dim(outer)[1]>1) {
        for(j in which(lnp_copy$Confirmed==FALSE)) {
          uniqueinds <- unique(c(lnp_copy$C1[[j]],lnp_copy$C2[[j]],lnp_copy$C3[[j]]))
          uniqueinds <- uniqueinds[uniqueinds>0]
          lnp_copy$SumArea[j] = sum(outer[outer$image%in%uniqueinds,"TotalArea"])
        }
      }
      maxArea = max(lnp_copy$SumArea[lnp_copy$Confirmed==FALSE])
      test = which(lnp_copy$Confirmed==FALSE & lnp_copy$SumArea==maxArea)[1]
    }
    
    
    #solve IP
    inds = which(NDF$y1<LNP$y1[test] & NDF$y2<LNP$y2[test] & NDF$y3<LNP$y3[test])
    IPsolved=IPsolved+1
    if(length(inds)==0) { #LNP is weakly ND
      LNP$Confirmed[test]=TRUE
      LNP$StepConfirmed[test]=IPsolved
    } else { #LNP is dominated
      if(length(inds)==1) {
        #return first id
        foundInd = inds[1]
      } else {
        #return minimal weighted sum
        feas <- NDF[inds,]
        feas$sum = feas$y1+feas$y2+feas$y3
        foundInd = inds[which.min(feas$sum)]
      }
      
      NDF$StepFound[foundInd] = IPsolved
      #previnds <- c(1,2,3,which())
      LNPnew <- KlamrothUpdateMaximal(LNP, NDF[NDF$StepFound>0 & NDF$StepFound<IPsolved,], 
                                      NDF[foundInd,])
      LNPnew <- ComputeKway(LNPnew)
      LNP=LNPnew
    }
    LNPlist[[paste0("step",IPsolved)]] <- LNP
  }
  dataout = data.frame(IPsolved=IPsolved)
  return(list(LNPlist,NDF,dataout))
}

subsetLNP <- function(LNP,im) {
  if(dim(LNP)[1]==0) return(LNP)
  LNPsubset <- LNP[1,]
  for(i in 1:dim(LNP)[1]) {
    if(im%in%LNP$C1[[i]] || im%in%LNP$C2[[i]] || im%in%LNP$C3[[i]]) {
      LNPsubset<-rbind(LNPsubset,LNP[i,])
    }
  }
  LNPsubset<-LNPsubset[-1,]
  return(LNPsubset)
}
kernelWeight <- function(im) {
  dif<-1/(im)
  dif[which(dif==Inf)] <- 0
  sumdif<-sum(dif)
  kw <- as.numeric(c(dif[1]/sumdif, dif[2]/sumdif, dif[3]/sumdif))
  return(kw)
}
pointDistance <- function(p1,p2) {
  dist <- sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2)
  return(as.numeric(dist))
}
triangleArea <- function(p1,p2,p3) {
  a<-pointDistance(p1,p2)
  b<-pointDistance(p1,p3)
  c<-pointDistance(p2,p3)
  if(min(a,b,c)<0.00001) return(0)
  s<-0.5*(a+b+c)
  return(sqrt(max(0,s*(s-a)*(s-b)*(s-c))))
}
polarTransform <- function(center,point) {
  radius=pointDistance(center,point)
  if(radius==0) angle=0
  else {
    deltaX=point[1]-center[1]
    deltaY=point[2]-center[2]
    #if(deltaX==0 && deltaY>0) angle=pi/2
    #if(deltaX==0 && deltaY<0) angle=-pi/2
    #else {
    angle=atan2(deltaY,deltaX)
    #}  
  }
  return(c(radius,angle))
}
computeLNP <- function(NDF,ids) {
  new <- data.frame(y1=max(NDF$y1[NDF$id%in%ids]),
                    y2=max(NDF$y2[NDF$id%in%ids]),
                    y3=max(NDF$y3[NDF$id%in%ids]))
  c1<-NDF$id[which(NDF$id%in%ids & NDF$y1==new$y1)]
  c2<-NDF$id[which(NDF$id%in%ids & NDF$y2==new$y2)]
  c3<-NDF$id[which(NDF$id%in%ids & NDF$y3==new$y3)]
  new$C1 <- I(list(c1))
  new$C2 <- I(list(c2))
  new$C3 <- I(list(c3))
  new$Kway <- length(unique(c(c1,c2,c3)))
  new$Confirmed <- FALSE
  new$StepConfirmed <- 0
  return(new)
}
recursiveLNPcheck <- function(NDF,LNP,maxlambda) {
  kw1=kernelWeight(LNP[1,c("y1","y2","y3")])
  if(max(kw1)<=maxlambda) return(TRUE)
  imgs=unique(c(LNP$C1[[1]],LNP$C2[[1]],LNP$C3[[1]]))
  if(length(imgs)<=2) return(FALSE)
  for(i in imgs) {
    subi = setdiff(imgs,i)
    if(sum(subi>0)==0) next
    lnp=computeLNP(NDF,subi)
    kw2=kernelWeight(lnp[1,c("y1","y2","y3")])
    #edge detection
    for(j in 1:3) {
      if((kw1[j]<maxlambda) != (kw2[j]<maxlambda)) { #one below and one above
         #(LNP[1,paste0("w",j)]>maxlambda && lnp[1,paste0("w",j)]<maxlambda)) {
          mu=(maxlambda-kw1[j]) / (kw2[j]-kw1[j])
          midpoint = (1-mu)*kw1 + (mu)*kw2
          if(max(midpoint)<=maxlambda) return(TRUE)
      }
    }
    
    if(recursiveLNPcheck(NDF,lnp,maxlambda)) return(TRUE)
  }
  return(FALSE)
}
minTchebAugmented <- function(NDF,lambda) {
  NDF<-NDF[-(1:3),]
  NDF$weightedy1 = NDF$y1*lambda[1]
  NDF$weightedy2 = NDF$y2*lambda[2]
  NDF$weightedy3 = NDF$y3*lambda[3]
  NDF$weightedT = apply(NDF[,c("weightedy1","weightedy2","weightedy3")],1,max)
  NDF$weightedT = NDF$weightedT+0.0001*(NDF$y1+NDF$y2+NDF$y3)
  inds=which.min(NDF$weightedT)
  return(NDF$id[inds[1]])
}

TriangleUnique <- function(tridf,newtri) {
  inds1=tridf$tri[which(tridf$lambda1==newtri$lambda1[1] && tridf$lambda2==newtri$lambda2[1])]
  inds2=tridf$tri[which(tridf$lambda1==newtri$lambda1[2] && tridf$lambda2==newtri$lambda2[2])]
  inds3=tridf$tri[which(tridf$lambda1==newtri$lambda1[3] && tridf$lambda2==newtri$lambda2[3])]
  return(length(intersect(intersect(inds1,inds2),inds3))==0)
}

uniqueC <- function(LNP,newlnp) {
  sub=subset(LNP,Kway==newlnp$Kway & y1==newlnp$y1 & y2==newlnp$y2 & y3==newlnp$y3)
  if(dim(sub)[1]==0) return(TRUE)
  for(i in 1:dim(sub)[1]) {
    if(setequal(LNP$C1[[i]],newlnp$C1[[1]]) &&
       setequal(LNP$C2[[i]],newlnp$C2[[1]]) &&
       setequal(LNP$C3[[i]],newlnp$C3[[1]])) return(FALSE)
  }
  return(TRUE)
}

ImpliedLNPs <- function(LNP,NDF,im) {
  n=1
  outList=LNP
  while(n<=dim(LNP)[1]) {
    lnp=LNP[n,]
    C=unique(c(unlist(lnp$C1),unlist(lnp$C2),unlist(lnp$C3)))
    removable=setdiff(C,im)
    for(r in removable) {
      newC = setdiff(removable,r)
      newlnp<-computeLNP(NDF,c(im,setdiff(removable,r)))
      newlnp[,c("w1","w2","w3")] <- kernelWeight(newlnp[,c("y1","y2","y3")])
      output=FALSE
      inds=which(outList$y1==newlnp$y1[1] & outList$y2==newlnp$y2[1] & outList$y3==newlnp$y3[1])
      if(length(inds)==0) {
        output=TRUE
      } else if(newlnp$Kway>max(outList$Kway[inds])) {
        outList <- outList[-inds,]
        output=TRUE
      }
      if(output) outList = rbind(outList,newlnp)
      furtherdecompose=TRUE
      if(newlnp$Kway==2) {
        furtherdecompose=FALSE
      } else if(!uniqueC(LNP,newlnp)) {
        furtherdecompose=FALSE
      }
      if(furtherdecompose) LNP = rbind(LNP,newlnp)
    }
    n=n+1
  }
  return(outList)
}

removeDominating <- function(imgs,ys) {
  i=1
  while(i<=dim(imgs)[1]) {
    candidate<-imgs[i,]
    boolvec=imgs$y1==imgs$y1
    for(j in ys) {
      boolvec=boolvec & (candidate[1,paste0("y",j)]<imgs[,paste0("y",j)])
    }
    dominated<-which(boolvec)
    if(length(dominated)>0) {
      imgs=imgs[-i,]
    } else {
      i=i+1
    }
  }
  return(imgs)
}

CompromisePrimal <- function(baseNDF,maxlambda,initLambda) {
  NDF <- InitializeNDF(baseNDF)
  NDF$StepFound = 0
  NDF$StepFound[1:3] = 1
  
  LNP <- InitializeLNP(NDF)
  LNP <- ComputeKway(LNP)
  LNPlist <- list("step0"=LNP)
  IPsolved = 1
  for(i in 1:dim(initLambda)[1]) {
    lambda=as.numeric(initLambda[i,c("w1","w2","w3")])
    id = minTchebAugmented(NDF,lambda)
    IPsolved=IPsolved+1
    if(NDF$StepFound[NDF$id==id]==0) {
      NDF$StepFound[NDF$id==id] = IPsolved
      LNPnew <- KlamrothUpdateMaximal(LNP, NDF[NDF$StepFound>0 & NDF$StepFound<IPsolved,], 
                                      NDF[NDF$id==id,])
      LNPnew <- ComputeKway(LNPnew)
      LNP=LNPnew 
    }
    LNPlist[[paste0("step",IPsolved)]] <- LNP
  }
  Skipped=0
  while(sum(LNP$Confirmed==FALSE)>0) {
    #policy for choosing LNP
    test = which(LNP$Confirmed==FALSE)[1] #first
    #test if LNP within compromise region
    if(!recursiveLNPcheck(NDF,LNP[test,],maxlambda)) {
      LNP$Confirmed[test]=TRUE
      Skipped=Skipped+1
    } else {
      #solve IP
      inds = which(NDF$y1<LNP$y1[test] & NDF$y2<LNP$y2[test] & NDF$y3<LNP$y3[test])
      IPsolved=IPsolved+1
      if(length(inds)==0) { #LNP is weakly ND
        LNP$Confirmed[test]=TRUE
        LNP$StepConfirmed[test]=IPsolved
      } else { #LNP is dominated
        foundInd = inds[1]
        NDF$StepFound[foundInd] = IPsolved
        #previnds <- c(1,2,3,which())
        LNPnew <- KlamrothUpdateMaximal(LNP, NDF[NDF$StepFound>0 & NDF$StepFound<IPsolved,], 
                                        NDF[foundInd,])
        LNPnew <- ComputeKway(LNPnew)
        LNP=LNPnew
      }
    }
    
    LNPlist[[paste0("step",IPsolved)]] <- LNP
  }
  dataout = data.frame(IPsolved=IPsolved,Skipped=Skipped)
  return(list(LNPlist,NDF,dataout))
}

ComponentPerimeterByC <- function(LNP,NDF,im) {
  y<-as.numeric(NDF[NDF$id==im,c("y1","y2","y3")])
  sub<-subsetLNP(LNP,im) #subset of maximal LNPs
  if(dim(sub)[1]==0) return(NULL)
  sub[,c("w1","w2","w3")] <- t(apply(sub[,c("y1","y2","y3")], 1, kernelWeight))
  
  #generate all LNPs
  all<-ImpliedLNPs(sub,NDF,im)
  #subset for C1 and order
  sub1=subset(all,y1==y[1])
  sub1=removeDominating(sub1,2:3)
  sub1=sub1[order(-sub1$y2,sub1$y3),]
  #subset for C2 and order
  sub2=subset(all,y2==y[2])
  sub2=removeDominating(sub2,c(3,1))
  sub2=sub2[order(-sub2$y3,sub2$y1),]
  #subset for C3 and order
  sub3=subset(all,y3==y[3])
  sub3=removeDominating(sub3,1:2)
  sub3=sub3[order(-sub3$y1,sub3$y2),]
  #return final list
  p<-rbind(sub1,sub2[-1,],sub3[-1,])
  return(p)
}

ComponentInnerApproxByC <- function(LNP,NDF,im) {
  y<-as.numeric(NDF[NDF$id==im,c("y1","y2","y3")])
  sub<-subsetLNP(LNP,im) #subset of maximal LNPs
  if(dim(sub)[1]==0) return(NULL)
  sub[,c("w1","w2","w3")] <- t(apply(sub[,c("y1","y2","y3")], 1, kernelWeight))
  
  #limit to confirmed LNPs
  sub<-subset(sub,Confirmed==TRUE)
  all=ImpliedLNPs(sub,NDF,im)
  #subset for C1 and order
  sub1=subset(all,y1==y[1])
  #sub1=ImpliedLNPs(sub1,NDF,im)
  sub1=removeDominating(sub1,2:3)
  sub1=sub1[order(-sub1$y2,sub1$y3),]
  #subset for C2 and order
  sub2=subset(all,y2==y[2])
  #sub2=ImpliedLNPs(sub2,NDF,im)
  sub2=removeDominating(sub2,c(3,1))
  sub2=sub2[order(-sub2$y3,sub2$y1),]
  #subset for C3 and order
  sub3=subset(all,y3==y[3])
  #sub3=ImpliedLNPs(sub3,NDF,im)
  sub3=removeDominating(sub3,1:2)
  sub3=sub3[order(-sub3$y1,sub3$y2),]
  #return final list
  p<-list(sub1,sub2,sub3)
  return(p)
}

ComputeAllPerimeters <- function(NDF,LNP) {
  ids <- NDF$id
  ids <- ids[ids>0]
  #p<-ComponentPerimeter(LNP,NDF,ids[1])
  p<-ComponentPerimeterByC(LNP,NDF,ids[1])
  if(!is.null(p)) {
    p$Perimeter <- ids[1]
  }
  perimeterset <- p
  if(length(ids)==1) return(perimeterset)
  for(i in ids[-1]) {
    #p<-ComponentPerimeter(LNP,NDF,i)
    p<-ComponentPerimeterByC(LNP,NDF,i)
    if(!is.null(p)) {
      p$Perimeter <- i
      perimeterset <- rbind(perimeterset,p)
    }
  }
  return(perimeterset)
}

TriangulateInner <- function(NDF,LNP) {
  ids <- NDF$id
  ids <- ids[ids>0]
  triangledf <- data.frame(tri=0,image=0,lambda1=0,lambda2=0)
  areadf <- data.frame(tri=0,image=0,tarea=0)
  tri=1
  p<-ComponentInnerApproxByC(LNP,NDF,ids[1])
  if(!is.null(p)) {
    center <- as.numeric(NDF[NDF$id==ids[1],c("w1","w2")])
    for(i in 1:3) {
      subdf <- p[[i]]
      if(dim(subdf)[1]<2) next
      for(j in 1:(dim(subdf)[1]-1)) {
        p1<-as.numeric(subdf[j,c("w1","w2")])
        p2<-as.numeric(subdf[j+1,c("w1","w2")])
        
        newtri <- data.frame(tri=tri, image=ids[1], 
                             lambda1=c(center[1],p1[1],p2[1]),
                             lambda2=c(center[2],p1[2],p2[2]))
        triangledf <- rbind(triangledf,newtri)
        areadf <- rbind(areadf,as.numeric(c(tri,ids[1],triangleArea(center,p1,p2))))
        tri=tri+1
      }
    }
  }
  if(length(ids)==1) {
    triangledf <- triangledf[-1,]
    areadf <- areadf[-1,]
    return(list(triangledf,areadf))
  }
  for(id in ids[-1]) {
    p<-ComponentInnerApproxByC(LNP,NDF,id)
    if(!is.null(p)) {
      center <- as.numeric(NDF[NDF$id==id,c("w1","w2")])
      for(i in 1:3) {
        subdf <- p[[i]]
        if(dim(subdf)[1]<2) next
        for(j in 1:(dim(subdf)[1]-1)) {
          p1<-as.numeric(subdf[j,c("w1","w2")])
          p2<-as.numeric(subdf[j+1,c("w1","w2")])
          
          newtri <- data.frame(tri=tri, image=id, 
                               lambda1=c(center[1],p1[1],p2[1]),
                               lambda2=c(center[2],p1[2],p2[2]))
          triangledf <- rbind(triangledf,newtri)
          areadf <- rbind(areadf,as.numeric(c(tri,id,triangleArea(center,p1,p2))))
          tri=tri+1
        }
      }
      
    }
  }
  triangledf <- triangledf[-1,]
  areadf <- areadf[-1,]
  return(list(triangledf,areadf))
}

TriangulatePerimeterSet <- function(NDF,perimeterset) {
  ##Triangulate data
  triangledf <- data.frame(tri=0,image=0,lambda1=0,lambda2=0)
  areadf <- data.frame(tri=0,image=0,tarea=0)
  linesegdf <- data.frame(startx=0.0,starty=0.0,endx=0.0,endy=0.0,image=0)
  
  tri<-1
  for(i in unique(perimeterset$Perimeter)) {
    subdf <- subset(perimeterset,Perimeter==i)
    center <- as.numeric(NDF[NDF$id==i,c("w1","w2")])
    for(j in 1:dim(subdf)[1]) {
      p1<-as.numeric(subdf[j,c("w1","w2")])
      if(j<dim(subdf)[1]) {
        p2<-as.numeric(subdf[j+1,c("w1","w2")])
      } else {
        p2<-as.numeric(subdf[1,c("w1","w2")])
      }
      
      newtri <- data.frame(tri=tri, image=i, 
                           lambda1=c(center[1],p1[1],p2[1]),
                           lambda2=c(center[2],p1[2],p2[2]))
      if(TriangleUnique(triangledf,newtri)) {
        triangledf <- rbind(triangledf,newtri)
        areadf <- rbind(areadf,as.numeric(c(tri,i,triangleArea(center,p1,p2))))
        linesegdf <- rbind(linesegdf,as.numeric(c(p1[1],p1[2],p2[1],p2[2],i)))
        tri=tri+1
      } 
      
    }
  }
  triangledf <- triangledf[-1,]
  areadf <- areadf[-1,]
  linesegdf <- linesegdf[-1,]
  return(list(triangledf,areadf,linesegdf))
}

SummarizeOuterApprox <- function(NDF,LNPlist) {
  NDF[,c("w1","w2","w3")] <- t(apply(NDF[,c("y1","y2","y3")], 1, kernelWeight))
  
  steps <- NDF$StepFound
  steps <- steps[steps>1]
  steps <- steps[order(steps)]
  
  OuterApprox <- data.frame(Step=0,image=0,TotalArea=0)
  for(s in steps) {
    LNP <- LNPlist[[paste0("step",s)]]
    LNP[,c("w1","w2","w3")] <- t(apply(LNP[,c("y1","y2","y3")], 1, kernelWeight))
    subNDF <- NDF[NDF$StepFound<=s,]
    perimset<-ComputeAllPerimeters(subNDF,LNP)
    trilist <- TriangulatePerimeterSet(subNDF,perimset)
    #triangledf <- trilist[[1]]
    areadf <- trilist[[2]]
    #linesegdf <- trilist[[3]]
    AA <- areadf %>% group_by(image) %>% summarise(TotalArea = sum(tarea))
    AA$Step=s
    OuterApprox = rbind(OuterApprox,AA)
  }
  OuterApprox = OuterApprox[-1,]
  OuterApprox$image <- as.factor(OuterApprox$image)
  return(OuterApprox)
}

SingleStepOuterApprox <- function(NDF,LNP,step) {
  NDF[,c("w1","w2","w3")] <- t(apply(NDF[,c("y1","y2","y3")], 1, kernelWeight))
  
  #OuterApprox <- data.frame(Step=0,image=0,TotalArea=0)
  #for(s in steps) {
    #LNP <- LNPlist[[paste0("step",s)]]
    LNP[,c("w1","w2","w3")] <- t(apply(LNP[,c("y1","y2","y3")], 1, kernelWeight))
    subNDF <- NDF[NDF$StepFound<=step,]
    perimset<-ComputeAllPerimeters(subNDF,LNP)
    trilist <- TriangulatePerimeterSet(subNDF,perimset)
    areadf <- trilist[[2]]
    AA <- areadf %>% group_by(image) %>% summarise(TotalArea = sum(tarea))
    AA$image <- as.factor(AA$image)
    trilist[[2]] <- AA
    #AA$Step=s
    #OuterApprox = AA
  #}
  #OuterApprox = OuterApprox[-1,]
  
  return(trilist)
}

SingleStepInnerApprox <- function(NDF,LNP,step) {
  NDF[,c("w1","w2","w3")] <- t(apply(NDF[,c("y1","y2","y3")], 1, kernelWeight))
  
  LNPconfirmed <- subset(LNP,Confirmed==TRUE)
  if(dim(LNPconfirmed)[1]==0) next
  trilist <- TriangulateInner(subNDF,LNPconfirmed) 
  areadf <- trilist[[2]]
  AA <- areadf %>% group_by(image) %>% summarise(TotalArea = sum(tarea))
  AA$image <- as.factor(AA$image)
  trilist[[2]] <- AA
  
  return(trilist)
}

OuterMinusInner <- function(InnerApprox,OuterApprox) {
  AllApprox <- data.frame(Image=0,Step=0,Outer=0,Inner=0)
  nsteps=max(InnerApprox$Step,OuterApprox$Step)
  for(i in levels(OuterApprox$image)) {
    firststep=min(OuterApprox$Step[OuterApprox$image==i])
    combined <- data.frame(Image=i,Step=firststep:nsteps,Outer=0,Inner=0)
    for(s in combined$Step) {
      ind=which(combined$Step==s)
      subdf <- subset(OuterApprox,Step<=s & image==i)
      if(dim(subdf)[1]>0) combined$Outer[ind]=subdf$TotalArea[subdf$Step==max(subdf$Step)]
      subdf <- subset(InnerApprox,Step<=s & image==i)
      if(dim(subdf)[1]>0) combined$Inner[ind]=subdf$TotalArea[subdf$Step==max(subdf$Step)]
    }
    AllApprox<-rbind(AllApprox,combined)
  }
  AllApprox<-AllApprox[-1,]
  AllApprox$OuterMinusInner <- AllApprox$Outer-AllApprox$Inner
  return(AllApprox)
}

SummarizeInnerOuterApprox <- function(NDF,LNPlist) {
  NDF[,c("w1","w2","w3")] <- t(apply(NDF[,c("y1","y2","y3")], 1, kernelWeight))
  
  #steps <- NDF$StepFound
  steps <- as.numeric(substr(names(LNPlist),5,10))
  steps <- steps[steps>1]
  steps <- steps[order(steps)]
  
  OuterApprox <- data.frame(Step=0,image=0,TotalArea=0)
  InnerApprox <- data.frame(Step=0,image=0,TotalArea=0)
  for(s in steps) {
    LNP <- LNPlist[[paste0("step",s)]]
    LNP[,c("w1","w2","w3")] <- t(apply(LNP[,c("y1","y2","y3")], 1, kernelWeight))
    #Outer
    subNDF <- NDF[NDF$StepFound<=s,]
    perimset<-ComputeAllPerimeters(subNDF,LNP)
    trilist <- TriangulatePerimeterSet(subNDF,perimset)
    outertriangledf <- trilist[[1]]
    areadf <- trilist[[2]]
    AA <- areadf %>% group_by(image) %>% summarise(TotalArea = sum(tarea))
    AA$Step=s
    OuterApprox = rbind(OuterApprox,AA)
    #Inner
    LNPconfirmed <- subset(LNP,Confirmed==TRUE)
    if(dim(LNPconfirmed)[1]==0) next
    trilist <- TriangulateInner(subNDF,LNPconfirmed) 
    #perimset<-ComputeAllPerimeters(subNDF,LNPconfirmed)
    #trilist <- TriangulatePerimeterSet(subNDF,perimset)
    innertriangledf <- trilist[[1]]
    if(dim(innertriangledf)[1]==0) next
    areadf <- trilist[[2]]
    AA <- areadf %>% group_by(image) %>% summarise(TotalArea = sum(tarea))
    AA$Step=s
    InnerApprox = rbind(InnerApprox,AA)
    
  }
  OuterApprox = OuterApprox[-1,]
  InnerApprox = InnerApprox[-1,]
  OuterApprox$image <- as.factor(OuterApprox$image)
  InnerApprox$image <- as.factor(InnerApprox$image)
  AllApprox <- OuterMinusInner(InnerApprox,OuterApprox)
  return(list(InnerApprox,OuterApprox,AllApprox))
}

representationMetrics <- function(NDF) {
  metrics <- data.frame(Step=sort(unique(NDF$StepFound)),
                        Dmean=0,
                        Dmax=0)
  metrics<-metrics[-1,]
  true <- NDF[NDF$id>0,1:5]
  w <- 1/sapply(true[,c("y1","y2","y3")],max)
  #w = 1/c(max(true$y1),max(true$y2),max(true$y3))
  for(s in sort(true$Step)) {
    rep <- subset(true,StepFound<=s)
    minwdist = c()
    for(id in true$id) {
      img <- true[true$id==id,c("y1","y2","y3")]
      rep$y1diff <- w[1]*abs(rep$y1-img$y1)
      rep$y2diff <- w[2]*abs(rep$y2-img$y2)
      rep$y3diff <- w[3]*abs(rep$y3-img$y3)
      rep$maxWdiff <- apply(rep[,c("y1diff","y2diff","y3diff")], 1, max)
      minwdist = c(minwdist,min(rep$maxWdiff))
    }
    metrics$Dmean[metrics$Step==s] = mean(minwdist)
    metrics$Dmax[metrics$Step==s] = max(minwdist)
  }
  return(metrics)
}




