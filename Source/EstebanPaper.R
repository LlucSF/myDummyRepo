library(ggplot2)
THSimg <- rMSI::LoadMsiData("/home/lluc/Documents/PhDprogram/Data/BrainTHS/20151001-Brain2_CS7_Au-proc.tar")
THSmtx <- rMSIproc::LoadPeakMatrix("/home/lluc/Documents/PhDprogram/Data/BrainTHS/postcode2/mergeddata-peaks-Au.zip")

THSmtx3 <- THSmtx
THSmtx3$intensity <- THSmtx3$intensity[(sum(THSmtx3$numPixels[1:2])+1):sum(THSmtx3$numPixels),]
THSmtx3$SNR <- THSmtx3$SNR[(sum(THSmtx3$numPixels[1:2])+1):sum(THSmtx3$numPixels),]
THSmtx3$area <- THSmtx3$area[(sum(THSmtx3$numPixels[1:2])+1):sum(THSmtx3$numPixels),]
THSmtx3$pos <- THSmtx3$pos[(sum(THSmtx3$numPixels[1:2])+1):sum(THSmtx3$numPixels),]
THSmtx3$posMotors <- THSmtx3$posMotors[(sum(THSmtx3$numPixels[1:2])+1):sum(THSmtx3$numPixels),]
THSmtx3$normalizations <- THSmtx3$normalizations[(sum(THSmtx3$numPixels[1:2])+1):sum(THSmtx3$numPixels),]
THSmtx3$names <- THSmtx3$names[-(1:2)] 
THSmtx3$uuid <- THSmtx3$uuid[-(1:2)] 
THSmtx3$numPixels <- THSmtx3$numPixels[-(1:2)] 

clus <- kmeans(THSmtx3$intensity,centers = 5,iter.max = 100)
rMSIproc::plotClusterImage(THSmtx3,clusters = clus$cluster)
THSmtx3$intensity  <- THSmtx3$intensity[!(clus$cluster==clus$cluster[1]),] 
THSmtx3$area  <- THSmtx3$area[!(clus$cluster==clus$cluster[1]),] 
THSmtx3$pos  <- THSmtx3$pos[!(clus$cluster==clus$cluster[1]),] 
THSmtx3$SNR  <- THSmtx3$SNR[!(clus$cluster==clus$cluster[1]),] 
THSmtx3$posMotors  <- THSmtx3$posMotors[!(clus$cluster==clus$cluster[1]),] 
THSmtx3$normalizations  <- THSmtx3$normalizations[!(clus$cluster==clus$cluster[1]),]
THSmtx3$numPixels <- nrow(THSmtx3$intensity)

clus2 <- kmeans(scale(THSmtx3$intensity),centers = 7,iter.max = 100)
rMSIproc::plotClusterImage(THSmtx3,clusters = clus2$cluster)

rMSI::plotMassImageByPeak(THSimg,mass.peak = c(336.769,330.761,298.811))

################################### FIG S1 ###################################################


ClustersZeroperCent <- data.frame(Cluster1 = rep(0, times = ncol(THSmtx3$intensity)),
                                  Cluster2 = rep(0, times = ncol(THSmtx3$intensity)),
                                  Cluster3 = rep(0, times = ncol(THSmtx3$intensity)),
                                  Cluster4 = rep(0, times = ncol(THSmtx3$intensity)),
                                  Cluster5 = rep(0, times = ncol(THSmtx3$intensity)),
                                  Cluster6 = rep(0, times = ncol(THSmtx3$intensity)),
                                  Cluster7 = rep(0, times = ncol(THSmtx3$intensity))
                                  )

for(j in 1:ncol(ClustersZeroperCent))
{
  for(i in 1:ncol(THSmtx3$intensity))
  {
    ClustersZeroperCent[i,j] <- (length((which(THSmtx3$intensity[(which(clus2$cluster==j)),i]/THSmtx3$normalizations$AcqTic[(which(clus2$cluster==j))]<(2.5*10^-5))))*100)/length((which(clus2$cluster==j)))
  }
}


tidyClus <- tidyr::gather(ClustersZeroperCent)
tidyClus$index <- rep(1:length(THSmtx3$mass),times = 7) 

ggplot(data = tidyClus) + geom_point(mapping = aes(y = value, x = index), pch = 18, color = "blue") + theme_bw() + facet_grid(.~key) +
  labs(y = "Percentage of pixels with zero values", x = "Ions index")

################################### FIG S2 ###################################################
dist1 <- rnorm(mean = 0.04,sd =  0.001,n = 500)
dist1 <- c(dist1,rep(0,times = 500))
dist2 <- rnorm(mean = 0.025,sd = 0.001,n = 900)
dist2 <- c(dist2,rep(0,times = 100))

FigS2 <- data.frame(A = dist1, B = dist2)
FigS2 <- tidyr::gather(FigS2)
FigS2$key[which(FigS2$key=="A")] <- "A: mean = 0.04" 
FigS2$key[which(FigS2$key=="B")] <- "B: mean = 0.025" 

ggplot(FigS2) + geom_histogram(mapping = aes(x = value,fill = key, y = stat(count/1000)),bins = 100,position ="identity") + 
  theme_bw() + 
  labs(x = "Magnitude", y = "Probability", fill = "Variable", subtitle = "")

t.test(c(dist1,dist2))
wilcox.test(c(dist1,dist2))
kruskal.test(c(dist1,dist2),g = rep(1:2,each =1000))

dist1 <- rnorm(mean = 0.04,sd =  0.001,n = 1000)
dist2 <- rnorm(mean = 0.025,sd = 0.001,n = 1000)

################################### FIG S3 ###################################################
THSmtx4 <- THSmtx3
THSmtx4$intensity <- (THSmtx4$intensity/THSmtx4$normalizations$AcqTic)
pvalue <- rep(NA,times = length(THSmtx3$mass))
log2FC <- rep(NA,times = length(THSmtx3$mass))
for(i in 1:ncol(THSmtx3$intensity))
{
  pvalue[i] <- wilcox.test(x = THSmtx4$intensity[clus2$cluster==5,i], y = THSmtx4$intensity[clus2$cluster==3,i])$p.value
  log2FC[i] <- log2(mean(THSmtx4$intensity[clus2$cluster==5,i])/mean(THSmtx4$intensity[clus2$cluster==3,i]))  
}

FigS3 <- data.frame(pvalue = -log(pvalue), log2FC = log2FC,mass = signif(THSmtx4$mass,6) )

ggplot(data = FigS3, mapping = aes(x = log2FC, y = pvalue)) +
  geom_point(color = "blue", alpha = 0.6) + 
  geom_text(data = subset(FigS3,(pvalue>400) & (abs(log2FC)>1)),mapping = aes(label = mass),nudge_y = 10,nudge_x = -0.05,size = 3,check_overlap = T) +
  geom_hline(yintercept = 400,color = "red") + geom_vline(xintercept = 1,color = "red") + geom_vline(xintercept = -1,color = "red") +
  theme_bw() + labs(y = "p-value", x = "Fold change")

################################### FIG S4 ###################################################

FigS4 <- data.frame(Int = THSmtx4$intensity[,which.min(abs(THSmtx4$mass-836.6))], clus = clus2$cluster)
  
ggplot(data = subset(FigS4, clus == 5 | clus == 3)) + 
  geom_histogram(mapping = aes(x = Int,
                               fill = as.factor(clus),
                               y = stat(count)),
                               bins = 75,
                               position ="identity",
                               color = "grey",
                               alpha = 0.6
                 ) + 
  theme_bw() + 
  scale_fill_discrete(labels =c("1","5")) +
  labs(x = "Magnitude", 
       y = "Count", 
       fill = "Cluster", 
       title = paste("Ion: ",signif(THSmtx4$mass[which.min(abs(THSmtx4$mass-836.6))],7)," m/z")
       )

################################### FIG S5 ###################################################

FigS4 <- data.frame(Int = THSmtx4$intensity[,which.min(abs(THSmtx4$mass-146.9796))], clus = clus2$cluster)

ggplot(data = subset(FigS4, clus == 5 | clus == 3)) + 
  geom_histogram(mapping = aes(x = Int,
                               fill = as.factor(clus),
                               y = stat(count)),
                 bins = 75,
                 position ="identity",
                 color = "grey",
                 alpha = 0.6
  ) + 
  theme_bw() + 
  scale_fill_discrete(labels =c("1","5")) +
  labs(x = "Magnitude", 
       y = "Count", 
       fill = "Cluster", 
       title = paste("Ion: ",signif(THSmtx4$mass[which.min(abs(THSmtx4$mass-146.9796))],7)," m/z")
  )

################################### FIG S6 ###################################################

FigS4 <- data.frame(Int = THSmtx4$intensity[,which.min(abs(THSmtx4$mass-824.6287))], clus = clus2$cluster)

ggplot(data = subset(FigS4, clus == 5 | clus == 3)) + 
  geom_histogram(mapping = aes(x = Int,
                               fill = as.factor(clus),
                               y = stat(count)),
                 bins = 75,
                 position ="identity",
                 color = "grey",
                 alpha = 0.6
  ) + 
  theme_bw() + 
  scale_fill_discrete(labels =c("1","5")) +
  labs(x = "Magnitude", 
       y = "Count", 
       fill = "Cluster", 
       title = paste("Ion: ",signif(THSmtx4$mass[which.min(abs(THSmtx4$mass-824.6287))],7)," m/z")
  )


################################### FIG S7 ###################################################

FigS4 <- data.frame(Int = THSmtx4$intensity[,which.min(abs(THSmtx4$mass-478.3167))], clus = clus2$cluster)

ggplot(data = subset(FigS4, clus == 5 | clus == 3)) + 
  geom_histogram(mapping = aes(x = Int,
                               fill = as.factor(clus),
                               y = stat(count)),
                 bins = 75,
                 position ="identity",
                 color = "grey",
                 alpha = 0.6
  ) + 
  theme_bw() + 
  scale_fill_discrete(labels =c("1","5")) +
  labs(x = "Magnitude", 
       y = "Count", 
       fill = "Cluster", 
       title = paste("Ion: ",signif(THSmtx4$mass[which.min(abs(THSmtx4$mass-478.3167))],7)," m/z")
  )

