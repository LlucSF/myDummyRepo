
###############################################################################################################################
plotMassScatter <- function(PeakMtx,colA, colB, slope = NA,Normalization = NA)
{
  
  if(is.na(Normalization))
  {
    df <- data.frame(A = PeakMtx$intensity[,colA], B = PeakMtx$intensity[,colB])
  }else
   {
    df <- data.frame(A = (PeakMtx$intensity/PeakMtx$normalizations[[Normalization]])[,colA], B = (PeakMtx$intensity/PeakMtx$normalizations[[Normalization]])[,colB])
   }
  
  ggplot2::ggplot(data = df) +
  ggplot2::geom_point(mapping = ggplot2::aes(x = A, y = B), size = 0.5) + 
  ggplot2::xlab(signif(PeakMtx$mass[colA], digits = 7)) +
  ggplot2::ylab(signif(PeakMtx$mass[colB], digits = 7)) +
  ggplot2::geom_abline(mapping = ggplot2::aes(slope = slope, intercept = 0), colour = "red") +
  ggplot2::theme_bw()
}
plotMassScatter(R7mtx,colA =  49,colB =  53,slope =  1.643)

###############################################################################################################################
NumberOfC <- function(str)
{
  tmp <- strsplit(str,"H")[[1]][1]
  tmp <- as.numeric(strsplit(tmp,"C")[[1]][2])
  return(tmp)
}

###############################################################################################################################
plotScoreSummary <- function(iso,scoreCutoff)   
{
  isoRF <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isoRM <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isoRI <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isotopeNames <- paste("M+",as.character(1:(length(iso)-1)),sep = "")
  colnames(isoRF) <- isotopeNames
  colnames(isoRM) <- isotopeNames
  colnames(isoRI) <- isotopeNames
  
  for( i in 1:(length(iso)-1) )
  {
    for( j in 1:length(iso[[i]]) )
    {
      if(length(iso[[i]][[j]])==1)
      {
        indx <- 1
      }
      else 
      {
        index <- which.min(abs(iso[[i]][[j]][4,which(iso[[i]][[j]][5,]>scoreCutoff)]))
      }
      
      indx <- which.min(iso[[i]][[j]][4,])
      isoRF[[i]][j] <- iso[[i]][[j]][5,indx]
      isoRM[[i]][j] <- iso[[i]][[j]][6,indx]
      isoRI[[i]][j] <- iso[[i]][[j]][7,indx]
    }
  }
  
  isoRF <- tidyr::gather(isoRF)
  isoRF$score <- rep("Final score",times = nrow(isoRF))
  obs <- c()
  for( i in 1:(length(iso)-1) )
  {
    obs <- c(obs,rep(sum(!is.na(isoRF$value[isoRF$key==isotopeNames[i]])),times = length(iso[[1]])))
  }
  isoRF$obs <- obs
  
  isoRM <- tidyr::gather(isoRM)
  isoRM$score <- rep("Morphology score",times = nrow(isoRM))
  isoRM$obs <- obs
  
  isoRI <- tidyr::gather(isoRI)
  isoRI$score <- rep("Intensity score",times = nrow(isoRI))
  isoRI$obs <- obs
  
  isoRe <- data.frame(isotope = c(isoRF$key,isoRM$key,isoRI$key),
                      scoreName = c(isoRF$score,isoRM$score,isoRI$score), 
                      scoreValue= as.numeric(c(isoRF$value,isoRM$value,isoRI$value)),
                      numberObservation = as.factor(c(isoRF$obs,isoRM$obs,isoRI$obs))
                      )
  
  isoRe$accepted <- rep(isoRF$value>scoreCutoff, times = 3)
  isoRe$accepted <- as.numeric(isoRe$accepted) 
  Acc <-length(which(isoRe$accepted[1:(nrow(isoRe)/3)]==1)) 
  Rej <-length(which(isoRe$accepted[1:(nrow(isoRe)/3)]==0)) 
  isoRe$accepted[is.na(isoRe$accepted)] <- 0


  
  isoRe$accepted <- as.factor(isoRe$accepted)
  levels(isoRe$accepted) <- c(paste("Rejected (",Rej,")",sep = ""),
                              paste("Accepted (",Acc,")",sep = ""))
  ##TODO el rejected esta malament! este que descomptar els que eren NA
  
  g <-  ggplot2::ggplot(data = isoRe,
                        mapping = ggplot2::aes(x = isotope,
                                               y = scoreValue, 
                                               colour = numberObservation)
                        ) + 
        ggplot2::geom_violin(na.rm = TRUE,
                             fill = "bisque1",
                             adjust = .8,
                             draw_quantiles = c(0.25, 0.5, 0.75),
                             scale="width"
                             ) + 
        ggplot2::geom_jitter(mapping = ggplot2::aes(shape = accepted),
                             size = 1.2,
                             na.rm = T
                             ) +
        ggplot2::scale_shape_manual(values = c(4,19)) +
        ggplot2::facet_grid(.~scoreName) +
        ggplot2::geom_hline(data = subset(isoRe, scoreName == "Final score"), 
                            mapping = ggplot2::aes(yintercept = scoreCutoff),
                            color = "gray30", 
                            linetype = 5,
                            size = 0.5
                            ) +
        ggplot2::labs(y = "Score value", 
                      x = "Isotopic ion", 
                      colour = "Tested ions",
                      shape = "Test results"
                      ) + 
        ggplot2::theme_bw() + 
        ggplot2::guides(colour = ggplot2::guide_legend(order = 1,reverse = T),
                        shape = ggplot2::guide_legend(order = 2, reverse = T)
                        ) +
        ggplot2::ggtitle(paste("Score cutoff:",scoreCutoff)) 

  return(g)
}
plotScoreSummary(iso,0.7)
###############################################################################################################################
mtx <- R7mtx
img <- R7img
scoreCutoff <- 0.7
RunIsoAndPlotScoreSummary <- function(mtx,img,scoreCutoff,normalization = "raw",IsoNumber,Scans)   
{
  iso <- rMSIproc::Deisotoping(PeakMtx = mtx,MassChannelsVec = img$mass,IsoNumber = IsoNumber,Scans = Scans,Score = scoreCutoff,Normalization = normalization)
  
  isoRF <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isoRM <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isoRI <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isoAI <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isoMS <- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isoPPM<- as.data.frame(matrix(nrow = length(iso[[1]]), ncol = length(iso)-1))
  isotopeNames <- paste("M+",as.character(1:(length(iso)-1)),sep = "")
  colnames(isoRF) <- isotopeNames
  colnames(isoRM) <- isotopeNames
  colnames(isoRI) <- isotopeNames
  colnames(isoAI) <- isotopeNames
  colnames(isoMS) <- isotopeNames
  colnames(isoPPM)<- isotopeNames
  
  for( i in 1:(length(iso)-1) )
  {
    for( j in 1:length(iso[[i]]) )
    {
      if(length(iso[[i]][[j]])==1)
      {
        indx <- 1
      }
      else 
      {
        index <- which.min(abs(iso[[i]][[j]][4,which(iso[[i]][[j]][5,]>scoreCutoff)]))
      }
      
      indx <- which.min(iso[[i]][[j]][4,])
      isoRF[[i]][j] <- iso[[i]][[j]][5,indx]
      isoRM[[i]][j] <- iso[[i]][[j]][6,indx]
      isoRI[[i]][j] <- iso[[i]][[j]][7,indx]
      isoAI[[i]][j] <- mean(mtx$intensity[,iso[[i]][[j]][3,indx]])
      isoMS[[i]][j] <- mtx$mass[iso[[i]][[j]][2,indx]]
      isoPPM[[i]][j]<- abs(iso[[i]][[j]][4,indx])
    }
  }
  
  isoRF <- tidyr::gather(isoRF)
  isoRF$score <- rep("Final score",times = nrow(isoRF))
  obs <- c()
  for(i in 1:(length(iso)-1))
  {
    obs <- c(obs,rep(sum(!is.na(isoRF$value[isoRF$key==isotopeNames[i]])),times = length(iso[[1]])))
  }
  isoRF$obs <- obs
  
  isoRM <- tidyr::gather(isoRM)
  isoRM$score <- rep("Morphology score",times = nrow(isoRM))
  isoRM$obs <- obs
  
  isoRI <- tidyr::gather(isoRI)
  isoRI$score <- rep("Intensity score",times = nrow(isoRI))
  isoRI$obs <- obs
  
  isoAI <- tidyr::gather(isoAI)
  isoMS <- tidyr::gather(isoMS)
  isoPPM<- tidyr::gather(isoPPM)
  
  isoRe <- data.frame(isotope = c(isoRF$key,isoRM$key,isoRI$key),
                      scoreName = c(isoRF$score,isoRM$score,isoRI$score), 
                      scoreValue= as.numeric(c(isoRF$value,isoRM$value,isoRI$value)),
                      numberObservation = as.factor(c(isoRF$obs,isoRM$obs,isoRI$obs)),
                      averageIntensity = rep(isoAI$value,times = 3),
                      monoisotopicmass = rep(isoMS$value,times = 3),
                      ppmError = rep(isoPPM$value,times = 3)
  )
  
  isoRe$accepted <- rep(isoRF$value>scoreCutoff, times = 3)
  isoRe$accepted <- as.numeric(isoRe$accepted) 
  Acc <-length(which(isoRe$accepted[1:(nrow(isoRe)/3)]==1)) 
  Rej <-length(which(isoRe$accepted[1:(nrow(isoRe)/3)]==0)) 
  isoRe$accepted[is.na(isoRe$accepted)] <- 0
  isoRe$accepted <- as.factor(isoRe$accepted)
  levels(isoRe$accepted) <- c(paste("Rejected (",Rej,")",sep = ""),
                              paste("Accepted (",Acc,")",sep = ""))

  
  g <-  ggplot2::ggplot() + 
    ggplot2::geom_violin(data = isoRe, 
                         mapping = ggplot2::aes(x = isotope,
                                                y = scoreValue,
                                                fill = numberObservation
                                                ),
                         na.rm = TRUE,
                         adjust = .8,
                         draw_quantiles = c(0.25, 0.5, 0.75),
                         scale="count"
    ) + 
    ggplot2::geom_jitter(data = isoRe,
                         mapping = ggplot2::aes(x = isotope,
                                                y = scoreValue,
                                                shape = accepted,
                                                color = ppmError
                                                ),
                         size = 1.2,
                         na.rm = T,
                         inherit.aes = F
    ) +
    ggplot2::scale_shape_manual(values = c(4,19)) +
    ggplot2::facet_grid(.~scoreName) +
    ggplot2::geom_hline(data = subset(isoRe, scoreName == "Final score"), 
                        mapping = ggplot2::aes(yintercept = scoreCutoff),
                        color = "gray30", 
                        linetype = 5,
                        size = 0.5
    ) +
    ggplot2::labs(y = "Score value", 
                  x = "Isotopic ion", 
                  color = "ppm error",
                  shape = "Test results",
                  fill = "Number of ions"
    ) + 
    ggplot2::theme_bw() + 
    ggplot2::guides(color = ggplot2::guide_colorbar(order = 1),
                    shape = ggplot2::guide_legend(order = 2, reverse = T),
                    fill = ggplot2::guide_legend(order = 3, reverse = T)
    ) +
    ggplot2::scale_color_gradientn(colours = rainbow(7)) +
    ggplot2::ggtitle(paste("Score cutoff:",scoreCutoff)) 
  
  print(g)
  return(iso)
}
RunIsoAndPlotScoreSummary(R7mtx, R7img, 0.9)

##############################################################################################################################
ionSlope <- function(PeakMtx,colA,colB)
{
  return(lm(PeakMtx$intensity[,colA]~PeakMtx$intensity[,colB])$coefficients[[2]])
}

